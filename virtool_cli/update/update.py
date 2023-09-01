import json
from pathlib import Path
import asyncio
import aiofiles
from typing import Optional
import logging
from structlog import BoundLogger
from urllib.error import HTTPError

from virtool_cli.utils.logging import base_logger
from virtool_cli.utils.ref import get_otu_paths, get_isolate_paths
from virtool_cli.utils.hashing import generate_hashes, get_unique_ids
from virtool_cli.utils.ncbi import fetch_upstream_accessions, fetch_upstream_records, NCBI_REQUEST_INTERVAL
from virtool_cli.utils.evaluate import evaluate_sequence, get_qualifiers
from virtool_cli.utils.format import format_sequence
from virtool_cli.accessions.helpers import search_by_otu_id

DEFAULT_INTERVAL = 0.001


def run(
    otu_path: Path, catalog: Path, debugging = False
):
    """
    """
    filter_class = logging.DEBUG if debugging else logging.INFO
    logging.basicConfig(
        format="%(message)s",
        level=filter_class,
    )
    
    [ _, otu_id ] = otu_path.name.split('--')
    listing_path = search_by_otu_id(otu_id, catalog)

    asyncio.run(update_otu(otu_path, listing_path))

async def update_otu(otu_path: Path, listing_path: Path):
    """
    """
    listing = json.loads(listing_path.read_text())

    # extract taxon ID and _id hash from listing filename
    [ taxid, otu_id ] = (listing_path.stem).split('--')
    
    logger = base_logger.bind(taxid=taxid, otu_id=otu_id)

    record_data = await request_new_records(listing)
    if not record_data:
        return
    
    otu_updates = await process_records(
        records=record_data, listing=listing, 
        auto_evaluate=True, logger=logger)
    if not otu_updates:
        return
    
    unique_iso, unique_seq = await get_unique_ids([otu_path])
    
    await write_records(
        otu_path, otu_updates, 
        unique_iso, unique_seq,
        logger=logger
    )
    
async def request_new_records(
    listing: dict, 
    logger: BoundLogger = base_logger
) -> list:
    """
    """
    try:
        new_accessions = await fetch_upstream_accessions(
            listing['taxid'], listing, logger=logger)
        
        await asyncio.sleep(NCBI_REQUEST_INTERVAL)

        record_data = await fetch_upstream_records(new_accessions, logger)
        
        await asyncio.sleep(NCBI_REQUEST_INTERVAL)
    
    except HTTPError as e:
        logger.error(e)
        return
    
    except Exception as e:
        logger.exception(e)
        return []
    
    return record_data

async def process_records(
    records, listing: dict,
    auto_evaluate: bool = True, 
    logger: BoundLogger = base_logger
):
    """
    Evaluates whether the record should be included
    """
    filter_set = set(listing['accessions']['included'])
    filter_set.update(listing['accessions']['excluded'])

    if auto_evaluate:
        
        if not listing.get('schema', []):
            logger.warning('Missing schema. moving on...')
            return False
        
        required_parts = dict()
        if len(listing['schema']) < 1:
            for part in listing['schema']:
                if part['required']:
                    part_length = part.get('length', -1)
                    if part_length < 0:
                        logger.warning(
                            f"{part['name']} lacks a length listing.")
                        
                    required_parts[part['name']] = part_length
                    
        else:
            part = listing['schema'][0]
            if (part_name := part.get('name', None)) is None:
                part_name = 'unlabelled'
            
            part_length = part.get('length', -1)
            if part_length < 0:
                logger.warning(
                    f"{part['name']} lacks a length listing.")
            
            required_parts[part_name] = part_length

        logger = logger.bind(required_parts=required_parts)
    
    auto_excluded = []
    otu_updates = []
    for seq_data in records:
        [ accession, _ ] = (seq_data.id).split('.')
        seq_qualifier_data = get_qualifiers(seq_data.features)

        if accession in filter_set:
            logger.debug('Accession already exists', accession=seq_data.id)
            return False
        
        if auto_evaluate:
            if '_' in seq_data.id:
                logger.warning(
                    'This is a RefSeq accession. This data may already exist under a different accession.', 
                    accession=seq_data.id)

            inclusion_passed = evaluate_sequence(
                seq_data=seq_data, 
                seq_qualifier_data=seq_qualifier_data,
                required_parts=required_parts, 
                logger=logger
            )
        
            if not inclusion_passed:
                logger.debug(f'GenBank #{accession} did not pass the curation test...', 
                    accession=accession, including=False
                )
                auto_excluded.append(accession)
                continue

        isolate_type = find_isolate(seq_qualifier_data)
        if isolate_type is None:
            continue
        isolate_name = seq_qualifier_data.get(isolate_type)[0]

        seq_dict = await format_sequence(seq_data, seq_qualifier_data)

        if 'segment' not in seq_dict:
            seq_dict['segment'] = listing.get("schema")[0]['name']
        
        isolate = { 
            'source_name': isolate_name, 
            'source_type': isolate_type
        }
        seq_dict['isolate'] = isolate
        otu_updates.append(seq_dict)
    
    if auto_excluded:
        logger.info(
            'Consider adding these accessions to the catalog exclusion list',
             auto_excluded=auto_excluded)

    return otu_updates

async def write_records(
    otu_path: Path, 
    new_sequences,
    unique_iso, unique_seq, 
    logger = base_logger
):
    """
    """
    ref_isolates = await label_isolates(otu_path)

    try:
        seq_hashes = generate_hashes(
            excluded=unique_seq, n=len(new_sequences))
    except Exception as e:
        logger.exception(e)
        return e
    
    logger.debug(f'Writing {len(new_sequences)} sequences...')

    for seq_data in new_sequences:
        
        isolate_data = seq_data.pop('isolate')
        isolate_name = isolate_data['source_name']
        isolate_type = isolate_data['source_type']

        if isolate_name in ref_isolates:
            iso_hash = ref_isolates[isolate_name]['id']
            logger.debug(
                'Existing isolate name found', 
                iso_name=isolate_name, 
                iso_hash=iso_hash
            )

        else:
            try:
                iso_hash = generate_hashes(excluded=unique_iso, n=1).pop()
            except Exception as e:
                logger.exception(e)
            
            logger.debug('Assigning new isolate hash', 
                iso_name=isolate_name,
                iso_hash=iso_hash)
            
            new_isolate = await store_isolate(
                isolate_name, isolate_type, iso_hash, otu_path)
            
            unique_iso.add(iso_hash)
            ref_isolates[isolate_name] = new_isolate
            
            logger.info('Created a new isolate directory', 
                path=str(otu_path / f'iso_hash'))
        
        iso_path = otu_path / iso_hash
        
        seq_hash = seq_hashes.pop()
        logger.debug('Assigning new sequence', 
            seq_hash=seq_hash,
        )
        
        await store_sequence(seq_data, seq_hash, iso_path)
        unique_seq.add(seq_hash)

        logger.info(f"Wrote new sequence '{seq_hash}'", 
            path=str(iso_path / f'{seq_hash}.json'))
        
async def label_isolates(otu_path: Path) -> dict:
    """
    Return all isolates present in an OTU directory

    :param otu_path: Path to an OTU directory
    :return: A dictionary of isolates indexed by source_name
    """
    if not otu_path.exists():
        return FileNotFoundError

    isolates = {}
    
    for iso_path in get_isolate_paths(otu_path):
        with open(iso_path / "isolate.json", "r") as f:
            isolate = json.load(f)
        isolates[isolate.get('source_name')] = isolate
    
    return isolates

def find_isolate(isolate_data: dict) -> Optional[str]:
    """
    Determine the source type in a Genbank record

    :param isolate_data: Dictionary containing qualifiers in a features section of a Genbank record
    :return:
    """
    for qualifier in ["isolate", "strain"]:
        if qualifier in isolate_data:
            return qualifier

    return None

async def store_isolate(
    source_name: str, source_type: str, 
    iso_hash: str, otu_path: Path
) -> str:
    """
    Creates a new isolate folder for an OTU

    :param source_name: Assigned source name for an accession
    :param source_type: Assigned source type for an accession
    :param iso_hash: Unique ID number for this new isolate
    :param otu_path: Path to the parent OTU
    :return: The unique isolate id (aka iso_hash)
    """
    iso_path = otu_path / iso_hash
    iso_path.mkdir()

    new_isolate = {
        "id": iso_hash,
        "source_type": source_type,
        "source_name": source_name,
        "default": False,
    }

    async with aiofiles.open(iso_path / "isolate.json", "w") as f:
        await f.write(json.dumps(new_isolate, indent=4))
    
    return new_isolate

async def store_sequence(
    sequence: dict, seq_hash: str, iso_path: Path
) -> str:
    """
    Write sequence to isolate directory within the src directory

    :param sequence: Dictionary containing formatted sequence data
    :param seq_hash: Unique ID number for this new sequence
    :param iso_path: Path to the parent isolate
    :return: The unique sequence id (aka seq_hash)
    """
    sequence['_id'] = seq_hash
    seq_path = iso_path / f'{seq_hash}.json'
    
    async with aiofiles.open(seq_path, "w") as f: 
        await f.write(
            json.dumps(sequence, indent=4, sort_keys=True)
        )

    return seq_hash