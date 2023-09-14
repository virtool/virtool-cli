import json
from pathlib import Path
import asyncio
import aiofiles
from typing import Optional
import logging
from structlog import BoundLogger
from urllib.error import HTTPError

from virtool_cli.utils.logging import base_logger
from virtool_cli.utils.reference import get_otu_paths, get_isolate_paths
from virtool_cli.utils.hashing import generate_hashes, get_unique_ids
from virtool_cli.utils.ncbi import fetch_upstream_accessions, fetch_upstream_records, NCBI_REQUEST_INTERVAL
from virtool_cli.utils.evaluate import evaluate_sequence, get_qualifiers
from virtool_cli.utils.format import format_sequence
from virtool_cli.accessions.helpers import search_by_otu_id

DEFAULT_INTERVAL = 0.001


def run(
    otu_path: Path, catalog: Path, 
    auto_evaluate: bool = False, 
    debugging: bool = False
):
    """
    CLI entry point for update.update.run()
    
    Requests updates for a single OTU directory
    Searches the catalog for a matching catalog listing and requests new updates if it finds one.

    :param otu_path: Path to a OTU directory
    :param catalog: Path to a catalog directory
    :param auto_evaluate: Auto-evaluation flag, enables automatic filtering for fetched results
    :param debugging: Enables verbose logs for debugging purposes
    """
    filter_class = logging.DEBUG if debugging else logging.INFO
    logging.basicConfig(
        format="%(message)s",
        level=filter_class,
    )

    if auto_evaluate:
        base_logger.warning(
            'Auto-evaluation is in active development and may produce false negatives.'
        )
    
    [ _, otu_id ] = otu_path.name.split('--')
    listing_path = search_by_otu_id(otu_id, catalog)
    
    if listing_path:
        asyncio.run(update_otu(otu_path, listing_path, auto_evaluate))

    else:
        base_logger.critical('Listing not found for this OTU.')
        base_logger.info('Use virtool acc init to create a new accession catalog.')

async def update_otu(
    otu_path: Path, listing_path: Path, 
    auto_evaluate: bool = False
):
    """
    Requests new records for a single taxon ID
    and writes new data under the corresponding path.

    :param otu_path: Path to a OTU directory
    :param listing_path: Path to a listing in an accession catalog directory
    """
    src_path = otu_path.parent
    listing = json.loads(listing_path.read_text())

    # extract taxon ID and _id hash from listing filename
    [ taxid, otu_id ] = (listing_path.stem).split('--')
    
    logger = base_logger.bind(taxid=taxid, otu_id=otu_id)

    record_data = await request_new_records(listing)
    if not record_data:
        return
    
    otu_updates = await process_records(
        records=record_data, listing=listing, 
        auto_evaluate=auto_evaluate, logger=logger)
    if not otu_updates:
        return
    
    # List all isolate and sequence IDs presently in src
    unique_iso, unique_seq = await get_unique_ids(get_otu_paths(src_path))
    
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
    :param listing_path: Path to a listing in an accession catalog directory
    :param logger: Optional entry point for a shared BoundLogger
    """
    try:
        new_accessions = await fetch_upstream_accessions(
            listing=listing, logger=logger)
        
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
    Takes sequence records and:
        1) Evaluates whether those records should be added to the database,
        2) Formats the records into a smaller dictionary
        3) Returns new formatted dicts in a list
    
    WARNING: Auto-evaluation is still under active development, especially multipartite filtering
    
    :param listing_path: Path to a listing in an accession catalog directory
    :param auto_evaluate: Boolean flag for whether automatic evaluation functions
        should be run
    :param logger: Optional entry point for a shared BoundLogger
    """
    filter_set = set(listing['accessions']['included'])
    filter_set.update(listing['accessions']['excluded'])

    if auto_evaluate:
        if not listing.get('schema', []):
            logger.warning('Missing schema. moving on...')
            return False
        
        if len(listing['schema']) < 1:
            required_parts = get_lengthdict_monopartite(listing['schema'], logger)
                    
        else:
            required_parts = get_lengthdict_monopartite(listing['schema'], logger)

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

        seq_dict = format_sequence(
            record=seq_data, qualifiers=seq_qualifier_data, 
            logger=logger
        )

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
    new_sequences: dict,
    unique_iso: set, unique_seq: set, 
    logger: BoundLogger = base_logger
):
    """
    :param unique_iso: Set of all unique isolate IDs present in the reference
    :param unique_seq: Set of all unique sequence IDs present in the reference
    :param otu_path: A path to an OTU directory under a src reference directory
    :param logger: Optional entry point for a shared BoundLogger
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
                path=str(otu_path / iso_hash))
        
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

def find_isolate(record_features: dict) -> Optional[str]:
    """
    Determine the source type in a Genbank record

    :param record_features: Dictionary containing qualifiers in a features section of a Genbank record
    :return: "isolate" or "strain" if listed as such in the qualifiers, else None
    """
    for qualifier in ["isolate", "strain"]:
        if qualifier in record_features:
            return qualifier

    return None

async def store_isolate(
    source_name: str, source_type: str, 
    iso_hash: str, otu_path: Path
) -> str:
    """
    Creates a new isolate directory and metadata file under an OTU directory,
    then returns the metadata in dict form

    :param source_name: Assigned source name for an accession
    :param source_type: Assigned source type for an accession
    :param iso_hash: Unique ID number for this new isolate
    :param otu_path: Path to the parent OTU
    :return: The unique isolate id
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

def get_lengthdict_multipartite(schema: list, logger: BoundLogger = base_logger):
    """
    """
    required_parts = {}

    for part in schema:
        
        if part['required']:
            part_length = part.get('length', -1)
            if part_length < 0:
                logger.warning(
                    f"{part['name']} lacks a length listing.")
                
            required_parts[part['name']] = part_length
    
    return required_parts

def get_lengthdict_monopartite(schema: list, logger: BoundLogger = base_logger):
    """
    """
    part = schema[0]

    if (part_name := part.get('name', None)) is None:
        part_name = 'unlabelled'
    
    part_length = part.get('length', -1)
    if part_length < 0:
        logger.warning(
            f"{part['name']} lacks a length listing.")
    
    return { part_name: part_length }