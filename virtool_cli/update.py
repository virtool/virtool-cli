import json
from pathlib import Path
import asyncio
import aiofiles
from typing import Optional
import structlog
import logging
from urllib.error import HTTPError
from Bio import Entrez, SeqIO

from virtool_cli.utils.logging import base_logger
from virtool_cli.utils.ref import get_otu_paths, get_isolate_paths, search_otu_by_id
from virtool_cli.utils.hashing import generate_hashes, get_unique_ids
from virtool_cli.utils.ncbi import NCBI_REQUEST_INTERVAL
from virtool_cli.accessions.helpers import search_by_id, filter_catalog
from virtool_cli.accessions.evaluate import evaluate

DEFAULT_INTERVAL = 0.001


def run(src: Path, catalog: Path, debugging: bool = False):
    """
    Calls the parent routine

    :param src: Path to a reference directory
    :param catalog: Path to a catalog directory
    :param debugging: Enables verbose logs for debugging purposes
    """
    filter_class = logging.DEBUG if debugging else logging.INFO
    logging.basicConfig(
        format="%(message)s",
        level=filter_class,
    )
    logger = structlog.get_logger()
    logger = logger.bind(src=str(src), catalog=str(catalog))

    logger.info('Updating src directory accessions using catalog listings...')

    asyncio.run(
        update(src_path=src, catalog_path=catalog)
    )

async def update(src_path: Path, catalog_path: Path):
    """
    Creates 2 queues:
        1) upstream: Holds raw NCBI GenBank data,
        2) write: Holds formatted sequence data and isolate data
    
    Monitors 3 asynchrononous processes:
        1) fetcher: 
            Requests and retrieves new accessions from NCBI GenBank
            and pushes results to upstream queue
        2) processor: 
            Pulls Genbank data from upstream queue, 
            formats into dict form and pushes to write queue
        3) writer: 
            Pulls formatted sequences from write queue
            and writes json to the correct location in the 
            src directory

    :param src_path: Path to a reference directory
    :param catalog_path: Path to a catalog directory
    """

    # Filter out cached listings that are not present in this src directory
    included_listings = filter_catalog(src_path, catalog_path)
    storage = {}

    # Holds raw NCBI GenBank data
    upstream_queue = asyncio.Queue()

    # Holds formatted sequence data and isolate data
    write_queue = asyncio.Queue()

    # Requests and retrieves new accessions from NCBI GenBank
    # and pushes results to upstream queue
    fetcher = asyncio.create_task(
        fetcher_loop(included_listings, queue=upstream_queue))
    
    # Pulls Genbank data from upstream queue, formats into dict form 
    # and pushes to write queue
    asyncio.create_task(
        processor_loop(upstream_queue, write_queue))
    
    # Pulls formatted sequences from write queue, checks isolate metadata
    # and writes json to the correct location in the src directory
    asyncio.create_task(
        writer_loop(src_path, storage, write_queue, dry_run=False))

    await asyncio.gather(*[fetcher], return_exceptions=True)

    await upstream_queue.join() # wait until the consumer has processed all items
    
    await write_queue.join()

    return

async def fetcher_loop(
    listing_paths: list,
    queue: asyncio.Queue 
) -> None:
    """
    Loops through selected OTU listings from accession catalogue, 
    indexed by NCBI taxon ID, and:
        1) requests NCBI Genbank for accession numbers not extant
            in catalog,
        2) loops through retrieved new accession numbers and 
            requests relevant record data from NCBI Genbank
        3) Pushes new records and corresponding OTU information
            to a queue for formatting
    
    :param listing_paths: A list of paths to listings from the accession catalog
    :param queue: Queue holding fetched NCBI GenBank data
    """
    base_logger.debug("Starting fetcher...")

    for path in listing_paths:
        acc_listing = json.loads(path.read_text())

        # extract taxon ID and _id hash from listing filename
        [ taxid, otu_id ] = (path.stem).split('--')
        
        logger = base_logger.bind(taxid=taxid, otu_id=otu_id)

        try:
            new_accessions = await fetch_upstream_accessions(
                taxid, acc_listing, reqseq_only=True)
        except Exception as e:
            logger.exception(e)
            continue
        
        logger.debug('New accessions', new=new_accessions)

        record_data = []
        for accession in new_accessions:
            new_data = await fetch_upstream_records(accession, logger)
            record_data.append(new_data)
            await asyncio.sleep(NCBI_REQUEST_INTERVAL)
        
        if not record_data:
            continue

        packet = { 
            'taxid': taxid, 'otu_id': otu_id, 
            'listing': acc_listing,
            'data': record_data, 
        }

        await queue.put(packet)
        logger.debug(
            f'Pushed {len(new_accessions)} requests to upstream queue', 
            uids=new_accessions, n_requests=len(new_accessions), taxid=taxid)
        await asyncio.sleep(DEFAULT_INTERVAL)

async def processor_loop(
    upstream_queue: asyncio.Queue, 
    downstream_queue: asyncio.Queue
):
    """
    Awaits fetched sequence data from the fetcher:
        1) Formats the sequence data into reference-compatible dictionaries,
        2) Checks for validity,
        3) Pushes the formatted data into the downstream queue to be dealt with by the writer

    :param upstream_queue: Queue holding NCBI GenBank data pushed by the fetcher,
    :param downstream_queue: Queue holding formatted sequence and isolate data processed by this loop
    """
    base_logger.debug("Starting processor...")

    while True:
        fetch_packet = await upstream_queue.get()
        
        listing = fetch_packet['listing']
        taxid = fetch_packet['taxid']
        otu_id = fetch_packet['otu_id']
        logger = base_logger.bind(taxid=taxid)

        filter_set = set(listing['accessions']['included'].keys())

        otu_updates = []
        for seq_list in fetch_packet['data']:
            for seq_data in seq_list:
                [ accession, version ] = (seq_data.id).split('.')
                
                if accession in filter_set:
                    logger.debug('Accession already exists', accession=seq_data.id)
                    continue
                
                seq_qualifier_data = await get_qualifiers(seq_data.features)

                if len(listing.get('schema')) > 1:
                    if 'segment' not in seq_qualifier_data:
                        continue

                isolate_type = await find_isolate(seq_qualifier_data)
                if isolate_type is None:
                    continue
                isolate_name = seq_qualifier_data.get(isolate_type)[0]

                seq_dict = await format_sequence(seq_data, seq_qualifier_data)
                if 'segment' not in seq_dict:
                    seq_dict['segment'] = listing.get("schema")[0]
                isolate = { 
                    'source_name': isolate_name, 
                    'source_type': isolate_type
                }
                seq_dict['isolate'] = isolate
                otu_updates.append(seq_dict)

        if not otu_updates:
            await asyncio.sleep(DEFAULT_INTERVAL)
            upstream_queue.task_done()
            continue

        processed_packet = { 'taxid': taxid, 'otu_id': otu_id, 'data': otu_updates }

        await downstream_queue.put(processed_packet)
        logger.debug(
            f'Pushed {len(otu_updates)} new accessions to downstream queue')
        await asyncio.sleep(DEFAULT_INTERVAL)
        upstream_queue.task_done()

async def writer_loop(
    src_path: Path, 
    storage: dict, 
    queue: asyncio.Queue, 
    dry_run: bool = False
):
    """
    TO-DO: write changes back to local catalog?

    :param src_path: Path to a reference directory
    """
    base_logger.debug("Starting writer...")

    unique_iso, unique_seq = await get_unique_ids(get_otu_paths(src_path))

    while True:
        packet = await queue.get()
        
        taxid = packet['taxid']
        otu_id = packet['otu_id']
        new_sequence_set = packet['data']

        log = base_logger.bind(otu_id=otu_id, taxid=taxid)

        otu_path = search_otu_by_id(src_path, otu_id)
        ref_isolates = await label_isolates(otu_path)

        try:
            seq_hashes = generate_hashes(
                excluded=unique_seq, n=len(new_sequence_set))
        except Exception as e:
            log.exception(e)
            return e
        
        log.debug(f'Writing {len(new_sequence_set)} sequences...')

        for seq_data in new_sequence_set:
            
            isolate_data = seq_data.pop('isolate')
            isolate_name = isolate_data['source_name']
            isolate_type = isolate_data['source_type']

            if isolate_name in ref_isolates:
                iso_hash = ref_isolates[isolate_name]['id']
                log.debug(
                    'Existing isolate name found', 
                    iso_name=isolate_name, 
                    iso_hash=iso_hash
                )

            else:
                try:
                    iso_hash = generate_hashes(excluded=unique_iso, n=1).pop()
                except Exception as e:
                    log.exception(e)
                
                log.debug('Assigning new isolate hash', 
                    iso_name=isolate_name,
                    iso_hash=iso_hash)
                
                new_isolate = await store_isolate(
                    isolate_name, isolate_type, iso_hash, otu_path)
                
                unique_iso.add(iso_hash)
                ref_isolates[isolate_name] = new_isolate
                
                log.info(
                    'Created a new isolate directory', 
                    path=str(otu_path.relative_to(src_path) / f'iso_hash'))
            
            iso_path = otu_path / iso_hash
            
            seq_hash = seq_hashes.pop()
            log.debug('Assigning new sequence', 
                seq_hash=seq_hash,
            )
            
            await store_sequence(seq_data, seq_hash, iso_path)
            unique_seq.add(seq_hash)

            log.info(f"Wrote new sequence '{seq_hash}'", 
                path=str(iso_path.relative_to(src_path) / f'{seq_hash}.json'))
    
        await asyncio.sleep(DEFAULT_INTERVAL)
        queue.task_done()

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

async def fetch_upstream_accessions(
    taxid: int, 
    listing: dict,
    reqseq_only: bool = True
) -> list:
    """
    :param taxid: OTU Taxon ID
    :param listing: Corresponding listing from the accession catalog for this OTU in dict form
    :return: A list of accessions from NCBI Genbank for the taxon ID, 
        sans included and excluded accessions
    """

    filter_set = set()
    filter_set.update(listing['accessions']['included'].values())
    filter_set.update(listing['accessions']['excluded'].values())

    logger = base_logger.bind(taxid=taxid)
    logger.debug('Exclude catalogued UIDs', catalogued=filter_set)

    upstream_uids = []

    if reqseq_only:
        otu_searchname = "+".join(listing.get("name").split(" "))
        entrez_accessions = Entrez.read(
            Entrez.esearch(
                db="nucleotide", 
                term=f"{otu_searchname}[Organism] AND RefSeq[keyword]"
            )
        )

        for entrez_id in entrez_accessions['IdList']:
            upstream_uids.append(entrez_id)
    
    else:
        entrez_accessions = Entrez.read(
            Entrez.elink(
                dbfrom="taxonomy", db="nucleotide", 
                id=str(taxid), idtype="acc")
        )

        for linksetdb in entrez_accessions[0]["LinkSetDb"][0]["Link"]:
            upstream_uids.append(linksetdb["Id"])

    upstream_set = set(upstream_uids)

    return list(upstream_set.difference(filter_set))

async def fetch_upstream_records(
    fetch_list: list, 
    logger: structlog.BoundLogger = base_logger
) -> list:
    """
    Take a list of accession numbers and request the records from NCBI GenBank
    
    :param fetch_list: List of accession numbers to fetch from GenBank
    :param logger: Structured logger

    :return: A list of GenBank data converted from XML to dicts if possible, 
        else an empty list
    """
    try:
        handle = Entrez.efetch(
            db="nucleotide", id=fetch_list, rettype="gb", retmode="text"
        )
    except HTTPError as e:
        logger.error(f'{e}, moving on...')
        # logger.exception(e)
        return []
    
    ncbi_records = SeqIO.to_dict(SeqIO.parse(handle, "gb"))

    if ncbi_records is None:
        return []
    
    try:
        accession_list = [record for record in ncbi_records.values() if record.seq]
        return accession_list
    except Exception as e:
        logger.exception(e)
        raise e

# def filter_catalog(
#     src_path, catalog_path
# ) -> list:
#     """
#     Return paths for cached accession catalogues that are included in the source reference

#     :param src_path: Path to a reference directory
#     :param catalog_path: Path to an accession record directory
#     :return: A list of paths to relevant listings in the accession catalog
#     """
#     otu_paths = get_otu_paths(src_path)
#     included_listings = []
    
#     for path in otu_paths:
#         otu_id = (path.name).split('--')[1]

#         included_listings.append(
#             search_by_id(otu_id, catalog_path)
#         )
    
#     return included_listings

async def get_qualifiers(seq: list) -> dict:
    """
    Get relevant qualifiers in a Genbank record

    :param seq: SeqIO features object for a particular accession
    :return: Dictionary containing all qualifiers in the source field of the features section of a Genbank record
    """
    features = [feature for feature in seq if feature.type == "source"]
    
    isolate_data = {}

    for feature in features:
        for qualifier in feature.qualifiers:
            isolate_data[qualifier] = feature.qualifiers.get(qualifier)

    return isolate_data

async def find_isolate(isolate_data: dict) -> Optional[str]:
    """
    Determine the source type in a Genbank record

    :param isolate_data: Dictionary containing qualifiers in a features section of a Genbank record
    :return:
    """
    for qualifier in ["isolate", "strain"]:
        if qualifier in isolate_data:
            return qualifier

    return None

async def format_sequence(record: SeqIO.SeqRecord, qualifiers: dict) -> dict:
    """
    Creates a new sequence file for a given isolate

    :param path: Path to a isolate folder
    :param record: Genbank record object for a given accession
    :param qualifiers: Dictionary containing all qualifiers in the source field of the features section of a Genbank record
    :return: A new sequence dictionary if possible, else an empty dict if not
    """
    logger = base_logger.bind(accession=record.id)
    # logger.debug(f"{qualifiers}")
    try:
        seq_dict = {
            "accession": record.id,
            "definition": record.description,
            "host": qualifiers.get("host")[0] if qualifiers.get("host") is not None else None,
            "sequence": str(record.seq),
        }
        if 'segment' in qualifiers:
            seq_dict['segment'] = qualifiers.get("segment")[0] if qualifiers.get("segment") is not None else None

        return seq_dict
    
    except Exception as e:
        logger.exception(e)
        return {}

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