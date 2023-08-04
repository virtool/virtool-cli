import json
from pathlib import Path
import asyncio
import aiofiles
from typing import Optional
import structlog
from logging import INFO, DEBUG
from urllib.error import HTTPError

from Bio import Entrez, SeqIO

from virtool_cli.utils.ref import get_otu_paths, get_isolate_paths, parse_otu, parse_isolates, map_otus, search_otu_by_id
from virtool_cli.utils.hashing import generate_hashes, get_unique_ids
from virtool_cli.utils.ncbi import NCBI_REQUEST_INTERVAL
from virtool_cli.accessions import search_accessions_by_id

logger = structlog.get_logger()


def run(src: Path, records: Path, debugging: bool = False):
    """
    Runs the asynchronous routines to find new isolates for all OTU in a reference

    :param src: Path to a reference directory
    """
    filter_class = DEBUG if debugging else INFO
    structlog.configure(
        wrapper_class=structlog.make_filtering_bound_logger(filter_class))

    asyncio.run(update(src, records))

async def update(src: Path, records: Path):
    """
    """
    upstream_queue = asyncio.Queue()
    write_queue = asyncio.Queue()

    included_records = filter_records(src, records)
    storage = {}

    fetcher = asyncio.create_task(fetcher_loop(upstream_queue, included_records))
    processor = asyncio.create_task(processor_loop(upstream_queue, write_queue))
    writer = asyncio.create_task(writer_loop(write_queue, src, storage))

    await asyncio.gather(*[fetcher], return_exceptions=True)

    await upstream_queue.join() # wait until the consumer has processed all items
    
    await write_queue.join()

    # print(storage)

async def fetcher_loop(queue: asyncio.Queue, record_paths: list):
    """
    """
    for record_path in record_paths:
        record = json.loads(record_path.read_text())

        taxid = (record_path.name).split('--')[0]
        taxid_log = logger.bind(taxid=taxid)

        new_accessions = await fetch_upstream_new(taxid, record)
        # await asyncio.sleep(NCBI_REQUEST_INTERVAL)
        
        await taxid_log.adebug('New accessions', new=new_accessions)

        data = []
        for accession in new_accessions:
            new_data = await fetch_accession_data(accession, taxid_log)
            data.append(new_data)
            # await asyncio.sleep(NCBI_REQUEST_INTERVAL)

        packet = { 'taxid': taxid, 'otu_id': record['_id'], 'data': data }

        await queue.put(packet)
        taxid_log.debug(
            'Pushed %d requests to upstream queue', 
            len(new_accessions), taxid=taxid)
        await asyncio.sleep(0.7)

async def processor_loop(
    upstream_queue: asyncio.Queue, 
    downstream_queue: asyncio.Queue):
    """
    """
    while True:
        fetch_packet = await upstream_queue.get()
        
        taxid = fetch_packet['taxid']
        otu_id = fetch_packet['otu_id']
        taxid_log = logger.bind(taxid=taxid)

        otu_updates = []
        for seq_list in fetch_packet['data']:
            for seq_data in seq_list:
                seq_qualifier_data = await get_qualifiers(seq_data.features)

                isolate_type = await find_isolate(seq_qualifier_data)
                if isolate_type is None:
                    continue
                isolate_name = seq_qualifier_data.get(isolate_type)[0]
                
                # if isolate_name not in otu_updates:
                #     otu_updates[isolate_name] = []

                seq_dict = await format_sequence(seq_data, seq_qualifier_data)
                isolate = { 
                    'source_name': isolate_name, 
                    'source_type': isolate_type
                }
                seq_dict['isolate'] = isolate
                otu_updates.append(seq_dict)
                # otu_updates[isolate_name].append(seq_dict)
            
                # taxid_log.debug(
                #     'Processed data', 
                #     accession=seq_dict['accession'], 
                #     isolate=isolate_name)
                
        processed_packet = { 'taxid': taxid, 'otu_id': otu_id, 'data': otu_updates }

        await downstream_queue.put(processed_packet)
        taxid_log.debug(
            f'Pushed new accessions to downstream queue')
        await asyncio.sleep(0.7)
        upstream_queue.task_done()

async def writer_loop(queue, src_path, storage):
    """
    """
    unique_iso, unique_seq = await get_unique_ids(get_otu_paths(src_path))
    # generate_hashes(n=1, length=8, mixed_case=False, excluded=unique_seq)

    while True:
        packet = await queue.get()
        
        taxid = packet['taxid']
        otu_id = packet['otu_id']
        new_sequence_set = packet['data']

        log = logger.bind(otu_id=otu_id, taxid=taxid)

        otu_path = search_otu_by_id(src_path, otu_id)
        ref_isolates = await label_isolates(otu_path)

        print_new(new_sequence_set)

        try:
            seq_hashes = generate_hashes(n=len(new_sequence_set), length=8, mixed_case=False, excluded=unique_seq)
        except Exception as e:
            log.exception(e)
            return e

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
                    iso_hash = generate_hashes(
                        n=1, length=8, mixed_case=False, excluded=unique_iso).pop()
                except Exception as e:
                    log.exception(e)
                
                log.debug('Assigning new isolate hash', 
                    iso_name=isolate_name,
                    iso_hash=iso_hash)
                await store_isolate(isolate_name, isolate_type, iso_hash, otu_path, indent=2)
                unique_iso.add(iso_hash)
            
            iso_path = otu_path / iso_hash
            
            seq_hash = seq_hashes.pop()
            log.debug('Assigning new sequence', 
                seq_name=isolate_name,
                seq_hash=seq_hash)
            
            await store_sequence(seq_data, seq_hash, iso_path, indent=2)
            unique_seq.add(seq_hash)
            
        #     try:
        #         seq_hashes = generate_hashes(n=4, length=8, mixed_case=False, excluded=unique_seq)
        #     except Exception as e:
        #         log.exception(e)
        #         return e

        #     for sequence in new_sequence_set.get(source_name):
        #         seq_hash = seq_hashes.pop()
        #         log.debug('Assigning new sequence', 
        #             iso_name=source_name,
        #             iso_hash=iso_hash)
                
        #         await store_sequence(sequence, seq_hash, iso_path, indent=2)
        #         unique_seq.add(seq_hash)

        await asyncio.sleep(0.1)
        queue.task_done()

async def store_isolate(
    source_name: str, source_type: str, 
    iso_hash: str, otu_path: Path, indent=0
):
    """
    Creates a new isolate folder for an OTU

    :param path: Path to an OTU
    :param source_name: Assigned source name for an accession
    :param source_type: Assigned source type for an accession
    :param unique_ids: Set containing unique Virtool ids for all isolates in a reference
    :return: A newly generated unique id
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
        await f.write(json.dumps(new_isolate, indent=indent))
    
    return iso_hash

async def store_sequence(
    sequence: dict, seq_hash: str, iso_path: Path, indent: int = 0):
    """
    """
    sequence['_id'] = seq_hash
    seq_path = iso_path / f'{seq_hash}.json'
    
    async with aiofiles.open(seq_path, "w") as f: 
        await f.write(json.dumps(sequence, indent=indent, sort_keys=True))

    return seq_hash

async def label_isolates(otu_path):
    """
    """
    isolates = {}
    
    for iso_path in get_isolate_paths(otu_path):
        with open(iso_path / "isolate.json", "r") as f:
            isolate = json.load(f)
        isolates[isolate.get('source_name')] = isolate
    
    return isolates

async def fetch_upstream_new(taxid: int, record: dict) -> list:
    """
    """
    upstream_accessions = []
    entrez_record = Entrez.read(
        Entrez.elink(
            dbfrom="taxonomy", db="nucleotide", 
            id=str(taxid), idtype="acc")
    )

    for linksetdb in entrez_record[0]["LinkSetDb"][0]["Link"]:
        upstream_accessions.append(linksetdb["Id"])

    filter_set = set(
        record['accessions']['included'] + \
        record['accessions']['excluded'])
    upstream_set = set(upstream_accessions)

    return list(upstream_set.difference(filter_set))

async def fetch_accession_data(fetch_list: list, log) -> list:
    """
    """
    try:
        handle = Entrez.efetch(
            db="nucleotide", id=fetch_list, rettype="gb", retmode="text"
        )
    except HTTPError as e:
        log.exception(e)
        return []
    
    ncbi_records = SeqIO.to_dict(SeqIO.parse(handle, "gb"))

    if ncbi_records is None:
        return []
    
    try:
        accession_list = [record for record in ncbi_records.values() if record.seq]
        return accession_list
    except Exception as e:
        log.exception(e)
        return []

def filter_records(src, records):
    """
    """
    otu_paths = get_otu_paths(src)
    included_records = []
    
    for path in otu_paths:
        otu_id = (path.name).split('--')[1]
        included_records.append(search_accessions_by_id(records, otu_id))
    
    return included_records

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

async def format_sequence(accession: SeqIO.SeqRecord, data: dict) -> dict:
    """
    Creates a new sequence file for a given isolate

    :param path: Path to a isolate folder
    :param accession: Genbank record object for a given accession
    :param data: Dictionary containing all qualifiers in the source field of the features section of a Genbank record
    :return: The newly generated unique id
    """
    try:
        seq_dict = {
            "accession": accession.id,
            "definition": accession.description,
            "host": data.get("host")[0] if data.get("host") is not None else None,
            "sequence": str(accession.seq),
        }
        return seq_dict
    
    except Exception as e:
        logger.exception(e)
        return {}

def print_new(otu_record: dict) -> None:
    for sequence in otu_record:
        print(f"     {sequence['accession']}:") 
        print(f"     {sequence['definition']}")
        print(f"     {sequence['isolate']}")

        print()
    return


if __name__ == '__main__':
    debug = True
    
    REPO_DIR = '/Users/sygao/Development/UVic/Virtool/Repositories'
    
    project_path = Path(REPO_DIR) / 'ref-mini-mini'
    src_path = project_path / 'src'
    cache_path = project_path / '.cache'

    run(src_path, cache_path, debugging=True)