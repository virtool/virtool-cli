import json
from pathlib import Path
import asyncio
from typing import Optional
import structlog
from logging import INFO, DEBUG
from urllib.error import HTTPError

from Bio import Entrez, SeqIO

from virtool_cli.utils.ref import get_otu_paths, map_otus
from virtool_cli.utils.hashing import get_unique_ids
from virtool_cli.utils.ncbi import NCBI_REQUEST_INTERVAL
from virtool_cli.accessions import get_accessions_by_id

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

    included_records = filter_records(src, records)
    storage = {}

    requester = asyncio.create_task(fetch_loop(upstream_queue, included_records))
    processor = asyncio.create_task(process_loop(upstream_queue, storage))

    await asyncio.gather(*[requester], return_exceptions=True)
    await upstream_queue.join() # wait until the consumer has processed all items
    processor.cancel()

    await print_new(storage)

async def fetch_loop(queue: asyncio.Queue, record_paths: list):
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

        entry = { 'taxid': taxid, 'data': data }

        await queue.put(entry)
        taxid_log.debug(
            'Pushed %d requests to queue', 
            len(new_accessions), taxid=taxid)
        await asyncio.sleep(0.7)

async def process_loop(queue: asyncio.Queue, storage: dict):
    """
    """
    while True:
        entry = await queue.get()
        
        taxid = entry['taxid']
        taxid_log = logger.bind(taxid=taxid)
        storage[taxid] = {}
        for seq_list in entry['data']:
            # taxid_log.debug("Checking accession...")
            for seq_data in seq_list:
                seq_qualifier_data = await get_qualifiers(seq_data.features)

                isolate_type = await find_isolate(seq_qualifier_data)
                if isolate_type is None:
                    continue

                isolate_name = seq_qualifier_data.get(isolate_type)[0]
                if isolate_name not in storage[taxid]:
                    storage[taxid][isolate_name] = []

                seq_dict = await format_sequence(seq_data, seq_qualifier_data)
                storage[taxid][isolate_name].append(seq_dict)
            
                taxid_log.debug('Processed data', accession=seq_dict['accession'], isolate=isolate_name)

        # isolate_type = await find_isolate(seq_data)
        #     if isolate_type is None:
        #         continue

        #     isolate_name = seq_data.get(isolate_type)[0]
        #     taxid_log.debug('Isolate name found', isolate=isolate_name)

        await asyncio.sleep(0.1)
        queue.task_done()

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

    filter_set = set(record['accessions']['included'] + record['accessions']['excluded'])
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
        included_records.append(get_accessions_by_id(records, otu_id))
    
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
        id_unversioned = accession.id.split(".")[0]
        seq_dict = {
            "_id": id_unversioned,
            "accession": accession.id,
            "definition": accession.description,
            "host": data.get("host")[0] if data.get("host") is not None else None,
            "sequence": str(accession.seq),
        }
        return seq_dict
    
    except Exception as e:
        logger.exception(e)
        return {}

async def print_new(storage: dict):
    for taxid in storage:
        print(f"{taxid}")
        for isolate in storage[taxid]:
            print(f"  {isolate}")
            for sequence in storage[taxid][isolate]:
                # print(sequence)
                print(f"     {sequence['accession']}:") 
                print(f"     {sequence['definition']}")

        print()

if __name__ == '__main__':
    debug = True
    
    REPO_DIR = '/Users/sygao/Development/UVic/Virtool/Repositories'
    
    project_path = Path(REPO_DIR) / 'ref-mini-mini'
    src_path = project_path / 'src'
    cache_path = project_path / '.cache'

    run(src_path, cache_path, debugging=True)