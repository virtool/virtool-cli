from pathlib import Path
import json
import structlog
from logging import INFO, DEBUG

from virtool_cli.utils.ref import (
    parse_otu,
    get_otu_paths,
    get_sequence_paths
)

logger = structlog.get_logger()


def run(src: Path, cache: Path, debugging: bool = False):
    """
    :param src: Path to a reference directory
    :param cache: Path to a cache directory
    """

    filter_class = DEBUG if debugging else INFO
    structlog.configure(
        wrapper_class=structlog.make_filtering_bound_logger(filter_class))
    
    if not cache.exists():
        logger.info('Initializing .cache directory...', cache=str(cache))
        cache.mkdir()

    get_all_accessions(src, cache)

def get_all_accessions(src: Path, cache: Path):
    """
    :param src: Path to a reference directory
    :param cache: Path to a cache directory
    """
    cache_paths = [record for record in cache.glob('*.json')]
    
    if cache_paths:
        logger.info(
            'Cache is already populated. \n Checking accessions for accuracy...', 
            cache=str(cache))

        updated_list = check_all_records(src=src, cache=cache)
        
        if updated_list:
            logger.info(f'Updated {len(updated_list)} cache entries', updated=updated_list)
        else:
            logger.info(f'Cache is up to date')

    else:
        logger.info('Empty cache, initializing records...', cache=str(cache))

        all_current_accessions = generate_accession_records(src)
    
        for taxid in all_current_accessions:
            logger.debug(
                f'Processing OTU with NCBI taxon id {taxid}', 
                taxid=taxid, 
                current_accessions=all_current_accessions.get(taxid))
            
            try:
                write_record(
                    taxid=taxid, 
                    record=all_current_accessions.get(taxid),
                    cache=cache
                )
            except Exception as e:
                logger.exception(e)
    
    
    # logger.info('Accessions written to cache', cache=str(cache))

def check_all_records(src: Path, cache: Path):
    """
    :param src: Path to a reference directory
    :param cache: Path to a cache directory
    """
    cache_change_list = []

    for otu_path in get_otu_paths(src):
        try:
            otu_data = parse_otu(otu_path)
            taxid = otu_data.get('taxid', None)
            if taxid is None:
                logger.debug('taxid is null', name=otu_data['name'])
                continue

        except Exception as e:
            logger.exception('taxid not found')
            continue

        otu_log = logger.bind(name=otu_data['name'], otu_id=otu_data['_id'], taxid=taxid)
        record_path = cache / f"{taxid}--{otu_data.get('_id')}.json"

        check_record(
            taxid=taxid,
            record_path=record_path,
            otu_path=otu_path, 
            ref_accessions=get_current_accessions(otu_path), 
            cache=cache, 
            otu_log=otu_log
        )
    
    return cache_change_list

def check_record(taxid, record_path, otu_path, ref_accessions, cache, otu_log):
    """
    """
    if record_path.exists():
        with open(record_path, "r") as f:
            cached_record = json.load(f)

        if set(ref_accessions) != set(cached_record.get('accessions')['included']):
            otu_log.info('Changes found')
            print(record_path.name)
            print(cached_record['accessions']['included'])
            print(ref_accessions)

            cached_record.get('accessions')['included'] = ref_accessions
            
            update_record(record_path, ref_accessions, cached_record)
            otu_log.debug('Wrote new list')
            return taxid

    else:
        otu_log.info(f'No accession record for f{taxid}. Creating {record_path.name}')
        
        new_record = generate_record(
            otu_data=parse_otu(otu_path), 
            accession_list=ref_accessions)
        write_record(taxid, new_record, cache)
        return taxid
    
    return

def generate_accession_records(src: Path):
    """
    Initialize an accession cache list by pulling accessions
    from all OTU directories

    :param src: Path to a reference directory
    """
    all_accessions = {}

    counter = 0
    for otu_path in get_otu_paths(src):
        try:
            otu_data = parse_otu(otu_path)
            taxid = otu_data.get('taxid')
            if taxid is None:
                taxid = f'none{counter}'
                counter += 1
        except Exception as e:
            logger.warning('taxid not found', path=src(otu_path))
            continue
        accessions = get_current_accessions(otu_path)

        all_accessions[taxid] = generate_record(otu_data, accessions)
    
    return all_accessions

def generate_record(otu_data: dict, accession_list: list):
    """
    """
    otu_fetch_data = {}
    otu_fetch_data['_id'] = otu_data.get('_id')
    otu_fetch_data['name'] = otu_data.get('name')
    otu_fetch_data['accessions'] = {}
    otu_fetch_data['accessions']['included'] = accession_list
    
    return otu_fetch_data

def get_current_accessions(otu_path: Path):
    """
    Gets all accessions from an OTU directory

    :param src: Path to an OTU directory
    """
    accessions = []
    
    for isolate_path in get_otu_paths(otu_path):
        for sequence_path in get_sequence_paths(isolate_path):
            with open(sequence_path, "r") as f:
                sequence = json.load(f)
            accessions.append(sequence['accession'])

    return accessions

def write_record(
        taxid: int, 
        record: dict, 
        cache: Path,
        indent: bool = True):
    """
    Write accession file to cache directory

    :param taxid: OTU taxon id
    :param accessions: List of OTU's accessions
    :param cache: Cache directory
    :param indent: Indent flag
    """
    taxid_log = logger.bind(taxid=taxid)

    output_path = cache / f"{taxid}--{record['_id']}.json"

    taxid_log.debug('Writing accession', accession_path=output_path)

    record['accessions']['excluded'] = []

    with open(output_path, "w") as f:
        json.dump(record, f, indent=2 if indent else None)
    
    return

def update_record(
        path: Path,
        accessions: list,
        current_record: dict, 
        indent: bool = True):
    """
    Write accession file to cache directory

    :param taxid: OTU taxon id
    :param accessions: List of OTU's accessions
    :param cache: Cache directory
    :param indent: Indent flag
    """
    taxid_log = logger.bind(taxid=path.stem)

    taxid_log.debug('Updating accession', accession_path=path)

    current_record['accessions']['included'] = accessions

    with open(path, "w") as f:
        json.dump(current_record, f, indent=2 if indent else None)
    
    return


if __name__ == '__main__':
    REPO_DIR = '/Users/sygao/Development/UVic/Virtool/Repositories'
    
    src_path = Path(REPO_DIR) / 'ref-mini/src'
    cache_path = src_path.parent / '.cache'

    if not cache_path.exists():
        cache_path.mkdir()

    accession_path = Path(REPO_DIR) / 'ref-fetched-accessions/src'

    run(src=src_path, cache=cache_path, debugging=False)