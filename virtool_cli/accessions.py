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


def run(src: Path, catalog: Path, debugging: bool = False):
    """
    :param src: Path to a reference directory
    :param catalog: Path to an accession catalog directory
    :param debugging: Debugging flag
    """

    filter_class = DEBUG if debugging else INFO
    structlog.configure(
        wrapper_class=structlog.make_filtering_bound_logger(filter_class))
    
    if not catalog.exists():
        logger.info('Initializing .cache directory...', catalog=str(catalog))
        catalog.mkdir()

    get_catalog(src, catalog)

def split_pathname(path):
    """
    """
    [ taxid, otu_id ] = (path.stem).split('--')

    return taxid, otu_id

def search_by_id(catalog: Path, otu_id: str):
    """
    Searches records for a matching id and returns the first matching path in the accession records


    """
    matches = [listing for listing in catalog.glob(f'*--{otu_id}.json')]
    if matches:
        return matches[0]
    else:
        return FileNotFoundError

def get_catalog(src: Path, catalog: Path):
    """
    :param src: Path to a reference directory
    :param catalog: Path to an accession catalog directory
    """
    catalog_paths = [record for record in catalog.glob('*.json')]
    
    if catalog_paths:
        logger.info(
            'Catalog is already populated. \n Checking listings for accuracy...', 
            catalog=str(catalog))

        updated_list = check_listing(src=src, catalog=catalog)
        
        if updated_list:
            logger.info(f'Updated {len(updated_list)} cache entries', 
                updated=updated_list)
        else:
            logger.info(f'Cache is up to date')

    else:
        logger.info('Empty cache, initializing listings...', catalog=str(catalog))

        all_current_accessions = generate_catalog(src)
    
        for taxid in all_current_accessions:
            logger.debug(
                f'Processing OTU with NCBI taxon id {taxid}', 
                taxid=taxid, 
                current_accessions=all_current_accessions.get(taxid))
            
            try:
                write_listing(
                    taxid=taxid, 
                    record=all_current_accessions.get(taxid),
                    catalog=catalog
                )
            except Exception as e:
                logger.exception(e)
    
    
    # logger.info('Accessions written to cache', cache=str(cache))

def check_listing(src: Path, catalog: Path):
    """
    :param src: Path to a reference directory
    :param catalog: Path to a catalog directory
    """
    changed_listings = []

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
        record_path = catalog / f"{taxid}--{otu_data.get('_id')}.json"

        check_listing(
            taxid=taxid,
            listing_path=record_path,
            otu_path=otu_path, 
            ref_accessions=catalog_otu(otu_path), 
            cache=catalog, 
            log=otu_log
        )
    
    return changed_listings

def check_listing(
    taxid, listing_path, otu_path, ref_accessions, cache, log
):
    """
    """
    if listing_path.exists():
        with open(listing_path, "r") as f:
            cached_record = json.load(f)

        if set(ref_accessions) != set(cached_record.get('accessions')['included']):
            log.info('Changes found')
            print(listing_path.name)
            print(cached_record['accessions']['included'])
            print(ref_accessions)

            cached_record.get('accessions')['included'] = ref_accessions
            
            update_listing(listing_path, ref_accessions, cached_record)
            log.debug('Wrote new list')
            return taxid

    else:
        log.info(f'No accession record for f{taxid}. Creating {listing_path.name}')
        
        new_record = generate_record(
            otu_data=parse_otu(otu_path), 
            accession_list=ref_accessions)
        write_listing(taxid, new_record, cache)
        return taxid
    
    return

def generate_catalog(src: Path):
    """
    Initialize an accession catalog by pulling accessions
    from all OTU directories

    :param src: Path to a reference directory
    """
    catalog = {}

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
        accessions = catalog_otu(otu_path)

        catalog[taxid] = generate_record(otu_data, accessions)
    
    return catalog

def generate_record(otu_data: dict, accession_list: list):
    """
    """
    otu_fetch_data = {}
    otu_fetch_data['_id'] = otu_data.get('_id')
    otu_fetch_data['name'] = otu_data.get('name')
    otu_fetch_data['accessions'] = {}
    otu_fetch_data['accessions']['included'] = accession_list
    
    return otu_fetch_data

def catalog_otu(otu_path: Path) -> list:
    """
    Gets all accessions from an OTU directory and returns a list

    :param src: Path to an OTU directory
    """
    accessions = []
    
    for isolate_path in get_otu_paths(otu_path):
        for sequence_path in get_sequence_paths(isolate_path):
            with open(sequence_path, "r") as f:
                sequence = json.load(f)
            accessions.append(sequence['accession'])

    return accessions

def write_listing(
        taxid: int, 
        record: dict, 
        catalog: Path,
        indent: bool = True):
    """
    Write accession file to cache directory

    :param taxid: OTU taxon id
    :param accessions: List of OTU's accessions
    :param cache: Cache directory
    :param indent: Indent flag
    """
    taxid_log = logger.bind(taxid=taxid)

    output_path = catalog / f"{taxid}--{record['_id']}.json"

    taxid_log.debug('Writing accession', accession_path=output_path)

    record['accessions']['excluded'] = []

    with open(output_path, "w") as f:
        json.dump(record, f, indent=2 if indent else None)
    
    return

def update_listing(
        path: Path,
        accessions: list,
        listing: dict, 
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

    listing['accessions']['included'] = accessions

    with open(path, "w") as f:
        json.dump(listing, f, indent=2 if indent else None)
    
    return