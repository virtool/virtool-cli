from pathlib import Path
import json
import structlog
from logging import INFO, DEBUG

from virtool_cli.utils.ref import (
    parse_otu,
    get_otu_paths,
    get_sequence_paths
)

base_logger = structlog.get_logger()


def run(src_path: Path, catalog_path: Path, debugging: bool = False):
    """
    :param src: Path to a reference directory
    :param catalog: Path to an accession catalog directory
    :param debugging: Debugging flag
    """

    filter_class = DEBUG if debugging else INFO
    structlog.configure(
        wrapper_class=structlog.make_filtering_bound_logger(filter_class))
    
    if not catalog_path.exists():
        base_logger.info('Initializing .cache/catalog directory...', catalog=str(catalog_path))
        catalog_path.mkdir()

    get_catalog(src_path, catalog_path)

def split_pathname(path):
    """
    Split a filename formatted as 'taxid--otu_id.json' into (taxid, otu_id)

    :param path: Path to a listing in an accession catalog
    """
    [ taxid, otu_id ] = (path.stem).split('--')

    return taxid, otu_id

def search_by_id(otu_id: str, catalog_path: Path):
    """
    Searches records for a matching id and returns the first matching path in the accession records

    :param otu_id: Unique OTU ID string
    :param catalog_path: Path to an accession catalog directory
    """
    matches = [listing for listing in catalog_path.glob(f'*--{otu_id}.json')]
    if matches:
        return matches[0]
    else:
        return FileNotFoundError

def get_catalog(src: Path, catalog: Path):
    """
    :param src: Path to a reference directory
    :param catalog: Path to an accession catalog directory
    """
    logger = base_logger.bind(catalog=str(catalog))

    catalog_paths = [record for record in catalog.glob('*.json')]
    
    if catalog_paths:
        logger.info(
            'Catalog is already populated. \n Checking listings for accuracy...')

        updated_list = check_catalog(src_path=src, catalog_path=catalog)
        
        if updated_list:
            logger.info(f'Updated {len(updated_list)} cache entries', 
                updated=updated_list)
        else:
            logger.info(f'Cache is up to date')

    else:
        logger.info(
            'Empty cache, initializing listings...', 
            catalog=str(catalog))

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
                    catalog_path=catalog
                )
            except Exception as e:
                logger.exception(e)
    
    
    # logger.info('Accessions written to cache', cache=str(cache))

def check_catalog(src_path: Path, catalog_path: Path):
    """
    :param src_path: Path to a reference directory
    :param catalog_path: Path to a catalog directory
    """
    changed_listings = []

    for otu_path in get_otu_paths(src_path):
        logger = base_logger.bind(
            otu_id=otu_path.name.split('--')[1],
            path=str(otu_path.relative_to(otu_path.parents[1])))

        try:
            otu_data = parse_otu(otu_path)
            taxid = otu_data.get('taxid', None)
            if taxid is None:
                logger.debug('taxid=null', name=otu_data['name'])
                continue

        except Exception as e:
            logger.exception('taxid not found')
            continue

        logger.bind(name=otu_data['name'], taxid=taxid)
        record_path = catalog_path / f"{taxid}--{otu_data.get('_id')}.json"

        check_listing(
            taxid=taxid,
            ref_accessions=catalog_otu(otu_path), 
            otu_path=otu_path, 
            listing_path=record_path,
            catalog_path=catalog_path, 
            logger=logger
        )
    
    return changed_listings

def check_listing(
    taxid: int, 
    ref_accessions: list, 
    otu_path: Path, 
    listing_path: Path,
    logger: structlog.BoundLogger
):
    """
    :param taxid: Path to a reference directory
    :param ref_accessions: List of included accessions
    :param otu_path: Path to a OTU directory in a src directory
    :param listing_path: Path to a listing file in an accession catalog directory
    :param logger: Structured logger
    """
    if listing_path.exists():
        with open(listing_path, "r") as f:
            cached_record = json.load(f)

        if set(ref_accessions) != set(cached_record.get('accessions')['included']):
            logger.info('Changes found')
            print(listing_path.name)
            print(cached_record['accessions']['included'])
            print(ref_accessions)

            cached_record.get('accessions')['included'] = ref_accessions
            
            update_listing(listing_path, ref_accessions, cached_record)
            logger.debug('Wrote new list')
            return taxid

    else:
        logger.info(f'No accession record for f{taxid}. Creating {listing_path.name}')
        
        new_record = generate_record(
            otu_data=parse_otu(otu_path), 
            accession_list=ref_accessions)
        
        write_listing(taxid, new_record, catalog_path=listing_path.parent)
        
        return taxid
    
    return ''

def generate_catalog(src: Path):
    """
    Initialize an accession catalog by pulling accessions from all OTU directories

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
            base_logger.warning('taxid not found', path=src(otu_path))
            continue
        accessions = catalog_otu(otu_path)

        catalog[taxid] = generate_record(otu_data, accessions)
    
    return catalog

def generate_record(otu_data: dict, accession_list: list):
    """
    :param otu_data: OTU data in dict form
    :param accession_list: list of included accesssions
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

    :param otu_path: Path to an OTU directory
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
        catalog_path: Path,
        indent: bool = True):
    """
    Write accession file to cache directory

    :param taxid: OTU taxon id
    :param accessions: List of OTU's accessions
    :param catalog: Path to an accession catalog
    :param indent: Indent flag
    """
    taxid_log = base_logger.bind(taxid=taxid)

    output_path = catalog_path / f"{taxid}--{record['_id']}.json"

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
    :param listing: 
    :param indent: Indent flag
    """
    taxid_log = base_logger.bind(taxid=path.stem)

    taxid_log.debug('Updating accession', accession_path=path)

    listing['accessions']['included'] = accessions

    with open(path, "w") as f:
        json.dump(listing, f, indent=2 if indent else None)
    
    return