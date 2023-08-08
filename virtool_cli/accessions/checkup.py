from pathlib import Path
import json
import structlog
from logging import INFO, DEBUG

from virtool_cli.accessions.helpers import get_listings, split_pathname

logger = structlog.get_logger()


def run(catalog: Path, debugging: bool = False):
    """
    Check catalog for outstanding issues

    :param catalog: Path to an accession catalog directory
    :param debugging: Debugging flag
    """

    filter_class = DEBUG if debugging else INFO
    structlog.configure(
        wrapper_class=structlog.make_filtering_bound_logger(filter_class))
    
    if not catalog.exists():
        logger.critical('Catalog directory not found at path', path=str(catalog))
        return
    
    find_duplicates(catalog)
    

def find_duplicates(catalog):
    logger = structlog.get_logger()

    unique_taxids = []

    for listing_path in get_listings(catalog):
        [ taxid, otu_id ] = split_pathname(listing_path)

        logger = logger.bind(
            path=str(listing_path.relative_to(catalog.parent)), 
            taxid=taxid, 
            otu_id=otu_id
        )

        if taxid in unique_taxids:
            logger.warning("Duplicate found!")
        else:
            unique_taxids.add(taxid)
    
    # logger.info(f'{unique_taxids}')

if __name__ == '__main__':
    debug = True
    
    REPO_DIR = '/Users/sygao/Development/UVic/Virtool/Repositories'
    
    project_path = Path(REPO_DIR) / 'ref-plant-viruses'
    src_path = project_path / 'src'
    # catalog_path = project_path / '.cache/catalog'
    catalog_path = Path(REPO_DIR) / 'ref-fetched-accessions/src'

    run(catalog=catalog_path, debugging=True)
