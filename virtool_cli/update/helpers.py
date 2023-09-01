from pathlib import Path
from structlog import BoundLogger

from virtool_cli.utils.ref import get_otu_paths
from virtool_cli.accessions.helpers import search_by_otu_id

def generate_branchname(dirname: str) -> str:
    """
    Creates an automatically-generated branchname
    """
    return 'auto__' + dirname

def filter_catalog(
    src_path: Path, catalog_path: Path,
    logger: BoundLogger
) -> list:
    """
    Return paths for cached accession catalogues that are included in the source reference

    :param src_path: Path to a reference directory
    :param catalog_path: Path to an accession record directory
    :return: A list of paths to relevant listings in the accession catalog
    """
    otu_paths = get_otu_paths(src_path)
    included_listings = []
    
    for path in otu_paths:
        try:
            [ _, otu_id ] = (path.name).split('--')
        except Exception as e:
            logger.exception(f'{e}')

        included_listings.append(
            search_by_otu_id(otu_id, catalog_path)
        )
    
    return included_listings
