from pathlib import Path
import asyncio
import structlog
from logging import INFO, DEBUG

from virtool_cli.accessions.initialize import initialize
from virtool_cli.accessions.update import update

base_logger = structlog.get_logger()


def run(src_path: Path, catalog_path: Path, debugging: bool = False):
    """
    :param src: Path to a reference directory
    :param catalog: Path to an accession catalog directory
    :param debugging: Enables verbose logs for debugging purposes
    """

    filter_class = DEBUG if debugging else INFO
    structlog.configure(
        wrapper_class=structlog.make_filtering_bound_logger(filter_class))
    
    if not catalog_path.exists():
        base_logger.info('Initializing .cache/catalog directory...', catalog=str(catalog_path))
        catalog_path.mkdir()

    get_catalog(src_path, catalog_path)

def get_catalog(src: Path, catalog: Path):
    """
    :param src: Path to a reference directory
    :param catalog: Path to an accession catalog directory
    """
    logger = base_logger.bind(catalog=str(catalog))

    catalog_paths = [record for record in catalog.glob('*.json')]
    
    if catalog_paths:
        logger.info(
            'Catalog is already populated.', catalogued=True)
        
        update_catalog(src, catalog, logger)
        
    else:
        logger.info(
            'Catalog is empty.', catalogued=False)

        generate_catalog(src, catalog, logger)

def generate_catalog(
    src_path: Path, catalog_path: Path, 
    logger: structlog.BoundLogger = base_logger
):
    """
    """
    logger.info("Initializing listings...", task="initialize")
    asyncio.run(initialize(src_path, catalog_path))

def update_catalog(
    src_path: Path, catalog_path: Path, 
    logger: structlog.BoundLogger = base_logger
):
    """
    """
    logger.info("Updating listings...", task="initialize")
    asyncio.run(update(src_path, catalog_path))