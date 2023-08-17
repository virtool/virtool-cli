from pathlib import Path
import asyncio
import structlog
import logging

from virtool_cli.utils.logging import base_logger
from virtool_cli.accessions.initialize import initialize
from virtool_cli.accessions.update import update


def run(src_path: Path, catalog_path: Path, debugging: bool = False):
    """
    :param src: Path to a reference directory
    :param catalog: Path to an accession catalog directory
    :param debugging: Enables verbose logs for debugging purposes
    """
    filter_class = logging.DEBUG if debugging else logging.INFO
    logging.basicConfig(
        format="%(message)s",
        level=filter_class,
    )
    
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
    :param catalog: Path to an accession catalog directory
    :param logger: Structlog logger
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