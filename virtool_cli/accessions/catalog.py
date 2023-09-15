from pathlib import Path
import asyncio
import logging

from virtool_cli.utils.logging import base_logger
from virtool_cli.accessions.initialize import initialize
from virtool_cli.accessions.update import update


def run(src_path: Path, catalog_path: Path, debugging: bool = False):
    """
    CLI entry point for accession.catalog.run()

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
        base_logger.info(
            'Initializing .cache/catalog directory...', 
            catalog=str(catalog_path)
        )
        catalog_path.mkdir()

    get_catalog(src_path, catalog_path)

def get_catalog(src: Path, catalog: Path):
    """
    Runs accession.catalog.update if a catalog directory is found, 
    else runs accession.catalog.update

    :param src: Path to a reference directory
    :param catalog: Path to an accession catalog directory
    """
    logger = base_logger.bind(catalog=str(catalog))

    catalog_paths = list(catalog.glob('*--*.json'))
    
    if catalog_paths:
        logger.info(
            'Catalog is already populated.', catalogued=True)
        
        logger.info("Updating listings...", task="update")
        asyncio.run(update(src, catalog))
        
    else:
        logger.info(
            'Catalog is empty.', catalogued=False)

        logger.info("Initializing listings...", task="initialize")
        asyncio.run(initialize(src, catalog))