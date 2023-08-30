from pathlib import Path
import structlog
import logging

from virtool_cli.utils.logging import base_logger
from virtool_cli.utils.ref import parse_otu, generate_otu_dirname

def run(src_path: Path, debugging: bool = False):
    """
    :param src_path: Path to a src database directory
    :param debugging: Enables verbose logs for debugging purposes
    """
    filter_class = logging.DEBUG if debugging else logging.INFO
    logging.basicConfig(
        format="%(message)s",
        level=filter_class,
    )

    logger = base_logger.bind(src_path=str(src_path))
    
    if not [alpha for alpha in src_path.glob('[a-z]')]:
        logger.info("Reference in src_path is already v2")
        return
    
    logger.info("Converting reference in src_path to v2...")
    try:
        flatten_src(src_path)
    except Exception as e:
        logger.error('Error occured during conversion')
        logger.exception(f'Error: {e}')

def flatten_src(src_path: Path):
    """
    Traverses through binning directories a to z and reassigns 
    all OTU directories to have src_path as a direct parent,
    then deletes the binning directory

    :param src_path: Path to a src database directory
    """

    for alpha in [alpha for alpha in src_path.glob('[a-z]')]:

        otu_paths = [otu for otu in alpha.iterdir() if otu.is_dir()]
        
        for otu_path in otu_paths:
            otu = parse_otu(otu_path)
            new_name = generate_otu_dirname(
                otu.get('name'), 
                otu.get('_id')
            )
            new_path = src_path / new_name
            otu_path.rename(new_path)
            print(new_path)

        # Delete alpha bin
        try:
            for chaff in alpha.iterdir():
                chaff.unlink()
        except Exception as e:
            return e
        alpha.rmdir()