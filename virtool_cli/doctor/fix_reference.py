from pathlib import Path
import logging

from virtool_cli.utils.logging import base_logger
from virtool_cli.utils.reference import get_otu_paths
from virtool_cli.doctor.fix_otu import repair_otu_data

def run(src_path: Path, debugging: bool = False):
    """
    CLI entry point for doctor.fix_reference.run()

    :param src_path: Path to a given reference directory
    :param debugging: Enables verbose logs for debugging purposes
    """
    filter_class = logging.DEBUG if debugging else logging.INFO
    logging.basicConfig(
        format="%(message)s",
        level=filter_class,
    )

    repair_data(src_path)

def repair_data(src_path: Path):
    """
    Fixes incorrect data in all OTUs under a reference

    :param src_path: Path to a given reference directory
    """
    otu_paths = get_otu_paths(src_path)
    for otu_path in otu_paths:
        logger = base_logger.bind(otu_path=str(otu_path.relative_to(src_path)))

        repair_otu_data(otu_path, logger)