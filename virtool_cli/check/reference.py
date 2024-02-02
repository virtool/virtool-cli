from pathlib import Path

import structlog
from structlog import get_logger

from virtool_cli.check.otu import check_otu
from virtool_cli.utils.reference import get_otu_paths

base_logger = structlog.get_logger()


def check_reference(path: Path) -> bool:
    """Inspects all files in the reference src path and logs a warning for each problem
    found.

    :param path: the path to the reference repository
    :return: whether the reference is valid
    """
    logger = get_logger("check_reference", path=path)
    logger.info("Running checks on reference directory...")

    ref_valid = True

    for otu_path in get_otu_paths(path / "src"):
        if not check_otu(otu_path):
            ref_valid = False
            logger.warning(f"OTU {otu_path.name} has issues", otu=otu_path.name)

    return ref_valid
