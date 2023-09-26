import json
from pathlib import Path
import re
import structlog
from structlog import get_logger, BoundLogger

from virtool_cli.utils.logging import DEFAULT_LOGGER, DEBUG_LOGGER
from virtool_cli.utils.reference import (
    get_otu_paths,
    get_isolate_paths,
    get_sequence_paths,
)


def run(src_path: Path, debugging: bool = False):
    """
    CLI entry point for doctor.checkup.run()

    :param src_path: Path to a given reference directory
    :param debugging: Enables verbose logs for debugging purposes
    """
    structlog.configure(wrapper_class=DEBUG_LOGGER if debugging else DEFAULT_LOGGER)
    logger = get_logger().bind(src=str(src_path), verbose=debugging)

    logger.info("Running checks on reference directory...")

    check_reference(src_path)


def check_reference(src_path: Path):
    """
    Inspects all files in the reference and logs a warning if a problem is found

    :param src_path: Path to a given reference directory
    """
    for otu_path in get_otu_paths(src_path):
        logger = base_logger.bind(otu_path=str(otu_path.relative_to(src_path)))

        check_otu(otu_path, logger)

        for isolate_path in get_isolate_paths(otu_path):
            for sequence_path in get_sequence_paths(isolate_path):
                check_sequence(sequence_path, logger)


def check_otu(otu_path: Path, logger: BoundLogger = get_logger()):
    """
    Checks an OTU file for bad content and logs a warning if problems are found

    :param otu_path: A path to an OTU directory under a src reference directory
    :param logger: Optional entry point for a shared BoundLogger
    """
    with open(otu_path / "otu.json", "r") as f:
        otu = json.load(f)

    if not otu.get("schema", []):
        logger.warning(f"{otu['_id']}: missing schema")


def check_sequence(sequence_path: Path, logger: BoundLogger = get_logger()):
    """
    Checks a sequence file for bad content and logs a warning if problems are found

    :param sequence_path: A path to a sequence file under an isolate directory
    :param logger: Optional entry point for a shared BoundLogger
    """
    logger = logger.bind(sequence_id=sequence_path.stem)

    sequence = json.loads(sequence_path.read_text())

    accession = sequence["accession"]

    if not verify_accession(accession):
        logger.error(
            f"Accession '{accession}' contains invalid characters or capitalization"
        )

    if "." not in accession:
        logger.warning(f"Version not included in accession={accession}")


def verify_accession(accession: str):
    """
    Returns True if the accession matches NCBI standards
    The accession must be free of characters other than digits,
    capital letters A-Z, '.', and '_'

    :param accession: An accession to be inspected
    :return: Boolean based on whether the accession is free of invalid characters
    """
    if re.search(r"([^A-Z_.0-9])", accession) is None:
        return True

    return False
