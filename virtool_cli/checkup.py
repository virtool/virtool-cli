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

base_logger = get_logger()


def run(src_path: Path, debugging: bool = False):
    """
    CLI entry point for doctor.checkup.run()

    :param src_path: Path to a given reference directory
    :param debugging: Enables verbose logs for debugging purposes
    """
    structlog.configure(wrapper_class=DEBUG_LOGGER if debugging else DEFAULT_LOGGER)
    logger = base_logger.bind(src=str(src_path), verbose=debugging)

    logger.info("Running checks on reference directory...")

    logger.debug("Debug flag is enabled")

    if check_reference(src_path):
        logger.info("Reference folder is valid.", valid_result=True)
    else:
        logger.warning("Reference folder is invalid.", valid_result=False)


def check_reference(src_path: Path) -> bool:
    """
    Inspects all files in the reference directory and logs a warning if a problem is found.

    :param src_path: Path to a given reference directory
    """
    logger = base_logger

    is_valid = True

    for otu_path in get_otu_paths(src_path):
        otu_logger = logger.bind(otu=otu_path.name)

        if not check_otu(otu_path, otu_logger):
            is_valid = False

        isolate_paths = get_isolate_paths(otu_path)
        if not isolate_paths:
            otu_logger.error("No accession data in OTU", otu=otu_path.name)
            is_valid = False

        else:
            for isolate_path in isolate_paths:
                if not check_isolate(isolate_path, otu_logger):
                    is_valid = False

                sequence_paths = get_sequence_paths(isolate_path)
                if not sequence_paths:
                    otu_logger.error(
                        "No sequences in isolate", isolate=isolate_path.name
                    )
                    is_valid = False

                for sequence_path in sequence_paths:
                    if not check_sequence(sequence_path, otu_logger):
                        is_valid = False

    if is_valid:
        logger.info("No issues found in src", src=str(src_path))
    else:
        logger.error(
            "Found issues in the reference. Refer to the reference for suggested corrections."
        )

    return is_valid


def check_otu(otu_path: Path, logger: BoundLogger = base_logger) -> bool:
    """
    Checks OTU metadata for bad content and logs a warning if problems are found.

    :param otu_path: A path to an OTU directory under a src reference directory
    :param logger: Optional entry point for a shared BoundLogger
    :return: True if an OTU's metadata is valid and there are sequences assigned under it,
        False if a problem is found
    """
    otu_metadata_path = otu_path / "otu.json"

    if not otu_metadata_path.exists():
        logger.error(f"{otu_path.name}: missing metadata")
        return False

    with open(otu_metadata_path, "r") as f:
        otu = json.load(f)

    if not otu.get("_id", ""):
        logger.error(f"OTU metadata lacks a unique Virtool ID")
        return False

    if not otu.get("schema", []):
        logger.error(f"OTU metadata lacks a schema")
        return False

    if otu.get("taxid", -1) < 0:
        logger.error(f"OTU metadata lacks a NCBI Taxonomy UID")
        return False

    return True


def check_isolate(isolate_path: Path, logger: BoundLogger = base_logger) -> bool:
    """
    Checks isolate metadata for bad content and logs a warning if problems are found.

    :param isolate_path: A path to an isolate directory under an OTU directory
    :param logger: Optional entry point for a shared BoundLogger
    :return: True if an isolate's metadata is valid, False if a problem is found
    """
    logger = logger.bind(isolate_id=isolate_path.stem)

    isolate_metadata_path = isolate_path / "isolate.json"

    if not isolate_metadata_path.exists():
        logger.error(f"{isolate_path.name}: missing metadata")
        return False

    with open(isolate_metadata_path, "r") as f:
        isolate = json.load(f)

    if not isolate.get("id", ""):
        logger.error(f"Isolate metadata lacks a unique Virtool ID")
        return False

    if not isolate.get("source_type", ""):
        logger.error(f"Isolate metadata lacks a type")
        return False

    if not isolate.get("source_name", []):
        logger.error(f"Isolate metadata lacks a name")
        return False

    if "default" in isolate:
        if type(isolate["default"]) != bool:
            logger.error(f"Isolate metadata lacks a valid default flag")
            return False
    else:
        logger.error(f"Isolate metadata lacks a valid default flag")
        return False

    return True


def check_sequence(sequence_path: Path, logger: BoundLogger = base_logger) -> bool:
    """
    Checks a sequence file for bad content and logs a warning if problems are found.

    :param sequence_path: A path to a sequence file under an isolate directory
    :param logger: Optional entry point for a shared BoundLogger
    :return: True if a sequence's metadata is valid, False if a problem is found
    """
    logger = logger.bind(sequence_id=sequence_path.stem)

    sequence = json.loads(sequence_path.read_text())

    accession = sequence["accession"]

    if not verify_accession(accession):
        logger.error(
            f"Accession '{accession}' contains invalid characters or capitalization"
        )
        return False

    if "." not in accession:
        logger.warning(f"Version not included in accession={accession}")
        return False

    return True


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
