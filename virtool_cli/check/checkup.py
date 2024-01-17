import json
import re
from pathlib import Path
from structlog import get_logger, BoundLogger

from virtool_cli.utils.reference import get_isolate_paths, get_sequence_paths

OTU_KEYS = ["_id", "name", "abbreviation", "schema", "taxid"]
ISOLATE_KEYS = ["id", "source_type", "source_name", "default"]
SEQUENCE_KEYS = ["_id", "accession", "definition", "host", "sequence"]


def check_otu(otu_path: Path, logger: BoundLogger = get_logger()) -> bool:
    """
    :param otu_path: A path to an OTU directory under a src reference directory
    :param logger: Optional entry point for a shared BoundLogger
    :return: True if an OTU's metadata is valid,
        False if a problem is found
    """
    is_valid = True

    try:
        [_, otu_id] = otu_path.name.split("--")
    except ValueError as e:
        logger.exception(e)
        return False

    otu_logger = logger.bind(otu_id=otu_id)

    otu_logger.debug(f"Running checks on OTU {otu_id}...")

    if not validate_otu(otu_path, otu_logger):
        is_valid = False

    isolate_paths = get_isolate_paths(otu_path)
    if not isolate_paths:
        otu_logger.error("No accession data in OTU", otu=otu_path.name)
        is_valid = False

    else:
        for isolate_path in isolate_paths:
            if not validate_isolate(isolate_path, logger):
                is_valid = False

            sequence_paths = get_sequence_paths(isolate_path)
            if not sequence_paths:
                otu_logger.error("No sequences in isolate", isolate=isolate_path.name)
                is_valid = False

            for sequence_path in sequence_paths:
                if not validate_sequence(sequence_path, logger):
                    is_valid = False

    return is_valid


def validate_otu(otu_path: Path, logger: BoundLogger = get_logger()) -> bool:
    """
    Checks OTU metadata for bad content and logs a warning if problems are found.

    :param otu_path: A path to an OTU directory under a src reference directory
    :param logger: Optional entry point for a shared BoundLogger
    :return: True if an OTU's metadata is valid and there are sequences assigned under it,
        False if a problem is found
    """
    otu_metadata_path = otu_path / "otu.json"

    try:
        with open(otu_metadata_path, "r") as f:
            otu = json.load(f)
    except FileNotFoundError:
        logger.error(f"Missing OTU metadata file")
        return False
    except json.JSONDecodeError:
        logger.error(f"OTU metadata is not JSON-parseable")
        return False

    is_valid = True

    if not otu.get("_id", ""):
        logger.error(f"OTU metadata lacks a unique Virtool ID")
        is_valid = False

    if not otu.get("schema", []):
        logger.warning(f"OTU metadata lacks a schema")

    if otu.get("taxid", None) is None:
        logger.warning(f"OTU metadata lacks a NCBI Taxonomy UID")

    return is_valid


def validate_isolate(isolate_path: Path, logger: BoundLogger = get_logger()) -> bool:
    """
    Checks isolate metadata for bad content and logs a warning if problems are found.

    :param isolate_path: A path to an isolate directory under an OTU directory
    :param logger: Optional entry point for a shared BoundLogger
    :return: True if an isolate's metadata is valid, False if a problem is found
    """
    isolate_id = isolate_path.stem
    logger = logger.bind(isolate_id=isolate_id)
    logger.debug(f"Running checks on isolate {isolate_id}...")

    isolate_metadata_path = isolate_path / "isolate.json"

    logger = logger.bind(path=str(isolate_metadata_path))

    # Run file checks
    try:
        with open(isolate_metadata_path, "r") as f:
            isolate = json.load(f)
    except FileNotFoundError:
        logger.error(f"Missing isolate data")
        return False
    except json.JSONDecodeError:
        logger.error(f"Isolate metadata is not JSON-parseable")
        return False

    is_valid = True

    # Run isolate metadata checks
    if not isolate.get("id", ""):
        logger.error(f"Isolate metadata lacks a unique Virtool ID in the JSON")
        is_valid = False

    if not isolate.get("source_type", ""):
        logger.warning(f"Isolate metadata lacks a source type")

    if not isolate.get("source_name", ""):
        logger.warning(f"Isolate metadata lacks a source name")

    if "default" in isolate:
        if type(isolate["default"]) is not bool:
            logger.error(f"Isolate metadata lacks a valid default flag")
            is_valid = False
    else:
        logger.error(f"Isolate metadata lacks a valid default flag")
        is_valid = False

    return is_valid


def validate_sequence(sequence_path: Path, logger: BoundLogger = get_logger()) -> bool:
    """
    Checks a sequence file for bad content and logs a warning if problems are found.

    :param sequence_path: A path to a sequence file under an isolate directory
    :param logger: Optional entry point for a shared BoundLogger
    :return: True if a sequence's metadata is valid, False if a problem is found
    """
    is_valid = True
    sequence_location = str(sequence_path.relative_to(sequence_path.parents[3]))

    sequence_id = sequence_path.stem
    logger = logger.bind(
        sequence_id=sequence_id,
        path=sequence_location,
    )
    logger.debug(f"Running checks on sequence {sequence_id}...")

    # Run file checks
    try:
        sequence_data = json.loads(sequence_path.read_text())
    except FileNotFoundError:
        logger.error(f"Missing sequence data")
        return False
    except json.JSONDecodeError:
        logger.error(f"Sequence metadata is not JSON-parseable")
        return False

    accession = sequence_data["accession"]

    # Run metadata checks
    if not sequence_data.get("_id", ""):
        logger.error(f"Sequence metadata lacks a unique Virtool ID in the JSON")
        is_valid = False

    if not verify_accession(accession):
        logger.error(
            f"Accession '{accession}' contains invalid characters or capitalization"
        )
        is_valid = False

    # if "." not in accession:
    #     logger.warning(f"Version not included in accession={accession}")
    #     return False

    logger.debug(f"Sequence {sequence_id} is valid")

    return is_valid


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
