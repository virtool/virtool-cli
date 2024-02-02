import json
import re
from pathlib import Path

from structlog import get_logger

from virtool_cli.utils.reference import get_isolate_paths, get_sequence_paths

logger = get_logger("check_otu")


def check_otu(otu_path: Path) -> bool:
    """:param otu_path: A path to an OTU directory under a src reference directory
    :return: True if an OTU's metadata is valid,
        False if a problem is found
    """
    is_valid = True

    try:
        [_, otu_id] = otu_path.name.split("--")
    except ValueError as e:
        logger.exception(e)
        return False

    log = logger.bind(otu_id=otu_id)

    log.debug(f"Running checks on OTU {otu_id}...")

    if not validate_otu(otu_path):
        is_valid = False

    isolate_paths = get_isolate_paths(otu_path)

    if not isolate_paths:
        log.error("No accession data in OTU", otu=otu_path.name)
        is_valid = False
    else:
        for isolate_path in isolate_paths:
            if not validate_isolate(isolate_path):
                is_valid = False

            sequence_paths = get_sequence_paths(isolate_path)

            if sequence_paths:
                if any(
                    not validate_sequence(sequence_path)
                    for sequence_path in sequence_paths
                ):
                    is_valid = False

            else:
                log.error("No sequences in isolate", isolate=isolate_path.name)
                is_valid = False

    if is_valid:
        log.info("OTU is valid")

    return is_valid


def validate_otu(otu_path: Path) -> bool:
    """Checks OTU metadata for bad content and logs a warning if problems are found.

    :param otu_path: the path to an otu directory to validate
    :return: whether the otu is valid
    """
    otu_metadata_path = otu_path / "otu.json"

    try:
        with open(otu_metadata_path) as f:
            otu = json.load(f)
    except FileNotFoundError:
        logger.error("Could not find otu.json")
        return False
    except json.JSONDecodeError:
        logger.error("Could not decode otu.json")
        return False

    is_valid = True

    if "_id" not in otu:
        logger.error("No _id in otu.json")
        is_valid = False

    if "schema" not in otu:
        logger.warning("No schema in otu.json")
        is_valid = False

    if not otu.get("taxid"):
        logger.warning("No taxid in otu.json")
        is_valid = False

    if is_valid:
        logger.info("OTU metadata is valid")

    return is_valid


def validate_isolate(isolate_path: Path) -> bool:
    """Checks isolate metadata for bad content and logs a warning if problems are found.

    :param isolate_path: A path to an isolate directory under an OTU directory
    :return: whether the isolate metadata is valid
    """
    isolate_id = isolate_path.stem
    isolate_json_path = isolate_path / "isolate.json"

    log = logger.bind(id=isolate_id, json_path=str(isolate_json_path))

    try:
        with open(isolate_json_path) as f:
            isolate = json.load(f)
    except FileNotFoundError:
        logger.warning("Could not find isolate.json")
        return False
    except json.JSONDecodeError:
        log.warning("Could not parse isolate.json")
        return False

    is_valid = True

    # Run isolate metadata checks
    if not isolate.get("id", ""):
        log.warning("Invalid or missing id field in isolate.json")
        is_valid = False

    if not isolate.get("source_type", ""):
        log.warning("Invalid or missing source_type field in isolate.json")
        is_valid = False

    if not isolate.get("source_name"):
        log.warning("Invalid or missing source_name field in isolate.json")
        is_valid = False

    if "default" not in isolate or not isinstance(isolate["default"], bool):
        log.error("Invalid or missing default field in isolate.json")
        is_valid = False

    if is_valid:
        log.info("Isolate metadata is valid")

    return is_valid


def validate_sequence(sequence_path: Path) -> bool:
    """Checks a sequence file for bad content and logs a warning if problems are found.

    :param sequence_path: A path to a sequence file under an isolate directory
    :return: whether the sequence metadata is valid
    """
    is_valid = True

    sequence_id = sequence_path.stem

    log = logger.bind(
        sequence_id=sequence_id,
        path=str(sequence_path.relative_to(sequence_path.parents[3])),
    )

    try:
        sequence_data = json.loads(sequence_path.read_text())
    except FileNotFoundError:
        log.error("Missing sequence data")
        return False
    except json.JSONDecodeError:
        log.error("Sequence metadata is not JSON-parseable")
        return False

    accession = sequence_data["accession"]

    # Run metadata checks
    if not sequence_data.get("_id", ""):
        log.error("Sequence metadata lacks a unique Virtool ID in the JSON")
        is_valid = False

    if not verify_accession(accession):
        log.error(
            f"Accession '{accession}' contains invalid characters or capitalization",
        )
        is_valid = False

    if is_valid:
        log.info("Sequence metadata is valid")

    return is_valid


def verify_accession(accession: str) -> bool:
    """Returns True if the accession matches NCBI standards
    The accession must be free of characters other than digits,
    capital letters A-Z, '.', and '_'

    :param accession: An accession to be inspected
    :return: Boolean based on whether the accession is free of invalid characters
    """
    return re.search(r"([^A-Z_.0-9])", accession) is None
