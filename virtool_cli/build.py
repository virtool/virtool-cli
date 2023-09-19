import json
from pathlib import Path
import arrow
import logging
from structlog import BoundLogger

from virtool_cli.utils.logging import base_logger
from virtool_cli.utils.reference import (
    get_otu_paths,
    get_isolate_paths,
    get_sequence_paths,
    is_v1,
)

OTU_KEYS = ["_id", "name", "abbreviation", "schema", "taxid"]

ISOLATE_KEYS = ["id", "source_type", "source_name", "default"]

SEQUENCE_KEYS = ["_id", "accession", "definition", "host", "sequence"]


def run(
    src_path: Path, output: Path, indent: bool, version: str, debugging: bool = False
):
    """
    Build a Virtool reference JSON file from a data directory.

    :param src_path: Path to database src directory
    :param output: The output path for the reference.json file
    :param indent: A flag to indicate whether the output file should be indented
    :param version: The version string to include in the reference.json file
    :param debugging: Enables verbose logs for debugging purposes
    """
    filter_class = logging.DEBUG if debugging else logging.INFO
    logging.basicConfig(
        format="%(message)s",
        level=filter_class,
    )

    if is_v1(src_path):
        base_logger.critical(
            "This reference database is a deprecated v1 reference."
            + 'Run "virtool ref migrate" before trying again.',
            src=str(src_path),
        )
        return

    try:
        build_from_src(src_path, output, indent, version)
    except FileNotFoundError as e:
        base_logger.exception(e)
    except Exception as e:
        base_logger.exception(e)


def build_from_src(src_path, output, indent, version):
    """
    Build a Virtool reference JSON file from a data directory.

    :param src_path: Path to database src directory
    :param output: The output path for the reference.json file
    :param indent: A flag to indicate whether the output file should be indented
    :param version: The version string to include in the reference.json file
    """
    logger = base_logger.bind(src=str(src_path), output=str(output))

    data = {"data_type": "genome", "organism": ""}

    if meta := parse_meta(src_path):
        logger.debug(
            "Metadata parsed", meta=meta, metadata_path=str(src_path / "meta.json")
        )

        data["data_type"] = meta.get("data_type", "genome")
        data["organism"] = (meta.get("organism", ""),)

    else:
        logger.warning(
            f'Metadata file not found at {src_path / "meta.json"}',
            metadata_path=str(src_path / "meta.json"),
        )

    otus = []

    for otu_path in get_otu_paths(src_path):
        try:
            otu = parse_otu_contents(otu_path)
        except (FileNotFoundError, json.JSONDecodeError) as e:
            logger.critical("Reference data at src_path is invalid.")
            logger.exception(e)
            return
        except Exception as e:
            logger.exception(e)
            return

        otus.append(otu)

        logger.debug(
            f"Added {otu['_id']} to reference data",
            path=str(otu_path.relative_to(src_path)),
        )

    data.update(
        {"otus": otus, "name": version, "created_at": arrow.utcnow().isoformat()}
    )

    with open(output, "w") as f:
        json.dump(data, f, indent=4 if indent else None, sort_keys=True)

    logger.info("Reference file built at output")


def parse_meta(src_path: Path) -> dict:
    """
    Deserializes and returns meta.json if found, else returns an empty dictionary.

    :param src_path: Path to database src directory
    :return: The deserialized meta.json object or an empty dictionary
    """
    try:
        with open(src_path / "meta.json", "r") as f:
            return json.load(f)
    except FileNotFoundError:
        return {}


def parse_alpha(alpha: Path) -> list:
    """
    Generates and returns a list with every OTU in the directory of the given alpha.
    Deprecated as of v2.

    :param alpha: Path to a given alpha directory in a reference
    :return: A list containing all the OTU in the given directory
    """
    return [otu for otu in alpha.iterdir() if otu.is_dir()]


def parse_otu_contents(otu_path: Path, logger: BoundLogger = base_logger) -> dict:
    """
    Traverses, deserializes and returns all data under an OTU directory.

    :param otu_path: Path to a OTU directory
    :param logger: Optional entry point for a shared BoundLogger
    :return: All isolate and sequence data under an OTU,
        deserialized and compiled in a dict
    """
    with open(otu_path / "otu.json", "r") as f:
        otu = json.load(f)

    logger = logger.bind(path=otu_path, otu_id=otu["_id"])

    isolates = []
    for isolate_path in get_isolate_paths(otu_path):
        with open(isolate_path / "isolate.json", "r") as f:
            isolate = json.load(f)

        logger = logger.bind(isolate_id=isolate["id"])

        sequences = []
        for sequence_path in get_sequence_paths(isolate_path):
            sequence = json.loads(sequence_path.read_text())

            sequences.append(sequence)
            logger.debug(
                f"Added sequence {sequence.get('accession')} under id={sequence.get('_id')}"
            )

        isolate["sequences"] = sequences

        isolates.append(isolate)

    otu["isolates"] = isolates

    return otu
