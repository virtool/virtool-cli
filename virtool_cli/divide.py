import json
from pathlib import Path
import shutil
import structlog
from structlog import get_logger

from virtool_cli.utils.logging import DEBUG_LOGGER, DEFAULT_LOGGER
from virtool_cli.utils.reference import generate_otu_dirname

OTU_KEYS = ["_id", "name", "abbreviation", "schema", "taxid"]

ISOLATE_KEYS = ["id", "source_type", "source_name", "default"]

SEQUENCE_KEYS = [
    "_id",
    "accession",
    "definition",
    "host",
    "segment",
    "sequence",
    "target",
]


def run(reference_path: Path, output_path: Path, debugging: bool = False):
    """
    Divide a reference.json file from Virtool into a src tree structure.

    :param reference_path: Path to a reference.json file
    :param output_path: Path to the where the src tree should be generated
    :param debugging: Enables verbose logs for debugging purposes
    """
    structlog.configure(wrapper_class=DEBUG_LOGGER if debugging else DEFAULT_LOGGER)
    logger = get_logger().bind(reference=str(reference_path), output=str(output_path))

    logger.info(f"Dividing {output_path.name} into {reference_path.name}...")

    shutil.rmtree(output_path, ignore_errors=True)
    output_path.mkdir()

    with open(reference_path, "r") as export_handle:
        data = json.load(export_handle)

        for otu in data.get("otus"):
            otu_path = build_otu(output_path, otu)

            logger.debug(
                f"Built {otu['_id']}",
                otu={"_id": otu["_id"], "path": str(otu_path.relative_to(output_path))},
            )

            isolates = otu.pop("isolates")

            for isolate in isolates:
                isolate_path = build_isolate(otu_path, isolate)

                sequences = isolate.pop("sequences")

                for sequence in sequences:
                    build_sequence(isolate_path, sequence)

        with open(output_path / "meta.json", "w") as f:
            json.dump(
                {"data_type": data["data_type"], "organism": data["organism"]},
                f,
                sort_keys=True,
            )


def build_otu(src_path: Path, otu: dict) -> Path:
    """
    Creates a directory for all OTUs that begin with a particular
    letter if it doesn't already exist. Generates a directory for a
    given OTU and copies key information about it to a otu.json file.

    :param src_path: Path to the where the src tree should be generated
    :param otu: Dictionary of an OTU
    :return: Path to a newly generated OTU directory
    """
    if "schema" not in otu:
        otu["schema"] = []

    otu_dirname = generate_otu_dirname(otu.get("name"), otu.get("_id"))
    otu_path = src_path / otu_dirname
    otu_path.mkdir()

    with open(otu_path / "otu.json", "w") as f:
        json.dump({key: otu.get(key) for key in OTU_KEYS}, f, indent=4, sort_keys=True)

    return otu_path


def build_isolate(otu_path: Path, isolate: dict) -> Path:
    """
    Creates a directory for a given isolate and generates
    a isolate.json with key information about it.

    :param otu_path: A path to an OTU directory under a src reference directory
    :param isolate: A dictionary containing isolate information
    :return: A path to a newly generated isolate directory
    """
    isolate_path = otu_path / isolate.get("id")
    isolate_path.mkdir()

    with open(isolate_path / "isolate.json", "w") as f:
        json.dump(
            {key: isolate[key] for key in ISOLATE_KEYS}, f, indent=4, sort_keys=True
        )

    return isolate_path


def build_sequence(isolate_path: Path, sequence: dict):
    """
    Generates a JSON file for a sequence under an isolate directory

    :param isolate_path: A path to a specified isolate directory under an OTU directory
    :param sequence: A dictionary containing sequence information
    """
    with open(isolate_path / f"{sequence.get('_id')}.json", "w") as f:
        json.dump(
            {key: sequence[key] for key in SEQUENCE_KEYS if key in sequence},
            f,
            indent=4,
            sort_keys=True,
        )
