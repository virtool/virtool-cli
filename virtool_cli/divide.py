import json
from pathlib import Path
import shutil
import logging

from virtool_cli.utils.logging import base_logger
from virtool_cli.utils.ref import generate_otu_dirname

OTU_KEYS = ["_id", "name", "abbreviation", "schema", "taxid"]

ISOLATE_KEYS = ["id", "source_type", "source_name", "default"]

SEQUENCE_KEYS = (
    "_id",
    "accession",
    "definition",
    "host",
    "segment",
    "sequence",
    "target",
)


def run(src_path: Path, output: Path, debugging: bool = False):
    """
    Divide a reference.json file from Virtool into a src tree.

    :param src_path: Path to a reference.json file
    :param output: Path to the where the src tree should be generated
    :param debugging: Enables verbose logs for debugging purposes
    """
    filter_class = logging.DEBUG if debugging else logging.INFO
    logging.basicConfig(
        format="%(message)s",
        level=filter_class,
    )
    
    shutil.rmtree(output, ignore_errors=True)
    output.mkdir()

    with open(src_path, "r") as export_handle:
        data = json.load(export_handle)

        for otu in data.get("otus"):
            otu_path = build_otu(output, otu)

            isolates = otu.pop("isolates")

            for isolate in isolates:

                isolate_path = build_isolate(otu_path, isolate)

                sequences = isolate.pop("sequences")

                for sequence in sequences:
                    build_sequence(isolate_path, sequence)

        with open(output / "meta.json", "w") as f:
            json.dump({"data_type": data["data_type"], "organism": data["organism"]}, f)


def build_otu(output: Path, otu: dict) -> str:
    """
    Creates a directory for all OTUs that begin with a particular
    letter if it doesn't already exist. Generates a directory for a
    given OTU and copies key information about it to a otu.json file.

    :param output: Path to the where the src tree should be generated
    :param otu: Dictionary of an OTU
    :return: Path to a newly generated OTU directory
    """
    # lower_name = otu["name"].lower()

    otu_path = output / generate_otu_dirname(otu.get('name'), otu.get('_id'))
    otu_path.mkdir()

    with open(otu_path / "otu.json", "w") as f:
        if "schema" not in otu:
            otu["schema"] = list()

        json.dump({key: otu.get(key) for key in OTU_KEYS}, f, indent=4)

    return otu_path


def build_isolate(otu_path: Path, isolate: dict) -> str:
    """
    Creates a directory for a given isolate and generates
    a isolate.json with key information about it.

    :param otu_path: A path to a specified OTU
    :param isolate: A dictionary containing isolate information
    :return: A path to a newly generated isolate directory
    """
    isolate_path = otu_path / isolate.get("id")
    isolate_path.mkdir()

    with open(isolate_path / "isolate.json", "w") as f:
        json.dump({key: isolate[key] for key in ISOLATE_KEYS}, f, indent=4)

    return isolate_path


def build_sequence(isolate_path: Path, sequence: dict):
    """
    Generates a JSON file for one of the isolate's sequences

    :param isolate_path: A path to a specified isolate
    :param sequence: A dictionary containing information on one of the isolates' sequences
    """
    with open(isolate_path / "{}.json".format(sequence["_id"]), "w") as f:
        json.dump(
            {key: sequence[key] for key in SEQUENCE_KEYS if key in sequence},
            f,
            indent=4,
        )
