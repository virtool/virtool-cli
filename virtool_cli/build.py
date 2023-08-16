import json
from pathlib import Path
from typing import Tuple
import arrow
import structlog
import logging

import virtool_cli.utils.logging
from virtool_cli.utils.ref import get_otu_paths

OTU_KEYS = ["_id", "name", "abbreviation", "schema", "taxid"]

ISOLATE_KEYS = ["id", "source_type", "source_name", "default"]

SEQUENCE_KEYS = ["_id", "accession", "definition", "host", "sequence"]


def run(src_path: Path, output: Path, 
    indent: bool, version: str, 
    debugging: bool = False):
    """
    Build a Virtool reference JSON file from a data directory.

    :param src_path: Path to database src directory
    :param output: The output path for the reference.json file
    :param indent: A flag to indicate whether the output file should be indented
    :param version: The version string to include in the reference.json file
    """
    filter_class = logging.DEBUG if debugging else logging.INFO
    logging.basicConfig(
        format="%(message)s",
        level=filter_class,
    )

    logger = structlog.get_logger()
    logger = logger.bind(src=str(src_path), output=str(output))

    meta = parse_meta(src_path)

    data = {
        "data_type": meta.get("data_type", "genome"),
        "organism": meta.get("organism", ""),
    }

    otus = list()

    otu_paths = get_otu_paths(src_path)

    for otu_path in otu_paths:

        otu, isolate_ids = parse_otu(otu_path)

        for isolate_path in isolate_ids:

            isolate, sequence_ids = parse_isolate(isolate_path)

            for sequence_path in sequence_ids:
                with open(sequence_path, "r") as f:
                    sequence = json.load(f)

                isolate["sequences"].append(sequence)

            otu["isolates"].append(isolate)

        otus.append(otu)

    try:
        with open(output, "w") as f:
            data.update(
                {"otus": otus, "name": version, "created_at": arrow.utcnow().isoformat()}
            )

            json.dump(data, f, indent=4 if indent else None, sort_keys=True)
            logger.info('Reference file built')
    except:
        logger.critical(
            'Reference file build failed')


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
        return dict()


def parse_alpha(alpha: Path) -> list:
    """
    Generates and returns a list with every OTU in the directory of the given alpha.
    Deprecated as of v2.

    :param alpha: Path to a given alpha directory in a reference
    :return: A list containing all the OTU in the given directory
    """
    return [otu for otu in alpha.iterdir() if otu.is_dir()]


def parse_otu(otu_path: Path) -> Tuple[dict, list]:
    """
    Creates an list of isolate IDs for a given OTU

    :param otu_path: Path to a OTU directory
    :return: The dictionary and list of isolate ids for a given OTU
    """
    with open(otu_path / "otu.json", "r") as f:
        otu = json.load(f)

    otu["isolates"] = list()

    isolate_ids = [
        i for i in otu_path.iterdir() if i.name != "otu.json" and i.name[0] != "."
    ]

    return otu, isolate_ids


def parse_isolate(isolate_path: Path) -> Tuple[dict, list]:
    """
    Creates a list of sequence IDs for a given sequence

    :param isolate_path: Path to a isolate directory
    :return: The dictionary and list of sequenece ids for a given isolate
    """
    with open(isolate_path / "isolate.json", "r") as f:
        isolate = json.load(f)

    isolate["sequences"] = list()

    sequence_ids = [
        i
        for i in isolate_path.glob('*.json')
        if i.name != "isolate.json" and i.name[0] != "."
    ]

    return isolate, sequence_ids
