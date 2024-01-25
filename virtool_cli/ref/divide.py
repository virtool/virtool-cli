import asyncio
import json
import shutil
from pathlib import Path

import aiofiles
import structlog

from virtool_cli.utils.logging import configure_logger
from virtool_cli.utils.reference import generate_otu_dirname

base_logger = structlog.get_logger()

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
    """Divide a reference.json file from Virtool into a src tree structure.

    :param reference_path: Path to a reference.json file
    :param output_path: Path to the where the src tree should be generated
    :param debugging: Enables verbose logs for debugging purposes
    """
    configure_logger(debugging)
    logger = base_logger.bind(reference=str(reference_path), output=str(output_path))

    logger.info(f"Dividing {output_path.name} into {reference_path.name}...")

    asyncio.run(divide_reference_file(reference_path, output_path))


async def divide_reference_file(reference_path: Path, output_path: Path):
    logger = base_logger.bind(reference=str(reference_path), output=str(output_path))

    shutil.rmtree(output_path, ignore_errors=True)
    output_path.mkdir()

    async with aiofiles.open(reference_path, "r") as export_handle:
        contents = await export_handle.read()
        data = json.loads(contents)

    for otu in data.get("otus"):
        otu_path = await build_otu(output_path, otu)

        logger.debug(
            f"Built {otu['_id']}",
            otu={"_id": otu["_id"], "path": str(otu_path.relative_to(output_path))},
        )

        isolates = otu.pop("isolates")

        for isolate in isolates:
            isolate_path = await build_isolate(otu_path, isolate)

            sequences = isolate.pop("sequences")

            for sequence in sequences:
                await build_sequence(isolate_path, sequence)

        # Create empty exclusion list
        with open(otu_path / "exclusions.json", "w") as f:
            json.dump([], f, indent=4, sort_keys=True)

    meta = {"data_type": data["data_type"], "organism": data["organism"]}
    async with aiofiles.open(output_path / "meta.json", "w") as f:
        await f.write(json.dumps(meta, sort_keys=True))

    logger.info("Wrote input reference to output_path")


async def build_otu(src_path: Path, otu: dict) -> Path:
    """Creates a directory for all OTUs that begin with a particular
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

    otu = {key: otu.get(key) for key in OTU_KEYS}
    with open(otu_path / "otu.json", "w") as f:
        json.dump(otu, f, indent=4, sort_keys=True)

    return otu_path


async def build_isolate(otu_path: Path, isolate: dict) -> Path:
    """Creates a directory for a given isolate and generates
    a isolate.json with key information about it.

    :param otu_path: A path to an OTU directory under a src reference directory
    :param isolate: A dictionary containing isolate information
    :return: A path to a newly generated isolate directory
    """
    isolate_path = otu_path / isolate.get("id")
    isolate_path.mkdir()

    isolate = {key: isolate[key] for key in ISOLATE_KEYS}
    with open(isolate_path / "isolate.json", "w") as f:
        json.dump(isolate, f, indent=4, sort_keys=True)

    return isolate_path


async def build_sequence(isolate_path: Path, sequence: dict):
    """Generates a JSON file for a sequence under an isolate directory

    :param isolate_path: A path to a specified isolate directory under an OTU directory
    :param sequence: A dictionary containing sequence information
    """
    sequence = {key: sequence[key] for key in SEQUENCE_KEYS if key in sequence}
    sequence_id = sequence.get("_id")

    with open(isolate_path / f"{sequence_id}.json", "w") as f:
        json.dump(sequence, f, indent=4, sort_keys=True)
