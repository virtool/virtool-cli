import asyncio
import json
from pathlib import Path
import arrow
import aiofiles
import structlog

from virtool_cli.utils.logging import configure_logger
from virtool_cli.utils.reference import (
    get_otu_paths,
    get_isolate_paths,
    get_sequence_paths,
    is_v1,
)

base_logger = structlog.get_logger()

OTU_KEYS = ["_id", "name", "abbreviation", "schema", "taxid"]

ISOLATE_KEYS = ["id", "source_type", "source_name", "default"]

SEQUENCE_KEYS = ["_id", "accession", "definition", "host", "sequence"]


def run(
    src_path: Path,
    output_path: Path,
    indent: bool,
    version: str,
    debugging: bool = False,
):
    """
    Build a Virtool reference JSON file from a data directory.

    :param src_path: Path to database src directory
    :param output_path: The output path for the reference.json file
    :param indent: A flag to indicate whether the output file should be indented
    :param version: The version string to include in the reference.json file
    :param debugging: Enables verbose logs for debugging purposes
    """
    configure_logger(debugging)

    logger = base_logger.bind(src=str(src_path))

    if is_v1(src_path):
        logger.error(
            "This reference database is a deprecated v1 reference."
            + 'Run "virtool ref migrate" before trying again.'
        )
        return

    try:
        asyncio.run(build_from_src(src_path, output_path, indent, version))
    except FileNotFoundError as e:
        logger.exception(e)
    except Exception as e:
        logger.exception(e)


async def build_from_src(src_path: Path, output_path: Path, indent: bool, version: str):
    """
    Build a Virtool reference JSON file from a data directory.

    :param src_path: Path to database src directory
    :param output_path: The output path for the reference.json file
    :param indent: A flag to indicate whether the output file should be indented
    :param version: The version string to include in the reference.json file
    """
    logger = base_logger.bind(src=str(src_path), output=str(output_path))

    meta = await parse_meta(src_path)
    if meta:
        logger.debug(
            "Metadata parsed", meta=meta, metadata_path=str(src_path / "meta.json")
        )
        data = {
            "data_type": meta.get("data_type", "genome"),
            "organism": meta.get("organism", ""),
        }

    else:
        logger.warning(
            f'Metadata file not found at {src_path / "meta.json"}',
            metadata_path=str(src_path / "meta.json"),
        )

        data = {"data_type": "genome", "organism": ""}

    otus = []

    for otu_path in get_otu_paths(src_path):
        try:
            otu = await parse_otu_contents(otu_path)
        except (FileNotFoundError, json.JSONDecodeError) as e:
            logger.error("Reference data at src_path is invalid.")
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

    try:
        async with aiofiles.open(output_path, "w") as f:
            await f.write(
                json.dumps(data, indent=4 if indent else None, sort_keys=True)
            )

        logger.info("Reference file built at output")

        print(str(output_path.resolve()))

    except Exception as e:
        logger.exception(e)


async def parse_meta(src_path: Path) -> dict:
    """
    Deserializes and returns meta.json if found, else returns an empty dictionary.

    :param src_path: Path to database src directory
    :return: The deserialized meta.json object or an empty dictionary
    """
    try:
        async with aiofiles.open(src_path / "meta.json", "r") as f:
            contents = await f.read()
            return json.loads(contents)

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


async def parse_otu_contents(otu_path: Path) -> dict:
    """
    Traverses, deserializes and returns all data under an OTU directory.

    :param otu_path: Path to an OTU directory
    :return: All isolate and sequence data under an OTU,
        deserialized and compiled in a dict
    """
    async with aiofiles.open(otu_path / "otu.json", "r") as f:
        contents = await f.read()
        otu = json.loads(contents)

    logger = base_logger.bind(path=otu_path, otu_id=otu["_id"])

    isolates = []

    for isolate_path in get_isolate_paths(otu_path):
        async with aiofiles.open(isolate_path / "isolate.json", "r") as f:
            contents = await f.read()
            isolate = json.loads(contents)

        sequences = []
        for sequence_path in get_sequence_paths(isolate_path):
            sequence = await parse_sequence(sequence_path)

            sequences.append(sequence)
            logger.debug(
                f"Added sequence {sequence.get('accession')} under id={sequence.get('_id')}",
                sequence_id=sequence,
                isolate_id=isolate["id"],
            )

        isolate["sequences"] = sequences

        isolates.append(isolate)

    otu["isolates"] = isolates

    return otu


async def parse_sequence(path):
    """
    Asynchronously reads a sequence file and returns the contents as a dictionary.

    :param path:
    :return: A deserialized sequence
    """
    async with aiofiles.open(path, "r") as f:
        contents = await f.read()
        sequence = json.loads(contents)

    return sequence
