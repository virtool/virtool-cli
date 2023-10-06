import json
from pathlib import Path
import aiofiles
import structlog

from virtool_cli.utils.reference import get_isolate_paths
from virtool_cli.utils.id_generator import generate_unique_ids
from virtool_cli.utils.format import format_isolate


async def write_records(
    otu_path: Path,
    new_sequences: list,
    unique_iso: set,
    unique_seq: set,
    logger: structlog.BoundLogger = structlog.get_logger(),
):
    """
    :param otu_path: A path to an OTU directory under a src reference directory
    :param new_sequences: List of new sequences under the OTU
    :param unique_iso: Set of all unique isolate IDs present in the reference
    :param unique_seq: Set of all unique sequence IDs present in the reference
    :param logger: Optional entry point for a shared BoundLogger
    """
    ref_isolates = await label_isolates(otu_path)

    try:
        seq_hashes = generate_unique_ids(
            n=len(new_sequences), excluded=list(unique_seq)
        )
    except Exception as e:
        logger.exception(e)
        return e

    logger.debug(f"Writing {len(new_sequences)} sequences...")

    for seq_data in new_sequences:
        isolate_data = seq_data.pop("isolate")
        isolate_name = isolate_data["source_name"]
        isolate_type = isolate_data["source_type"]

        if isolate_name in ref_isolates:
            iso_id = ref_isolates[isolate_name]["id"]
            logger.debug(
                "Existing isolate name found", iso_name=isolate_name, iso_hash=iso_id
            )

        else:
            try:
                iso_id = generate_unique_ids(n=1, excluded=list(unique_iso)).pop()
            except Exception as e:
                logger.exception(e)
                continue

            logger.debug(
                "Assigning new isolate hash", iso_name=isolate_name, iso_hash=iso_id
            )

            try:
                new_isolate = format_isolate(isolate_name, isolate_type, iso_id)
            except Exception as e:
                logger.exception(e)
                continue

            await store_isolate(new_isolate, iso_id, otu_path)

            unique_iso.add(iso_id)
            ref_isolates[isolate_name] = new_isolate

            logger.info("Created a new isolate directory", path=str(otu_path / iso_id))

        iso_path = otu_path / iso_id

        seq_hash = seq_hashes.pop()
        logger.debug(
            "Assigning new sequence",
            seq_hash=seq_hash,
        )
        try:
            await store_sequence(seq_data, seq_hash, iso_path)
        except Exception as e:
            logger.exception(e)
        unique_seq.add(seq_hash)

        logger.info(
            f"Wrote new sequence '{seq_hash}'", path=str(iso_path / f"{seq_hash}.json")
        )


async def store_isolate(
    isolate: dict, isolate_id: str, otu_path: Path
):
    """
    Creates a new isolate directory and metadata file under an OTU directory,
    then returns the metadata in dict form

    :param isolate: Dictionary containing isolate metadata
    :param isolate_id: Unique ID number for this new isolate
    :param otu_path: Path to the parent OTU
    :return: The unique isolate id
    """
    iso_path = otu_path / isolate_id
    iso_path.mkdir()

    async with aiofiles.open(iso_path / "isolate.json", "w") as f:
        await f.write(json.dumps(isolate, indent=4))


async def store_sequence(sequence: dict, sequence_id: str, iso_path: Path):
    """
    Write sequence to isolate directory within the src directory

    :param sequence: Dictionary containing formatted sequence data
    :param sequence_id: Unique ID number for this new sequence
    :param iso_path: Path to the parent isolate
    :return: The unique sequence id (aka seq_hash)
    """
    sequence["_id"] = sequence_id
    seq_path = iso_path / f"{sequence_id}.json"

    async with aiofiles.open(seq_path, "w") as f:
        await f.write(json.dumps(sequence, indent=4, sort_keys=True))


async def label_isolates(otu_path: Path) -> dict:
    """
    Return all isolates present in an OTU directory

    :param otu_path: Path to an OTU directory
    :return: A dictionary of isolates indexed by source_name
    """
    if not otu_path.exists():
        raise FileNotFoundError

    isolates = {}

    for iso_path in get_isolate_paths(otu_path):
        async with aiofiles.open(iso_path / "isolate.json", "r") as f:
            contents = await f.read()
            isolate = json.loads(contents)
        isolates[isolate.get("source_name")] = isolate

    return isolates


async def read_otu(path: Path) -> dict:
    """
    Returns a json file in dict form

    :param path: Path to an OTU directory under a reference source
    :return: Deserialized OTU data in dict form
    """
    async with aiofiles.open(path / "otu.json", "r") as f:
        contents = await f.read()
        otu = json.loads(contents)

    return otu
