import json
from pathlib import Path
import aiofiles


async def store_isolate(
    source_name: str, source_type: str, iso_hash: str, otu_path: Path
) -> str:
    """
    Creates a new isolate directory and metadata file under an OTU directory,
    then returns the metadata in dict form

    :param source_name: Assigned source name for an accession
    :param source_type: Assigned source type for an accession
    :param iso_hash: Unique ID number for this new isolate
    :param otu_path: Path to the parent OTU
    :return: The unique isolate id
    """
    iso_path = otu_path / iso_hash
    iso_path.mkdir()

    new_isolate = {
        "id": iso_hash,
        "source_type": source_type,
        "source_name": source_name,
        "default": False,
    }

    async with aiofiles.open(iso_path / "isolate.json", "w") as f:
        await f.write(json.dumps(new_isolate, indent=4))

    return new_isolate


async def store_sequence(sequence: dict, seq_hash: str, iso_path: Path) -> str:
    """
    Write sequence to isolate directory within the src directory

    :param sequence: Dictionary containing formatted sequence data
    :param seq_hash: Unique ID number for this new sequence
    :param iso_path: Path to the parent isolate
    :return: The unique sequence id (aka seq_hash)
    """
    sequence["_id"] = seq_hash
    seq_path = iso_path / f"{seq_hash}.json"

    async with aiofiles.open(seq_path, "w") as f:
        await f.write(json.dumps(sequence, indent=4, sort_keys=True))

    return seq_hash
