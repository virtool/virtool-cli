import json
from pathlib import Path
import aiofiles


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
