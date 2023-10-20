import asyncio
import json
from pathlib import Path
from concurrent.futures.thread import ThreadPoolExecutor
from random import choice
from string import ascii_letters, ascii_lowercase, digits
from typing import Iterable, Tuple, Optional
from urllib.error import HTTPError

import aiofiles
import aiojobs
from Bio import Entrez, SeqIO
import structlog

from virtool_cli.utils.reference import get_otu_paths, get_unique_ids
from virtool_cli.utils.ncbi import NCBI_REQUEST_INTERVAL

logger = structlog.get_logger()


async def isolate(src_path: Path):
    """
    Runs routines to find new isolates for OTU in a reference directory and writes newfound accessions to local cache

    :param src_path: Path to a given reference directory
    """
    scheduler = aiojobs.Scheduler()
    asyncio.get_event_loop().set_default_executor(ThreadPoolExecutor())

    coros = []
    q = asyncio.Queue()

    existing_accessions = get_cache()

    # get paths for all OTU in the directory
    paths = get_otu_paths(src_path)

    # get mapping of all OTU paths to its OTU JSON file
    otu_path_map = map_otus(paths)

    for path in paths:
        taxid = otu_path_map[path].get("taxid", None)
        name = otu_path_map[path].get("name", "")

        if taxid is None:
            logger.info("No taxid assigned to OTU", taxid=taxid, name=name)
            continue

        accessions = existing_accessions.get(str(taxid))

        coros.append(fetch_otu_isolates(str(taxid), name, path, accessions, paths, q))

    while len(coros) != 0:
        try:
            if scheduler.active_count < scheduler.limit:
                await scheduler.spawn(coros.pop())
                await asyncio.sleep(NCBI_REQUEST_INTERVAL)
        except KeyboardInterrupt:
            await scheduler.close()
            break

    results_to_cache = list()
    while True:
        try:
            results_to_cache.append(q.get_nowait())
        except asyncio.QueueEmpty:
            if not scheduler.active_count:
                break

        await asyncio.sleep(0.01)

    await scheduler.close()

    updated_cache = update_cache(existing_accessions, results_to_cache)

    write_cache(updated_cache)


async def fetch_otu_isolates(
    taxid: str,
    name: str,
    path: Path,
    accessions: list,
    paths: list,
    queue: asyncio.Queue,
):
    """
    Fetch new isolates for a given OTU, creating a directory for each

    :param taxid: The taxon id for a given OTU
    :param name: Name of a given OTU
    :param path: Path to a given OTU in a reference directory
    :param accessions: List of accessions for the given taxon id from local cache
    :param paths: List of all paths to every OTU in a reference
    :param queue: Asyncio queue to put newfound accessions onto for caching purposes
    """

    otu_log = logger.bind(name=name, taxid=taxid)

    isolates = await get_isolates(path)

    records, new_accessions = await asyncio.get_event_loop().run_in_executor(
        None, get_records, accessions, taxid
    )

    isolate_ids, sequence_ids = await get_unique_ids(paths)

    if records is None:
        otu_log.info("Found 0 isolates, could not link taxid to nucleotide records")
        return

    new_isolates = dict()

    for accession in [accession for accession in records.values() if accession.seq]:
        accession_data = await get_qualifiers(accession.features)

        # try to find isolate type and name automatically
        isolate_type = await find_isolate(accession_data)

        if isolate_type is None:
            continue

        isolate_name = accession_data.get(isolate_type)[0]

        # see if isolate directory already exists and create it if needed
        if isolate_name not in isolates:
            new_isolates[isolate_name] = isolate_type

            new_id = await store_isolate(path, isolate_name, isolate_type, isolate_ids)
            isolate_ids.add(new_id)
            isolates[isolate_name] = new_id

        # create new sequence file

        isolate_path = path / isolates.get(isolate_name)

        new_id = await store_sequence(
            isolate_path, accession, accession_data, sequence_ids
        )
        sequence_ids.add(new_id)

        await queue.put((taxid, new_accessions))

    await log_results(name, taxid, new_isolates, otu_log)


def get_records(accessions: list, taxid: str) -> Tuple[Optional[dict], list]:
    """
    Uses an OTU's taxid to find accessions numbers if needed, then fetches Genbank records for each accession number

    :param accessions: A list of accession ids to retrieve or None to retrieve accession numbers for the given taxid
    :param taxid: The taxid to find accessions for
    :return: A dictionary with all the accession records and a list containing all accession ids for caching purposes
    """

    if accessions is None:
        accessions = get_accession_numbers(taxid)

    return fetch_records(accessions), accessions


def get_accession_numbers(taxid: str) -> list:
    """
    Links an OTU's taxid to all associated accession numbers in the Nucleotide database

    :param taxid: The taxid to find accessions for
    :return: List containing accession numbers to fetch records for
    """
    accessions = []

    record = Entrez.read(
        Entrez.elink(dbfrom="taxonomy", db="nucleotide", id=taxid, idtype="acc")
    )

    for linksetdb in record[0]["LinkSetDb"][0]["Link"]:
        accessions.append(linksetdb["Id"])

    return accessions


def fetch_records(accessions: list) -> Optional[dict]:
    """
    Fetches a Genbank record for each accession number in a given list, then parses all of the records into a
    SeqIO dictionary that can be used to access specific fields in each record

    :param accessions: List of accession numbers to be fetched
    :return: A SeqIO dictionary object
    """
    try:
        handle = Entrez.efetch(
            db="nucleotide", id=accessions, rettype="gb", retmode="text"
        )
    except HTTPError:
        return None

    return SeqIO.to_dict(SeqIO.parse(handle, "gb"))


async def log_results(
    name: str, taxid: str, new_isolates: dict, console: structlog.BoundLogger
):
    """
    Log isolate discovery results to console

    :param name: Name of the OTU to log
    :param taxid: Taxid of the OTU to log
    :param new_isolates: Dictionary containing newly found isolates mapped to their type
    :param console: Rich console object used to log results
    """
    new_isolate_count = len(new_isolates)

    console.info(
        f"Found {new_isolate_count} new isolates",
        name=name,
        taxon_id=taxid,
        n_isolates=new_isolate_count,
    )


async def get_qualifiers(seq: list) -> dict:
    """
    Get relevant qualifiers in a Genbank record

    :param seq: SeqIO features object for a particular accession
    :return: Dictionary containing all qualifiers in the source field of the features section of a Genbank record
    """
    features = [feature for feature in seq if feature.type == "source"]
    isolate_data = {}

    for feature in features:
        for qualifier in feature.qualifiers:
            isolate_data[qualifier] = feature.qualifiers.get(qualifier)
    return isolate_data


async def find_isolate(isolate_data: dict) -> Optional[str]:
    """
    Determine the source type in a Genbank record

    :param isolate_data: Dictionary containing qualifiers in a features section of a Genbank record
    :return:
    """
    for qualifier in ["isolate", "strain"]:
        if qualifier in isolate_data:
            return qualifier

    return None


async def store_sequence(
    path: Path, accession: SeqIO.SeqRecord, data: dict, unique_ids: set
) -> Optional[str]:
    """
    Creates a new sequence file for a given isolate

    :param path: Path to a isolate folder
    :param accession: Genbank record object for a given accession
    :param data: Dictionary containing all qualifiers in the source field of the features section of a Genbank record
    :param unique_ids: Set containing unique Virtool ids for all sequences in a reference
    :return: The newly generated unique id
    """
    sequences = await get_sequences(path)

    # strip accession version
    accession_id = accession.id.split(".")[0]

    # check if accession doesn't already exist
    if accession_id in sequences:
        return

    # generate new sequence id
    new_id = random_alphanumeric(8, False, unique_ids)

    seq_file = {
        "_id": new_id,
        "accession": accession_id,
        "definition": accession.description,
        "host": data.get("host")[0] if data.get("host") is not None else None,
        "sequence": str(accession.seq),
    }

    async with aiofiles.open(path / (new_id + ".json"), "w") as f:
        await f.write(json.dumps(seq_file, indent=4))

    return new_id


async def store_isolate(
    path: Path, source_name: str, source_type: str, unique_ids: set
) -> str:
    """
    Creates a new isolate folder for an OTU

    :param path: Path to an OTU
    :param source_name: Assigned source name for an accession
    :param source_type: Assigned source type for an accession
    :param unique_ids: Set containing unique Virtool ids for all isolates in a reference
    :return: A newly generated unique id
    """

    new_id = random_alphanumeric(8, False, unique_ids)
    (path / new_id).mkdir()

    new_isolate = {
        "id": new_id,
        "source_type": source_type,
        "source_name": source_name,
        "default": False,
    }

    async with aiofiles.open(path / new_id / "isolate.json", "w") as f:
        await f.write(json.dumps(new_isolate, indent=4))

    return new_id


def get_cache() -> dict:
    """
    Fetches the local cache containing taxon ids mapped to their accessions

    :return: A dictionary that maps taxon ids to found accessions
    """
    try:
        with open(".cli/accessions.json", "r") as f:
            return json.load(f)
    except (NotADirectoryError, FileNotFoundError):
        return {}


def update_cache(existing_cache: dict, results: list) -> dict:
    """
    Updates the existing accessions cache with fetch results

    :param existing_cache: The existing cache from the .cli folder
    :param results: List of tuples pulled out of the asyncio Queue containing each taxid and their accessions
    """
    cache_update = dict(results)

    return {**existing_cache, **cache_update}


def write_cache(accessions: dict):
    """
    Cache a mapping of taxon ids to all accessions found.

    :param accessions: Dictionary containing an updated mapping of taxids to their accessions
    """
    try:
        (Path.cwd() / ".cli").mkdir()
    except FileExistsError:
        pass

    with open(".cli/accessions.json", "w") as f:
        json.dump(accessions, f, indent=4)


def random_alphanumeric(
    length: int = 6, mixed_case: bool = False, excluded: Optional[Iterable[str]] = None
) -> str:
    """
    Generates a random string composed of letters and numbers.

    :param length: the length of the string.
    :param mixed_case: included alpha characters will be mixed case instead of lowercase
    :param excluded: strings that may not be returned.
    :return: a random alphanumeric string.
    """
    excluded = set(excluded or list())

    characters = digits + (ascii_letters if mixed_case else ascii_lowercase)

    candidate = "".join([choice(characters) for _ in range(length)])

    if candidate not in excluded:
        return candidate

    return random_alphanumeric(length=length, excluded=excluded)


async def get_isolates(path: Path) -> dict:
    """
    Returns a mapping to every isolate and their folder name.

    :param path: A path to a OTU directory in a reference
    :return: A mapping of all of an OTU's isolates to their folder name (id)
    """
    isolates = dict()

    for folder in path.iterdir():
        # ignore the otu.json file in the OTU folder and parse through isolate folders
        if folder.is_dir() and folder / "isolate.json" in folder.iterdir():
            async with aiofiles.open(folder / "isolate.json", "r") as f:
                isolate = json.loads(await f.read())

            isolates[isolate["source_name"]] = folder.name

    return isolates


async def get_sequences(path: Path) -> dict:
    """
    Returns a mapping of sequence accessions to their file name in a isolate directory.

    :param path: A path to an isolate directory in a reference
    :return: A mapping of all accessions in an isolate to their file name (id)
    """
    sequences = dict()

    for sequence_id in path.glob("*.json"):
        if sequence_id.name != "isolate.json":
            async with aiofiles.open(sequence_id, "r") as f:
                sequence = json.loads(await f.read())
                sequences[sequence.get("accession")] = sequence_id.name

    return sequences


def map_otus(paths: list) -> dict:
    """
    Returns a mapping of every OTU path to their deserialized OTU dictionary.

    :param paths: List of paths to all OTU in a reference
    :return: A mapping of every OTU path to its OTU dictionary
    """
    path_taxid_map = {}

    for path in paths:
        path_taxid_map[path] = json.loads((path / "otu.json").read_text())

    return path_taxid_map


def run(src_path: Path):
    """
    Runs the asynchronous routines to find new isolates for all OTU in a reference

    :param src_path: Path to a reference directory
    """
    asyncio.run(isolate(src_path))
