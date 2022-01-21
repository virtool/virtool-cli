import asyncio
import csv
import json
import shutil
from collections import defaultdict
from dataclasses import dataclass, asdict
from hashlib import md5
from json import JSONDecodeError
from pathlib import Path
from typing import List

import aiofiles
from rich.console import Console
from rich.prompt import Prompt

from virtool_cli.utils import get_otu_paths


@dataclass
class OTUInfo:
    ref: str
    name: str
    directory: str


@dataclass
class IsolateInfo:
    ref: str
    otu_name: str
    id: str
    isolate: dict
    sequences: list
    taxid: str


@dataclass
class UnknownIsolate:
    otu_name: str
    isolate_id: str
    sequence_accessions: list


def run(
    source_src_path: str,
    target_src_path: str,
    resume: bool,
    in_place: bool,
    output: str,
    output_unknown_isolates: bool,
    unknown_output: str,
):
    asyncio.run(
        merge_refs(
            source_src_path,
            target_src_path,
            resume,
            in_place,
            output,
            output_unknown_isolates,
            unknown_output,
        )
    )


async def merge_refs(
    source_src_path: str,
    target_src_path: str,
    resume: bool,
    in_place: bool,
    output: str,
    output_unknown_isolates: bool,
    unknown_output: str,
):
    """
    Runs routines to merge two references.

    :param source_src_path: path to src for source reference
    :param target_src_path: path to src for target reference
    :param resume: flag to determine whether in-progress merge should be continued using cache
    :param in_place: flag to determine whether output should be written directly to directory of target reference
    :param output_unknown_isolates: flag to determine whether to write CSV file for isolates with unknown source name
    or source type
    :param unknown_output: name of output CSV file for isolates with unknown source name or source type
    :param output: name of output directory
    """
    source_src_path = Path(source_src_path)
    target_src_path = Path(target_src_path)

    otu_infos_by_taxid, no_taxid_found = await organize_otus_by_taxid(
        source_src_path, target_src_path
    )

    if in_place:
        output = target_src_path.parent.name

    if resume:
        existing_cache = get_cache()
        for taxid in existing_cache.keys():
            otu_infos_by_taxid.pop(taxid, None)

        if output_unknown_isolates:
            try:
                shutil.copyfile(
                    Path("./.cli/unknown_isolates.txt"),
                    source_src_path.parent.parent / unknown_output,
                )
            except FileNotFoundError:
                (Path.cwd() / ".cli").mkdir()
                Path("./.cli/unknown_isolates.txt").touch()

            with open(Path("./.cli/unknown_isolates.txt"), "r+") as f:
                f.truncate()
    else:
        # empty cache
        write_cache({})

        # empty unknown isolates cache
        try:
            with open(source_src_path.parent.parent / unknown_output, "r+") as f:
                f.truncate()

        except (NotADirectoryError, FileNotFoundError):
            if output_unknown_isolates:
                Path(source_src_path.parent.parent / unknown_output).touch()

        try:
            with open(Path("./.cli/unknown_isolates.txt"), "r+") as f:
                f.truncate()

        except (NotADirectoryError, FileNotFoundError):
            if output_unknown_isolates:
                (Path.cwd() / ".cli").mkdir(exist_ok=True)
                Path("./.cli/unknown_isolates.txt").touch()

    for taxid, otu_infos in otu_infos_by_taxid.items():
        (
            isolate_infos_by_source_name,
            unknown_isolates,
        ) = await organize_isolates_by_source_name(
            otu_infos, source_src_path, target_src_path, taxid, output_unknown_isolates
        )

        if output_unknown_isolates:
            write_unknown_isolates(unknown_isolates, unknown_output, source_src_path)

        isolate_infos_to_add = await compare_isolates(isolate_infos_by_source_name)

        otu_results_path = (
            source_src_path.parent.parent
            / output
            / "src"
            / otu_infos[0].directory[0]
            / otu_infos[0].directory
        )

        result = {
            taxid: {
                "otu_path": str(otu_results_path),
                "isolate_infos": isolate_infos_to_add,
            }
        }

        update_cache(result)

    verified_isolates = get_cache()
    await write_isolates(verified_isolates)


async def organize_otus_by_taxid(
    source_src_path: Path, target_src_path: Path
) -> [defaultdict, List]:
    """
    Organize OTUInfo objects by taxid.

    Iterate through all OTUs in source and target src directories and create OTUInfo objects for each OTU with a taxid.

    If any OTUs are found without corresponding taxid, add them to no_taxid_found.

    :param source_src_path: path to src for source ref
    :param target_src_path: path to src for target ref

    :return: OTUInfo organized by taxid, list of OTUs for which no taxid was found
    """
    otu_paths = get_otu_paths(target_src_path) + get_otu_paths(source_src_path)

    otus_by_taxid = defaultdict(list)
    no_taxid_found = list()

    for otu_path in otu_paths:
        with open(otu_path / "otu.json", "r") as f:
            otu = json.loads(f.read())
            ref = "source" if str(source_src_path) in str(otu_path) else "target"

            if otu.get("taxid", None):
                otus_by_taxid[str(otu["taxid"])].append(
                    OTUInfo(ref, otu["name"], directory=otu_path.name)
                )
            else:
                no_taxid_found.append(otu)

    return otus_by_taxid, no_taxid_found


async def organize_isolates_by_source_name(
    otu_infos: list,
    source_src_path: Path,
    target_src_path: Path,
    taxid: str,
    output_unknown_isolates: bool,
) -> [defaultdict, list]:
    """
    Iterate through OTUInfo objects with the same taxid and generate dictionary containing IsolateInfo objects
    organized by the source name of each isolate.

    :param otu_infos: list of OTUInfo objects for each taxid
    :param source_src_path: path to src for source reference
    :param target_src_path: path to src for target reference
    :param taxid: taxid for OTU to which isolates belong
    :param output_unknown_isolates: flag to determine whether to write CSV file for isolates with unknown source name
    or source type

    :return: Dictionary of isolate paths organized by source name
    """
    isolate_infos_by_source_name = defaultdict(list)
    unknown_isolates = list()

    for otu_info in otu_infos:
        otu_path = generate_otu_path(otu_info, source_src_path, target_src_path)

        for folder in otu_path.iterdir():
            # ignore the otu.json file in the OTU folder and parse through isolate folders

            if (
                folder.name != "otu.json"
                and folder / "isolate.json" in folder.iterdir()
            ):
                sequences = await get_sequences(folder)

                async with aiofiles.open(folder / "isolate.json", "r") as f:
                    isolate = json.loads(await f.read())

                    if (
                        isolate["source_name"].lower() == "unknown"
                        or isolate["source_type"].lower == "unknown"
                    ):
                        if output_unknown_isolates:
                            accessions = await get_accessions_for_unknown_isolate(
                                folder
                            )
                            unknown_isolates.append(
                                UnknownIsolate(otu_info.name, isolate["id"], accessions)
                            )

                    isolate_infos_by_source_name[isolate["source_name"].lower()].append(
                        IsolateInfo(
                            ref=otu_info.ref,
                            otu_name=otu_info.name,
                            id=isolate["id"],
                            isolate=isolate,
                            sequences=sequences,
                            taxid=taxid,
                        )
                    )

    return isolate_infos_by_source_name, unknown_isolates


async def get_accessions_for_unknown_isolate(isolate_path: Path) -> list:
    """
    Gather accessions for all sequences of isolate with unknown source type or source name.

    :param isolate_path: path to directory containing isolate.json file and sequences
    """
    accessions = list()

    for path in isolate_path.iterdir():
        if path.name != "isolate.json":
            async with aiofiles.open(path, "r") as f:
                sequence = json.loads(await f.read())
                accessions.append(sequence["accession"])

    return sorted(accessions)


def write_unknown_isolates(
    unknown_isolates: list, unknown_output: str, source_src_path: Path
):
    """
    Write CSV file containing information about isolates with unknown source type or source name

    :param unknown_isolates: list of UnknownIsolate objects
    :param unknown_output: name for CSV file
    :param source_src_path: path to source ref directory
    """
    output_path = source_src_path.parent.parent / unknown_output
    output_path.touch(exist_ok=True)

    with open(output_path, "a", newline="") as handle:
        writer = csv.writer(
            handle, delimiter="\t", escapechar=" ", quoting=csv.QUOTE_NONE
        )

        for unknown_isolate in unknown_isolates:
            writer.writerow(
                [
                    unknown_isolate.otu_name,
                    unknown_isolate.isolate_id,
                    sorted(unknown_isolate.sequence_accessions),
                ]
            )

    unknown_isolates_cache = Path("./.cli/unknown_isolates.txt")

    shutil.copyfile(output_path, unknown_isolates_cache)


def generate_otu_path(
    otu_info: OTUInfo, source_src_path: Path, target_src_path: Path
) -> Path:
    """
    Generate path to OTU directory from OTUInfo object

    :param otu_info: OTUInfo object
    :param source_src_path: path to src for source reference
    :param target_src_path: path to src for target reference

    :return: path to OTU directory
    """
    if otu_info.ref == "source":
        otu_path = source_src_path / otu_info.directory[0] / otu_info.directory
    else:
        otu_path = target_src_path / otu_info.directory[0] / otu_info.directory

    return otu_path


async def get_sequences(isolate_path: Path) -> list:
    """
    Returns a list of sequences for a given OTU

    :param isolate_path: A path to an otu directory in a reference
    :return: A list of all sequence dictionaries found for given OTU
    """
    sequences = list()
    file_paths = list(isolate_path.iterdir())

    for file in sorted(file_paths):
        if file.name != "isolate.json":
            async with aiofiles.open(file, "r") as f:
                sequence = json.loads(await f.read())
                sequences.append(sequence)

    return sequences


async def compare_isolates(isolates_to_compare: dict) -> List:
    """
    Prompt user input to determine which isolates should be included in output.

    If only one IsolateInfo object is stored under a source name  in isolates_to_compare and the isolate is from the
    source ref, the user is prompted to confirm that it is a new isolate. If the isolate is from the target ref, it is
    automatically added to verified_isolates without prompting for user input.

    If two IsolateInfo objects are stored under a source name in isolates_to_compare, then the user is prompted to
    confirm that they are the same. If the isolates are the same, the target isolate is added to verified_isolates.

    If the isolates aren't the same, the user is allowed to select one, both or neither of the isolates to add to
    verified_isolates.

    :param isolates_to_compare: dictionary of IsolateInfo objects organized by source_name

    :return: verified isolates to add for each OTU
    """
    verified_isolates = list()

    add_source_isolate = False
    add_target_isolate = False

    console = Console()

    for _, isolate_infos in isolates_to_compare.items():
        prompt = prepare_prompt(isolate_infos)

        if len(isolate_infos) == 1:
            if isolate_infos[0].ref == "source":
                console.print(prompt, style="cyan")
                new_isolate = Prompt.ask("", choices=["Y", "N"], default="Y")

                if new_isolate == "Y":
                    add_source_isolate = True
            else:
                add_target_isolate = True

        else:
            console.print(prompt, style="cyan")
            same = Prompt.ask("", choices=["Y", "N"], default="Y")

            if same == "Y":
                add_target_isolate = True
            else:
                prompt = """Which isolate should be added to result?\n
        1) Ⓢ Source
        2) Ⓣ Target
        3) Both
        4) Neither\n
                """
                console.print(prompt, style="cyan")
                to_add = Prompt.ask("", choices=["1", "2", "3", "4"])
                if to_add == "1":
                    add_source_isolate = True
                elif to_add == "2":
                    add_target_isolate = True
                elif to_add == "3":
                    add_source_isolate = True
                    add_target_isolate = True
                else:
                    return verified_isolates

        for isolate_info in isolate_infos:
            # convert to dict for caching
            if isolate_info.ref == "source":
                source_isolate_info = asdict(isolate_info)
            else:
                target_isolate_info = asdict(isolate_info)

        if add_target_isolate:
            verified_isolates.append(target_isolate_info)
        if add_source_isolate:
            verified_isolates.append(source_isolate_info)

    return verified_isolates


def prepare_prompt(isolate_infos: list) -> str:
    """
    Prepare prompt for user to compare isolates.

    :param isolate_infos: IsolateInfo objects for isolates to be compared

    :return: prompt for user to be asked
    """
    if len(isolate_infos) == 1:
        prompt = "Is this a new isolate?\n"
    else:
        prompt = "Are these isolates the same?\n"

    for isolate_info in isolate_infos:
        if isolate_info.ref == "source":
            symbol = "Ⓢ"
        else:
            symbol = "Ⓣ"

        source_type_caps = isolate_info.isolate["source_type"].capitalize()
        source_name = isolate_info.isolate["source_name"]

        isolate_name = f"{source_type_caps} {source_name}"

        prompt += f"\n\t• {symbol}  {isolate_info.ref.capitalize()} - {isolate_name} ({isolate_info.taxid})\n"

        for sequence in isolate_info.sequences:
            accession = sequence["accession"]
            definition = sequence["definition"]
            length = len(sequence["sequence"])
            sequence_hash = md5(sequence["sequence"].encode())

            prompt += f"\t\t• {accession} | {definition} | Length={length} | Hash={sequence_hash.hexdigest()[0:8]}\n"

    return prompt


async def write_isolates(verified_isolates: dict):
    """
    Write isolate.json files for each isolate in to their respective folders

    :param verified_isolates: dictionary containing information about isolates
    """
    for _, result in verified_isolates.items():

        for isolate_info in result["isolate_infos"]:
            # revert to IsolateInfo object from dict in cache
            isolate_info = IsolateInfo(**isolate_info)

            sequences = isolate_info.sequences
            isolate = isolate_info.isolate

            isolate_path = Path(result["otu_path"]) / isolate["id"]

            isolate_path.mkdir(parents=True, exist_ok=True)

            async with aiofiles.open(isolate_path / "isolate.json", "w") as f:
                await f.write(json.dumps(isolate, indent=4))

            await write_sequences(isolate_path, sequences)


async def write_sequences(isolate_path: Path, sequences: list):
    """
    Write sequences to their respective sequence_id.json files

    :param isolate_path: path to directory with isolate.json file
    :param sequences: sequences to be written to output reference
    """
    for sequence in sequences:
        sequence_id = sequence["_id"]
        async with aiofiles.open(isolate_path / f"{sequence_id}.json", "w") as f:
            await f.write(json.dumps(sequence, indent=4))


def write_cache(isolates: dict):
    """
    Cache a mapping of taxon ids to all isolates found.

    :param isolates: Dictionary containing an updated mapping of taxids to their accessions
    """
    (Path.cwd() / ".cli").mkdir(exist_ok=True)
    (Path.cwd() / ".cli" / "isolates.json").touch(exist_ok=True)

    with open(".cli/isolates.json", "w") as f:
        json.dump(isolates, f, indent=4)


def get_cache() -> dict:
    """
    Fetches the local cache containing taxon ids mapped to their isolates

    :return: A dictionary that maps taxon ids to found accessions
    """
    try:
        with open(".cli/isolates.json", "r") as f:
            return json.load(f)
    except (NotADirectoryError, FileNotFoundError, JSONDecodeError):
        return {}


def update_cache(result: dict) -> dict:
    """
    Updates the existing isolate cache with fetch results

    :param result: List of tuples pulled out of the asyncio Queue containing each taxid and their accessions
    """
    existing_cache = get_cache()

    updated_cache = {**existing_cache, **result}

    write_cache(updated_cache)

    return updated_cache
