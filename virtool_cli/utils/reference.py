import re
from pathlib import Path
from typing import Optional, Tuple


def is_v1(src_path: Path) -> bool:
    """
    Returns True if the given reference directory is formatted under version 1 guidelines.
    Determines if virtool ref migrate should be run before using this version of virtool-cli
    
    :param src_path: Path to a src database directory
    :return: Boolean value depending on whether alphabetized bins are found.
    """
    alpha_bins = list(src_path.glob('[a-z]'))

    return bool(alpha_bins)


def get_otu_paths(src_path: Path) -> list:
    """
    Generates a list of paths to all OTUs in a src directory.

    :param src_path: Path to a src database directory
    :return: List of paths to all OTU in a src directory
    """
    return [otu for otu in src_path.iterdir() if otu.is_dir()]


def get_isolate_paths(otu_path: Path) -> list:
    """
    Generates a list of paths to all isolates in an OTU directory.

    :param otu_path: Path to an OTU directory under a reference directory
    :return: List of paths to all OTU in a src directory
    """
    return [iso_path for iso_path in otu_path.iterdir() if iso_path.is_dir()]


def get_sequence_paths(isolate_path: Path) -> list:
    """
    Generates a list of paths to all sequences in an isolate directory.

    :param isolate_path: Path to an isolate directory under an OTU directory
    :return: A list of paths to all sequence files in an isolate directory
    """
    sequence_ids = [
        i for i in isolate_path.glob('*.json') if i.stem != "isolate"
    ]

    return sequence_ids


def search_otu_by_id(otu_id: str, src_path: Path) -> Optional[Path]:
    """
    Searches filenames in the database by unique id and
    returns the path if it finds a match

    :param src_path: Path to a reference database directory
    :param otu_id: Virtool OTU unique ID
    """
    for path in src_path.glob(f'*--{otu_id}'):
        if path.is_dir():
            return path
    
    return None


def generate_otu_dirname(name: str, otu_id: str = '') -> str:
    """
    Takes in a human-readable string, replaces whitespace and symbols
    and adds the Virtool hash id as a suffix
    
    :param name: Human-readable, searchable name of the OTU
    :param otu_id: ID hash of OTU 
    :return: A directory name in the form of 'converted_otu_name--taxid'
    """
    no_plus = name.replace('+', 'plus ')
    no_symbols = re.split(r'[():/-]+', no_plus)
    joined = ' '.join(no_symbols)
    no_whitespace = re.sub(r'[\s]+', "_", joined)

    dirname = no_whitespace.lower()
    dirname += '--' + otu_id

    return dirname


async def get_unique_ids(otu_paths: list) -> Tuple[set, set]:
    """
    Returns sets containing unique random alphanumeric ids for both the isolates and the sequences

    :param otu_paths: List of paths to all OTU in a reference
    :return: Sets containing unique ids for both isolates and sequences
    """
    isolate_ids = set()
    sequence_ids = set()

    for otu_path in otu_paths:

        for isolate_path in get_isolate_paths(otu_path):
            isolate_ids.add(isolate_path.name)

            for seq_path in get_sequence_paths(isolate_path):
                sequence_ids.add(seq_path.stem)

    return isolate_ids, sequence_ids


async def get_unique_otu_ids(src_path: Path) -> list:
    """
    Returns a list of all unique OTU ids included in the reference database

    :param src_path: Path to a Virtool reference database directory
    :return: List of all unique OTU ids
    """
    unique_otus = [
        otu_path.stem.split("--")[1]
        for otu_path in src_path.iterdir()
        if otu_path.is_dir()
    ]
    return unique_otus
