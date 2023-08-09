from random import choice
from string import ascii_letters, ascii_lowercase, digits
from typing import Iterable, Optional, Tuple

from virtool_cli.utils.ref import get_isolate_paths, get_sequence_paths

def generate_random_alphanumeric(
    length: int = 8,
    mixed_case: bool = False,
    excluded: Optional[Iterable[str]] = None
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

    return generate_random_alphanumeric(length=length, excluded=excluded)

def generate_hashes(
    excluded: list,
    n: int = 1,
    length: int = 8,
    mixed_case: bool = False
):
    """
    """
    new_uniques = set()
    while len(new_uniques) < n:
        new_uniques.add(generate_random_alphanumeric(length, mixed_case, excluded))
    
    return new_uniques
    

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