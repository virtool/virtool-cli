from random import choice
from string import ascii_letters, ascii_lowercase, digits
from typing import Iterable, Optional, Tuple

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

async def get_unique_ids(paths: list) -> Tuple[set, set]:
    """
    Returns sets containing unique random alphanumeric ids for both the isolates and the sequences

    :param paths: List of paths to all OTU in a reference
    :return: Sets containing unique ids for both isolates and sequences
    """
    isolate_ids = set()
    sequence_ids = set()

    for path in paths:
        for isolate_id in path.iterdir():
            if isolate_id.is_dir():
                isolate_ids.add(isolate_id)
                for seq_id in isolate_id.iterdir():
                    if seq_id.name != "isolate.json" and seq_id.name != ".DS_Store":
                        sequence_ids.add(seq_id.name.rstrip(".json"))

    return isolate_ids, sequence_ids