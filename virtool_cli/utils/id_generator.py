from random import choice
from string import ascii_letters, ascii_lowercase, digits


def generate_random_alphanumeric(
    length: int = 8,
    mixed_case: bool = False,
    excluded: list | set | None = None
) -> str:
    """
    Generates a random string composed of letters and numbers.

    :param length: the length of the string.
    :param mixed_case: included alpha characters will be mixed case instead of lowercase
    :param excluded: strings that may not be returned.
    :return: a random alphanumeric string.
    """
    if excluded is None:
        excluded_set = set()
    else:
        excluded_set = set(excluded)

    characters = digits + (ascii_letters if mixed_case else ascii_lowercase)

    candidate = "".join([choice(characters) for _ in range(length)])

    if candidate not in excluded_set:
        return candidate

    return generate_random_alphanumeric(length=length, excluded=excluded_set)


def generate_unique_ids(
    n: int = 1,
    length: int = 8,
    mixed_case: bool = False,
    excluded: list | set | None = None
) -> set:
    """
    Generates an n-length list of unique alphanumeric IDs.

    :param n: The number of strings to be generated
    :param mixed_case: included alpha characters will be mixed case instead of lowercase
    :param length: The length of each string
    :param excluded: List of alphanumeric strings that should be excluded from generation
    """
    if excluded is None:
        excluded_set = set()
    else:
        excluded_set = set(excluded)

    new_uniques = set()
    while len(new_uniques) < n:
        new_uniques.add(
            generate_random_alphanumeric(length, mixed_case, excluded_set)
        )
    
    return new_uniques
