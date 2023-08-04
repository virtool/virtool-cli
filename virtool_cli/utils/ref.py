import json
import re
from pathlib import Path


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

    :param src_path: Path to a src database directory
    :return: List of paths to all OTU in a src directory
    """
    return [iso_path for iso_path in otu_path.iterdir() if iso_path.is_dir()]

def get_sequence_paths(isolate_path: Path) -> list:
    """
    Generates a list of paths to all sequences in an isolate directory.

    :param src_path: Path to a src database directory
    :return: List of paths to all OTU in a src directory
    """
    sequence_ids = [
        i
        for i in isolate_path.glob('*.json')
        if i.name != "isolate.json" and i.name[0] != "."
    ]

    return sequence_ids

def generate_otu_dirname(name: str, id: str = ''):
    """
    Takes in a human-readable string, replaces whitespace and symbols
    and adds the hash id
    
    :param name: Human-readable, searchable name of the OTU
    :param id: ID hash of OTU 
    :return: A directory name in the form of 'converted_otu_name--taxid'
    """
    no_plus = name.replace('+', 'plus ')
    no_symbols = re.split(r'[():/-]+', no_plus)
    joined = ' '.join(no_symbols)
    no_whitespace = re.sub(r'[\s]+', "_", joined)

    dirname = no_whitespace.lower()
    dirname += '--' + id

    return dirname

def search_otu_by_id(src_path: Path, otu_id: str):
    matches = [otu for otu in src_path.glob(f'*--{otu_id}')]
    if matches:
        return matches[0]
    else:
        return FileNotFoundError

def parse_otu(path: Path) -> dict:
    """
    Returns a json file in dict form

    :param paths: Path to an OTU in a reference
    :return: OTU data in dict form
    """
    with open(path / "otu.json", "r") as f:
        otu = json.load(f)

    return otu

def parse_isolates(otu_path: Path) -> dict:
    """
    Returns a json file in dict form

    :param paths: Path to an OTU in a reference
    :return: OTU data in dict form
    """
    isolates = {}
    for iso_path in get_isolate_paths(otu_path):
        with open(iso_path / "isolate.json", "r") as f:
            isolate = json.load(f)
        isolates[iso_path.name] = isolate

    return isolates

def map_otus(paths: list) -> dict:
    """
    Returns a mapping of every OTU path to their deserialized OTU dictionary.

    :param paths: List of paths to all OTU in a reference
    :return: A mapping of every OTU path to its OTU dictionary
    """
    path_taxid_map = dict()

    for path in paths:
        path_taxid_map[path] = parse_otu(path)

    return path_taxid_map