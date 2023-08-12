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
    Generates a list of paths to all OTUs in a src directory.

    :param src_path: Path to a src database directory
    :return: List of paths to all OTU in a src directory
    """
    return [iso_path for iso_path in otu_path.iterdir() if iso_path.is_dir()]

def get_sequence_paths(isolate_path: Path) -> list:
    """
    Generates a list of paths to all OTUs in a src directory.

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
    no_whitespace = name.replace(" ", "_")
    no_plus = no_whitespace.replace('+', 'and')
    no_symbols = re.split(r'[()/-]+', no_plus)

    dirname = no_symbols
    dirname += '--' + id

    return dirname

def parse_otu(path: Path) -> dict:
    """
    Returns a json file in dict form

    :param paths: Path to an OTU in a reference
    :return: OTU data in dict form
    """
    with open(path / "otu.json", "r") as f:
        otu = json.load(f)

    return otu