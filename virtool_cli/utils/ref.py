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

def generate_otu_dirname(otu_name: str):
    """
    """
    no_whitespace = otu_name.replace(" ", "_")
    no_plus = no_whitespace.replace('+', 'and')
    no_symbols = re.split(r'[()/-]+', no_plus)

    return no_symbols

def parse_otu(path: Path) -> dict:
    """
    Returns a json file in dict form

    :param paths: Path to an OTU in a reference
    :return: OTU data in dict form
    """
    with open(path / "otu.json", "r") as f:
        otu = json.load(f)

    return otu