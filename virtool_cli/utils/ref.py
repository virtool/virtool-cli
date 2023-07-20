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
    # otu_name.replace(" ", "_").replace("/", "_").lower()

    return no_symbols