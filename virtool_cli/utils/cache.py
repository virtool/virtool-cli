from pathlib import Path
import json


def write_to_cache(data: dict, cache: Path):
    with open(cache / "taxid_table.json", "w") as f:
        json.dump(data, f, indent=2)


def generate_taxid_table(src_path: Path):
    """
    :param src_path:
    """
    taxid_table = {}
    for otu_path in src_path.iterdir():
        if not is_otu_directory(otu_path):
            continue

        with open(otu_path / "otu.json", "r") as f:
            otu = json.load(f)

        taxid = int(otu.get('taxid', -1))
        if taxid > 0:
            taxid_table[taxid] = otu_path.name

    return taxid_table


def is_otu_directory(path):
    if path.is_dir:
        if (path / "otu.json").exists():
            return True

    return False
