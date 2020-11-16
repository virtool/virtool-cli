import json
import os

from virtool_cli.utils import get_otu_paths, get_otus


def repair(src):
    paths = get_otu_paths(src)
    otus = get_otus(paths)
    otus_to_update = {}

    for path in paths:
        fix_folder_name(path, otus[path])
        otu = fix_taxid(otus[path])

        # if a new otu is returned then it should be updated
        if otu:
            otus_to_update[path] = otu

    write_otus(otus_to_update)


def fix_folder_name(path, otu):
    try:
        folder_name = path.rsplit("/", 1)[1]
        new_folder_name = otu.get("name").replace(" ", "_").replace("/", "_").lower()

        if folder_name != new_folder_name:
            new_path = path.replace(folder_name, new_folder_name)
            os.rename(path, new_path)
    except AttributeError:
        # OTU must not have name
        pass


def fix_taxid(otu):
    try:
        taxid = otu.get("taxid")
        if type(taxid) is not int and not None:
            otu["taxid"] = int(taxid)
            return otu
    except TypeError:
        # OTU must not have a taxid
        return None

    return None


def write_otus(otus):
    for path in otus.keys():
        with open(os.path.join(path, "otu.json"), "w") as f:
            json.dump(otus[path], f, indent=4)


def run(src):
    repair(src)
