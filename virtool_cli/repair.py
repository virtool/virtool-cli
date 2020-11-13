from virtool_cli.utils import get_otu_paths, get_otus
import os


def repair(src):
    paths = get_otu_paths(src)
    otus = get_otus(paths)

    for path in paths:
        fix_folder_name(path, otus[path])


def fix_folder_name(path, otu):
    try:
        folder_name = path.rsplit("/", 1)[1]
        new_folder_name = otu.get("name").replace(" ", "_").lower()

        if folder_name != new_folder_name:
            new_path = path.replace(folder_name, new_folder_name)
            os.rename(path, new_path)
    except AttributeError:
        # OTU doesn't have a name
        pass


def run(src):
    repair(src)


run("tests/files/src")
