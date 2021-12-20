import json
import pathlib
import shutil

from virtool_cli.utils import create_otu_path

OTU_KEYS = ["_id", "name", "abbreviation", "schema", "taxid"]

ISOLATE_KEYS = ["id", "source_type", "source_name", "default"]

SEQUENCE_KEYS = ["_id", "accession", "definition", "host", "sequence"]


def run(src_path: pathlib.Path, output: pathlib.Path):
    """
    Divide a reference.json file from Virtool into a src tree.

    :param src_path: Path to a reference.json file
    :param output: Path to the where the src tree should be generated
    """

    shutil.rmtree(output, ignore_errors=True)
    output.mkdir()

    with open(src_path, "r") as export_handle:
        data = json.load(export_handle)

        for otu in data["otus"]:

            otu_path = build_otu(output, otu)

            isolates = otu.pop("isolates")

            for isolate in isolates:

                isolate_path = build_isolate(otu_path, isolate)

                sequences = isolate.pop("sequences")

                for sequence in sequences:
                    build_sequence(isolate_path, sequence)

        with open(output / "meta.json", "w") as f:
            json.dump({"data_type": data["data_type"], "organism": data["organism"]}, f)


def build_otu(output: pathlib.Path, otu: dict) -> str:
    """
    Creates a directory for all OTUs that begin with a particular
    letter if it doesn't already exist. Generates a directory for a
    given OTU and copies key information about it to a otu.json file.

    :param output: Path to the where the src tree should be generated
    :param otu: Dictionary of an OTU
    :return: Path to a newly generated OTU directory
    """
    lower_name = otu["name"].lower()
    first_letter = lower_name[0]

    try:
        (output / first_letter).mkdir()
    except FileExistsError:
        pass

    otu_path = create_otu_path(lower_name, output, first_letter)
    otu_path.mkdir()

    with open(otu_path / "otu.json", "w") as f:
        if "schema" not in otu:
            otu["schema"] = list()

        json.dump({key: otu.get(key) for key in OTU_KEYS}, f, indent=4)

    return otu_path


def build_isolate(otu_path: pathlib.Path, isolate: dict) -> str:
    """
    Creates a directory for a given isolate and generates
    a isolate.json with key information about it.

    :param otu_path: A path to a specified OTU
    :param isolate: A dictionary containing isolate information
    :return: A path to a newly generated isolate directory
    """
    isolate_path = otu_path / isolate["id"]
    isolate_path.mkdir()

    with open(isolate_path / "isolate.json", "w") as f:
        json.dump({key: isolate[key] for key in ISOLATE_KEYS}, f, indent=4)

    return isolate_path


def build_sequence(isolate_path: pathlib.Path, sequence: dict):
    """
    Generates a JSON file for one of the isolate's sequences

    :param isolate_path: A path to a specified isolate
    :param sequence: A dictionary containing information on one of the isolates' sequences
    """
    with open(isolate_path / "{}.json".format(sequence["_id"]), "w") as f:
        json.dump({key: sequence[key] for key in SEQUENCE_KEYS}, f, indent=4)
