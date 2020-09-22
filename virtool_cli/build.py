import os
import json
import arrow

OTU_KEYS = [
    "_id",
    "name",
    "abbreviation",
    "schema"
]

ISOLATE_KEYS = [
    "id",
    "source_type",
    "source_name",
    "default"
]

SEQUENCE_KEYS = [
    "_id",
    "accession",
    "definition",
    "host",
    "sequence"
]


def build(src_path: str, output: str, indent: bool, version: str):
    """
    Build a Virtool reference JSON file from a data directory.

    Parameters:
        src_path (str): Path to database src directory
        output (str): The output path for the reference.json file
        indent (bool): A flag to indicate whether the output file should be indented
        version (str): The version string to include in the reference.json file

    """
    meta = parse_meta(src_path)

    data = {
        "data_type": meta.get("data_type", "genome"),
        "organism": meta.get("organism", ""),
    }

    otus = list()

    alpha_paths = os.listdir(src_path)

    try:
        alpha_paths.remove("meta.json")
    except ValueError:
        pass

    for alpha in alpha_paths:
        otu_paths = parse_alpha(src_path, alpha)

        for otu_path in otu_paths:

            otu, isolate_ids = parse_otu(otu_path)

            for isolate_path in [os.path.join(otu_path, i) for i in isolate_ids]:

                isolate, sequence_ids = parse_isolate(isolate_path)

                for sequence_path in [os.path.join(isolate_path, i) for i in sequence_ids]:
                    with open(sequence_path, "r") as f:
                        sequence = json.load(f)

                    isolate["sequences"].append(sequence)

                otu["isolates"].append(isolate)

            otus.append(otu)

    with open(output, "w") as f:
        data.update({
            "otus": otus,
            "name": version,
            "created_at": arrow.utcnow().isoformat()
        })

        json.dump(data, f, indent=4 if indent else None, sort_keys=True)


def parse_meta(src_path: str) -> list:
    """
    Deserializes and returns meta.json if found, else returns an empty dictionary.

    Parameters:
        src_path (str): Path to database src directory

    Returns:
        The deserialized meta.json object or an empty dictionary.

    """
    try:
        with open(os.path.join(src_path, "meta.json"), "r") as f:
            return json.load(f)
    except FileNotFoundError:
        return dict()


def parse_alpha(src_path: str, alpha: str) -> list:
    """
    Generates and returns a list with every OTU in the directory of the given alpha

    Parameters:
        src_path (str): Path to database src directory
        alpha (str): Alphabetical character that denotes the name of the directory to search

    Returns:
        A list containing all the OTU in the given directory

    """
    return [os.path.join(src_path, alpha, otu) for otu in os.listdir(os.path.join(src_path, alpha))]


def parse_otu(otu_path: str) -> (dict, list):
    """
    Creates an list of isolate IDs for a given OTU

    Parameters:
        otu_path (str): Path to a OTU directory

    Returns:
        otu (dict): The dictionary for a given OTU
        isolate_ids (list): List of all isolate IDs for a given OTU

    """
    with open(os.path.join(otu_path, "otu.json"), "r") as f:
        otu = json.load(f)

    otu["isolates"] = list()

    isolate_ids = [i for i in os.listdir(otu_path) if i != "otu.json" and i[0] != "."]

    return otu, isolate_ids


def parse_isolate(isolate_path: str) -> (dict, list):
    """
    Creates a list of sequence IDs for a given sequence

    Parameters:
        isolate_path (str): Path to a isolate directory

    Returns:
        isolate (dict): The dictionary for a given OTU
        sequence_ids (list): List of all sequence IDs for a given isolate

    """
    with open(os.path.join(isolate_path, "isolate.json"), "r") as f:
        isolate = json.load(f)

    isolate["sequences"] = list()

    sequence_ids = [i for i in os.listdir(isolate_path) if i != "isolate.json" and i[0] != "."]

    return isolate, sequence_ids
