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

def build(src_path, output, indent, version):

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

            otu, isolate_ids = parse_otu(src_path, otu_path)

            for isolate_path in [os.path.join(otu_path, i) for i in isolate_ids]:
                
                isolate, sequence_ids = parse_isolate(src_path, isolate_path)

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


def parse_meta(src_path):
    try:
        with open(os.path.join(src_path, "meta.json"), "r") as f:
            return json.load(f)
    except FileNotFoundError:
        return dict()


def parse_alpha(src_path, alpha):
    return [os.path.join(src_path, alpha, otu) for otu in os.listdir(os.path.join(src_path, alpha))]

    
def parse_otu(src_path, otu_path):
    with open(os.path.join(otu_path, "otu.json"), "r") as f:
        otu = json.load(f)

    otu["isolates"] = list()

    isolate_ids = [i for i in os.listdir(otu_path) if i != "otu.json" and i[0] != "."]

    return otu, isolate_ids


def parse_isolate(src_path, isolate_path):
    with open(os.path.join(isolate_path, "isolate.json"), "r") as f:
        isolate = json.load(f)

    sequence_ids = [i for i in os.listdir(isolate_path) if i != "isolate.json" and i[0] != "."]

    isolate["sequences"] = list()

    return isolate, sequence_ids
