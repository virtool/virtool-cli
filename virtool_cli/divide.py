import os
import json
import shutil


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

def divide(src_path, output):
    shutil.rmtree(output, ignore_errors=True)
    os.mkdir(output)

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

        with open(os.path.join(output, "meta.json"), "w") as f:
            json.dump({
                "data_type": data["data_type"],
                "organism": data["organism"]
            }, f)


def build_otu(output, otu):
    lower_name = otu["name"].lower()
    first_letter = lower_name[0]

    try:
        os.mkdir(os.path.join(output, first_letter))
    except FileExistsError:
        pass

    otu_path = os.path.join(output, first_letter, lower_name.replace(" ", "_").replace("/", "_"))
    os.mkdir(otu_path)

    with open(os.path.join(otu_path, "otu.json"), "w") as f:
        if "schema" not in otu:
            otu["schema"] = list()
                
        json.dump({key: otu[key] for key in OTU_KEYS}, f, indent=4)
    
    return otu_path

def build_isolate(otu_path, isolate):
    isolate_path = os.path.join(otu_path, isolate["id"])
    os.mkdir(isolate_path)

    with open(os.path.join(isolate_path, "isolate.json"), "w") as f:
        json.dump({key: isolate[key] for key in ISOLATE_KEYS}, f, indent=4)

    return isolate_path

def build_sequence(isolate_path, sequence):
    with open(os.path.join(isolate_path, "{}.json".format(sequence["_id"])), "w") as f:
        json.dump({key: sequence[key] for key in SEQUENCE_KEYS}, f, indent=4)
