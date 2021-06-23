import collections
import json
import operator
import os.path
import subprocess
from pathlib import Path

from Bio import GenBank, SeqIO
from Bio import Entrez
from typing import Tuple, List


def get_taxonomy(seq_ids: List[str]) -> Tuple[str, str]:
    """
    Takes in list of sequence IDs and makes calls to NCBI database to gather family and genus information for each ID

    :param seq_ids: list of sequence IDs
    :return: string representations of family and genus dictionaries containing occurrences of each
    """
    families = {}
    genera = {}

    for seq_id in seq_ids:
        handle = Entrez.efetch(db="protein", id=seq_id, rettype="gb", retmode="text")

        for record in GenBank.parse(handle):
            family = record.taxonomy[-2]
            if family in families:
                families[family] += 1
            else:
                families[family] = 1

            genus = record.taxonomy[-1]
            if genus in genera:
                genera[genus] += 1
            else:
                genera[genus] = 1

    return str(families), str(genera)


def parse_log(cluster_name: str, output: Path) -> dict:
    """
    Parses log files and gathers count, alen, length, eff_nseqs, mean entropy and total_entropy

    :param cluster_name: cluster name from which to gather log file
    :param output: path to output folder where .log files are stored
    :return: log_data, dictionary containing log_data
    """
    log_path = output / "intermediate_files" / Path(f"{cluster_name}.log")
    log_data = {}

    with log_path.open("r") as handle:
        for line in handle:
            line_data = line.strip().split()
            if len(line_data) > 0 and "#" not in line_data[0]:
                log_data["count"] = int(line_data[2])
                log_data["alen"] = int(line_data[3])
                log_data["length"] = int(line_data[4])
                log_data["eff_nseq"] = float(line_data[5])
                log_data["mean_entropy"] = float(line_data[6])
                log_data["total_entropy"] = log_data["mean_entropy"] * float(log_data["count"])
                break

    return log_data


def parse_stat(cluster_name: str, output: Path) -> Tuple[float, float]:
    """
    Calls hmmstat and captures compKL and p rel E data from CompletedProcess output

    :param cluster_name: name of cluster file
    :param output: path to output directory
    :return: p rel E and compKL values
    """
    hmm_file = output / "intermediate_files" / f"{cluster_name}.hmm"
    hmmstat_cmd = ["hmmstat", hmm_file]
    stat_data = str(subprocess.run(hmmstat_cmd, capture_output=True)).split("\\n")

    prelE = str
    compKL = str
    for line in stat_data:
        if not line.startswith("#") and not line.startswith("CompletedProcess"):
            line = line.strip().split()
            prelE = float(line[-2])
            compKL = float(line[-1])
            break

    return float(prelE), float(compKL)


def get_gi(seq_id):
    """
    Makes calls to NCBI database to gather gi information for each sequence

    TODO: figure out why these calls aren't returning anything
    :param seq_id: sequence id for sequence, used to access correct .bgk file
    :return: record.gi, GenInfo identifier
    """
    handle = Entrez.efetch(db="protein", id=seq_id, rettype="gb", retmode="text")
    for record in GenBank.parse(handle):
        if record.gi:
            return record.gi
    return None


def get_names(annotation):
    """
    Returns three most common names in "entries" list and returns them to be output in json file

    :param annotation: dictionary containing all information to be output to json file for one cluster
    :return: top three names
    """
    names = [entry["name"] for entry in annotation["entries"]]
    top_three = collections.Counter(names).most_common(3)
    return [entry[0] for entry in top_three]


def clusters_to_json(cluster_files: List[Path], output: Path):
    """
    parses all filtered fasta cluster files and creates annotation dictionaries

    all data from each annotation is written to master.json file

    :param cluster_files: clustered, filtered fasta files from vfam pipeline
    :param output: Path to output directory containing intermediate files from vfam pipeline
    """
    output_path = output / "master.json"
    annotations = list()

    for cluster_file in cluster_files:
        annotation = {"cluster": os.path.basename(cluster_file).split("_")[1],
                      "names": list(),
                      "entries": list()
                      }

        log_data = parse_log(os.path.basename(cluster_file), output)
        for key in log_data:
            annotation[key] = log_data[key]

        stat_data = parse_stat(os.path.basename(cluster_file), output)
        annotation["prelE"] = stat_data[0]
        annotation["compKL"] = stat_data[1]

        with cluster_file.open("r") as handle:
            seq_ids = []
            for record in SeqIO.parse(handle, "fasta"):
                seq_ids.append(record.id)
                name = record.description.split("[")[0].split()
                name = " ".join(name[1:])

                annotation["entries"].append({
                    "gi": get_gi(record.id),
                    "accession": record.id,
                    "name": name,
                    "organism": record.description.split("[")[1].replace("]", "").strip()
                })

            taxonomy = get_taxonomy(seq_ids)
            annotation["families"] = json.loads(taxonomy[0].replace("'", '"'))
            annotation["genera"] = json.loads(taxonomy[1].replace("'", '"'))

        annotation["names"] = get_names(annotation)
        annotation["entries"] = sorted(annotation["entries"], key=operator.itemgetter("accession"))

        annotations.append(annotation)

    annotations = sorted(annotations, key=operator.itemgetter("cluster"))
    with open(output_path, "wt") as f:
        json.dump(annotations, f, indent=4)

