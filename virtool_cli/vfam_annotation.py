import csv
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
    Entrez.email = "eroberts9789@gmail.com"
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
    Parses log files and gathers nseqs, alen, mlen, eff_nseqs and relent information

    :param cluster_name: cluster name to gather log file for
    :param output: path to output folder where .log files are stored
    :return: log_data, dictionary containing log log_data
    """
    log_path = output / "intermediate_files" / Path(f"{cluster_name}.log")

    log_data = {
        "nseqs": int,
        "alen": int,
        "mlen": int,
        "eff_nseq": float,
        "relent": float,
    }

    with log_path.open("r") as handle:
        for line in handle:
            line_data = line.strip().split()
            if len(line_data) > 0 and "#" not in line_data[0]:
                log_data["nseq"] = line_data[2]
                log_data["alen"] = line_data[3]
                log_data["mlen"] = line_data[4]
                log_data["eff_nseq"] = line_data[5]
                log_data["relent"] = line_data[6]
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

    return prelE, compKL


def parse_clusters(cluster_files: List[Path], output: Path):
    """
    parses all filtered fasta cluster files and creates annotation dictionaries

    makes calls to get_taxononmy(), parse_log(), and parse_stat() to gather all info to write in write_annotation()

    :param cluster_files: clustered, filtered fasta files from vfam pipeline
    :param output: Path to output directory containing intermediate files from vfam pipeline
    """
    for cluster_file in cluster_files:
        annotation = {
            "cluster_name": os.path.basename(cluster_file),
            "nseq": int,
            "alen": int,
            "mlen": int,
            "eff_nseq": float,
            "relent": float,
            "prelE": float,
            "compKL": float,
            "families": str,
            "genera": str,
            "headers": []
        }
        log_data = parse_log(annotation["cluster_name"], output)
        for key in log_data:
            annotation[key] = log_data[key]

        stat_data = parse_stat(annotation["cluster_name"], output)
        annotation["prelE"] = stat_data[0]
        annotation["compKL"] = stat_data[1]

        with cluster_file.open("r") as handle:
            seq_ids = []
            for record in SeqIO.parse(handle, "fasta"):
                seq_ids.append(record.id)
                annotation["headers"].append(record.description)

            taxonomy = get_taxonomy(seq_ids)
            annotation["families"] = taxonomy[0]
            annotation["genera"] = taxonomy[1]

        write_annotation(annotation, output)


def write_annotation(annotation: dict, output: Path):
    """
    Writes annotation files using csv

    :param annotation: dictionary containing all annotation information to be output
    :param output: path to output directory where annotation file will be written
    """
    output_name = str(annotation["cluster_name"]) + ".annotation"
    output_path = output / "intermediate_files" / output_name

    with open(output_path, "w") as handle:
        writer = csv.writer(handle, delimiter="\t", escapechar=" ", quoting=csv.QUOTE_NONE)

        writer.writerow(["CLUSTER",
                         annotation["cluster_name"].split("_")[1]])
        writer.writerow(["NUM_SEQ",
                         annotation["nseq"]])
        writer.writerow(["EFFECTIVE_NUM_SEQS",
                         annotation["eff_nseq"]])
        writer.writerow(["LENGTH",
                         annotation["mlen"]])
        writer.writerow(["RELATIVE_ENTROPY_PER_POSITION",
                         annotation["relent"]])
        writer.writerow(["TOTAL_RELATIVE_ENTROPY",
                         "{:.2f}".format((float(annotation["relent"]) * float(annotation["mlen"])))])
        writer.writerow(["ALIGNED_COLUMNS",
                         annotation["alen"]])
        writer.writerow(["MEAN_RELATIVE_ENTROPY_PER_POSITION",
                         annotation["prelE"]])
        writer.writerow(["KULLBACK_LEIBLER_DIVERGENCE",
                         annotation["compKL"]])
        writer.writerow(["FAMILIES",
                         annotation["families"]])
        writer.writerow(["GENERA",
                         annotation["genera"]])

        writer.writerow(["FASTA SEQUENCE TITLES:"])
        writer = csv.writer(handle, delimiter="\n", escapechar=" ", quoting=csv.QUOTE_NONE)
        writer.writerow(annotation["headers"])
