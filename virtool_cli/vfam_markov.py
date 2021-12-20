import subprocess
import sys

from Bio import SeqIO, SearchIO
from pathlib import Path
from typing import List
from virtool_cli.vfam_console import console


def write_abc(
    blast_results_path: Path, polyprotein_ids: List[str], prefix=None
) -> Path:
    """
    Takes in BLAST results file path and list of polyprotein sequence IDs to not include in output.

    Gathers query ID, sequence ID, and evalue for each QueryResult object in BLAST results file and writes to abc_path.

    :param blast_results_path: path to BLAST file produced in blast_all_by_all step
    :param polyprotein_ids: list of sequence IDs of polyprotein-like records to not be included in output
    :param prefix: Prefix for intermediate and result files
    :return: path to .abc file
    """
    abc_name = "blast.abc"
    if prefix:
        abc_name = f"{prefix}_{abc_name}"

    abc_path = blast_results_path.parent / Path(abc_name)

    with abc_path.open("w") as abc_file:

        for query_result in SearchIO.parse(blast_results_path, "blast-tab"):
            for hit in query_result.hits:
                q_id = query_result.id
                s_id = hit.id
                for hsp in hit.hsps:
                    e_value = hsp.evalue

                    if q_id not in polyprotein_ids and s_id not in polyprotein_ids:
                        abc_line = "\t".join([q_id, s_id, str(e_value)]) + "\n"
                        abc_file.write(abc_line)

    return abc_path


def blast_to_mcl(
    blast_results_path: Path,
    polyprotein_ids: List[str],
    inflation_num=None,
    prefix=None,
) -> Path:
    """
    Converts sequences in blast_results_path to a .abc file.

    Calls mcxload on .abc file to generate a .mci and .tab file.

    Calls mcl on .tab file to generate newline-separated clusters.

    :param blast_results_path: path to BLAST file produced in all_by_all BLAST step
    :param polyprotein_ids: list of sequence ids of polyprotein like sequences to not be included in output
    :param inflation_num: Inflation number to be used in mcl call
    :param prefix: Prefix for intermediate and result files
    :return: mcl_path to file containing newline-separated clusters
    """
    abc_path = write_abc(blast_results_path, polyprotein_ids, prefix)

    mci_name = "blast.mci"
    tab_name = "blast.tab"
    mcl_name = "blast.mcl"

    if prefix:
        mci_name = f"{prefix}_{mci_name}"
        tab_name = f"{prefix}_{tab_name}"
        mcl_name = f"{prefix}_{mcl_name}"

    mci_path = blast_results_path.parent / Path(mci_name)
    tab_path = blast_results_path.parent / Path(tab_name)
    mcl_path = blast_results_path.parent / Path(mcl_name)

    mcxload_cmd = [
        "mcxload",
        "-abc",
        abc_path,
        "--stream-mirror",
        "--stream-neg-log10",
        "-stream-tf",
        "" "ceil(200)" "",
        "-o",
        mci_path,
        "-write-tab",
        tab_path,
    ]
    try:
        subprocess.run(mcxload_cmd, check=True, shell=False)
    except FileNotFoundError:
        console.print("Dependency mcxload not found in path", style="red")
        sys.exit(1)

    mcl_cmd = ["mcl", mci_path, "-use-tab", tab_path, "-o", mcl_path]

    if inflation_num:
        mcl_cmd += ["-I", inflation_num]

    try:
        subprocess.run(mcl_cmd, check=True, shell=False)
    except FileNotFoundError:
        console.print("Dependency mcl not found in path.", style="red")
        sys.exit(1)

    return mcl_path


def mcl_to_fasta(mcl_path: Path, clustered_fasta_path: Path, prefix=None) -> List[Path]:
    """
    Takes MCL clusters and a path to a clustered FASTA file, creates numbered FASTA files for each MCL cluster.

    :param mcl_path: path to mcl results file from blast_to_mcl step
    :param clustered_fasta_path: path to clustered FASTA from cd-hit step
    :param prefix: Prefix for intermediate and result files
    :return: list of paths to seperated  files
    """
    fasta_path = clustered_fasta_path.parent
    mcl_path_dict = {}

    with mcl_path.open("r") as handle:
        for line_num, line in enumerate(handle):
            fasta_name = f"cluster_{line_num + 1}"
            if prefix:
                fasta_name = f"{prefix}_{fasta_name}"

            for record_id in line.strip().split("\t"):
                mcl_path_dict[record_id] = fasta_path / Path(fasta_name)
                open(mcl_path_dict[record_id], "w").close()

    for record in SeqIO.parse(clustered_fasta_path, "fasta"):

        if record.id in mcl_path_dict:
            with mcl_path_dict[record.id].open("a") as fasta_path:
                SeqIO.write(record, fasta_path, "fasta")

    console.print(
        f"✔ Produced {len(list(set(mcl_path_dict.values())))} FASTA cluster files from MCL clusters.",
        style="green",
    )
    return list(set(mcl_path_dict.values()))
