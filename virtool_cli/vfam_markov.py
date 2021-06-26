import subprocess

from Bio import SeqIO
from pathlib import Path
from typing import List
from virtool_cli.vfam_polyprotein import Alignment


def write_abc(blast_results_path: Path, polyprotein_ids: List[str], prefix=None) -> Path:
    """
    Takes in BLAST results path and list of polyprotein_ids to not include, writes a .abc file with desired alignments.

    :param blast_results_path: path to blast file produced in all_by_all blast step
    :param polyprotein_ids: list of sequence ids of polyprotein-like records to not be included in output
    :param prefix: Prefix for intermediate and result files
    :return: path to abc file
    """
    abc_name = "blast.abc"
    if prefix:
        abc_name = f"{prefix}_{abc_name}"

    abc_path = blast_results_path.parent / Path(abc_name)

    with blast_results_path.open("r") as blast_file:
        with abc_path.open("w") as abc_file:
            for line in blast_file:

                alignment = Alignment(line)
                if alignment.query not in polyprotein_ids and alignment.subject not in polyprotein_ids:
                    abc_line = "\t".join([alignment.query, alignment.subject, alignment.evalue]) + "\n"
                    abc_file.write(abc_line)

    return abc_path


def blast_to_mcl(blast_results_path: Path, polyprotein_ids: List[str], inflation_num=None, prefix=None) -> Path:
    """
    Converts sequences in blast_results_path to a .abc file.

    Calls mcxload on .abc file to generate a .mci and .tab file.

    Calls mcl on .tab file to generate newline-separated clusters.

    :param blast_results_path: path to blast file produced in all_by_all blast step
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
        "-abc", abc_path,
        "--stream-mirror",
        "--stream-neg-log10",
        "-stream-tf", ""'ceil(200)'"",
        "-o", mci_path,
        "-write-tab", tab_path
    ]
    subprocess.run(mcxload_cmd, check=True, shell=False)

    mcl_cmd = [
        "mcl", mci_path,
        "-use-tab", tab_path,
        "-o", mcl_path
    ]

    if inflation_num:
        mcl_cmd += ["-I", inflation_num]

    subprocess.run(mcl_cmd, check=True, shell=False)

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

    for record in SeqIO.parse(clustered_fasta_path, "fasta"):

        if record.id in mcl_path_dict:

            with mcl_path_dict[record.id].open("a") as fasta_path:
                SeqIO.write(record, fasta_path, "fasta")

    return list(set(mcl_path_dict.values()))
