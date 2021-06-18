import subprocess

from Bio import SeqIO
from pathlib import Path
from typing import List
from virtool_cli.vfam_polyprotein import Alignment


def write_abc(blast_results: Path, polyproteins: list, prefix) -> Path:
    """
    Takes in blast results file and list of polyproteins to not include, writes a .abc file with desired alignments

    :param blast_results: blast file produced in all_by_all blast step
    :param polyproteins: list of polyprotein like sequences to not be included in output
    :param prefix: Prefix for intermediate and result files
    :return: path to abc file
    """
    abc_name = "blast.abc"
    if prefix:
        abc_name = f"{prefix}_{abc_name}"

    abc_path = blast_results.parent / Path(abc_name)

    with blast_results.open("r") as blast_file:
        with abc_path.open("w") as abc_file:
            alignments = (Alignment(line) for line in blast_file)
            while True:
                try:
                    alignment = next(alignments)
                    if alignment.query not in polyproteins and alignment.subject not in polyproteins:
                        abc_line = "\t".join([alignment.query, alignment.subject, alignment.evalue]) + "\n"
                        abc_file.write(abc_line)
                except StopIteration:
                    break
    return abc_path


def blast_to_mcl(blast_results: Path, polyproteins: list, inflation_num, prefix) -> Path:
    """
    Converts sequences not included in polyprotein_sequences to a .abc file

    calls mcxload on .abc file to generate a .mci and .tab file

    calls mcl on .tab file to generate newline-separated clusters

    :param blast_results:blast file produced in all_by_all blast step
    :param polyproteins: list of polyprotein like sequences to not be included in output
    :param inflation_num: Inflation number to be used in mcl call
    :param prefix: Prefix for intermediate and result files
    :return: mcl_path to file containing newline-separated clusters
    """
    abc_path = write_abc(blast_results, polyproteins, prefix)

    mci_name = "blast.mci"
    tab_name = "blast.tab"
    mcl_name = "blast.mcl"

    if prefix:
        mci_name = f"{prefix}_{mci_name}"
        tab_name = f"{prefix}_{tab_name}"
        mcl_name = f"{prefix}_{mcl_name}"

    mci_path = blast_results.parent / Path(mci_name)
    tab_path = blast_results.parent / Path(tab_name)
    mcl_path = blast_results.parent / Path(mcl_name)

    mcxload_cmd = ["mcxload", "-abc", abc_path, "--stream-mirror", "--stream-neg-log10", "-stream-tf", ""'ceil(200)'"",
                   "-o", mci_path, "-write-tab", tab_path]
    subprocess.run(mcxload_cmd, check=True, shell=False)

    mcl_cmd = ["mcl", mci_path, "-use-tab", tab_path, "-o", mcl_path]
    if inflation_num:
        mcl_cmd += ["-I", inflation_num]

    subprocess.run(mcl_cmd, check=True, shell=False)

    return mcl_path


def mcl_to_fasta(mcl_path: Path, clustered_fasta: Path, prefix) -> List[Path]:
    """
    Takes mcl clusters and a clustered fasta file, creates numbered fasta files for each mcl cluster

    :param mcl_path: path to mcl results file from blast_to_mcl step
    :param clustered_fasta: path to clustered fasta file from cd-hit step
    :param prefix: Prefix for intermediate and result files
    :return: list of paths to seperated fasta files
    """
    fasta_path = clustered_fasta.parent
    mcl_dict = {}
    line_num = 0

    with mcl_path.open('r') as handle:
        lines = (line for line in handle)
        while True:
            try:
                line = next(lines)
                line_num += 1

                fasta_name = f"cluster_{line_num}"
                if prefix:
                    fasta_name = f"{prefix}_{fasta_name}"

                for record_id in line.strip().split("\t"):
                    mcl_dict[record_id] = fasta_path / Path(fasta_name)

            except StopIteration:
                break

    for record in SeqIO.parse(clustered_fasta, "fasta"):
        if record.id in mcl_dict:
            with mcl_dict[record.id].open('a') as fasta_path:
                SeqIO.write(record, fasta_path, "fasta")

    return list(set(mcl_dict.values()))
