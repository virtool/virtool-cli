import os
import subprocess


from Bio import SeqIO
from pathlib import Path
from virtool_cli.vfam_polyprotein import Alignment


def write_abc(blast_results: Path, polyproteins: list) -> Path:
    """
    Takes in blast results file and list of polyproteins to not include, writes a .abc file with desired alignments

    :param blast_results: blast file produced in all_by_all blast step
    :param polyproteins: list of polyprotein like sequences to not be included in output
    """
    abc_path = Path(blast_results).parent / "all_by_all.abc"

    with blast_results.open('r') as blast_file:
        with abc_path.open('w') as abc_file:
            for line in blast_file:

                alignment = Alignment(line)
                if alignment.query not in polyproteins and alignment.subject not in polyproteins:
                    abc_line = '\t'.join([alignment.query, alignment.subject, alignment.evalue]) + "\n"
                    abc_file.write(abc_line)

    return abc_path


def blast_to_mcl(blast_results, polyproteins, inflation_num):
    """
    Converts sequences not included in polyprotein_sequences to a .abc file

    calls mcxload on .abc file to generate a .mci and .tab file

    calls mcl on .tab file to generate newline-separated clusters

    :param blast_results:blast file produced in all_by_all blast step
    :param polyproteins: list of polyprotein like sequences to not be included in output
    :param inflation_num: Inflation number to be used in mcl call
    :return: mcl_path to file containing newline-separated clusters
    """
    abc_file = write_abc(blast_results, polyproteins)

    mci_path = Path(blast_results).parent / "all_by_all.mci"
    tab_path = Path(blast_results).parent / "all_by_all.tab"
    mcl_path = Path(blast_results).parent / "all_by_all.mcl"

    mcxload_cmd = ["mcxload", "--stream-mirror", "-abc", abc_file, "-o", mci_path, "-write-tab", tab_path]
    subprocess.run(mcxload_cmd)

    if not inflation_num:
        mcl_cmd = ["mcl", mci_path, "-use-tab", tab_path, "-o", mcl_path]
    else:
        mcl_cmd = ["mcl", mci_path, "-use-tab", tab_path, "-I", inflation_num, "-o", mcl_path]

    subprocess.run(mcl_cmd)

    return mcl_path


def mcl_to_fasta(mcl_path, collapsed_fasta):
    """
    Takes mcl clusters and fasta file, outputs a list of numbered fasta files to be used later
    """
    mcl_dict = {}
    line_num = 0
    fasta_path = Path(collapsed_fasta).parent .parent / "fasta_files"

    if not fasta_path.exists():
        fasta_path.mkdir()


    with mcl_path.open('r') as handle:
        for line in handle:
            line_num += 1
            fasta_name = f"cluster_{line_num}"

            for record_id in line.rstrip().split("\t"):
                mcl_dict[record_id] = fasta_path / fasta_name

    for record in SeqIO.parse(collapsed_fasta, "fasta"):
        if record.id in mcl_dict:
            with mcl_dict[record.id].open('a') as fasta_path:
                SeqIO.write(record, fasta_path, "fasta")

    return list(mcl_dict.values())




