import os.path
import subprocess

from pathlib import Path
from typing import List

NUM_CORES = 16


def batch_muscle_call(fasta_files: list) -> List[Path]:
    """
    This function takes in a list of fasta files and makes msas using MUSCLE

    :param fasta_files: list of fasta files to produce MSAs
    :return: list of alignment files generated in fasta format
    """
    msa_files = []
    for fasta_file in fasta_files:

        msa_file = Path(f"{fasta_file}.msa")
        log_file = Path(f"{fasta_file}.log")
        msa_files.append(msa_file)

        muscle_cmd = ["muscle", "-in", fasta_file, "-out", msa_file, "-log", log_file, "-quiet"]
        subprocess.run(muscle_cmd, check=True, shell=False)

    return msa_files


def batch_hmm_call(msa_files: list) -> List[Path]:
    """
    Takes in a  list of fasta msa files, and builds HMMs for each of the using HMMer

    :param msa_files: list of msa files from batch_muscle_call step
    :return: hmm_files, a list of hmm files produced
    """
    hmm_files = []
    for msa_file in msa_files:

        hmm_file = Path(f"{os.path.splitext(msa_file)[0]}.hmm")
        log_file = Path(f"{os.path.splitext(msa_file)[0]}.log")
        hmm_files.append(hmm_file)

        hmm_build_cmd = ["hmmbuild", "--informat", "afa", "-o", log_file, "--cpu", str(NUM_CORES), hmm_file, msa_file]
        subprocess.run(hmm_build_cmd, check=True, shell=False)

    return hmm_files


def concatenate_hmms(hmm_files: list, output: Path, prefix) -> Path:
    """
    Takes in a list of hmm files containing individual profiles and writes them all to a master results file

    :param hmm_files: list of individual hmm files
    :param output: Path to output directory
    :param prefix: Prefix for intermediate and result files
    """
    output_name = "master.hmm"
    if prefix:
        output_name = f"{prefix}_{output_name}"

    output_file = output / Path(output_name)

    with output_file.open("w") as o_handle:
        for hmm_file in hmm_files:
            with open(hmm_file) as h_handle:
                lines = (line for line in h_handle)
                while True:
                    try:
                        o_handle.write(next(lines))
                    except StopIteration:
                        break
    return output_file


def organize_intermediates(output: Path):
    """
    Organizes intermediate files by type and sorts them into different folders in the output directory

    :param output: Path to output directory containing intermediate and master file
    """
    intermediate_path = output / Path("intermediate_files")
    file_names = os.listdir(intermediate_path)

    hmm_files = output / "hmm_files"
    if not hmm_files.exists():
        hmm_files.mkdir()

    blast_files = output / "blast_files"
    if not blast_files.exists():
        blast_files.mkdir()

    fasta_files = output / "fasta_files"
    if not fasta_files.exists():
        fasta_files.mkdir()

    log_files = output / "log_files"
    if not log_files.exists():
        log_files.mkdir()

    msa_files = output / "msa_files"
    if not msa_files.exists():
        msa_files.mkdir()

    for file_name in file_names:
        current_path = output / "intermediate_files" / file_name

        if ".hmm" in file_name:
            new_path = hmm_files / file_name
            os.rename(current_path, new_path)

        elif "blast" in file_name:
            new_path = blast_files / file_name
            os.rename(current_path, new_path)

        elif ".msa" in file_name:
            new_path = msa_files / file_name
            os.rename(current_path, new_path)

        elif ".log" in file_name:
            new_path = log_files / file_name
            os.rename(current_path, new_path)

        else:
            new_path = fasta_files / file_name
            os.rename(current_path, new_path)

    os.rmdir(intermediate_path)
