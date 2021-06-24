import os.path
import subprocess

from pathlib import Path
from typing import List, Optional

NUM_CORES = 16


def batch_muscle_call(fasta_paths: List[Path]) -> List[Path]:
    """
    Takes in a list of paths to FASTA files and makes MSAs using MUSCLE.

    :param fasta_paths: list of FASTA paths to produce MSAs
    :return: list of alignment file paths generated in fasta format
    """
    msa_paths = []
    for fasta_path in fasta_paths:
        msa_path = Path(f"{fasta_path}.msa")
        log_path = Path(f"{fasta_path}.log")
        msa_paths.append(msa_path)

        muscle_cmd = [
            "muscle",
            "-in", fasta_path,
            "-out", msa_path,
            "-log", log_path,
            "-quiet"
        ]

        subprocess.run(muscle_cmd, check=True, shell=False)

    return msa_paths


def batch_hmm_call(msa_paths: list) -> List[Path]:
    """
    Takes in a list of FASTA MSA paths, and builds HMMs for each of the using HMMer.

    :param msa_paths: list of paths to msa files from batch_muscle_call step
    :return: hmm_paths, a list of paths to hmm files produced
    """
    hmm_paths = []
    for msa_path in msa_paths:
        hmm_path = Path(f"{os.path.splitext(msa_path)[0]}.hmm")
        log_path = Path(f"{os.path.splitext(msa_path)[0]}.log")
        hmm_paths.append(hmm_path)

        hmm_build_cmd = ["hmmbuild", "--informat", "afa", "-o", log_file, "--cpu", str(NUM_CORES), hmm_file, msa_file]
        subprocess.run(hmm_build_cmd, check=True, shell=False)

    return hmm_paths


def concatenate_hmms(hmm_paths: list, output: Path, prefix: Optional[str]) -> Path:
    """
    Takes in a list of paths to hmm files containing individual profiles and writes them all to a master results file.

    :param hmm_paths: list of paths to individual hmm files
    :param output: Path to output directory
    :param prefix: Prefix for intermediate and result files
    """
    output_name = "master.hmm"
    if prefix:
        output_name = f"{prefix}_{output_name}"

    output_path = output / Path(output_name)

    with output_path.open("w") as o_handle:
        for hmm_path in hmm_paths:
            with hmm_path.open("r") as h_handle:
                for line in h_handle:
                    o_handle.write(line)

    return output_path
