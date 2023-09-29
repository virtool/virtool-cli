import os.path
import subprocess
import sys

from pathlib import Path
from typing import List, Optional
from structlog import get_logger

NUM_CORES = 16


def batch_muscle_call(fasta_paths: List[Path]) -> List[Path]:
    """
    Takes in a list of paths to FASTA files and makes MSAs using MUSCLE.

    :param fasta_paths: list of FASTA paths to produce MSAs
    :return: list of alignment file paths generated in fasta format
    """
    logger = get_logger()

    msa_paths = []
    for fasta_path in fasta_paths:
        msa_path = Path(f"{fasta_path}.msa")
        log_path = Path(f"{fasta_path}.log")
        msa_paths.append(msa_path)

        muscle_cmd = [
            "muscle",
            "-align",
            fasta_path,
            "-output",
            msa_path,
            "-log",
            log_path,
            "-quiet",
        ]
        try:
            subprocess.run(muscle_cmd, check=True, shell=False)
        except FileNotFoundError:
            logger.error("Dependency 'muscle' not found in path")
            sys.exit(1)

    num_msa = len(msa_paths)

    logger.info(
        f"Produced {num_msa} MSAs from {len(fasta_paths)} FASTA cluster files.",
        count=num_msa,
    )

    return msa_paths


def batch_hmm_call(msa_paths: List[Path]) -> List[Path]:
    """
    Takes in a list of FASTA MSA paths, and builds HMMs for each of the using HMMer.

    :param msa_paths: list of paths to msa files from batch_muscle_call step
    :return: hmm_paths, a list of paths to hmm files produced
    """
    logger = get_logger()

    hmm_paths = []
    for msa_path in msa_paths:
        hmm_path = Path(f"{os.path.splitext(msa_path)[0]}.hmm")
        log_path = Path(f"{os.path.splitext(msa_path)[0]}.log")
        hmm_paths.append(hmm_path)

        hmmer_cmd = [
            "hmmbuild",
            "--informat",
            "afa",
            "-o",
            log_path,
            "--cpu",
            str(NUM_CORES),
            hmm_path,
            msa_path,
        ]
        try:
            subprocess.run(hmmer_cmd, check=True, shell=False)
        except FileNotFoundError:
            logger.error("Dependency hmmbuild not found in path.")
            sys.exit(1)

    num_profiles = len(hmm_paths)
    logger.info(
        f"✔ Collected {num_profiles} HMM profiles produced by hmmbuild.",
        count=num_profiles,
    )

    return hmm_paths


def concatenate_hmms(
    hmm_paths: List[Path], output: Path, prefix: Optional[str]
) -> Path:
    """
    Takes in a list of paths to HMM files containing individual profiles and writes them all to a master results file.

    :param hmm_paths: list of paths to individual hmm files
    :param output: Path to output directory
    :param prefix: Prefix for intermediate and result files
    """
    logger = get_logger()

    output_name = "vfam.hmm"
    if prefix:
        output_name = f"{prefix}_{output_name}"

    output_path = output / Path(output_name)

    with output_path.open("w") as o_handle:
        for hmm_path in hmm_paths:
            with hmm_path.open("r") as h_handle:
                for line in h_handle:
                    o_handle.write(line)

    logger.info(f"Master HMM profile built in output_path", output_path=output_path)

    return output_path
