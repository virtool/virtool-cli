import os.path
import subprocess


from pathlib import Path
NUM_CORES = 16


def batch_muscle_call(fasta_files):
    """
    This function takes in a list of fasta files and makes msas using MUSCLE

    :param fasta_files: list of fasta files to pruduce msas
    :return: list of alignment files generated in fasta format
    """
    msa_files = []
    for fasta_file in fasta_files:

        msa_file = f"{fasta_file}.msa"
        log_file = f"{fasta_file}.log"
        msa_files.append(msa_file)

        muscle_cmd = ["muscle", "-in", fasta_file, "-out", msa_file, "-log", log_file, "-quiet"]
        subprocess.call(muscle_cmd)

    return msa_files


def batch_hmm_call(msa_files):
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
        subprocess.call(hmm_build_cmd)

    return hmm_files


def concatenate_hmms(hmm_files, output, prefix):
    if prefix:
        output_name = f"{prefix}_master.hmm"
    else:
        output_name = "master.hmm"
    output_file = output / Path(output_name)

    with output_file.open('w') as o_handle:
        for hmm_file in hmm_files:
            with hmm_file.open('r') as h_handle:
                for line in h_handle:
                    o_handle.write(line)

    return output_file
