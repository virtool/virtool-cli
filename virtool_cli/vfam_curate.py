import sys
from urllib.error import HTTPError

from Bio import SeqIO, Entrez
from pathlib import Path
from typing import List


def get_input_paths(src_path: Path) -> List[Path]:
    """
    Takes in a path to directory containing input FASTA files, returns a list of paths to the FASTA files.

    :param src_path: Path to input source directory containing unfiltered FASTA files
    :return: input_paths, list of paths to input files if any files are found
    """
    input_paths = list(src_path.iterdir())

    if input_paths:
        return input_paths

    print("No files found in input directory")
    sys.exit(1)


def group_input_paths(input_paths: List[Path], no_named_phages: bool) -> list:
    """
    Takes input_paths from get_input_paths() as input and yields records.

    Filters out records with "phage" in their description if no_named_phages is True.

    :param input_paths: list of paths to input FASTA files
    :param no_named_phages: bool that dictates whether phage records are filtered out by name or not
    :return: list of all records found in input paths
    """
    for input_path in input_paths:
        for record in SeqIO.parse(input_path, "fasta"):
            if no_named_phages:
                if "phage" not in record.description:
                    yield record
            else:
                yield record


def remove_dupes(records: iter, sequence_min_length: int) -> list:
    """
    Iterates through records and filters out repeated records and sequences shorter than sequence_min_length.

    :param records: iterable of records gathered in group_input_paths()
    :param sequence_min_length: minimum length of sequence to be included in output
    """
    record_seqs = list()
    no_dupes = list()
    for record in records:
        if record.seq not in record_seqs and len(record.seq) > sequence_min_length:
            record_seqs.append(record.seq)
            no_dupes.append(record)
    return no_dupes


def get_taxonomy(records: iter):
    """
   Makes calls to NCBI database using Bio.Entrez to gather taxonomic information for each record.

   If record isn't found in NCBI, record isn't included in output.

   :param records: list of records from group_input_paths() step.
   :return: record_taxonomy, a dictionary containing {record ID: taxonomy} pairs
   """
    record_taxonomy = dict()
    for record in records:
        try:
            handle = Entrez.efetch(db="protein", id=record.id, rettype="gb", retmode="text")
            for seq_record in SeqIO.parse(handle, "genbank"):
                record_taxonomy[seq_record.id] = seq_record.annotations["taxonomy"][-2:]

        except (HTTPError, AttributeError):
            continue
    return record_taxonomy


def write_curated_recs(records: iter, output: Path, taxonomy_ids: List[str], prefix=None) -> Path:
    """
    Writes filtered records to output if taxonomic information was found for the record in get_taxonomy().

    :param records: list of records from all protein files without keyword "phage"
    :param output: Path to output directory for profile HMMs and intermediate files
    :param taxonomy_ids: list of record IDs for which taxonomy was found
    :param prefix: Prefix for intermediate and result files
    :return: Path to curated FASTA file without repeats or phages
    """
    output_dir = output / Path("intermediate_files")

    if not output_dir.exists():
        output_dir.mkdir()

    output_name = "curated_records.faa"
    if prefix:
        output_name = f"{prefix}_{output_name}"
    output_path = output_dir / Path(output_name)

    SeqIO.write((record for record in records if record.id in taxonomy_ids), Path(output_path), "fasta")

    return output_path
