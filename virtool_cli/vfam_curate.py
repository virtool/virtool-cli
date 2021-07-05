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
    :return: records found in input paths
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
    Iterates through records, yields record if sequence isn't a duplicate and is longer than sequence_min_length.

    :param records: iterable of records gathered in group_input_paths()
    :param sequence_min_length: minimum length of sequence to be included in output
    :return: records with sequences longer than sequence_min_length, without duplicates
    """
    record_seqs = list()

    for record in records:
        if record.seq not in record_seqs and len(record.seq) > sequence_min_length:
            record_seqs.append(record.seq)
            yield record


def write_no_dupes(no_dupes: iter, output: Path, prefix: str) -> Path:
    """
    Writes records from no_dupes() generator function to no_duplicate_records.faa.

    :param no_dupes: iterable of records gathered in remove_dupes()
    :param output: path to output directory
    :param prefix:  Prefix for intermediate and result files
    :return: path to no_duplicate_records.faa in output directory
    """
    output_dir = output / Path("intermediate_files")

    if not output_dir.exists():
        output_dir.mkdir()

    output_name = "no_duplicate_records.faa"

    if prefix:
        output_name = f"{prefix}_{output_name}"

    output_path = output_dir / Path(output_name)

    SeqIO.write(no_dupes, Path(output_path), "fasta")

    return output_path


def get_taxonomy(no_dupes_path: Path) -> dict:
    """
    Iterates through records in no_dupes_path, Makes calls to NCBI database using Bio.Entrez.

    Gathers taxonomic information for each record, if record isn't found in NCBI, record is filtered out.

    :param no_dupes_path: path to file containing all records from write_no_dupes()
    :return: record_taxonomy, a dictionary containing {record ID: [family, genus]} key-value pairs
    """
    record_taxonomy = dict()

    with no_dupes_path.open("r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            try:
                handle = Entrez.efetch(db="protein", id=record.id, rettype="gb", retmode="text")

                for seq_record in SeqIO.parse(handle, "genbank"):
                    record_taxonomy[seq_record.id] = seq_record.annotations["taxonomy"][-2:]

            except (HTTPError, AttributeError):
                continue

    return record_taxonomy


def write_curated_recs(no_dupes_path, taxonomy_ids: List[str], prefix=None) -> Path:
    """
    Iterates through records in no_dupes_path, writes records to output if taxonomy was found in get_taxonomy().

    :param no_dupes_path: path to file containing all records from write_no_dupes()
    :param taxonomy_ids: list of record IDs for which taxonomy was found
    :param prefix: Prefix for intermediate and result files
    :return: Path to curated FASTA file without repeats or phages
    """
    output_name = "curated_records.faa"

    if prefix:
        output_name = f"{prefix}_{output_name}"

    output_path = no_dupes_path.parent / Path(output_name)

    with no_dupes_path.open("r") as handle:
        SeqIO.write(
            (record for record in SeqIO.parse(handle, "fasta") if record.id in taxonomy_ids),
            Path(output_path),
            "fasta"
        )

    return output_path
