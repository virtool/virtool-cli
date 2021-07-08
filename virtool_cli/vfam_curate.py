import gzip
import subprocess
import sys
from abc import ABC
from html.parser import HTMLParser
from urllib.error import HTTPError

from Bio import SeqIO
from subprocess import SubprocessError
from pathlib import Path
from typing import List
from virtool_cli.vfam_console import console


def get_genbank_files(output: Path) -> List[Path]:
    """
    Gathers .html file from https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/ using wget.

    Parses .html file to find .gpff filenames. Gathers .gpff files using wget.

    :param output: path to output directory
    :return: genbank_file_paths, a list of paths to the .gpff files in project directory
    """
    output_path = Path(output / "genbank_input.html")
    viral_release_url = "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/"

    wget_cmd = [
        "wget", viral_release_url,
        "-O", output_path,
    ]

    try:
        subprocess.run(wget_cmd)
    except (SubprocessError, HTTPError):
        console.print(f"Error gathering .html file from {viral_release_url}.", style="red")
        sys.exit(1)

    parser = ViralProteinParser()

    with output_path.open("r") as handle:
        parser.feed(handle.read())
    file_names = parser.close()

    genbank_file_paths = list()
    for file_name in file_names:
        output_path = output / file_name
        wget_cmd = [
            "wget", f"{viral_release_url}/{file_name}",
            "-O", output_path
        ]

        genbank_file_paths.append(output_path)

        try:
            subprocess.run(wget_cmd)
        except (SubprocessError, HTTPError):
            genbank_file_paths.remove(file_name)
            console.print(f"✘ Error gathering {file_name} from {viral_release_url}.", style="red")
            console.print(f"Record data from {file_name} will not be included in output.", style="red")
            continue

    console.print(f"✔ Retrieved {len(genbank_file_paths)} .gpff files from {viral_release_url}.", style="green")

    return genbank_file_paths


def get_input_paths(src_path: Path) -> List[Path]:
    """
    Takes in a path to directory containing input FASTA files, returns a list of paths to the FASTA files.

    :param src_path: Path to input source directory containing unfiltered FASTA files
    :return: input_paths, list of paths to input files if any files are found
    """
    input_paths = list(src_path.iterdir())

    if input_paths:
        console.print(f"✔ Retrieved {len(input_paths)} files from input directory.", style="green")
        return input_paths

    console.print("No files found in input directory.", style="red")
    sys.exit(1)


def group_input_paths(input_paths: List[Path], no_named_phages: bool) -> list:
    """
    Takes in paths to genbank files as input and yields records.

    Filters out records with "phage" in their description if no_named_phages is True.

    :param input_paths: list of paths to input FASTA files
    :param no_named_phages: bool that dictates whether phage records are filtered out by name or not
    :return: records found in input paths
    """
    phage_count = 0
    record_count = 0

    for input_path in input_paths:

        if str(input_path).endswith(".gz"):
            handle = gzip.open(input_path, "rt")
        else:
            handle = open(input_path, "r")

        for record in SeqIO.parse(handle, "genbank"):
            record_count += 1
            if no_named_phages:
                if "phage" in record.description:
                    phage_count += 1
                else:
                    yield record
            else:
                yield record

        handle.close()

    console.print(f"✔ Retrieved {record_count} records from {len(input_paths)} input files.", style="green")

    if no_named_phages:
        console.print(f"✔ Filtered out {phage_count} phage records by name.", style="green")


def remove_dupes(records: iter, sequence_min_length: int) -> list:
    """
    Iterates through records, yields record if sequence isn't a duplicate and is longer than sequence_min_length.

    :param records: iterable of records gathered in group_input_paths()
    :param sequence_min_length: minimum length of sequence to be included in output
    :return: records with sequences longer than sequence_min_length, without duplicates
    """
    record_seqs = list()
    dupes_count = 0

    for record in records:
        if record.seq not in record_seqs and len(record.seq) > sequence_min_length:
            record_seqs.append(record.seq)
            yield record
        else:
            dupes_count += 1

    console.print(f"✔ Filtered out {dupes_count} duplicate records.", style="green")


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

    output_name = "no_duplicate_records.gpff"

    if prefix:
        output_name = f"{prefix}_{output_name}"

    output_path = output_dir / Path(output_name)

    SeqIO.write(no_dupes, Path(output_path), "genbank")

    return output_path


def get_taxonomy(no_dupes_path: Path) -> dict:
    """
    Iterates through genbank records in no_dupes_path using SeqIO.

    Gathers taxonomic information for each record, stores in record_taxonomy dictionary.

    :param no_dupes_path: path to file containing all records from write_no_dupes()
    :return: record_taxonomy, a dictionary containing {record ID: [family, genus]} key-value pairs
    """
    record_taxonomy = dict()

    with no_dupes_path.open("r") as handle:
        for record in SeqIO.parse(handle, "genbank"):

            record_taxonomy[record.id] = record.annotations["taxonomy"][-2:]

    return record_taxonomy


def write_curated_recs(no_dupes_path, prefix=None) -> Path:
    """
    Iterates through records in no_dupes_path, writes records to output if taxonomy was found in get_taxonomy().

    :param no_dupes_path: path to file containing all records from write_no_dupes()
    :param prefix: Prefix for intermediate and result files
    :return: Path to curated FASTA file without repeats or phages
    """
    output_name = "curated_records.faa"

    if prefix:
        output_name = f"{prefix}_{output_name}"

    output_path = no_dupes_path.parent / Path(output_name)

    with no_dupes_path.open("r") as handle:
        SeqIO.write(
            (record for record in SeqIO.parse(handle, "genbank")),
            Path(output_path),
            "fasta"
        )

    return output_path


class ViralProteinParser(HTMLParser, ABC):
    """Parser used to gather viral protein data"""
    file_names = list()

    def handle_data(self, data):
        if data.startswith("viral") and ".gpff" in data:
            self.file_names.append(data)

    def close(self):
        return self.file_names
