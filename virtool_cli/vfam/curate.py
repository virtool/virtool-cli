import sys
from pathlib import Path
from typing import List
import urllib.request
from urllib.error import HTTPError, URLError
from html.parser import HTMLParser
import gzip
from abc import ABC
from Bio import SeqIO
from structlog import get_logger

base_logger = get_logger()


VIRAL_RELEASE_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/"


def get_genbank_files(output: Path) -> List[Path]:
    """
    Gathers .html file from https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/ using a urllib request.

    Parses .html file with ViralProteinParser, an HTMLParser subclass to find .gpff filenames.

    Gathers .gpff files using urllib requests.

    :param output: path to output directory
    :return: genbank_file_paths, a list of paths to the .gpff files in project directory gathered from NCBI
    """
    logger = base_logger
    viral_release_url = VIRAL_RELEASE_URL

    if not viral_release_url.startswith("http"):
        logger.error(f"{viral_release_url} is not a valid HTTP URL.")
        sys.exit(1)

    parser = ViralProteinParser()
    try:
        with urllib.request.urlopen(viral_release_url) as html_file:
            parser.feed(html_file.read().decode("utf-8"))
    except (HTTPError, URLError):
        logger.error(f"Error fetching .html file from {viral_release_url}")
        sys.exit(1)

    file_names = parser.close()

    genbank_file_paths = list()

    for file_name in file_names:
        try:
            output_path = output / file_name

            with urllib.request.urlopen(
                f"{viral_release_url}/{file_name}"
            ) as gpff_file:
                with open(output_path, "wb") as f:
                    f.write(gpff_file.read())

                genbank_file_paths.append(output_path)

        except (HTTPError, URLError):
            logger.error(f"Error retrieving {file_name} from {viral_release_url}")
            logger.error(f"Record data from {file_name} will not be included in output")

            continue

    if genbank_file_paths:
        logger.info(
            f"Retrieved {len(genbank_file_paths)} .gpff files from {viral_release_url}."
        )
        return genbank_file_paths

    logger.error(f"Retrieved 0 .gpff files from {viral_release_url}")
    sys.exit(1)


def get_input_paths(src_path: Path) -> List[Path]:
    """
    Takes in a path to directory containing input FASTA files, returns a list of paths to the FASTA files.

    :param src_path: Path to input source directory containing unfiltered FASTA files
    :return: input_paths, list of paths to input files if any files are found
    """
    logger = base_logger

    input_paths = list(src_path.iterdir())

    if input_paths:
        logger.info(f"Retrieved {len(input_paths)} files from input directory.")
        return input_paths

    logger.error("No files found in input directory.")
    sys.exit(1)


def group_input_paths(
    input_paths: List[Path], no_named_phages: bool, sequence_min_length: int
) -> list:
    """
    Takes in paths to genbank files as input and yields records.

    Filters out duplicate records and records with sequences shorter than sequence_min_length.

    Filters out records with "phage" in their description if no_named_phages is True.

    :param input_paths: list of paths to input genbank files
    :param no_named_phages: bool that dictates whether phage records are filtered out by name or not
    :param sequence_min_length: minimum length of sequence to be included in output
    :return: records found in input paths
    """
    logger = base_logger

    record_seqs = list()
    record_count = 0
    dupes_count = 0
    phage_count = 0

    for input_path in input_paths:
        if str(input_path).endswith(".gz"):
            handle = gzip.open(input_path, "rt")
        else:
            handle = open(input_path, "r")

        for record in SeqIO.parse(handle, "genbank"):
            record_count += 1

            if record.seq not in record_seqs and len(record.seq) > sequence_min_length:
                if no_named_phages:
                    if "phage" in record.description:
                        phage_count += 1
                    else:
                        record_seqs.append(record.seq)
                        yield record
                else:
                    record_seqs.append(record.seq)
                    yield record
            else:
                dupes_count += 1

        handle.close()

    logger.info(
        f"Retrieved {record_count} records from {len(input_paths)} input files.",
        count=record_count,
    )

    if no_named_phages:
        logger.info(
            f"Filtered out {phage_count} phage records by name.", count=phage_count
        )

    logger.info(f"Filtered out {dupes_count} duplicate records.", count=dupes_count)


def write_curated_records(curated_records: iter, output: Path, prefix: str) -> Path:
    """
    Writes records from no_dupes() generator function to no_duplicate_records.faa.

    :param curated_records: iterable of records gathered in group_input_paths()
    :param output: path to output directory
    :param prefix:  Prefix for intermediate and result files
    :return: path to curated_records.gpff in output directory
    """
    output_dir = output / Path("intermediate_files")

    if not output_dir.exists():
        output_dir.mkdir()

    output_name = "curated_records.gpff"

    if prefix:
        output_name = f"{prefix}_{output_name}"

    output_path = output_dir / Path(output_name)

    SeqIO.write(curated_records, Path(output_path), "genbank")

    return output_path


def get_taxonomy(curated_records_path: Path) -> dict:
    """
    Iterates through genbank records in curated_records_path using SeqIO.

    Gathers taxonomic information for each record, stores in record_taxonomy dictionary.

    :param curated_records_path: path to genbank file containing all records from write_no_dupes()
    :return: record_taxonomy, a dictionary containing {record ID: [family, genus]} key-value pairs
    """
    record_taxonomy = dict()

    with curated_records_path.open("r") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            record_taxonomy[record.id] = record.annotations["taxonomy"][-2:]

    return record_taxonomy


def genbank_to_fasta(curated_records_path, prefix=None) -> Path:
    """
    Iterates through genbank records in curated_records_path, writes all records to output in FASTA format.

    :param curated_records_path: path to file containing all genbank records from write_curated_records()
    :param prefix: Prefix for intermediate and result files
    :return: Path to curated FASTA file
    """
    output_name = "curated_records.faa"

    if prefix:
        output_name = f"{prefix}_{output_name}"

    output_path = curated_records_path.parent / Path(output_name)

    with curated_records_path.open("r") as handle:
        SeqIO.write(
            (record for record in SeqIO.parse(handle, "genbank")),
            Path(output_path),
            "fasta",
        )

    return output_path


class ViralProteinParser(HTMLParser, ABC):
    """Parser used to gather .gpff file names from NCBI viral release .html file"""

    def __init__(self):
        super().__init__()
        self.file_names = []

    def handle_data(self, data):
        if data.startswith("viral") and ".gpff" in data:
            self.file_names.append(data)

    def close(self):
        return self.file_names
