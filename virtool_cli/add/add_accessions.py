from pathlib import Path
import asyncio
import structlog

from virtool_cli.utils.logging import configure_logger
from virtool_cli.utils.ncbi import request_from_nucleotide
from virtool_cli.add.helpers import is_addable, get_no_fetch_lists
from virtool_cli.add.format import format_record
from virtool_cli.add.write import write_sequences_to_src

base_logger = structlog.get_logger()


def run(
    accessions_string: str, otu_path: Path, debugging: bool = False
):
    """
    CLI entry point for virtool_cli.add.accession.add_accessions()

    :param accessions_string: NCBI Taxonomy accessions to be added to the reference
    :param otu_path: Path to a reference directory
    :param debugging: Enables verbose logs for debugging purposes
    """
    configure_logger(debugging)

    logger = base_logger.bind(src=str(otu_path))

    accession_list = split_clean_csv_string(accessions_string, delimiter=",")

    logger.debug("Debug flag is enabled")

    if otu_path.exists():
        asyncio.run(add_accessions(accession_list, otu_path))

    else:
        logger.error("Given OTU path does not exist", otu=otu_path)


async def add_accessions(accessions: list, otu_path: Path):
    """
    Add a list of accessions to an OTU. Appropriate if you know the OTU path already.

    :param accessions: A list of accession numbers to query from NCBI
    :param otu_path: Path to an OTU directory
    """
    logger = base_logger.bind(accessions=accessions)

    record_list = await request_from_nucleotide(accessions)
    extant_list, exclusion_list = await get_no_fetch_lists(otu_path)

    # Process all records and add to new sequence list
    new_sequences = []
    for record in record_list:
        accession = record.id
        logger = logger.bind(accession=accession)

        addable = await is_addable(
            accession,
            otu_path=otu_path,
            exclusion_list=exclusion_list,
            logger=logger
        )
        if not addable:
            logger.warning("This accession will not be added.")
            continue

        new_sequence = await format_record(record, logger)
        new_sequences.append(new_sequence)

    new_sequence_paths = await write_sequences_to_src(
        sequences=new_sequences,
        otu_path=otu_path,
        src_path=otu_path.parent,
        logger=logger,
    )
    for path in new_sequence_paths:
        print(str(path))


def split_clean_csv_string(input_string: str, delimiter: str = ",") -> list[str]:
    """
    Splits comma-separated values into a list of strings.

    :param input_string: A raw comma-delineated string
    :param delimiter: A delimiter string for use with split(). Defaults to a comma (",")
    :return: A list of accession strings
    """
    raw_list = input_string.split(delimiter)
    stripped_list = []
    for item in raw_list:
        stripped_list.append(item.strip())

    return stripped_list
