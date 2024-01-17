from pathlib import Path
import asyncio
import structlog

from virtool_cli.utils.logging import configure_logger
from virtool_cli.utils.reference import is_v1
from virtool_cli.utils.ncbi import request_from_nucleotide
from virtool_cli.utils.cache import generate_taxid_table
from virtool_cli.add.helpers import is_addable, search_otu_path
from virtool_cli.add.format import format_record
from virtool_cli.add.write import write_sequences_to_src

base_logger = structlog.get_logger()


def run(
    accession: str,
    src_path: Path,
    debugging: bool = False,
):
    """
    CLI entry point for virtool_cli.add.accession.add_accession_auto()

    :param accession: NCBI Taxonomy accession to be added to the reference
    :param src_path: Path to a reference directory
    :param debugging: Enables verbose logs for debugging purposes
    """
    configure_logger(debugging)

    logger = base_logger.bind(src=str(src_path))

    if is_v1(src_path):
        logger.critical(
            'Reference folder "src" is a deprecated v1 reference.'
            + 'Run "virtool ref migrate" before trying again.'
        )
        return

    logger.debug("Debug flag is enabled")

    asyncio.run(
        add_accession(accession=accession, src_path=src_path)
    )


async def add_accession(accession: str, src_path: Path):
    """
    Takes a specified accession, fetches the corresponding record from NCBI Nucleotide.
    Finds a matching OTU directory in the reference and writes the new accession data
    under the existing OTU directory.

    :param accession: NCBI Taxonomy accession to be added to the reference
    :param src_path: Path to a reference directory
    """
    logger = base_logger.bind(accession=accession)

    record_list = await request_from_nucleotide([accession])
    seq_data = record_list.pop()

    otu_path = await search_otu_path(
        seq_data,
        src_path,
        taxid_table=generate_taxid_table(src_path),
        logger=logger
    )
    if not otu_path:
        logger.error("No matching OTU found.")
        return

    addable = await is_addable(accession, otu_path, logger=logger)
    if not addable:
        logger.warning("This accession will not be added.")
        return

    new_sequence = await format_record(seq_data, logger)

    new_sequence_paths = await write_sequences_to_src(
        sequences=[new_sequence],
        otu_path=otu_path,
        src_path=otu_path.parent,
        logger=logger,
    )

    for path in new_sequence_paths:
        print(str(path))


def find_taxon_id(db_xref: list[str]) -> int | None:
    """
    Searches the database cross-reference data for the associated NCBI taxonomy UID.

    :param db_xref: List of NCBI cross-reference information taken from NCBI taxonomy record
    :return: NCBI Taxonomy UID as an integer if found, None if not found
    """
    for xref in db_xref:
        [key, value] = xref.split(":")
        if key == "taxon":
            return int(value)

    return None
