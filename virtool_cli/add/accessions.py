from pathlib import Path

import click
from structlog import get_logger

from virtool_cli.add.format import format_record
from virtool_cli.add.helpers import (
    get_no_fetch_lists,
    is_addable,
    search_otu_path,
    write_sequences_to_src,
)
from virtool_cli.check.checkup import verify_accession
from virtool_cli.utils.cache import generate_taxid_table
from virtool_cli.utils.ncbi import request_from_nucleotide
from virtool_cli.utils.reference import is_v1


async def add_accession(accession: str, path: Path):
    """Takes a specified accession, fetches the corresponding record from NCBI Nucleotide.
    Finds a matching OTU directory in the reference and writes the new accession data
    under the existing OTU directory.

    :param accession: NCBI Taxonomy accession to be added to the reference
    :param path: the path to a reference repository
    """
    logger = get_logger(accession=accession)

    src_path = path / "src"

    if is_v1(src_path):
        click.echo(
            "Repository is a deprecated v1 reference.",
            err=True,
        )
        return

    if not verify_accession(accession):
        click.echo("Invalid accession", err=True)
        return

    record_list = await request_from_nucleotide([accession])
    record_list = record_list.pop()

    otu_path = await search_otu_path(
        record_list,
        src_path,
        taxid_table=generate_taxid_table(src_path),
        logger=logger,
    )

    if not otu_path:
        click.echo("No matching OTU found.", err=True)

    elif not await is_addable(accession, otu_path, logger=logger):
        click.echo("This accession will not be added.")

    else:
        new_sequence = await format_record(record_list)

        await write_sequences_to_src(
            sequences=[new_sequence],
            otu_path=otu_path,
            src_path=otu_path.parent,
            logger=logger,
        )


async def add_accessions(accessions: list, otu_path: Path):
    """Add a list of accessions to an OTU. Appropriate if you know the OTU path already.

    :param accessions: A list of accession numbers to query from NCBI
    :param otu_path: Path to an OTU directory
    """
    logger = get_logger("add_accessions", accessions=accessions)

    record_list = await request_from_nucleotide(accessions)
    extant_list, exclusion_list = await get_no_fetch_lists(otu_path)

    new_sequences = []

    for record in record_list:
        addable = await is_addable(
            record.id,
            otu_path=otu_path,
            exclusion_list=exclusion_list,
            logger=logger.bind(accession=record.id),
        )

        if not addable:
            click.echo(f"Accession {record.id} will not be added.")
            continue

        new_sequences.append(await format_record(record))

    await write_sequences_to_src(
        sequences=new_sequences,
        otu_path=otu_path,
        src_path=otu_path.parent,
        logger=logger,
    )
