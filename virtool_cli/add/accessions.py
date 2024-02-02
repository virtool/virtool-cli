from pathlib import Path

import click
from structlog import get_logger

from virtool_cli.add.format import format_record
from virtool_cli.add.helpers import (
    find_taxon_id,
    get_no_fetch_lists,
    is_addable,
    write_sequences_to_src,
)
from virtool_cli.check.otu import verify_accession
from virtool_cli.repo.cls import Repo, RepoSequence
from virtool_cli.utils.format import get_qualifiers
from virtool_cli.utils.ncbi import request_from_nucleotide


async def add_accession(accession: str, path: Path) -> RepoSequence | None:
    """Attempts to add the NCBI record associated with the ``accession`` to the
    reference.

    :param accession: the NCBI accession to be added to the reference
    :param path: the path to a reference repository
    """
    logger = get_logger(accession=accession)

    repo = Repo(path)

    logger.info("Adding accession")

    if not verify_accession(accession):
        click.echo("Invalid accession", err=True)
        return None

    record_list = await request_from_nucleotide([accession])
    record = record_list.pop()

    taxid = find_taxon_id(get_qualifiers(record.features)["db_xref"])

    if not taxid:
        click.echo("No taxon id found in record.", err=True)
        return None

    try:
        otu = repo.get_otu_by_taxid(taxid)
    except ValueError:
        logger.warning(
            "No matching OTU found for taxonomy ID.",
            accession=accession,
            taxid=taxid,
        )
        return None

    click.echo(f"Found matching OTU '{otu.name}' for accession '{accession}'.")

    if accession in otu.exclusions:
        logger.warning(
            "This accession has been previously excluded.",
            otu_name=otu.name,
            accession=accession,
        )
        return None

    sequence = await format_record(record)

    source_type = sequence["isolate"]["source_type"]
    source_name = sequence["isolate"]["source_name"]

    isolate = otu.get_isolate_by_name(source_type, source_name)

    if isolate is None:
        isolate = otu.add_isolate(source_type, source_name)
        click.echo(f"Created new isolate '{isolate.name}' with ID '{isolate.id}'.")
    else:
        click.echo(f"Found matching isolate '{isolate.name}' with ID '{isolate.id}'.")

    sequence = isolate.add_sequence(
        sequence["accession"],
        sequence["definition"],
        sequence["host"],
        sequence.get("segment", ""),
        sequence["sequence"],
    )

    if sequence:
        click.echo(
            f"Added sequence '{sequence.accession}' to isolate '{isolate.name}'.",
        )

    return sequence


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
