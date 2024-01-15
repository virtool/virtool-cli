from pathlib import Path
import asyncio
import structlog

from virtool_cli.utils.logging import DEFAULT_LOGGER, DEBUG_LOGGER
from virtool_cli.utils.reference import get_otu_paths, is_v1, get_unique_ids
from virtool_cli.utils.ncbi import request_from_nucleotide, fetch_isolate_metadata
from virtool_cli.utils.format import format_sequence, get_qualifiers, check_source_type
from virtool_cli.utils.storage import write_records, get_otu_accessions, fetch_exclusions
from virtool_cli.utils.cache import generate_taxid_table
from virtool_cli.add.helpers import is_accession_extant, find_taxon_id

base_logger = structlog.get_logger()


def run(
    accession: str,
    src_path: Path,
    debugging: bool = False,
):
    """
    CLI entry point for virtool_cli.add.accession.add_accession()

    :param accession: NCBI Taxonomy accession to be added to the reference
    :param src_path: Path to a reference directory
    :param debugging: Enables verbose logs for debugging purposes
    """
    structlog.configure(wrapper_class=DEBUG_LOGGER if debugging else DEFAULT_LOGGER)
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

    taxid_table = generate_taxid_table(src_path)

    record_list = await request_from_nucleotide([accession])
    seq_data = record_list.pop()

    otu_path = await get_otu_path(seq_data, src_path, taxid_table, logger=logger)
    if not otu_path:
        logger.error("No matching OTU found.")
        return

    try:
        extant_list = await get_otu_accessions(otu_path)
    except Exception as e:
        logger.exception(e)
        raise e

    try:
        accession_collision = await is_accession_extant(accession, extant_list)
        if accession_collision:
            logger.warning(
                "This accession already exists in the reference.",
                accession=accession,
            )
            return
    except Exception as e:
        logger.exception(e)
        raise e

    logger.debug("Accession not found in exclusion list")

    seq_qualifiers = get_qualifiers(seq_data.features)

    if isolate_type := check_source_type(seq_qualifiers):
        # Isolate metadata contained in qualifiers
        isolate = {
            "source_name": seq_qualifiers.get(isolate_type)[0],
            "source_type": isolate_type,
        }

    else:
        # Extract isolate metadata from NCBI Taxonomy docsum
        isolate = await fetch_isolate_metadata(find_taxon_id(seq_qualifiers["db_xref"]))

    new_sequence = format_sequence(record=seq_data, qualifiers=seq_qualifiers)
    new_sequence["isolate"] = isolate

    logger.debug(new_sequence)

    isolate_uids, sequence_uids = await get_unique_ids(get_otu_paths(src_path))

    try:
        new_sequence_paths = await write_records(
            otu_path,
            new_sequences=[new_sequence],
            unique_iso=isolate_uids,
            unique_seq=sequence_uids,
            logger=logger,
        )
    except Exception as e:
        logger.exception(e)
        new_sequence_paths = []

    for path in new_sequence_paths:
        print(str(path))


async def get_otu_path(
    seq_data, src_path: Path, taxid_table: dict, logger = base_logger
) -> Path | None:
    """
    Find an OTU directory in the reference using metadata
    from NCBI Nucleotide sequence records (such as Taxonomy UID and name).

    :param seq_data: Sequence data retrieved from NCBI Nucleotide and parsed as SeqRecord
    :param src_path: Path to a reference directory
    :param taxid_table:
    :param logger: Optional entry point for an existing BoundLogger
    """
    # Get taxon id and OTU id
    seq_qualifiers = get_qualifiers(seq_data.features)
    logger.debug(seq_qualifiers)

    # Generate a potential name for the directory and match it against OTU
    if taxid := find_taxon_id(seq_qualifiers["db_xref"]):
        logger.info("Taxon ID found", taxid=taxid)

        # Search for a matching OTU path using the taxid
        if taxid in taxid_table:
            otu_path = src_path / taxid_table[taxid]

            logger.info("OTU path found", otu_path=otu_path.name)

            return otu_path
        else:
            logger.warning("No preexisting OTU. Run add otu first.")

    else:
        logger.error("No taxon id found in metadata.")

    logger.error("No matching OTU found in src directory.")
    return None
