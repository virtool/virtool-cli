import json
from pathlib import Path
import asyncio

import structlog

from virtool_cli.utils.logging import DEFAULT_LOGGER, DEBUG_LOGGER
from virtool_cli.utils.reference import (
    get_otu_paths,
    is_v1,
    search_otu_by_id,
    generate_otu_dirname,
)
from virtool_cli.utils.ncbi import request_from_nucleotide
from virtool_cli.utils.id_generator import get_unique_ids
from virtool_cli.utils.format import format_sequence, get_qualifiers, check_source_type
from virtool_cli.utils.storage import write_records
from virtool_cli.catalog.helpers import get_otu_accessions


base_logger = structlog.get_logger()


def run_single(
    accession: str,
    src_path: Path,
    catalog_path: Path,
    debugging: bool = False,
):
    """
    CLI entry point for virtool_cli.add.accession.add_accession()

    :param accession: NCBI Taxonomy accession to be added to the reference
    :param src_path: Path to a reference directory
    :param catalog_path: Path to a catalog directory
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
        add_accession(accession=accession, src_path=src_path, catalog_path=catalog_path)
    )


def run_multiple(accessions: str, otu_path: Path, debugging: bool = False):
    """
    CLI entry point for virtool_cli.add.accession.add_accessions()

    :param accessions: NCBI Taxonomy accessions to be added to the reference
    :param otu_path: Path to a reference directory
    :param debugging: Enables verbose logs for debugging purposes
    """
    structlog.configure(wrapper_class=DEBUG_LOGGER if debugging else DEFAULT_LOGGER)
    logger = base_logger.bind(src=str(otu_path))

    accession_list = accessions.split(",")

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

    otu_accession_list = get_otu_accessions(otu_path)

    record_list = await request_from_nucleotide(accessions)

    new_sequences = []

    for record in record_list:
        accession = record.id
        logger = logger.bind(accession=accession)

        accession_collision = await check_accession_collision(
            accession, otu_accession_list, logger
        )
        if accession_collision:
            logger.warning(
                f"{accession} already in OTU, moving on...", accession=accession
            )
            continue

        seq_qualifiers = get_qualifiers(record.features)

        isolate_type = check_source_type(seq_qualifiers)
        if not isolate_type:
            return
        isolate_name = seq_qualifiers.get(isolate_type)[0]
        isolate = {"source_name": isolate_name, "source_type": isolate_type}

        new_sequence = format_sequence(record=record, qualifiers=seq_qualifiers)
        new_sequence["isolate"] = isolate

        new_sequences.append(new_sequence)

    isolate_uids, sequence_uids = await get_unique_ids(get_otu_paths(otu_path.parent))

    try:
        await write_records(
            otu_path, new_sequences, isolate_uids, sequence_uids, logger=logger
        )
        logger.info(
            f'OTU written to {otu_path.name}. Use "virtool ref add accessions" to add accessions.'
        )
    except Exception as e:
        logger.exception(e)


async def add_accession(accession: str, src_path: Path, catalog_path: Path):
    """
    Takes a specified accession, fetches the corresponding record from NCBI Taxonomy
    and writes it to the reference directory.

    :param accession: NCBI Taxonomy accession to be added to the reference
    :param src_path: Path to a reference directory
    :param catalog_path: Path to an accession catalog directory
    """
    logger = base_logger.bind(accession=accession)

    record_list = await request_from_nucleotide([accession])
    seq_data = record_list.pop()

    otu_path = await get_otu_path(seq_data, src_path, catalog_path, logger)

    otu_accession_list = get_otu_accessions(otu_path)
    accession_collision = await check_accession_collision(
        accession, otu_accession_list, logger
    )
    if accession_collision:
        logger.warning(f"{accession} already in OTU, moving on...", accession=accession)
        return

    seq_qualifiers = get_qualifiers(seq_data.features)

    isolate_type = check_source_type(seq_qualifiers)
    if not isolate_type:
        return
    isolate_name = seq_qualifiers.get(isolate_type)[0]
    isolate = {"source_name": isolate_name, "source_type": isolate_type}

    new_sequence = format_sequence(record=seq_data, qualifiers=seq_qualifiers)
    new_sequence["isolate"] = isolate

    logger.debug(new_sequence)

    isolate_uids, sequence_uids = await get_unique_ids(get_otu_paths(src_path))

    try:
        await write_records(
            otu_path, [new_sequence], isolate_uids, sequence_uids, logger=logger
        )
    except Exception as e:
        logger.exception(e)


async def get_otu_path(seq_data, src_path: Path, catalog_path: Path, logger):
    """
    Find a OTU directory in the reference using metadata
    from NCBI Nucleotide sequence records (such as Taxonomy UID and name).

    :param seq_data: Sequence data retrieved from NCBI Nucleotide and parsed as SeqRecord
    :param src_path: Path to a reference directory
    :param catalog_path: Path to an accession catalog directory
    :param logger: Optional entry point for an existing BoundLogger
    """
    # Get taxon id and OTU id
    seq_qualifiers = get_qualifiers(seq_data.features)
    logger.debug(seq_qualifiers)

    taxid = find_taxon_id(seq_qualifiers["db_xref"])
    if taxid < 0:
        logger.error("No taxon id found!")

    taxid_matches = list(catalog_path.glob(f"{taxid}--*.json"))
    if taxid_matches:
        # Get the OTU id directly from the listing
        listing_path = taxid_matches.pop()

        logger.info(
            "Found matching listing in catalog.",
            listing=listing_path.name,
        )

        matching_listing = json.loads(listing_path.read_text())

        otu_id = matching_listing["_id"]

        if otu_path := search_otu_by_id(otu_id, src_path):
            logger.debug(
                "Matching OTU directory found",
                listing=listing_path.name,
                otu=otu_path.name,
            )

            return otu_path

        logger.error("No matching OTU found in src directory.")

        return None

    else:
        # Generate a potential name for the directory and match it against OTU
        dummy_name = generate_otu_dirname(name=seq_data["ScientificName"], otu_id="")

        otu_matches = list(src_path.glob(f"{dummy_name}*"))
        if otu_matches:
            otu_path = otu_matches.pop()

            logger.debug("Unlisted matching OTU found", otu=otu_path.name)

            return otu_path

        else:
            logger.error("No matching OTU found in src directory.")
            return None


async def check_accession_collision(
    new_accession: str, accession_list: list, logger
) -> bool:
    """
    Check if a new accession already exists in a list of already-assessed accessions.

    :param new_accession: A new accession
    :param accession_list: A list of accessions that should not be added anew
    :param logger: Optional entry point for an existing BoundLogger
    :return: True if the accession collides with the accession list, False if not
    """
    for existing_accession in accession_list:
        if new_accession.split(".")[0] == existing_accession.split(".")[0]:
            logger.error(
                "This accession already exists in the reference. Consider editing the existing sequence."
            )
            return True

    return False


def find_taxon_id(db_xref: list[str]) -> int:
    """
    Searches the database cross-reference data for the associated taxon ID.

    :param db_xref: List of NCBI cross-reference information taken from NCBI taxonomy record
    :return: NCBI taxon ID as an integer if found, -1 if not found.
    """
    print(f"xrefs: {db_xref}")

    for xref in db_xref:
        [key, value] = xref.split(":")
        if key == "taxon":
            return int(value)

    return -1
