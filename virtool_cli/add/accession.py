from pathlib import Path
import asyncio

import structlog

from virtool_cli.utils.logging import DEFAULT_LOGGER, DEBUG_LOGGER
from virtool_cli.utils.reference import (
    get_otu_paths,
    is_v1,
    generate_otu_dirname,
    search_otu_by_id,
    get_unique_ids,
)
from virtool_cli.utils.ncbi import request_from_nucleotide, fetch_isolate_metadata
from virtool_cli.utils.format import format_sequence, get_qualifiers, check_source_type
from virtool_cli.utils.storage import write_records
from virtool_cli.catalog.listings import parse_listing
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


def run_multiple(
    accessions_string: str, otu_path: Path, catalog_path: Path, debugging: bool = False
):
    """
    CLI entry point for virtool_cli.add.accession.add_accessions()

    :param accessions_string: NCBI Taxonomy accessions to be added to the reference
    :param otu_path: Path to a reference directory
    :param catalog_path: Path to a catalog directory
    :param debugging: Enables verbose logs for debugging purposes
    """
    structlog.configure(wrapper_class=DEBUG_LOGGER if debugging else DEFAULT_LOGGER)
    logger = base_logger.bind(src=str(otu_path))

    accession_list = split_clean_csv_string(accessions_string, delimiter=",")

    logger.debug("Debug flag is enabled")

    if otu_path.exists():
        asyncio.run(add_accessions(accession_list, otu_path))

    else:
        logger.error("Given OTU path does not exist", otu=otu_path)

    logger.debug(
        "Future placeholder for listing update run", catalog_path=str(catalog_path)
    )


async def add_accessions(accessions: list, otu_path: Path):
    """
    Add a list of accessions to an OTU. Appropriate if you know the OTU path already.

    :param accessions: A list of accession numbers to query from NCBI
    :param otu_path: Path to an OTU directory
    """
    logger = base_logger.bind(accessions=accessions)

    otu_accession_list = get_otu_accessions(otu_path)

    logger.debug(otu_accession_list)

    logger.debug(accessions)

    record_list = await request_from_nucleotide(accessions)

    new_sequences = []

    for record in record_list:
        accession = record.id
        logger = logger.bind(accession=accession)

        accession_collision = await check_accession_collision(
            accession, otu_accession_list
        )
        if accession_collision:
            logger.warning(
                f"{accession} already in OTU, moving on...", accession=accession
            )
            continue

        seq_qualifiers = get_qualifiers(record.features)

        if isolate_type := check_source_type(seq_qualifiers):
            # Isolate metadata contained in qualifiers
            isolate = {
                "source_name": seq_qualifiers.get(isolate_type)[0],
                "source_type": isolate_type,
            }

        else:
            # Extract isolate metadata from NCBI Taxonomy docsum
            isolate = await fetch_isolate_metadata(
                find_taxon_id(seq_qualifiers["db_xref"])
            )

        new_sequence = format_sequence(record=record, qualifiers=seq_qualifiers)
        new_sequence["isolate"] = isolate

        new_sequences.append(new_sequence)

    isolate_uids, sequence_uids = await get_unique_ids(get_otu_paths(otu_path.parent))

    try:
        new_sequence_paths = await write_records(
            otu_path, new_sequences, isolate_uids, sequence_uids, logger=logger
        )
        logger.info(f"Accessions written to {otu_path.name}.")

        for path in new_sequence_paths:
            print(str(path))

    except Exception as e:
        logger.exception(e)


async def add_accession(accession: str, src_path: Path, catalog_path: Path):
    """
    Takes a specified accession, fetches the corresponding record from NCBI Nucleotide.
    Finds a matching OTU directory in the reference and writes the new accession data
    under the existing OTU directory.

    :param accession: NCBI Taxonomy accession to be added to the reference
    :param src_path: Path to a reference directory
    :param catalog_path: Path to an accession catalog directory
    """
    logger = base_logger.bind(accession=accession)

    record_list = await request_from_nucleotide([accession])
    seq_data = record_list.pop()

    otu_path = await get_otu_path(seq_data, src_path, catalog_path, logger)
    if not otu_path:
        logger.error("No matching OTU found.")
        return

    otu_accession_list = get_otu_accessions(otu_path)
    accession_collision = await check_accession_collision(accession, otu_accession_list)
    if accession_collision:
        logger.warning(
            "This accession already exists in the reference. Consider editing the existing sequence.",
            accession=accession,
        )
        return

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
    seq_data, src_path: Path, catalog_path: Path, logger
) -> Path | None:
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

    # Generate a potential name for the directory and match it against OTU
    if "organism" not in seq_qualifiers:
        logger.error("Could not find taxon identifier")
        raise ValueError
    otu_name = seq_qualifiers["organism"][0]
    dummy_name = generate_otu_dirname(name=otu_name, otu_id="")

    otu_matches = list(src_path.glob(f"{dummy_name}*"))
    if otu_matches:
        otu_path = otu_matches.pop()

        logger.debug("Matching OTU found", otu_path=str(otu_path))

        return otu_path

    logger.warning(
        "OTU could not be found by name. Your OTU naming scheme may not be up to date."
    )

    logger.info(
        "Attempting catalog match...",
        catalog_path=str(catalog_path),
    )
    if taxid := find_taxon_id(seq_qualifiers["db_xref"]):
        logger.info("Taxon ID found", taxid=taxid)

        # Search for a matching OTU path using the taxid
        if otu_path := await get_otu_path_from_taxid(taxid, src_path, catalog_path):
            logger.info("OTU path found", otu_path=otu_path.name)
            return otu_path

    else:
        logger.error("No taxon id found in metadata.")

    logger.error("No matching OTU found in src directory.")
    return None


async def get_otu_path_from_taxid(
    taxid: int, src_path: Path, catalog_path: Path
) -> Path | None:
    """
    Use NCBI Taxonomy UID to retrieve a matching OTU ID from catalog if possible.

    :param taxid: NCBI Taxonomy UID
    :param src_path: Path to a reference directory
    :param catalog_path: Path to an accession catalog directory
    :return: Matching OTU directory path if found, None if not found
    """
    # Attempt a direct catalog match
    taxid_matches = list(catalog_path.glob(f"{taxid}--*.json"))
    if not taxid_matches:
        return None

    # Get the OTU id directly from the listing
    listing_path = taxid_matches.pop()

    matching_listing = await parse_listing(listing_path)

    otu_id = matching_listing["_id"]

    if otu_path := search_otu_by_id(otu_id, src_path):
        return otu_path

    return None


async def check_accession_collision(new_accession: str, accession_list: list) -> bool:
    """
    Check if a new accession already exists in a list of already-assessed accessions.

    :param new_accession: A new accession
    :param accession_list: A list of accessions that should not be added anew
    :return: True if the accession collides with the accession list, False if not
    """
    return any(
        new_accession.split(".")[0] == existing_accession.split(".")[0]
        for existing_accession in accession_list
    )


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
