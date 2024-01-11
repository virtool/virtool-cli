from pathlib import Path
import asyncio
import structlog

from virtool_cli.utils.logging import DEFAULT_LOGGER, DEBUG_LOGGER
from virtool_cli.utils.reference import is_v1, get_otu_paths, get_unique_ids
from virtool_cli.utils.ncbi import request_from_nucleotide, fetch_isolate_metadata
from virtool_cli.utils.format import format_sequence, get_qualifiers, check_source_type
from virtool_cli.utils.storage import write_records, get_otu_accessions

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
    structlog.configure(wrapper_class=DEBUG_LOGGER if debugging else DEFAULT_LOGGER)
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

    otu_accession_list = await get_otu_accessions(otu_path)

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
