import json
from pathlib import Path
import asyncio

import structlog
from structlog import get_logger, BoundLogger

# from typing import Optional, Tuple
# from urllib.error import HTTPError

from virtool_cli.utils.logging import DEFAULT_LOGGER, DEBUG_LOGGER
from virtool_cli.utils.reference import is_v1, get_otu_paths, get_isolate_paths
from virtool_cli.utils.ncbi import request_from_nucleotide, NCBI_REQUEST_INTERVAL

# from virtool_cli.utils.id_generator import generate_unique_ids, get_unique_ids
# from virtool_cli.utils.storage import store_isolate, store_sequence
from virtool_cli.utils.format import format_sequence, get_qualifiers, find_isolate

DEFAULT_INTERVAL = 0.001


def run(
    accession: str,
    src_path: Path,
    catalog_path: Path,
    debugging: bool = False,
):
    """
    CLI entry point for update.update_ref.run()

    Requests updates for all OTU directories under a source reference

    :param accession:
    :param src_path: Path to a reference directory
    :param catalog_path: Path to a catalog directory
    :param debugging: Enables verbose logs for debugging purposes
    """
    structlog.configure(wrapper_class=DEBUG_LOGGER if debugging else DEFAULT_LOGGER)
    logger = get_logger()

    if is_v1(src_path):
        logger.critical(
            'reference folder "src" is a deprecated v1 reference.'
            + 'Run "virtool ref migrate" before trying again.'
        )
        return

    logger.debug("Debug flag is enabled")

    asyncio.run(add_accession(accession=accession, src_path=src_path))


async def add_accession(accession: str, src_path: Path) -> dict:
    logger = get_logger().bind(accession=accession)

    record_list = await request_from_nucleotide([accession])
    seq_data = record_list.pop()

    seq_qualifiers = get_qualifiers(seq_data.features)
    logger.debug(seq_qualifiers)

    taxid = find_taxon_id(seq_qualifiers["db_xref"])
    if taxid < 0:
        logger.error("No taxon ID found!")

    logger.debug(f"taxon_id: {taxid}")

    if not (isolate_type := find_isolate(seq_qualifiers)):
        return {}

    isolate_name = seq_qualifiers.get(isolate_type)[0]

    logger.debug(f"isolate={isolate_name}")

    new_sequence = format_sequence(
        record=seq_data, qualifiers=seq_qualifiers, logger=logger
    )

    logger.debug(new_sequence)

    pass


def find_taxon_id(db_xref: list[str]) -> int:
    """
    Searches the database cross-reference data for the associated taxon ID.
    :param db_xref:
    :return: NCBI taxon ID as an integer if found, -1 if not found.
    """
    print(f"xrefs: {db_xref}")

    for xref in db_xref:
        [key, value] = xref.split(":")
        if key == "taxon":
            return int(value)

    return -1


if __name__ == "__main__":
    debug = True

    REPO_DIR = "/Users/sygao/Development/UVic/Virtool/Repositories"

    project_path = Path(REPO_DIR) / "ref-plant-viruses"
    src_path = project_path / "src"
    catalog_path = Path(REPO_DIR) / "ref-accession-catalog/catalog"

    ACCESSION = "NC_017829"

    asyncio.run(add_accession(accession=ACCESSION, src_path=src_path))
