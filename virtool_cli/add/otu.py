import json
from pathlib import Path
import asyncio
import structlog
from urllib.error import HTTPError

from virtool_cli.utils.logging import DEFAULT_LOGGER, DEBUG_LOGGER
from virtool_cli.utils.reference import is_v1, generate_otu_dirname
from virtool_cli.utils.id_generator import generate_unique_ids, get_unique_otu_ids
from virtool_cli.utils.ncbi import fetch_taxonomy_record

OTU_KEYS = ["_id", "name", "abbreviation", "schema", "taxid"]

base_logger = structlog.get_logger()


def run(
    taxid: int,
    src_path: Path,
    catalog_path: Path,
    debugging: bool = False,
):
    """
    CLI entry point for virtool_cli.add.otu module

    Requests updates for a single OTU directory
    Searches the catalog for a matching catalog listing and requests new updates if it finds one.

    :param taxid:
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

    logger.info(f"Creating OTU directory for NCBI Taxonomy uid {taxid}")

    asyncio.run(add_otu(taxid=taxid, src_path=src_path, catalog_path=catalog_path))


async def add_otu(taxid: int, src_path: Path, catalog_path: Path):
    logger = structlog.get_logger().bind(taxon_id=taxid)

    try:
        taxonomy_data = await fetch_taxonomy_record(taxon_id=taxid)
        if not taxonomy_data:
            logger.error("Could not find a record under this taxon ID on NCBI Taxonomy")
            return

    except HTTPError:
        logger.error(
            "Could not download this data. Please wait a few seconds before trying again"
        )
        return

    taxid_matches = list(catalog_path.glob(f"{taxid}--*.json"))
    if taxid_matches:
        match = taxid_matches.pop()
        logger.error(
            "An OTU with this taxon ID already exists in the reference database.",
            listing=match.name,
        )
        return

    unique_otu_ids = await get_unique_otu_ids(src_path)
    new_id = generate_unique_ids(1, excluded=unique_otu_ids).pop()

    new_otu = generate_otu(taxonomy_data, new_id)
    if not new_otu:
        logger.error("Invalid data from NCBI Taxonomy")

    logger.debug(new_otu)

    await write_otu(new_otu, src_path)

    logger.info("Wrote OTU data")


def generate_otu(taxonomy_data: dict, new_id: str) -> dict:
    """
    Create a new OTU from a given Id from NCBI Taxonomy

    :param taxonomy_data:
    :param new_id:
    :return:
    """
    for key in ["ScientificName", "Id"]:
        if key not in taxonomy_data:
            return {}

    otu = {
        "_id": new_id,
        "name": taxonomy_data.get("ScientificName", ""),
        "abbreviation": "",
        "schema": [],
        "taxid": int(taxonomy_data.get("Id", "-1")),
    }

    return otu


async def write_otu(otu: dict, src_path: Path):
    """
    Generate a new directory for given OTU metadata and
    store the metadata under in otu.json

    :param otu: Dict of OTU metadata
    :param src_path: Path to a reference directory
    """
    dirname = generate_otu_dirname(name=otu["name"], otu_id=otu["_id"])

    otu_path = src_path / dirname
    otu_path.mkdir()

    with open(otu_path / "otu.json", "w") as f:
        json.dump(otu, f, indent=4, sort_keys=True)
