from pathlib import Path
import structlog
from structlog import BoundLogger, get_logger

from virtool_cli.utils.storage import get_otu_accessions, fetch_exclusions
from virtool_cli.utils.format import get_qualifiers


async def get_no_fetch_lists(otu_path):
    extant_list = await get_otu_accessions(otu_path)
    exclusion_list = await fetch_exclusions(otu_path)

    return extant_list, exclusion_list


async def is_addable(
    accession: str,
    otu_path: Path,
    extant_list: list = None,
    exclusion_list: list = None,
    force: bool = False,
    logger: BoundLogger = structlog.get_logger(),
):
    if extant_list is None:
        extant_list = await get_otu_accessions(otu_path)
    if exclusion_list is None:
        exclusion_list = await fetch_exclusions(otu_path)

    accession_collision = await is_accession_extant(accession, exclusion_list)
    if accession_collision:
        logger.warning(
            "This accession is on the exclusion list.",
            accession=accession,
        )
        if not force:
            return False

    accession_collision = await is_accession_extant(accession, extant_list)
    if accession_collision:
        logger.warning(
            "This accession already exists in the reference.",
            accession=accession,
        )
        return False

    return True


async def is_accession_extant(new_accession: str, excluded_accessions: list) -> bool:
    """
    Check if a new accession already exists in a list of already-assessed accessions.

    :param new_accession: A new accession
    :param excluded_accessions: A list of accessions that should not be added anew
    :return: True if the accession collides with the accession list, False if not
    """
    for extant_accession in excluded_accessions:
        new_accession_stripped = new_accession.split(".")

        if new_accession_stripped[0] == extant_accession:
            return True

    return False


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


async def search_otu_path(
    seq_data, src_path: Path, taxid_table: dict, logger: BoundLogger = get_logger()
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
