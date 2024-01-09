from pathlib import Path
import asyncio
import structlog
from structlog import BoundLogger
from urllib.error import HTTPError

from virtool_cli.utils.ncbi import (
    request_linked_accessions,
    request_from_nucleotide,
    NCBI_REQUEST_INTERVAL,
)
from virtool_cli.utils.storage import get_otu_accessions, fetch_exclusions

base_logger = structlog.get_logger()


async def get_no_fetch_set(otu_path):
    included = await get_otu_accessions(otu_path)
    excluded = await fetch_exclusions(otu_path)

    return set(included + excluded)


async def request_new_records(
    taxid: int, no_fetch_set: set, logger: BoundLogger = base_logger
) -> list:
    """
    :param taxid:
    :param no_fetch_set:
    :param logger: Optional entry point for a shared BoundLogger
    """
    try:
        upstream_accessions = await request_linked_accessions(taxon_id=taxid)
        await asyncio.sleep(NCBI_REQUEST_INTERVAL)
    except HTTPError as e:
        logger.error(e)
        return []
    except Exception as e:
        logger.exception(e)
        return []

    if upstream_accessions:
        upstream_set = set(upstream_accessions)
        new_accessions = list(upstream_set.difference(no_fetch_set))
    else:
        return []

    if not new_accessions:
        return []

    try:
        record_data = await request_from_nucleotide(new_accessions)
        await asyncio.sleep(NCBI_REQUEST_INTERVAL)

    except HTTPError as e:
        logger.error(e)
        return []

    except Exception as e:
        logger.exception(e)
        return []

    return record_data
#
#
# async def fetch_upstream_accessions(
#     taxid: int, included: list, excluded: list, logger: BoundLogger = base_logger
# ) -> list:
#     """
#     Requests a list of all uninspected accessions associated with an OTU's taxon ID
#
#     :param listing: Corresponding catalog listing for this OTU
#     :param logger: Optional entry point for an existing BoundLogger
#     :return: A list of accessions from NCBI Genbank for the taxon ID,
#         sans included and excluded accessions
#     """
#     logger = logger.bind(taxid=taxid)
#     included_set = set(included)
#
#     try:
#         upstream_accessions = await request_linked_accessions(taxon_id=taxid)
#     except HTTPError:
#         logger.error("Could not retrieve new accessions from NCBI.")
#         return []
#
#     upstream_set = set(upstream_accessions)
#
#     return list(upstream_set.difference(included_set))
