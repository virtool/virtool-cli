from pathlib import Path
import asyncio
from structlog import BoundLogger, get_logger
from urllib.error import HTTPError

from virtool_cli.utils.format import process_default
from virtool_cli.utils.ncbi import (
    request_linked_accessions,
    request_from_nucleotide,
    NCBI_REQUEST_INTERVAL,
)
from virtool_cli.utils.storage import get_otu_accessions, fetch_exclusions


async def get_no_fetch_set(otu_path: Path):
    """ """
    included = get_otu_accessions(otu_path)
    excluded = await fetch_exclusions(otu_path)

    return set(included + excluded)


async def request_new_records(
    taxid: int, no_fetch_set: set, logger: BoundLogger = get_logger()
) -> list:
    """
    :param taxid:
    :param no_fetch_set:
    :param logger: Optional entry point for a shared BoundLogger
    """
    try:
        upstream_accessions = await request_linked_accessions(taxon_id=taxid)
    except HTTPError as e:
        logger.exception(e)
        raise e
    except Exception as e:
        logger.exception(e)
        return []

    await asyncio.sleep(NCBI_REQUEST_INTERVAL)

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


async def process_records(
    records: list,
    metadata: dict,
    no_fetch_set: set,
    auto_evaluate: bool = True,
    logger: BoundLogger = get_logger(),
) -> list:
    """
    Takes sequence records and:
        1) Evaluates whether those records should be added to the database,
        2) Formats the records into a smaller dictionary
        3) Returns new formatted dicts in a list

    WARNING: Auto-evaluation is still under active development,
        especially multipartite filtering

    :param records: SeqRecords retrieved from the NCBI Nucleotide database
    :param metadata:
    :param no_fetch_set:
    :param auto_evaluate: Boolean flag for whether automatic evaluation functions
        should be run
    :param logger: Optional entry point for a shared BoundLogger
    :return: A list of valid sequences formatted for the Virtool reference database
    """
    try:
        otu_updates, auto_excluded = await process_default(
            records, metadata, no_fetch_set, logger
        )
    except Exception as e:
        logger.exception(e)
        raise e

    if auto_excluded:
        logger.info(
            "Consider adding these accessions to the exclusion list",
            auto_excluded=auto_excluded,
        )

    if otu_updates:
        return otu_updates

    return []
