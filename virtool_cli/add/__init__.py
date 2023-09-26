import json
from pathlib import Path
import asyncio

# from typing import Optional, Tuple
import logging
from structlog import BoundLogger
from urllib.error import HTTPError

from virtool_cli.utils.logging import base_logger

# from virtool_cli.utils.reference import get_otu_paths, get_isolate_paths
# from virtool_cli.utils.id_generator import generate_unique_ids, get_unique_ids
# from virtool_cli.utils.ncbi import (
#     request_linked_accessions,
#     request_from_nucleotide,
#     NCBI_REQUEST_INTERVAL,
# )
# from virtool_cli.update.evaluate import (
#     evaluate_sequence,
#     get_qualifiers,
#     get_lengthdict_monopartite,
#     get_lengthdict_multipartite,
# )
# from virtool_cli.utils.storage import store_isolate, store_sequence
# from virtool_cli.utils.format import format_sequence
from virtool_cli.catalog.helpers import search_by_otu_id

DEFAULT_INTERVAL = 0.001


def run(
    otu_path: Path,
    catalog_path: Path,
    auto_evaluate: bool = False,
    debugging: bool = False,
):
    """
    CLI entry point for update.update.run()

    Requests updates for a single OTU directory
    Searches the catalog for a matching catalog listing and requests new updates if it finds one.

    :param otu_path: Path to a OTU directory
    :param catalog_path: Path to a catalog directory
    :param auto_evaluate: Auto-evaluation flag, enables automatic filtering for fetched results
    :param debugging: Enables verbose logs for debugging purposes
    """
    filter_class = logging.DEBUG if debugging else logging.INFO
    logging.basicConfig(
        format="%(message)s",
        level=filter_class,
    )

    if auto_evaluate:
        base_logger.warning(
            "Auto-evaluation is in active development and may produce false negatives."
        )

    [_, otu_id] = otu_path.name.split("--")
    listing_path = search_by_otu_id(otu_id, catalog_path)

    if listing_path:
        asyncio.run(update_otu(otu_path, listing_path, auto_evaluate))

    else:
        base_logger.critical("Listing not found for this OTU.")
        base_logger.info("Use virtool acc init to create a new accession catalog.")
