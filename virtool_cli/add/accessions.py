import re
from pathlib import Path

import click
from structlog import get_logger

from virtool_cli.add.otu import add_otu
from virtool_cli.ncbi.client import NCBIClient
from virtool_cli.ref.repo import EventSourcedRepo as Repo
from virtool_cli.ref.otu import group_genbank_records_by_isolate


def add_accessions(
    repo: Repo, taxid: int, accessions: list, ignore_cache: bool = False
):
    """Add a list of accessions to an OTU. Appropriate if you know the OTU path already.

    :param accessions: A list of accession numbers to query from NCBI
    :param otu_path: Path to an OTU directory
    """
    logger = get_logger("add_accessions", accessions=accessions)

    client = NCBIClient.from_repo(repo.path, ignore_cache=ignore_cache)

    fetch_list = list(set(accessions).difference(set(otu.nofetch)))

    records = client.fetch_genbank_records(fetch_list)
    if records:
        grouped_records = group_genbank_records_by_isolate(records)


def verify_accession(accession: str) -> bool:
    """Returns True if the accession matches NCBI standards
    The accession must be free of characters other than digits,
    capital letters A-Z, '.', and '_'

    :param accession: An accession to be inspected
    :return: Boolean based on whether the accession is free of invalid characters
    """
    return re.search(r"([^A-Z_.0-9])", accession) is None
