import sys

from structlog import get_logger

from virtool_cli.ncbi.client import NCBIClient
from virtool_cli.ref.repo import EventSourcedRepo as Repo

OTU_KEYS = ["_id", "name", "abbreviation", "schema", "taxid"]


logger = get_logger()
