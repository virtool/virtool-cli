from git import Repo
from structlog import BoundLogger

from virtool_cli.utils.logging import base_logger

def generate_branchname(dirname: str) -> str:
    """
    Creates an automatically-generated branchname
    """
    return 'auto__' + dirname

def remove_branches(
    repo: Repo, branchname: list, 
    logger: BoundLogger = base_logger
):
    """
    Deletes a list of branches from the repo
    """
    for otu_branchname in branchname:
        if otu_branchname in repo.heads:
            try:
                repo.delete_head(otu_branchname, '-D')
                logger.debug('Deleted branch', branch=otu_branchname)
            except Exception as e:
                logger.error(f'Branch {otu_branchname} could not be deleted: {e}')