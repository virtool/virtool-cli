from pathlib import Path
from git import Repo
import logging
from structlog import BoundLogger

from virtool_cli.utils.logging import base_logger
from virtool_cli.utils.ref import get_otu_paths
from virtool_cli.git.helpers import generate_branchname

def run_cleanup(
    repo_path: Path,
    debugging: bool = False
):
    """
    """
    filter_class = logging.DEBUG if debugging else logging.INFO
    logging.basicConfig(
        format="%(message)s",
        level=filter_class,
    )

    logger = base_logger.bind(repository=str(repo_path))
    
    repo = Repo(repo_path)

    cleanup_repo(repo, src_path=repo_path/'src')

    logger.info('Removed all auto-- branches')

def cleanup_repo(
    repo: Repo, src_path: Path, 
    logger: BoundLogger = base_logger
):
    """
    Deletes all auto-made branches from the repo
    """
    all_otu_paths = get_otu_paths(src_path)

    for otu_path in all_otu_paths:
        otu_branchname = generate_branchname(otu_path.name)
        if otu_branchname in repo.heads:
            try:
                repo.delete_head(otu_branchname, '-D')
                logger.debug('Deleted branch', branch=otu_branchname)
            except Exception as e:
                logger.error(f'Branch {otu_branchname} could not be deleted: {e}')