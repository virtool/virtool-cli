from pathlib import Path
from git import Repo
import logging
from structlog import BoundLogger

from virtool_cli.utils.logging import base_logger
from virtool_cli.utils.ref import get_otu_paths
from virtool_cli.git.clean import cleanup_repo
from virtool_cli.git.helpers import generate_branchname

def run_split(
    repo_path: Path, 
    branch: str, main_branch: str, 
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

    logger.debug("Printing branches", branches=repo.heads)

    new_branches = split_branches(
        repo_path, repo,
        branch_name=branch, 
        main_branch_name=main_branch,
        logger=logger)
    
    if new_branches:
        logger.info(
            f"Updated OTUs from '{branch}' split into individual branches", 
        n_autobranches=len(new_branches))

def split_branches(
    repo_path: Path,
    repo: Repo,
    branch_name: str, 
    main_branch_name='main',
    logger: BoundLogger = base_logger
):
    """
    """
    if repo.bare:
        raise Exception
    git = repo.git
    
    src_path = repo_path / 'src'
    
    heads = repo.heads

    main_branch = heads[main_branch_name]
    main_tree = main_branch.commit.tree

    updated_branch = heads[branch_name]
    update_tree = updated_branch.commit.tree

    logger.debug(
        f'Checking updated branch against base branch', 
        base_branch=main_branch_name, updated_branch=branch_name
    )
    
    all_otu_paths = get_otu_paths(src_path)
    changed_otu_paths = []

    for otu_path in all_otu_paths:
        gitpath = str(otu_path.relative_to(repo_path))
        
        otu_diffs = update_tree[gitpath].diff(main_tree[gitpath])
        n_changes = len(otu_diffs)

        if n_changes > 0:
            logger.debug(f'{n_changes} changes found', path=gitpath)
            changed_otu_paths.append(otu_path)
        else:
            logger.debug(f'No changes found', path=gitpath)
    
    logger.info(f"{len(changed_otu_paths)} out of {len(all_otu_paths)} OTUs have pending changes")

    new_branches = []
    for otu_path in changed_otu_paths:
        gitpath = str(otu_path.relative_to(repo_path))
        otu_log = logger.bind(path=gitpath)

        new_branchname = generate_branchname(otu_path.name)
        
        if new_branchname in repo.heads:
            logger.warning(
                'Found a previous auto-fetched branch of this name, ' + \
                'removing all auto-branches first...')
            repo.head.reference = main_branch
            repo.head.reset(index=True, working_tree=True)
            cleanup_repo(repo, src_path)
        
        try:
            otu_log.debug('Creating new branch...')
            new_branch = repo.create_head(new_branchname)

            repo.head.reference = new_branch
            repo.head.reset(index=True, working_tree=True)
            git.checkout(branch_name, gitpath)
            git.commit(m=f'Auto-fetched updates for {otu_path.name}')
            
            otu_log.debug('OTU changes committed to branch', branch=new_branchname)
            new_branches.append(new_branchname)

        except Exception as e:
            logger.critical(f'Branch creation and checkout failed, deleting all auto-branches...')
            repo.head.reference = main_branch
            repo.head.reset(index=True, working_tree=True)
            cleanup_repo(repo, src_path)
            return []

        repo.head.reference = main_branch
        repo.head.reset(index=True, working_tree=True)
    
    return new_branches