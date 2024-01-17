import json
from pathlib import Path
from shutil import copytree
import structlog
from virtool_cli.utils.logging import DEFAULT_LOGGER, DEBUG_LOGGER

base_logger = structlog.get_logger()


def run(repo_path: Path, debugging: bool = False):
    """
    Creates an empty repository for Virtool reference data

    :param repo_path: Path to repository src directory
    :param debugging: Enables verbose logs for debugging purposes
    """
    structlog.configure(wrapper_class=DEBUG_LOGGER if debugging else DEFAULT_LOGGER)
    logger = base_logger.bind(repo=str(repo_path))

    try:
        logger.info("Creating new reference database at path...")
        initialize_reference(repo_path)
    except FileNotFoundError as e:
        logger.exception(e)
    except Exception as e:
        logger.exception(e)


def initialize_reference(repo_path: Path):
    """
    Creates an empty directory structure for Virtool reference data,
    complete with hidden .cache/catalog subdirectories

    :param repo_path: Path to repository src directory
    """
    logger = base_logger.bind(repo=str(repo_path))

    if not repo_path.exists():
        repo_path.mkdir()
        logger.info("Created a new directory at repo_path")

    src_path = repo_path / "src"
    if not src_path.exists():
        src_path.mkdir()
        logger.debug("Created a new src directory under repo_path")

    with open(src_path / "meta.json", "w") as f:
        json.dump({"data_type": "genome"}, f)

    package_path = Path(__file__).parents[1]
    model_path = package_path / "assets/github"
    logger.debug(model_path)

    github_path = repo_path / ".github"
    copytree(model_path, github_path)

    cache_path = repo_path / ".cache"
    if not cache_path.exists():
        cache_path.mkdir()
        logger.debug("Created a new hidden cache directory under repo_path")

    logger.info("Reference repository at repo_path is now complete.")

    # Sent new repo path to stdout
    print(str(repo_path))
