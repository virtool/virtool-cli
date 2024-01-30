"""Initialize new reference repositories.

TODO: Add command line options for setting the data_type and organism.

"""
import json
from pathlib import Path
from shutil import copytree

from structlog import get_logger


def init_reference(path: Path):
    """Create an empty reference repository.

    :param path: the path to initialize the repository at
    """
    logger = get_logger("Initialize reference", path=str(path))

    logger.info("Initializing new reference repository")

    if path.is_file():
        logger.error("The target path is a file")
        return

    path.mkdir(parents=True, exist_ok=True)

    if any(path.iterdir()):
        logger.error("The directory is not empty")
        return

    with open(path / ".gitignore", "w") as f:
        f.write(".cache\n")

    src_path = path / "src"
    src_path.mkdir()

    with open(src_path / "meta.json", "w") as f:
        json.dump({"data_type": "genome", "organism": ""}, f)

    copytree(Path(__file__).parent.parent.parent / "assets/github", path / ".github")

    (path / ".cache").mkdir()

    logger.info("Complete")
