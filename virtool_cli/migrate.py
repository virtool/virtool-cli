from pathlib import Path
import asyncio
import structlog
import json

from virtool_cli.utils.logging import configure_logger
from virtool_cli.utils.reference import generate_otu_dirname, is_v1
from virtool_cli.utils.storage import read_otu

base_logger = structlog.get_logger()


def run(src_path: Path, debugging: bool = False):
    """
    :param src_path: Path to a src database directory
    :param debugging: Enables verbose logs for debugging purposes
    """
    configure_logger(debugging)
    logger = base_logger.bind(src_path=str(src_path))

    if not is_v1(src_path):
        logger.info("Reference in src_path is already v2")

        return

    # Runs flatten once confirmed as a v1 reference
    try:
        asyncio.run(flatten_src(src_path))
        logger.info("Converted reference in src_path to v2")
    except Exception as e:
        logger.error("Error occured during conversion")
        logger.exception(f"Error: {e}")


async def flatten_src(src_path: Path):
    """
    Traverses through binning directories a to z and:
        1) Reassigns all OTU directories directly under src_path,
        2) Deletes the binning directory

    :param src_path: Path to a src database directory
    """
    for alpha in [alpha for alpha in src_path.glob("[a-z]") if alpha.is_dir()]:
        alpha_logger = base_logger.bind(bin=alpha.name)

        otu_paths = [otu for otu in alpha.iterdir() if otu.is_dir()]

        alpha_logger.debug(
            f"Flattening bin '{alpha.name}'", otus=[otu.name for otu in otu_paths]
        )

        for otu_path in otu_paths:
            init_exclusions(otu_path)

            otu = await read_otu(otu_path)
            new_name = generate_otu_dirname(otu.get("name"), otu.get("_id"))
            new_path = src_path / new_name

            alpha_logger.debug(
                f"Renaming '{otu_path.relative_to(src_path)}'"
                + f"to '{new_path.relative_to(src_path)}'..."
            )

            otu_path.rename(new_path)

        # Delete alpha bin
        try:
            for chaff in alpha.iterdir():
                chaff.unlink()
            alpha.rmdir()
        except Exception as e:
            alpha_logger.error("Bin deletion failed")
            alpha_logger.exception(e)


def init_exclusions(otu_path):
    with open(otu_path / "exclusions.json", "w") as f:
        json.dump([], f, indent=4, sort_keys=True)
