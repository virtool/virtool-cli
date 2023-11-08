from pathlib import Path

import structlog
from structlog import BoundLogger

from virtool_cli.check.checkup import (
    validate_otu,
    validate_isolate,
    validate_sequence,
)
from virtool_cli.utils.logging import DEBUG_LOGGER, DEFAULT_LOGGER
from virtool_cli.utils.reference import (
    get_otu_paths,
    get_isolate_paths,
    get_sequence_paths,
)

base_logger = structlog.get_logger()


def run(src_path: Path, debugging: bool = False):
    """
    CLI entry point for doctor.checkup.run()

    :param src_path: Path to a given reference directory
    :param debugging: Enables verbose logs for debugging purposes
    """
    structlog.configure(wrapper_class=DEBUG_LOGGER if debugging else DEFAULT_LOGGER)
    logger = base_logger.bind(src=str(src_path), verbose=debugging)

    logger.debug("Debug flag is enabled")

    if check_reference(src_path):
        logger.info("Reference folder is valid.", valid_result=True)
    else:
        logger.error("Reference folder is invalid.", valid_result=False)


def check_reference(src_path: Path, logger: BoundLogger = base_logger) -> bool:
    """
    Inspects all files in the reference directory and logs a warning if a problem is found.

    :param src_path: Path to a given reference directory
    :param logger: Optional entry point for a shared BoundLogger
    """
    ref_valid = True

    logger.info("Running checks on reference directory...")

    watch_list = []

    for otu_path in get_otu_paths(src_path):
        otu_valid = True

        try:
            [_, otu_id] = otu_path.name.split("--")
        except ValueError as e:
            logger.exception(e)
            ref_valid = False
            continue

        otu_logger = logger.bind(otu_id=otu_id)

        otu_logger.debug(f"Running checks on OTU {otu_id}...")

        if not validate_otu(otu_path, otu_logger):
            otu_valid = False

        isolate_paths = get_isolate_paths(otu_path)
        if not isolate_paths:
            otu_logger.error("No accession data in OTU", otu=otu_path.name)
            otu_valid = False

        else:
            for isolate_path in isolate_paths:
                if not validate_isolate(isolate_path, logger):
                    otu_valid = False

                sequence_paths = get_sequence_paths(isolate_path)
                if not sequence_paths:
                    otu_logger.error(
                        "No sequences in isolate", isolate=isolate_path.name
                    )
                    otu_valid = False

                for sequence_path in sequence_paths:
                    if not validate_sequence(sequence_path, logger):
                        otu_valid = False

        if not otu_valid:
            ref_valid = False
            watch_list.append(otu_path.name)

    # Print outstanding issues to stdout
    print(watch_list)

    return ref_valid
