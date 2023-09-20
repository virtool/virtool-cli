import json
from typing import Optional
from pathlib import Path
from structlog import BoundLogger
import logging
import re

from virtool_cli.utils.logging import base_logger
from virtool_cli.utils.reference import get_isolate_paths, get_sequence_paths


def run(otu_path: Path, src_path: Path, debugging: bool = False):
    """
    CLI entry point for doctor.fix_otu.run()

    :param otu_path: Path to a OTU directory under a reference directory
    :param src_path: Path to a given reference directory
    :param debugging: Enables verbose logs for debugging purposes
    """
    filter_class = logging.DEBUG if debugging else logging.INFO
    logging.basicConfig(
        format="%(message)s",
        level=filter_class,
    )

    logger = base_logger.bind(otu_path=otu_path.name)
    logger.info("Inspecting OTU for repairs...", src=str(src_path))

    repair_otu(otu_path, logger)


def repair_otu(otu_path: Path, logger: BoundLogger = base_logger):
    """
    Inspects an OTU and its contents for repairable issues and corrects them if found

    :param otu_path: Path to a OTU directory under a reference directory
    :param logger: Optional entry point for an existing BoundLogger
    """
    logger.debug("Getting OTU data...")
    try:
        with open(otu_path / "otu.json", "r") as f:
            otu = json.load(f)
    except Exception as e:
        logger.exception(e)

    checked_otu = repair_otu_data(otu)
    if otu != checked_otu:
        with open(otu_path / "otu.json", "w") as f:
            json.dump(otu, f, indent=4, sort_keys=True)

    for isolate_path in get_isolate_paths(otu_path):
        for sequence_path in get_sequence_paths(isolate_path):
            logger = logger.bind(sequence_path=str(sequence_path.relative_to(otu_path)))

            logger.debug(f"Inspecting sequence '{sequence_path.stem}'...")

            repair_sequence(sequence_path=sequence_path, logger=logger)


def repair_otu_data(otu: dict):
    """
    Deserializes an OTU's otu.json and updates the dictionary if issues are found
    and returns the dictionary

    :param otu_path: Path to a OTU directory under a reference directory
    """
    new_otu = otu.copy()
    if type(otu.get("taxid")) != int:
        new_otu["taxid"] = None

    if "schema" not in otu:
        new_otu["schema"] = []

    return otu


def repair_sequence(sequence_path: Path, logger: BoundLogger = base_logger):
    """
    :param logger: Optional entry point for an existing BoundLogger
    """
    sequence = json.loads(sequence_path.read_text())

    # Automatically repair misspelled accessions where possible
    verified_accession = correct_accession(sequence["accession"])
    if "." not in verified_accession:
        # assume this is version 1 of the accession
        verified_accession += ".1"

    if sequence["accession"] != verified_accession:
        logger.debug(f"Changes made, writing changes to {sequence_path.name}...")

        sequence["accession"] = verified_accession
        with open(sequence_path, "w") as f:
            json.dump(sequence, f, indent=4, sort_keys=True)


def correct_accession(accession: str) -> str:
    """
    :param accession: The accession to be corrected
    :return: Corrected accession
    """
    # Automatically repair misspelled accessions where possible
    if re.search(r"([^A-Z_.0-9])", accession) is None:
        return accession

    formatted_accession = accession
    formatted_accession = formatted_accession.strip()
    formatted_accession = formatted_accession.upper()
    formatted_accession = re.sub(r"-", r"_", formatted_accession)

    corrected_accession = formatted_accession

    return corrected_accession


def fix_taxid(otu: dict) -> Optional[dict]:
    """
    Ensures that each taxid inside every OTU's otu.json is of type int

    :param otu: A deserialized otu.json OTU
    :return: The modified otu parameter if it needs to be updated, else None
    """
    try:
        taxid = otu.get("taxid", None)
        if isinstance(taxid, str):
            return {**otu, "taxid": int(taxid)}
    except KeyError:
        # assure that taxid field is set to None
        return {**otu, "taxid": None}

    return None
