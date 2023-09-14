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
    Fixes incorrect reference data

    :param src_path: Path to a given reference directory
    """
    filter_class = logging.DEBUG if debugging else logging.INFO
    logging.basicConfig(
        format="%(message)s",
        level=filter_class,
    )
    
    logger = base_logger.bind(otu_path=otu_path.name)
    logger.info('Inspecting OTU for repairs...', src=str(src_path))

    repair_otu_data(otu_path, logger)

def repair_otu_data(otu_path, logger: BoundLogger = base_logger):
    """
    """
    repair_otu(otu_path, logger)

    for isolate_path in get_isolate_paths(otu_path):
        for sequence_path in get_sequence_paths(isolate_path):
            logger = logger.bind(sequence_path=str(sequence_path.relative_to(otu_path)))
            logger.debug(f"Inspecting sequence '{sequence_path.stem}'...")
            
            repair_sequence(sequence_path=sequence_path, logger=logger)

def repair_otu(otu_path, logger: BoundLogger = base_logger):
    """
    """
    with open(otu_path / 'otu.json', "r") as f:
        otu = json.load(f)

    if type(otu.get('taxid')) != int:
        otu['taxid'] = None
    
    if 'schema' not in otu:
        otu['schema'] = []
        
    return otu

def repair_sequence(
    sequence_path, sequence: dict = {},
    logger = base_logger
):
    """
    """
    if not sequence:
        with open(sequence_path, "r") as f:
            sequence = json.load(f)
        
    # logger = logger.bind(accession=f{sequence['accession']}")

    # Automatically repair misspelled accessions where possible
    verified_accession = verify_accession(sequence['accession'])
    if '.' not in verified_accession:
        # assume this is version 1 of the accession
        verified_accession += '.1'
    
    if sequence['accession'] != verified_accession:
        sequence['accession'] = verified_accession
        with open(sequence_path, "w") as f:
            json.dump(sequence, f, indent=4)

def verify_accession(original):
    """
    """
    # Automatically repair misspelled accessions where possible
    if re.search(r'([^A-Z_.0-9])', original) is None:
        return original
    
    formatted_accession = format_accession(original)

    return formatted_accession


def format_accession(original):
    """
    """
    formatted_accession = original
    formatted_accession = formatted_accession.strip()
    formatted_accession = formatted_accession.upper()
    formatted_accession = re.sub(r'-', r'_', formatted_accession)
    
    return formatted_accession

def fix_taxid(otu: dict) -> Optional[dict]:
    """
    Ensures that each taxid inside every OTU's otu.json is of type int

    :param otu: A deserialized otu.json
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