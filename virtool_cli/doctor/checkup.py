import json
from pathlib import Path
import re
import logging

from virtool_cli.utils.logging import base_logger
from virtool_cli.utils.reference import get_otu_paths, get_isolate_paths, get_sequence_paths

def run(src_path: Path, debugging: bool = False):
    """
    :param src_path: Path to a given reference directory
    """
    filter_class = logging.DEBUG if debugging else logging.INFO
    logging.basicConfig(
        format="%(message)s",
        level=filter_class,
    )

    check_taxonomy(src_path)

    
def check_taxonomy(src_path: Path):
    """

    """
    for otu_path in get_otu_paths(src_path):
        logger = base_logger.bind(
            otu_path=str(otu_path.relative_to(src_path))
        )
        check_otu(otu_path, logger)

        for isolate_path in get_isolate_paths(otu_path):
            for sequence_path in get_sequence_paths(isolate_path):
                check_sequence(sequence_path, logger)

def check_otu(otu_path: Path, logger):
    """
    """
    with open(otu_path / 'otu.json', "r") as f:
        otu = json.load(f)
    
    if not otu.get('schema', []):
        logger.warning(f"{otu['_id']}: missing schema")
                
def check_sequence(sequence_path: Path, logger = base_logger):
    """

    """
    logger = logger.bind(sequence_id=sequence_path.stem)
    sequence = json.loads(sequence_path.read_text())
    
    accession = sequence['accession']

    if not verify_accession(accession):
        logger.error(f"Accession '{accession}' contains invalid characters or capitalization")

    if '.' not in accession:
        logger.warning(f"Version not included in accession={accession}")

def verify_accession(accession: str):
    """
    """
    if re.search(r'([^A-Z_.0-9])', accession) is None:
        return True
    
    return False