import json
from typing import Optional
from pathlib import Path
from structlog import BoundLogger
import logging
import re

from virtool_cli.utils.logging import base_logger
from virtool_cli.utils.reference import (
    get_otu_paths, get_isolate_paths, get_sequence_paths,
)
from virtool_cli.doctor.fix_otu import repair_otu_data

def run(src_path: Path, debugging: bool = False):
    """
    Fixes incorrect reference data

    :param src_path: Path to a given reference directory
    """
    filter_class = logging.DEBUG if debugging else logging.INFO
    logging.basicConfig(
        format="%(message)s",
        level=filter_class,
    )

    repair_data(src_path)

def repair_data(src_path):
    """
    """
    otu_paths = get_otu_paths(src_path)
    for otu_path in otu_paths:
        logger = base_logger.bind(otu_path=str(otu_path.relative_to(src_path)))

        repair_otu_data(otu_path, logger)

        # repair_otu(otu_path, logger)

        # for isolate_path in get_isolate_paths(otu_path):
        #     for sequence_path in get_sequence_paths(isolate_path):
        #         repair_sequence(
        #             sequence_path,
        #             logger.bind(seq_path=str(sequence_path.relative_to(src_path)))
        #         )

# def repair_otu(otu_path, logger = base_logger):
#     """
#     """
#     with open(otu_path / 'otu.json', "r") as f:
#         otu = json.load(f)

#     if type(otu.get('taxid')) != int:
#         otu['taxid'] = None
    
#     return otu

# def repair_sequence(sequence_path, sequence = {}, logger: BoundLogger = base_logger):
#     """
#     """
#     if not sequence:
#         with open(sequence_path, "r") as f:
#             sequence = json.load(f)
        
#     logger = logger.bind(accession=f"'{sequence['accession']}'")

#     # Automatically repair misspelled accessions where possible
#     verified_accession = verify_accession(sequence['accession'])
#     if '.' not in verified_accession:
#         # assume this is version 1 of the accession
#         verified_accession += '.1'
    
#     if sequence['accession'] != verified_accession:
#         sequence['accession'] = verified_accession
#         with open(sequence_path, "w") as f:
#             json.dump(sequence, f, indent=4)

# def verify_accession(original):
#     # Automatically repair misspelled accessions where possible
#     if re.search(r'([^A-Z_.0-9])', original) is None:
#         return original
    
#     formatted_accession = format_accession(original)

#     return formatted_accession


# def format_accession(original):
#     """
#     """
#     formatted_accession = original
#     formatted_accession = formatted_accession.strip()
#     formatted_accession = formatted_accession.upper()
#     formatted_accession = re.sub(r'-', r'_', formatted_accession)
    
#     return formatted_accession

# def fix_taxid(otu: dict) -> Optional[dict]:
#     """
#     Ensures that each taxid inside every OTU's otu.json is of type int

#     :param otu: A deserialized otu.json
#     :return: The modified otu parameter if it needs to be updated, else None
#     """
#     try:
#         taxid = otu.get("taxid", None)
#         if isinstance(taxid, str):
#             return {**otu, "taxid": int(taxid)}
#     except KeyError:
#         # assure that taxid field is set to None
#         return {**otu, "taxid": None}

#     return None