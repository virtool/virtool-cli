import json
from pathlib import Path
import logging

from virtool_cli.utils.logging import base_logger
from virtool_cli.utils.ref import (
    get_otu_paths
)

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

    
def check_taxonomy(src_path):
    """
    """
    for otu_path in get_otu_paths(src_path):
        logger = base_logger.bind(
            otu_path=str(otu_path.relative_to(src_path))
        )
        check_otu(otu_path, logger)

        # for isolate_path in get_isolate_paths(otu_path):
        #     for sequence_path in get_sequence_paths(isolate_path):
        #         check_sequence(sequence_path, logger)

def check_otu(otu_path, logger):
    with open(otu_path / 'otu.json', "r") as f:
        otu = json.load(f)
    
    if 'schema' not in otu:
        logger.warning('missing schema')
    if not otu['schema']:
        logger.warning('missing schema')
                
def check_sequence(sequence_path, logger = base_logger):
    with open(sequence_path, "r") as f:
        sequence = json.load(f)
    
    if '.' not in sequence['accession']:
        logger.warning('Version not in accession')


if __name__ == '__main__':
    debug = True
    
    REPO_DIR = '/Users/sygao/Development/UVic/Virtool/Repositories'
    
    project_path = Path(REPO_DIR) / 'ref-plant-viruses'
    src_path = project_path / 'src'
    catalog_path = Path(REPO_DIR) / 'ref-accession-catalog/catalog'

    run(src_path, debug)