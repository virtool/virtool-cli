import json
from pathlib import Path
from virtool_cli.utils.ref import get_otu_paths, parse_otu
import structlog

base_logger = structlog.get_logger()

def evaluate(otu_path):
    """
    Check schema for whether the OTU is monopartite or multipartite, deal as necessary
    """
    otu_data = parse_otu(otu_path)

    logger = base_logger.bind(
        name=otu_data.get('name', ''),
        otu_id=otu_data.get('_id', ''),
    )
    schema = otu_data.get('schema', [])

    if len(schema) > 1:
        logger.debug('Multipartite', schema=schema)
        
        for part in schema:
            print(part)

    else:
        logger.debug('Monopartite', schema=schema)

if __name__ == '__main__':
    debug = True
    
    REPO_DIR = '/Users/sygao/Development/UVic/Virtool/Repositories'
    
    project_path = Path(REPO_DIR) / 'ref-mini'
    src_path = project_path / 'src'
    # catalog_path = project_path / '.cache/catalog'
    catalog_path = Path(REPO_DIR) / 'ref-fetched-accessions/src'

    for otu_path in get_otu_paths(src_path):
        evaluate(otu_path)