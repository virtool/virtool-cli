from pathlib import Path
import structlog

from virtool_cli.utils.ref import parse_otu, generate_otu_dirname

logger = structlog.get_logger()

def run(src_path):
    """
    :param src_path: Path to a src database directory
    """
    log = logger.bind(src_path=str(src_path))
    
    if not [alpha for alpha in src_path.glob('[a-z]')]:
        log.info("Reference in src_path is already v2")
        return
    
    log.info("Converting reference in src_path to v2...")
    try:
        flatten_src(src_path)
    except Exception as e:
        log.error('Error occured during conversion')
        log.exception(f'Error: {e}')

def flatten_src(src_path: Path):
    """
    Traverses through binning directories a to z and reassigns 
    all OTU directories to have src_path as a direct parent,
    then deletes the binning directory

    :param src_path: Path to a src database directory
    """

    for alpha in [alpha for alpha in src_path.glob('[a-z]')]:

        otu_paths = [otu for otu in alpha.iterdir() if otu.is_dir()]
        
        for otu_path in otu_paths:
            otu = parse_otu(otu_path)
            new_name = generate_otu_dirname(
                otu.get('name'), 
                otu.get('_id')
            )
            new_path = src_path / new_name
            otu_path.rename(new_path)
            print(new_path)

        # Delete alpha bin
        try:
            for chaff in alpha.iterdir():
                chaff.unlink()
        except Exception as e:
            return e
        alpha.rmdir()