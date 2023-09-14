import json
import aiofiles
from pathlib import Path
from structlog import BoundLogger
from virtool_cli.utils.ncbi import fetch_taxid
from virtool_cli.utils.logging import base_logger


async def generate_listing(
    otu_data: dict, 
    accession_list: list, 
    sequence_metadata: dict,
    logger: BoundLogger = base_logger
) -> dict:
    """
    Generates a new listing for a given OTU and returns it as a dict

    :param otu_data: OTU data in dict form
    :param accession_list: list of included accesssions
    :param sequence_metadata: dict of sequence metadata, including average length
    :param logger: Optional entry point for an existing BoundLogger
    """
    catalog_listing = {}
    
    otu_id = otu_data.get('_id')
    taxid = otu_data.get('taxid', None)

    catalog_listing['_id'] = otu_id

    # Attempts to fetch the taxon id if none is found in the OTU metadata
    if taxid is None:
        logger.info('Taxon ID not found. Attempting to fetch from NCBI Taxonomy...')
        
        try:
            taxid = await fetch_taxid(otu_data.get('name', None))
        except Exception as e:
            logger.exception(e)

        if taxid is None:
            catalog_listing['taxid'] = 'none'
            logger.debug(f'Matching ID not found in NCBI Taxonomy. Setting taxid={taxid}')    
        else:
            catalog_listing['taxid'] = int(taxid)
            logger.debug(f'Matching ID found in NCBI Taxonomy. Setting taxid={taxid}')

    else:
        catalog_listing['taxid'] = int(taxid)

    catalog_listing['name'] = otu_data['name']

    catalog_listing['accessions'] = { 'included': accession_list }

    schema = otu_data.get('schema', [])
    logger.debug(f'{schema}')
    
    if schema:
        if len(schema) > 1:
            required_parts = get_required_parts(schema)
            try:
                length_dict = measure_multipartite(sequence_metadata, required_parts)
            except Exception as e:
                logger.exception(e)
            
            for part in schema:
                if part['name'] in length_dict.keys():
                    part['length'] = length_dict[part['name']]

        else:
            try:
                average_length = measure_monopartite(sequence_metadata)
                logger.debug(f'Average length is {average_length}')
                schema[0]['length'] = average_length
            except Exception as e:
                logger.exception(e)
    
    catalog_listing['schema'] = schema
    
    return catalog_listing

async def update_listing(data: dict, path: Path):
    """
    Updates the listing file at path with new data

    :param data: Updated listing data
    :param path: Path to the listing file
    """
    async with aiofiles.open(path, "w") as f: 
        await f.write(json.dumps(data, indent=2, sort_keys=True))

async def write_listing(
    taxid: int, 
    listing: dict, 
    catalog_path: Path,
    indent: bool = True,
    logger: BoundLogger = base_logger
):
    """
    Writes prepared listing data to a listing file under the catalog directory

    :param taxid: OTU taxon id
    :param accessions: List of OTU's accessions
    :param catalog: Path to an accession catalog
    :param indent: Indent flag
    :param logger: Optional entry point for an existing BoundLogger
    """

    output_path = catalog_path / f"{taxid}--{listing['_id']}.json"
    logger = logger.bind(
        taxid=taxid, 
        listing_path=str(output_path.relative_to(output_path.parent)))

    logger.debug('Writing accession listing...')

    listing['accessions']['excluded'] = {}

    with open(output_path, "w") as f:
        json.dump(listing, f, indent=2 if indent else None, sort_keys=True)

def get_required_parts(schema: list) -> dict:
    """
    """
    required_parts = []
    for part in schema:
        if part['required']:
            required_parts.append(part['name'])
    
    return required_parts

def measure_monopartite(sequence_metadata: dict) -> int:
    """
    """
    sequence_lengths = []
    print(sequence_metadata)
    for accession in sequence_metadata:
        metadata = sequence_metadata[accession]
        seq_length = metadata['length']
        sequence_lengths.append(seq_length)
    
    average_length = sum(sequence_lengths)/len(sequence_lengths)

    return int(average_length)

def measure_multipartite(sequence_metadata, part_list) -> dict:
    """
    """
    part_total_dict = { element: [] for index, element in enumerate(part_list) }
            
    for accession in sequence_metadata:
        metadata = sequence_metadata[accession]
        seq_length = metadata['length']
        
        segment_name = metadata.get('segment', None)
        if segment_name is None:
            continue
        if segment_name in part_total_dict.keys():
            part_total_dict[segment_name].append(seq_length)
    
    length_dict = {}
    for segment in part_total_dict:
        if not part_total_dict[segment]:
            length_dict[segment] = 0
            continue
        average_length = sum(part_total_dict[segment]) / len(part_total_dict[segment])
        length_dict[segment] = int(average_length)

    # return int(average_length)
    return length_dict