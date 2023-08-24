from pathlib import Path
import json
import asyncio
import aiofiles
from Bio import Entrez, SeqIO
import structlog
import logging
from urllib.error import HTTPError

from virtool_cli.utils.logging import base_logger
from virtool_cli.utils.ref import parse_otu, get_otu_paths, search_otu_by_id
from virtool_cli.utils.ncbi import fetch_accession_uids
from virtool_cli.accessions.initialize import generate_listing, write_listing
from virtool_cli.accessions.helpers import (
    parse_listing, split_pathname, get_otu_accessions, 
    get_catalog_paths, filter_catalog
)

def run(src: Path, catalog: Path, debugging: bool = False):
    """
    """
    filter_class = logging.DEBUG if debugging else logging.INFO
    logging.basicConfig(
        format="%(message)s",
        level=filter_class,
    )
    logger = base_logger.bind(catalog=str(catalog))
    
    logger.info(f"Starting RefSeq filter...")
    
    asyncio.run(filter_refseq_accessions(catalog))

async def filter_refseq_accessions(catalog: Path):
    """
    """
    for listing_path in catalog.glob('*.json'):
        listing = parse_listing(listing_path)
        
        logger = base_logger.bind(otu_id=listing['_id'], taxid=listing['taxid'], name=listing['name'])

        if listing['accessions']['excluded']:
            continue
        
        # print(listing['accessions']['included'])
        refseq_accessions = []
        for accession in listing['accessions']['included']:
            # Only RefSeq accessions include underscores
            if '_' in accession:
                refseq_accessions.append(accession)

        if not refseq_accessions:
            continue
        
        logger.debug('RefSeq accessions found', refseq=refseq_accessions)
        try:
            included_records = await fetch_upstream_records(refseq_accessions)
        except Exception as e:
            logger.exception(e)
            continue
        
        excluded_records = []
        for record in included_records:
            # if 'RefSeq' in record.annotations['keywords']:
            comment_list = record.annotations['comment'].split('.')

            # print(comment_list)
            for part in comment_list:
                if 'identical' in part or 'derived from' in part:
                    accession = part.strip().split()[-1:][0]
                    if accession not in listing['accessions']['excluded']:
                        excluded_records.append(accession)
                if 'derived from' in part:
                    accession = part.strip().split()[-1:][0]
        
        if excluded_records:
            try:
                new_excluded_accessions = await fetch_accession_uids(excluded_records)
            except Exception as e:
                logger.exception(e)

            logger.info('Adding accessions to exclusion list...', new=new_excluded_accessions)

            listing['accessions']['excluded'].update(new_excluded_accessions)

            async with aiofiles.open(listing_path, "w") as f: 
                await f.write(
                    json.dumps(listing, indent=2, sort_keys=True)
                )

async def fetch_upstream_accessions(
    listing: dict,
) -> list:
    """
    :param taxid: OTU Taxon ID
    :param listing: Corresponding listing from the accession catalog for this OTU in dict form
    :return: A list of accessions from NCBI Genbank for the taxon ID, 
        sans included and excluded accessions
    """
    taxid = listing['taxid']

    filter_set = set()
    filter_set.update(listing['accessions']['included'].values())
    filter_set.update(listing['accessions']['excluded'].values())

    logger = base_logger.bind(
        taxid=taxid, otu_id=listing['_id'])
    logger.debug('Exclude catalogued UIDs', catalogued=filter_set)

    upstream_uids = []
    
    entrez_accessions = Entrez.read(
        Entrez.elink(
            dbfrom="taxonomy", db="nucleotide", 
            id=str(taxid), idtype="acc")
    )

    for linksetdb in entrez_accessions[0]["LinkSetDb"][0]["Link"]:
        upstream_uids.append(linksetdb["Id"])

    upstream_set = set(upstream_uids)

    return list(upstream_set.difference(filter_set))

async def fetch_upstream_records(fetch_list: list, logger = base_logger) -> list:
    """
    Take a list of accession numbers and request the records from NCBI GenBank
    
    :param fetch_list: List of accession numbers to fetch from GenBank

    :return: A list of GenBank data converted from XML to dicts if possible, 
        else an empty list
    """

    try:
        handle = Entrez.efetch(
            db="nucleotide", id=fetch_list, rettype="gb"
        )
    except HTTPError as e:
        return []
    
    parsed = SeqIO.parse(handle, "gb")
    
    try:
        record_dict = SeqIO.to_dict(parsed)
    except ValueError:
        base_logger.error(
            'Found two copies of the same record. Remove duplicate from this.'
        )
        return []
    except:
        base_logger.exception(e)
        return []

    if record_dict is None:
        return []
    
    try:
        accession_list = [record for record in record_dict.values() if record.seq]
        return accession_list
    except Exception as e:
        base_logger.exception(e)
        return []

if __name__ == '__main__':
    debug = True
    
    REPO_DIR = '/Users/sygao/Development/UVic/Virtool/Repositories'
    
    project_path = Path(REPO_DIR) / 'ref-plant-viruses'
    src_path = project_path / 'src'
    # catalog_path = project_path / '.cache/catalog'
    catalog_path = Path(REPO_DIR) / 'ref-accession-catalog/catalog'
    # catalog_path = Path('/Users/sygao/Development/UVic/Virtool/TestSets/cotton_dupe/.cache/catalog')

    run(src_path, catalog_path, debug)
