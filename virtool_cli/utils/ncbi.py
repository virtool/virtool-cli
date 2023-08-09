import os
import asyncio
from urllib.error import HTTPError
from Bio import Entrez

Entrez.email = os.environ.get("NCBI_EMAIL")
Entrez.api_key = os.environ.get("NCBI_API_KEY")

NCBI_REQUEST_INTERVAL = 0.3 if Entrez.email and Entrez.api_key else 0.8


async def fetch_taxid(name: str) -> int:
    """
    Searches the NCBI taxonomy database for a given OTU name, 
    then fetches and returns its taxon id. 
    Returns None if no taxon id is found.

    :param name: Name of a given OTU
    :return: Taxonomy id for the given OTU
    """
    try:
        handle = Entrez.esearch(db="taxonomy", term=name)
        record = Entrez.read(handle)
        handle.close()
    except Exception as e:
        return e

    try:
        taxid = int(record["IdList"][0])
    except IndexError:
        taxid = None

    return taxid

async def fetch_accession_gids(accessions: list):
    """
    """
    acc_stringed = []
    for acc in accessions:
        acc_stringed.append(f'{acc}[Accession]')
    
    deline_list = ' OR '.join(acc_stringed)
    
    handle = Entrez.esearch(db="nucleotide", term=deline_list)
    record = Entrez.read(handle)
    handle.close()

    gids = []
    for entrez_id in record['IdList']:
        gids.append(entrez_id)
    
    return gids

async def fetch_primary_ids(accessions: list) -> list:
    """
    """
    try:
        gids = await fetch_accession_gids(accessions)
    except Exception as e:
        raise e

    await asyncio.sleep(NCBI_REQUEST_INTERVAL)

    try:
        handle = Entrez.esummary(db="nucleotide", id=gids)
        record = Entrez.read(handle)
        handle.close()
    except Exception as e:
        raise e

    indexed_accessions = {}
    for r in record:
        indexed_accessions[r.get('Caption')] = int(r.get('Id'))
        
        # if 'AccessionVersion' in r:
        #     indexed_accessions[gid] = r.get('AccessionVersion')
        # else:
    
    return indexed_accessions