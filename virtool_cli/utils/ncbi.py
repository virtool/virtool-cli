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

async def esearch_uids(accessions: list) -> list:
    """
    Queries Entrez esearch for a list of accessions

    :param name: List of accession numbers
    """
    acc_stringed = []
    for acc in accessions:
        acc_stringed.append(f'{acc}[Accession]')
    
    deline_list = ' OR '.join(acc_stringed)
    
    try:
        handle = Entrez.esearch(db="nucleotide", term=deline_list)
        record = Entrez.read(handle)
        handle.close()
    except Exception as e:
        return []
    
    if len(record) < 1:
        return []

    gids = []
    for entrez_id in record['IdList']:
        gids.append(entrez_id)
    
    return gids

async def fetch_accession_uids(accessions: list) -> dict:
    """
    Creates a dictionary keyed by the input accessions,
    fetches the corresponding UIDs from NCBI Genbank and 
    returns as many key-value pairs as possible.
    If a UID cannot be found for an accession, value defaults to None.add()

    :param name: List of accession numbers
    :return: Dictionary of Genbank UIDs keyed by corresponding accession number
    """
    indexed_accessions = { accession: None for index, accession in enumerate(accessions) }
    
    try:
        uids = await esearch_uids(accessions)
    except Exception as e:
        raise e

    await asyncio.sleep(NCBI_REQUEST_INTERVAL)

    if not uids:
        return indexed_accessions

    try:
        handle = Entrez.esummary(db="nucleotide", id=uids)
        record = Entrez.read(handle)
        handle.close()
    except Exception as e:
        raise e

    for r in record:
        if r.get('Caption') in indexed_accessions.keys():
            indexed_accessions[r.get('Caption')] = int(r.get('Id'))
        elif r.get('AccessionVersion') in indexed_accessions.keys():
            indexed_accessions[r.get('AccessionVersion')] = int(r.get('Id'))
    
    return indexed_accessions