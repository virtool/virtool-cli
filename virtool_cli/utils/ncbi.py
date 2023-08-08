import os
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
    handle = Entrez.esearch(db="taxonomy", term=name)
    record = Entrez.read(handle)
    handle.close()

    try:
        taxid = int(record["IdList"][0])
    except IndexError:
        taxid = None

    return taxid