import os
import asyncio
from urllib.error import HTTPError
from Bio import Entrez, SeqIO

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

async def fetch_docsums(fetch_list: list) -> list:
    """
    Take a list of accession numbers and requests documents from NCBI GenBank
    
    :param fetch_list: List of accession numbers to fetch from GenBank

    :return: A list of GenBank data converted from XML to dicts if possible, 
        else an empty list
    """
    try:
        handle = Entrez.efetch(db="nucleotide", id=fetch_list, rettype="docsum")
        record = Entrez.read(handle)
        handle.close()
    except Exception as e:
        return []
    
    await asyncio.sleep(NCBI_REQUEST_INTERVAL)

    return record

async def fetch_taxonomy_species(taxon_id):
    """
    """
    try:
        handle = Entrez.efetch(db="taxonomy", id=taxon_id, rettype="null")
        record = Entrez.read(handle)
        handle.close()
    except Exception as e:
        return []
    
    taxid = ''
    for r in record:
        
        if r['Rank'] == 'species':
            taxid = r['TaxId']
        
        else:
            lineage_list = r['LineageEx']

            while lineage_list:
                line_entry = lineage_list.pop()
                if line_entry['Rank'] == 'species':
                    taxid = line_entry['TaxId']
                    break
        
    if taxid:
        return taxid

async def fetch_taxonomy_rank(taxon_id):
    """
    """
    try:
        handle = Entrez.efetch(db="taxonomy", id=taxon_id, rettype="null")
        record = Entrez.read(handle)
        handle.close()
    except Exception as e:
        return []
    
    taxid = ''
    for r in record:
        taxid = r['Rank']
    return taxid

async def fetch_upstream_record_taxids(fetch_list: list) -> list:
    """
    Take a list of accession numbers and request the records from NCBI GenBank
    
    :param fetch_list: List of accession numbers to fetch from GenBank
    :param logger: Structured logger

    :return: A list of GenBank data converted from XML to dicts if possible, 
        else an empty list
    """
    try:
        handle = Entrez.efetch(db="nucleotide", id=fetch_list, rettype="docsum")
        record = Entrez.read(handle)
        handle.close()
    except Exception as e:
        return []
    
    taxids = set()

    for r in record:
        taxids.add(int(r.get('TaxId')))
    
    await asyncio.sleep(NCBI_REQUEST_INTERVAL)
    
    return list(taxids)

async def fetch_records(fetch_list: list) -> list:
    """
    Take a list of accession numbers and request the records from NCBI GenBank
    
    :param fetch_list: List of accession numbers to fetch from GenBank
    :param logger: Structured logger

    :return: A list of GenBank data converted from XML to dicts if possible, 
        else an empty list
    """
    try:
        handle = Entrez.efetch(
            db="nucleotide", id=fetch_list, rettype="gb", retmode="xml"
        )
        ncbi_records = Entrez.parse(handle)
        handle.close()
    except HTTPError as e:
        return []
    
    # ncbi_records = SeqIO.to_dict(SeqIO.parse(handle, "gb"))
    # handle.close()

    if ncbi_records is None:
        return []
    
    try:
        record_list = [record for record in ncbi_records.values() if record.seq]
        return record_list
    except Exception as e:
        return []

async def get_spelling(name, db='taxonomy'):
    urlsafe_name = name.replace(' ', '+')
    try:
        handle = Entrez.espell(db=db, term=urlsafe_name)
        record = Entrez.read(handle)
        handle.close()
    except Exception as e:
        return e
    
    return record['CorrectedQuery']

# async def process_records(records, listing: dict):
#     """
#     """
#     otu_updates = []
#     for seq_list in records:
#         for seq_data in seq_list:
#             [ accession, version ] = (seq_data.id).split('.')
            
#             if accession in filter_set:
#                 continue

#             otu_updates.

#     return otu_updates