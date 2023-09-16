import os
import asyncio
from urllib.error import HTTPError
from Bio import Entrez, SeqIO
from structlog import BoundLogger

Entrez.email = os.environ.get("NCBI_EMAIL")
Entrez.api_key = os.environ.get("NCBI_API_KEY")

NCBI_REQUEST_INTERVAL = 0.3 if Entrez.email and Entrez.api_key else 0.8


async def fetch_upstream_accessions(
    listing: dict,
    logger: BoundLogger
) -> list:
    """
    Requests a list of all uninspected accessions associated with an OTU's taxon ID
    
    :param listing: Corresponding catalog listing for this OTU
    :return: A list of accessions from NCBI Genbank for the taxon ID, 
        sans included and excluded accessions
    """
    taxid = listing.get('taxid')
    included_set = set(listing['accessions']['included'])
    excluded_set = set(listing['accessions']['excluded'])

    logger = logger.bind(taxid=taxid)
    logger.debug(
        'Exclude catalogued accessions', 
        included=included_set, excluded=excluded_set)

    upstream_accessions = []

    # Request results as accessions, not UIDs
    entrez_acclist = Entrez.read(
        Entrez.elink(
            dbfrom="taxonomy", db="nucleotide", 
            id=str(taxid), idtype="acc")
        )

    for linksetdb in entrez_acclist[0]["LinkSetDb"][0]["Link"]:
        accession = linksetdb["Id"]
        if accession.split('.')[0] not in excluded_set:
            upstream_accessions.append(accession)

    upstream_set = set(upstream_accessions)

    return list(upstream_set.difference(included_set))

async def fetch_upstream_records(
    fetch_list: list, 
    logger: BoundLogger
) -> list:
    """
    Take a list of accession numbers and request the records from NCBI GenBank
    
    :param fetch_list: List of accession numbers to fetch from GenBank
    :param logger: Structured logger
    :return: A list of GenBank data converted from XML to dicts if possible, 
        else an empty list
    """
    try:
        handle = Entrez.efetch(
            db="nucleotide", id=fetch_list, rettype="gb", retmode="text"
        )
    except HTTPError as e:
        logger.error(f'{e}, moving on...')
        return []
    
    ncbi_records = SeqIO.to_dict(SeqIO.parse(handle, "gb"))
    handle.close()

    if ncbi_records is None:
        return []
    
    try:
        accession_list = [record for record in ncbi_records.values() if record.seq]
        return accession_list
    except Exception as e:
        logger.exception(e)
        raise e

async def fetch_nuccore(fetch_list: list) -> list:
    """
    Take a list of accession numbers and request the corresponding records from NCBI Nucleotide
    
    :param fetch_list: List of accession numbers to fetch from GenBank

    :return: A list of GenBank data converted from XML to dicts if possible, 
        else an empty list
    """
    try:
        handle = Entrez.efetch(
            db="nucleotide", id=fetch_list, 
            rettype="gb", retmode="text"
        )
        ncbi_records = SeqIO.to_dict(SeqIO.parse(handle, "gb"))
        handle.close()
    except HTTPError as e:
        raise e

    if ncbi_records is None:
        return []
    
    try:
        accession_list = [record for record in ncbi_records.values() if record.seq]
        return accession_list
    except Exception as e:
        raise e

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
        raise e

    try:
        taxid = int(record["IdList"][0])
    except IndexError:
        taxid = None

    return taxid

async def fetch_taxonomy_species(taxon_id: str):
    """
    Fetches a record from NCBI Taxonomy and extracts the species taxid
    from its lineage list.

    :param taxon_id: NCBI Taxonomy UID
    :return: The NCBI Taxonomy ID of the OTU's species
    """
    try:
        handle = Entrez.efetch(db="taxonomy", id=taxon_id, rettype="null")
        record = Entrez.read(handle)
        handle.close()
    except HTTPError:
        return -1
    
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
        return int(taxid)

async def fetch_taxonomy_rank(taxon_id) -> str:
    """
    Fetches a record from NCBI Taxonomy and extracts its taxonomic rank

    :param taxon_id: NCBI Taxonomy UID
    :return: The taxonomic rank of this Taxonomy record as a string
    """
    try:
        handle = Entrez.efetch(db="taxonomy", id=taxon_id, rettype="null")
        record = Entrez.read(handle)
        handle.close()
    except Exception:
        return []
    
    rank = ''
    for r in record:
        rank = r['Rank']
    return rank

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
    except Exception:
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
    except HTTPError:
        return []

    if ncbi_records is None:
        return []
    
    try:
        record_list = [record for record in ncbi_records.values() if record.seq]
        return record_list
    except Exception:
        return []


async def get_spelling(name: str, db: str = 'taxonomy') -> str:
    """
    Takes the name of an OTU, requests an alternative spelling
    from the Entrez ESpell utility and returns the suggestion

    :param name: The OTU name that requires correcting
    :param db: NCBI Database to check against. Defaults to 'taxonomy'.
    :return: String containing NCBI-suggested spelling changes
    """
    urlsafe_name = name.replace(' ', '+')
    try:
        handle = Entrez.espell(db=db, term=urlsafe_name)
        record = Entrez.read(handle)
        handle.close()
    except Exception as e:
        raise e
    
    return record['CorrectedQuery']