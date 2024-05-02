import os
from http.client import IncompleteRead
from urllib.error import HTTPError

from Bio import Entrez, SeqIO

Entrez.email = os.environ.get("NCBI_EMAIL")
Entrez.api_key = os.environ.get("NCBI_API_KEY")

NCBI_REQUEST_INTERVAL = 0.3 if Entrez.email and Entrez.api_key else 0.8


async def request_linked_accessions(taxon_id: int) -> list:
    """Take an NCBI Taxonomy UID and return a list of linked accessions from the Nucleotide database

    :param taxon_id: NCBI Taxonomy UID as an integer
    :return: List of accessions that are linked to the input taxon ID
    """
    upstream_accessions = []

    try:
        elink_results = Entrez.elink(
            dbfrom="taxonomy",
            db="nuccore",
            id=str(taxon_id),
            idtype="acc",
        )
    except IncompleteRead:
        raise HTTPError("IncompleteRead")
    except HTTPError as e:
        raise e

    try:
        entrez_acclist = Entrez.read(elink_results)
    except RuntimeError as e:
        raise e

    if not entrez_acclist:
        return []

    for linksetdb in entrez_acclist[0]["LinkSetDb"][0]["Link"]:
        upstream_accessions.append(str(linksetdb["Id"]))

    return sorted(upstream_accessions)


async def request_from_nucleotide(fetch_list: list) -> list:
    """Take a list of accession numbers and request the corresponding records from NCBI Nucleotide

    :param fetch_list: List of accession numbers to fetch from GenBank
    :return: A list of GenBank data converted from GenBank entries to dicts if possible,
        else an empty list
    """
    try:
        handle = Entrez.efetch(
            db="nuccore",
            id=fetch_list,
            rettype="gb",
            retmode="text",
        )
        ncbi_records = SeqIO.to_dict(SeqIO.parse(handle, "gb"))
        handle.close()
    except IncompleteRead:
        raise HTTPError("IncompleteRead")
    except HTTPError as e:
        raise e

    if ncbi_records is None:
        return []

    try:
        accession_list = [record for record in ncbi_records.values() if record.seq]
        return accession_list

    except Exception as e:
        raise e


async def fetch_taxonomy_record(taxon_id) -> dict:
    """Fetches a record from NCBI Taxonomy and returns the contents as a dictionary if found

    :param taxon_id: NCBI Taxonomy UID
    :return: The taxonomic rank of this Taxonomy record as a string
    """
    try:
        handle = Entrez.efetch(
            db="taxonomy",
            id=taxon_id,
            rettype="docsum",
            retmode="xml",
        )
        record = Entrez.read(handle)
        handle.close()
    except IncompleteRead:
        raise HTTPError("IncompleteRead")
    except HTTPError:
        raise HTTPError

    if record:
        return record.pop()

    return {}


async def fetch_isolate_metadata(taxid: int) -> dict:
    """:param taxid: NCBI taxon ID in integer form
    :return: Dict of isolate name and type
    """
    taxid_docsum = await fetch_taxonomy_record(str(taxid))

    isolate_type = taxid_docsum.get("Rank", "unknown")
    isolate_name = taxid_docsum.get("ScientificName")

    return {"source_name": isolate_name, "source_type": isolate_type}
