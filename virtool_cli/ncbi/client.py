import os
import asyncio
from Bio import Entrez

from structlog import get_logger, BoundLogger
from http.client import IncompleteRead
from urllib.error import HTTPError
from urllib.parse import quote_plus

DEFAULT_INTERVAL = 0.001


class NCBIClient:
    def __init__(
        self,
        email: str | None = os.environ.get("NCBI_EMAIL"),
        api_key: str | None = os.environ.get("NCBI_API_KEY"),
        logger: BoundLogger = get_logger(),
    ):
        Entrez.email = email
        Entrez.api_key = api_key

        self.pause = 0.3 if email and api_key else 0.8

        self.logger = logger

    async def link_accessions(self, taxon_id: int) -> list:
        """
        Requests a cross-reference for NCBI Taxonomy and Nucleotide via ELink
        and returns the results as a list.

        :param taxon_id: A NCBI Taxonomy ID
        :return: A list of accessions linked to the Taxon Id
        """
        try:
            elink_results = Entrez.read(
                Entrez.elink(
                    dbfrom="taxonomy",
                    db="nuccore",
                    id=str(taxon_id),
                    idtype="acc",
                )
            )
        except IncompleteRead:
            raise RuntimeError("IncompleteRead")
        except HTTPError as e:
            raise e
        await asyncio.sleep(self.pause)

        if not elink_results:
            return []

        # Discards unneeded tables and formats needed table as a list
        for link_set_db in elink_results[0]["LinkSetDb"]:
            if link_set_db["LinkName"] == "taxonomy_nuccore":
                id_table = link_set_db["Link"]

                return [keypair["Id"] for keypair in id_table]

            else:
                raise ValueError("Retrieved incompatible link data")

    async def fetch_accessions(self, accessions: list) -> list:
        """
        Take a list of accession numbers, download the corresponding records from GenBank as XML
        and return the parsed records

        :param accessions: List of accession numbers to fetch from GenBank
        :return: A list of deserialized sequence records from NCBI Nucleotide
        """
        if not accessions:
            return []

        try:
            ncbi_records = await self._fetch_serialized_records(accessions)
        except (HTTPError, IncompleteRead) as e:
            raise e

        await asyncio.sleep(self.pause)

        if ncbi_records is None:
            return []

        record_list = ncbi_records

        return record_list

    async def fetch_accession(self, accession: str) -> dict:
        """
        A wrapper for the fetching of a single accession
        """
        record = await self.fetch_accessions([accession])

        if record:
            return record.pop()

    async def _fetch_serialized_records(self, accessions: list) -> list[dict]:
        """
        Requests XML GenBank records for a list of accessions and returns serialized results

        :param accessions: A list of accessions
        :return: A dictionary of accession-record key-pairs
        """
        with Entrez.efetch(
            db="nuccore", id=accessions, rettype="gb", retmode="xml"
        ) as f:
            records = Entrez.read(f)

        await asyncio.sleep(self.pause)

        return records

    async def _fetch_raw_records(self, accessions: list) -> str:
        """
        Requests XML GenBank records for a list of accessions and returns results as unparsed XML

        :param accessions: A list of accessions
        :return: A dictionary of accession-record key-pairs
        """
        with Entrez.efetch(
            db="nuccore", id=accessions, rettype="gb", retmode="xml"
        ) as f:
            raw_records = f.read()

        await asyncio.sleep(self.pause)

        return raw_records

    async def fetch_taxonomy(self, taxon_id: int, long=False) -> dict:
        """Requests a taxonomy record from NCBI Taxonomy"""
        try:
            if long:
                taxonomy = await self._fetch_taxon_long(taxon_id)

            else:
                taxonomy = await self._fetch_taxon_docsum(taxon_id)

            await asyncio.sleep(self.pause)

        except IncompleteRead:
            raise RuntimeError("IncompleteRead")
        except HTTPError:
            raise HTTPError

        return taxonomy

    async def fetch_taxonomy_id_by_name(self, name: str) -> int | None:
        """Returns a best-guess taxon ID for a given OTU name.

        Searches the NCBI taxonomy database for a given OTU name, then fetches and
        returns its taxonomy.

        Returns ``None`` if no matching taxonomy is found.

        :param name: the name of an otu
        :return: The taxonomy id for the given otu name
        """
        with Entrez.esearch(db="taxonomy", term=name) as f:
            record = Entrez.read(f)

        try:
            taxid = int(record["IdList"][0])
        except IndexError:
            return None

        return taxid

    async def _fetch_taxon_docsum(self, taxon_id):
        record = Entrez.read(
            Entrez.efetch(
                db="taxonomy",
                id=taxon_id,
                rettype="docsum",
                retmode="xml",
            )
        )

        return record[0]

    async def _fetch_taxon_long(self, taxon_id) -> dict:
        with Entrez.efetch(db="taxonomy", id=taxon_id, rettype="null") as f:
            record = Entrez.read(f)

        return record[0]

    async def fetch_taxon_rank(self, taxon_id: int) -> str:
        with Entrez.efetch(
            db="taxonomy", id=taxon_id, rettype="docsum", retmode="xml"
        ) as f:
            for r in Entrez.parse(f):
                return r["Rank"]

    async def fetch_species_taxid(self, taxid: str) -> int | None:
        """Gets the species taxid for the given lower-rank taxid.

        :param taxid: NCBI Taxonomy UID
        :return: The NCBI Taxonomy ID of the OTU's species
        """
        with Entrez.efetch(db="taxonomy", id=taxid, rettype="null") as f:
            record = Entrez.parse(f)

            for r in record:
                if r["Rank"] == "species":
                    return int(r["TaxId"])

                for line in r["LineageEx"]:
                    if line["Rank"] == "species":
                        return int(line["TaxId"])

        return None

    @staticmethod
    async def check_spelling(name: str, db: str = "taxonomy") -> str:
        """Takes the name of an OTU, requests an alternative spelling
        from the Entrez ESpell utility and returns the suggestion

        :param name: The OTU name that requires correcting
        :param db: NCBI Database to check against. Defaults to 'taxonomy'.
        :return: String containing NCBI-suggested spelling changes
        """
        with Entrez.espell(db=db, term=quote_plus(name)) as f:
            record = Entrez.read(f)

        if record:
            return record["CorrectedQuery"]

        return name
