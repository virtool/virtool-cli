import os
from pathlib import Path
from Bio import Entrez

from structlog import get_logger
from urllib.parse import quote_plus
from urllib.error import HTTPError

from virtool_cli.ncbi.error import IncompleteRecordsError, NCBIParseError
from virtool_cli.ncbi.utils import parse_nuccore
from virtool_cli.ncbi.cache import NCBICache
from virtool_cli.repo.cls import Repo, RepoOTU

Entrez.email = os.environ.get("NCBI_EMAIL")
Entrez.api_key = os.environ.get("NCBI_API_KEY")

DEFAULT_INTERVAL = 0.001

base_logger = get_logger()


class NCBIClient:
    def __init__(self, cache_path: Path):
        self.cache = NCBICache(cache_path)

    @classmethod
    def for_repo(cls, repo: Repo):
        """Initialize NCBIClient from a repo path"""
        return NCBIClient(repo.path / ".cache/ncbi")

    async def fetch_taxon_records(self, taxon_id: int) -> list[dict]:
        """Fetch all records linked to a taxonomy record.

        :param taxon_id: Taxonomy UID
        :return: A list of records
        """
        accessions = await self.link_accessions(taxon_id)

        base_logger.debug(
            "Fetching accessions...", taxid=taxon_id, accessions=accessions
        )

        records = await self.fetch_accessions(accessions)

        return records

    async def fetch_otu_updates(self, otu: RepoOTU, use_cached: bool = True):
        """Fetch updates for an extant OTU.
        Excludes blocked accessions automatically.

        :param otu: OTU data from repo
        :param use_cached: Cache use flag
        """
        logger = base_logger.bind(otu_id=otu.id)
        if use_cached:
            records = self.cache.load_records(otu.id)
            if records:
                logger.info("Cached records found", n_records=len(records))
                return records

        taxid_accessions = await self.link_accessions(otu.taxid)

        new_accessions = self.filter_accessions(otu, taxid_accessions)

        logger.debug("Fetching accessions...", new_accessions=new_accessions)

        if new_accessions:
            records = await self.fetch_accessions(list(new_accessions))

            return records

    @staticmethod
    def filter_accessions(otu: RepoOTU, accessions: list | set) -> list:
        """
        Takes a list of accessions and removes blocked accessions
        """
        accession_set = set(accessions)
        blocked_set = set(otu.blocked_accessions)

        return list(blocked_set.difference(accession_set))

    @staticmethod
    def process_record(record: dict):
        try:
            return parse_nuccore(record)
        except NCBIParseError as e:
            base_logger.error(f"Parse failure: {e}")

    @staticmethod
    def process_records(records: list[dict]):
        clean_records = []

        for record in records:
            try:
                clean_records.append(parse_nuccore(record))

            except NCBIParseError as e:
                base_logger.error(f"Parse failure: {e}")

        return clean_records

    @staticmethod
    async def link_accessions(taxon_id: int) -> list:
        """
        Requests a cross-reference for NCBI Taxonomy and Nucleotide via ELink
        and returns the results as a list.

        :param taxon_id: A NCBI Taxonomy ID
        :return: A list of accessions linked to the Taxon Id
        """
        elink_results = Entrez.read(
            Entrez.elink(
                dbfrom="taxonomy",
                db="nuccore",
                id=str(taxon_id),
                idtype="acc",
            )
        )

        if not elink_results:
            return []

        # Discards unneeded tables and formats needed table as a list
        for link_set_db in elink_results[0]["LinkSetDb"]:
            if link_set_db["LinkName"] == "taxonomy_nuccore":
                id_table = link_set_db["Link"]

                return [keypair["Id"] for keypair in id_table]

    @staticmethod
    async def fetch_accessions(accessions: list[str]) -> list[dict]:
        """
        Take a list of accession numbers, download the corresponding records
        from GenBank as XML and return the parsed records

        :param accessions: List of accession numbers to fetch from GenBank
        :return: A list of deserialized sequence records from NCBI Nucleotide
        """
        if not accessions:
            return []

        logger = base_logger.bind(accessions=accessions)

        try:
            records = NCBIClient._fetch_serialized_records(accessions)
            return records

        except IncompleteRecordsError as e:
            logger.error(e.message)

            if e.data:
                logger.debug("Partial results fetched, returning results...")
                return e.data

        except HTTPError as e:
            if e.code == 400:
                logger.error(f"{e}. Bad accessions?")
            else:
                logger.error(e)

        return []

    @staticmethod
    async def fetch_accession(accession: str) -> dict:
        """
        A wrapper for the fetching of a single accession
        """
        record = await NCBIClient.fetch_accessions([accession])

        if record:
            return record[0]

    @staticmethod
    def _fetch_serialized_records(accessions: list) -> list[dict] | None:
        """
        Requests XML GenBank records for a list of accessions
        and returns an equal-length list of serialized records.

        Raises an error if fewer records are fetched than accessions.

        :param accessions: A list of n accessions
        :return: A list of n deserialized records
        """
        with Entrez.efetch(
            db="nuccore", id=accessions, rettype="gb", retmode="xml"
        ) as f:
            records = Entrez.read(f)

        # Handle cases where not all accessions can be fetched
        if len(records) == len(accessions):
            return records

        raise IncompleteRecordsError("Bad accession in list", data=records)

    @staticmethod
    def _fetch_raw_records(accessions: list) -> str:
        """
        Requests XML GenBank records for a list of accessions
        and returns results as unparsed XML

        :param accessions: A list of accessions
        :return: XML data as an unparsed string
        """
        with Entrez.efetch(
            db="nuccore", id=accessions, rettype="gb", retmode="xml"
        ) as f:
            raw_records = f.read()

        return raw_records

    @staticmethod
    async def fetch_taxonomy(taxon_id: int) -> dict:
        """Requests a taxonomy record from NCBI Taxonomy"""
        return await NCBIClient._fetch_taxon_long(taxon_id)

    @staticmethod
    async def fetch_taxonomy_id_by_name(name: str) -> int | None:
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

    @staticmethod
    async def _fetch_taxon_docsum(taxon_id: int):
        record = Entrez.read(
            Entrez.efetch(
                db="taxonomy",
                id=taxon_id,
                rettype="docsum",
                retmode="xml",
            )
        )

        return record[0]

    @staticmethod
    async def _fetch_taxon_long(taxon_id: int) -> dict:
        with Entrez.efetch(db="taxonomy", id=taxon_id, rettype="null") as f:
            record = Entrez.read(f)

        return record[0]

    @staticmethod
    async def fetch_taxon_rank(taxon_id: int) -> str:
        taxonomy = await NCBIClient._fetch_taxon_docsum(taxon_id)

        return taxonomy["Rank"]

    async def fetch_species_taxid(self, taxid: int) -> int | None:
        """Gets the species taxid for the given lower-rank taxid.

        :param taxid: NCBI Taxonomy UID
        :return: The NCBI Taxonomy ID of the OTU's species
        """
        taxonomy = await self._fetch_taxon_long(taxid)

        if taxonomy["Rank"] == "species":
            return int(taxonomy["TaxId"])

        for line in taxonomy["LineageEx"]:
            if line["Rank"] == "species":
                return int(line["TaxId"])

        return None

    @staticmethod
    async def check_spelling(name: str, db: str = "taxonomy") -> str:
        """Takes the name of an OTU, requests an alternative spelling
        from the Entrez ESpell utility and returns the suggestion

        :param name: The OTU name that requires correcting
        :param db: Database to check against. Defaults to 'taxonomy'.
        :return: String containing NCBI-suggested spelling changes
        """
        with Entrez.espell(db=db, term=quote_plus(name)) as f:
            record = Entrez.read(f)

        if record:
            return record["CorrectedQuery"]

        return name
