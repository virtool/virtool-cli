import os
from pathlib import Path
from collections import namedtuple
from enum import StrEnum
from Bio import Entrez

from structlog import get_logger
from urllib.parse import quote_plus
from urllib.error import HTTPError

from virtool_cli.ncbi.error import IncompleteRecordsError, NCBIParseError
from virtool_cli.ncbi.model import NCBINuccore

from virtool_cli.ncbi.cache import NCBICache
from virtool_cli.repo.cls import Repo, RepoOTU

Entrez.email = os.environ.get("NCBI_EMAIL")
Entrez.api_key = os.environ.get("NCBI_API_KEY")

DEFAULT_INTERVAL = 0.001

base_logger = get_logger()


NuccorePacket = namedtuple("NuccorePacket", ["sequence", "source"])


class NCBIClient:
    def __init__(self, cache_path: Path):
        self.cache = NCBICache(cache_path)

    @classmethod
    def for_repo(cls, repo: Repo):
        """Initialize NCBIClient from a repo path"""
        return NCBIClient(repo.path / ".cache/ncbi")

    async def procure_from_taxid(
        self, taxid: int, use_cached: bool = True
    ) -> list[NCBINuccore]:
        """
        Fetch all linked accessions for an organism in NCBI Taxonomy
        and return results as a set of NCBIAccession and NCBISource

        :param taxid: NCBI Taxonomy UID as an integer
        :param use_cached: Cache use flag
        """

        logger = base_logger.bind(taxid=taxid)
        if use_cached:
            records = self.cache.load_nuccore(str(taxid))
            if records:
                logger.info("Cached records found", n_records=len(records))

        records = await NCBIClient.fetch_taxon_records(taxid)

        return NCBIClient.process_records(records)

    async def procure_updates(
        self, otu: RepoOTU, use_cached: bool = True
    ) -> list[NCBINuccore]:
        """
        Fetch updates for an extant OTU and return results
        as a set of NCBIAccession and NCBISource
        Excludes blocked accessions automatically.

        :param otu: OTU data from repo
        :param use_cached: Cache use flag
        """
        logger = base_logger.bind(otu_id=otu.id)

        if use_cached:
            records = self.cache.load_nuccore(otu.id)
            if records:
                logger.info("Cached records found", n_records=len(records))

                return NCBIClient.process_records(records)

        taxid_accessions = await NCBIClient.link_accessions(otu.taxid)

        new_accessions = NCBIClient.filter_accessions(otu, taxid_accessions)

        logger.debug("Fetching accessions...", new_accessions=new_accessions)

        if new_accessions:
            records = await NCBIClient.fetch_accessions(list(new_accessions))

            return NCBIClient.process_records(records)

    async def cache_from_taxid(self, taxid: int):
        records = await NCBIClient.fetch_taxon_records(taxid)

        self.cache.cache_nuccore(records, str(taxid))

    async def cache_updates(self, otu: RepoOTU, use_cached: bool = True):
        """Fetch and cache updates for an extant OTU.
        Excludes blocked accessions automatically.

        :param otu: OTU data from repo
        :param use_cached: Cache use flag
        """
        logger = base_logger.bind(otu_id=otu.id)

        if use_cached:
            records = self.cache.load_nuccore(otu.id)
            if records:
                logger.info("Cached records found", n_records=len(records))
                return records

        taxid_accessions = await NCBIClient.link_accessions(otu.taxid)

        new_accessions = NCBIClient.filter_accessions(otu, taxid_accessions)

        logger.debug("Fetching accessions...", new_accessions=new_accessions)

        if new_accessions:
            records = await NCBIClient.fetch_accessions(list(new_accessions))

            self.cache.cache_nuccore(records, otu.id)

            logger.debug("Cached records", n_records=len(records))

    @staticmethod
    async def fetch_taxon_records(taxid: int) -> list[dict]:
        """Fetch all records linked to a taxonomy record.

        :param taxid: Taxonomy UID
        :return: A list of records
        """
        accessions = await NCBIClient.link_accessions(taxid)

        base_logger.debug("Fetching accessions...", taxid=taxid, accessions=accessions)

        records = await NCBIClient.fetch_accessions(accessions)

        return records

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
    def process_records(records: list[dict]) -> list[NCBINuccore]:
        clean_records = []

        for record in records:
            try:
                clean_records.append(parse_nuccore(record))

            except NCBIParseError as e:
                base_logger.error(f"Parse failure: {e}")

        return clean_records

    @staticmethod
    def process_record(record: dict) -> NCBINuccore | None:
        try:
            return parse_nuccore(record)
        except NCBIParseError as e:
            base_logger.error(f"Parse failure: {e}")

            return None

    @staticmethod
    def filter_accessions(otu: RepoOTU, accessions: list | set) -> list:
        """
        Takes a list of accessions and removes blocked accessions
        """
        accession_set = set(accessions)
        blocked_set = set(otu.blocked_accessions)

        return list(blocked_set.difference(accession_set))

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

    @staticmethod
    async def fetch_species_taxid(taxid: int) -> int | None:
        """Gets the species taxid for the given lower-rank taxid.

        :param taxid: NCBI Taxonomy UID
        :return: The NCBI Taxonomy ID of the OTU's species
        """
        taxonomy = await NCBIClient._fetch_taxon_long(taxid)

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


class GBSeq(StrEnum):
    ACCESSION = "GBSeq_accession-version"
    DEFINITION = "GBSeq_definition"
    SEQUENCE = "GBSeq_sequence"
    LENGTH = "GBSeq_length"
    COMMENT = "GBSeq_comment"
    FEATURE_TABLE = "GBSeq_feature-table"


def parse_nuccore(raw: dict) -> NCBINuccore:
    source_dict = feature_table_to_dict(get_source_table(raw))

    record = NCBINuccore(
        accession=raw[GBSeq.ACCESSION],
        definition=raw[GBSeq.DEFINITION],
        sequence=raw[GBSeq.SEQUENCE],
        taxid=parse_taxid(source_dict),
    )

    if GBSeq.COMMENT in raw:
        record.comment = raw[GBSeq.COMMENT]

    if "isolate" in source_dict:
        record.isolate = source_dict["isolate"]

    if "strain" in source_dict:
        record.strain = source_dict["strain"]

    if "clone" in source_dict:
        record.strain = source_dict["clone"]

    if "host" in source_dict:
        record.isolate = source_dict["host"]

    if "segment" in source_dict:
        record.segment = source_dict["segment"]

    return record


def get_source_table(raw) -> dict | None:
    for feature in raw[GBSeq.FEATURE_TABLE]:
        if feature["GBFeature_key"] == "source":
            return feature

    raise NCBIParseError(
        keys=get_feature_table_keys(raw),
        message="Feature table does not contain source data",
    )


def get_feature_table_keys(raw):
    keys = []
    for feature in raw[GBSeq.FEATURE_TABLE]:
        keys.append(feature["GBFeature_key"])

    return keys


def feature_table_to_dict(feature: dict) -> dict:
    """Converts the feature table format to a dict"""
    qualifier_dict = {}
    for qualifier in feature["GBFeature_quals"]:
        qual_name = qualifier["GBQualifier_name"]
        qual_value = qualifier["GBQualifier_value"]
        qualifier_dict[qual_name] = qual_value

    return qualifier_dict


def parse_taxid(source_feature) -> int | None:
    value = source_feature["db_xref"]

    return value.split(":")[1]
