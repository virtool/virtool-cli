import os
from pathlib import Path
from enum import StrEnum
from Bio import Entrez

from structlog import get_logger
from urllib.parse import quote_plus
from urllib.error import HTTPError
from pydantic import ValidationError

from virtool_cli.ncbi.error import IncompleteRecordsError, NCBIParseError
from virtool_cli.ncbi.model import NCBINuccore, NCBISource
from virtool_cli.ncbi.cache import NCBICache

Entrez.email = os.environ.get("NCBI_EMAIL")
Entrez.api_key = os.environ.get("NCBI_API_KEY")

base_logger = get_logger()


class NCBIClient:
    def __init__(self, cache_path: Path):
        """
        :param cache_path: A path to a directory to be used as a cache
        """
        self.cache = NCBICache(cache_path)

    @classmethod
    def from_repo(cls, repo_path: Path) -> Repo:
        """Initializes the NCBI cache in the default subpath
        under a given repository

        :param repo_path: A path to a reference repository
        :return:
        """
        return NCBIClient(repo_path / ".cache/ncbi")

    @staticmethod
    async def procure_accessions(
        requested: list[str],
        blocked: list[str] | None = None,
    ) -> list[NCBINuccore]:
        """
        Filter an accession list, then fetch and validate NCBI Genbank records

        :TODO: Get rid of filtering. Just fetch records by a list of accessions.
        :TODO: Cache on accession per file for easier lookup in `NCBICache`.
        :TODO: Merge other accession fetching, caching, updating methods into this one and call this method
               `fetch_accessions`.

        :param requested:
        :param blocked:
        """
        if blocked is None:
            blocked = []
        logger = base_logger.bind(requested=requested, blocked=blocked)

        accessions = NCBIClient.filter_accessions(requested, blocked)
        if not accessions:
            logger.warning("No new accessions")
            return []

        records = await NCBIClient.fetch_by_accessions(accessions)

        if records:
            return NCBIClient.validate_records(records)
        else:
            return []

    async def procure_from_taxid(
        self,
        taxid: int,
        blocked_accessions: list[str] | None = None,
        use_cached: bool = True,
    ) -> list[NCBINuccore]:
        """
        Fetch all linked accessions for an organism in NCBI Taxonomy
        and return results as a set of NCBIAccession and NCBISource

        :TODO: Merge cache, fetch, procure into one method here. Rename.
        :TODO: Return a list of plain associated accessions. Leaving subsequent fetch of full Genbank records to the
               caller.

        :param taxid: NCBI Taxonomy UID as an integer
        :param blocked_accessions:
        :param use_cached: Cache use flag
        :return: A list of NCBINuccore parsed records
        """
        logger = base_logger.bind(taxid=taxid)
        if use_cached:
            records = self.cache.load_nuccore(str(taxid))
            if records:
                logger.info("Cached records found", n_records=len(records))

        all_accessions = await NCBIClient.link_accessions(taxid)

        return await NCBIClient.procure_accessions(
            requested=all_accessions, blocked=blocked_accessions
        )

    async def procure_updates(
        self,
        otu_id: str,
        taxid: int,
        blocked_accessions: list[str],
        use_cached: bool = True,
    ) -> list[NCBINuccore]:
        """
        Fetch updates for an extant OTU and return results as
        a list of NCBINuccore validated records
        Excludes blocked accessions automatically.

        :param otu_id: OTU data from repo
        :param taxid: NCBI Taxonomy UID as an integer
        :param blocked_accessions: A list of accessions to exclude from this operation
        :param use_cached: Cache use flag
        :return: A list of parsed records
        """
        logger = base_logger.bind(otu_id=otu_id)

        if use_cached:
            records = self.cache.load_nuccore(otu_id)
            if records:
                logger.info("Cached records found", n_records=len(records))

                return NCBIClient.validate_records(records)

        taxid_accessions = await NCBIClient.link_accessions(taxid)

        new_accessions = NCBIClient.filter_accessions(
            blocked_accessions, taxid_accessions
        )

        if new_accessions:
            logger.debug("Fetching accessions...", new_accessions=new_accessions)

            records = await NCBIClient.fetch_by_accessions(list(new_accessions))

            return NCBIClient.validate_records(records)

    async def cache_from_taxid(
        self, taxid: int, blocked_accessions: list[str] | None = None
    ):
        """Fetch all linked accessions for an organism in NCBI Taxonomy
        and cache the results

        :param taxid: NCBI Taxonomy UID as an integer
        :param blocked_accessions:
        """
        all_accessions = await NCBIClient.link_accessions(taxid)

        if blocked_accessions is None:
            blocked_accessions = []

        records = await NCBIClient.fetch_by_accessions(
            NCBIClient.filter_accessions(all_accessions, blocked_accessions)
        )

        self.cache.cache_nuccore(records, str(taxid))

    async def cache_updates(
        self,
        otu_id: str,
        taxid: int,
        blocked_accessions: list[str],
        use_cached: bool = True,
    ):
        """Fetch and cache updates for an extant OTU.
        Excludes blocked accessions automatically.

        :param otu_id: OTU data from repo
        :param taxid: NCBI Taxonomy UID as an integer
        :param blocked_accessions: A list of accessions to exclude from this operation
        :param use_cached: Cache use flag
        """
        logger = base_logger.bind(otu_id=otu_id)

        if use_cached:
            records = self.cache.load_nuccore(otu_id)
            if records:
                logger.info("Cached records found", n_records=len(records))
                return records

        taxid_accessions = await NCBIClient.link_accessions(taxid)

        new_accessions = NCBIClient.filter_accessions(
            blocked_accessions, taxid_accessions
        )

        logger.debug("Fetching accessions...", new_accessions=new_accessions)

        if new_accessions:
            records = await NCBIClient.fetch_by_accessions(list(new_accessions))

            self.cache.cache_nuccore(records, otu_id)

            logger.debug("Cached records", n_records=len(records))

    @staticmethod
    async def fetch_by_taxid(taxid: int) -> list[dict]:
        """Fetch all records linked to a taxonomy record.
        Usable when an OTU does not exist

        :param taxid: A Taxonomy UID
        :return: A list of Entrez-parsed records
        """
        accessions = await NCBIClient.link_accessions(taxid)

        base_logger.debug("Fetching accessions...", taxid=taxid, accessions=accessions)

        return await NCBIClient.fetch_by_accessions(accessions)

    @staticmethod
    async def fetch_by_accessions(accessions: list[str]) -> list[dict]:
        """
        Take a list of accession numbers, download the corresponding records
        from GenBank as XML and return Genbank XML-parsed records

        :param accessions: List of accession numbers to fetch from GenBank
        :return: A list of deserialized records from NCBI Nucleotide
        """
        if not accessions:
            return []

        logger = base_logger.bind(accessions=accessions)

        try:
            records = NCBIClient.__fetch_serialized_records(accessions)
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
    async def fetch_by_accession(accession: str) -> dict:
        """
        Wrapper that fetches a single accession and returns
        a single NCBI Genbank-parsed record

        :param accession: A single accession
        :return: A single NCBI Genbank-parsed record
        """
        record = await NCBIClient.fetch_by_accessions([accession])

        if record:
            return record[0]

    @staticmethod
    def validate_records(records: list[dict]) -> list[NCBINuccore]:
        """
        Process a list of raw Genbank dicts into validated NCBINuccore records.
        Logs an error if there is an issue with validation or parsing,
        but does not fail out.

        :param records: A list of NCBI Genbank dict records
        :return: A list of validated records as NCBINuccore
        """
        clean_records = []

        for record in records:
            try:
                clean_records.append(NCBIClient.validate_nuccore(record))

            except (ValidationError, NCBIParseError) as exc:
                accession = record.get(GBSeq.ACCESSION, "?")
                base_logger.error(f"{exc}", accession=accession)

        return clean_records

    @staticmethod
    def filter_accessions(new: list | set, blocked: list | set | None = None) -> list:
        """
        Takes a list of accessions and filters blocked accessions from it,
        leaving only new accessions

        :param new: A list of accessions
        :param blocked: A list of blocked accessions
        :return: All accessions in the "new" list that are
            not also in the "blocked" list
        """
        new_set = set(new)
        if blocked is None:
            blocked_set = set()
        else:
            blocked_set = set(blocked)

        return list(new_set.difference(blocked_set))

    @staticmethod
    def validate_nuccore(raw: dict) -> NCBINuccore:
        """
        Parses an NCBI Genbank record from a Genbank dict to a
        validated NCBINuccore

        :TODO: End function with `return NCBISource(**{...})`. Pydantic can take expanded dict and validate and get rid
               of unwanted fields.

        :param raw: A NCBI Genbank dict record, parsed by Bio.Entrez.Parser
        :return: A validated subset of Genbank record data
        """
        source_dict = NCBIClient.__feature_table_to_dict(
            NCBIClient.__get_source_table(raw)
        )

        source = NCBISource(taxid=NCBIClient.__parse_taxid(source_dict))

        for source_key in ("isolate", "strain", "clone", "host", "segment"):
            if source_key in source_dict:
                setattr(source, source_key, source_dict[source_key])

        record = NCBINuccore(
            accession=raw[GBSeq.ACCESSION],
            definition=raw[GBSeq.DEFINITION],
            sequence=raw[GBSeq.SEQUENCE].upper(),
            source=source,
        )

        if GBSeq.COMMENT in raw:
            record.comment = raw[GBSeq.COMMENT]

        return record

    @staticmethod
    async def link_accessions(taxon_id: int) -> list:
        """
        Requests a cross-reference for NCBI Taxonomy and Nucleotide via ELink
        and returns the results as a list.

        :TODO: Merge into caller

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
    def __fetch_serialized_records(accessions: list) -> list[dict] | None:
        """
        Requests XML GenBank records for a list of accessions
        and returns an equal-length list of serialized records.

        Raises an error if fewer records are fetched than accessions.

        :TODO: Merge up into caller. This will make exception handler a bit cleaner.
        :TODO: Change to single underscore for private methods.

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
    def __fetch_raw_records(accessions: list) -> str:
        """
        Requests XML GenBank records for a list of accessions
        and returns results as unparsed XML

        :TODO: Remove? It doesn't appear to be used anywhere.

        :param accessions: A list of accessions
        :return: XML data as an unparsed string
        """
        with Entrez.efetch(
            db="nuccore", id=accessions, rettype="gb", retmode="xml"
        ) as f:
            raw_records = f.read()

        return raw_records

    @staticmethod
    async def fetch_taxonomy_by_taxid(taxon_id: int) -> dict:
        """Requests a taxonomy record from NCBI Taxonomy"""
        return await NCBIClient.__fetch_taxon_long(taxon_id)

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
    async def __fetch_taxon_docsum(taxon_id: int):
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
    async def __fetch_taxon_long(taxon_id: int) -> dict:
        with Entrez.efetch(db="taxonomy", id=taxon_id, rettype="null") as f:
            record = Entrez.read(f)

        return record[0]

    @staticmethod
    async def fetch_taxon_rank(taxon_id: int) -> str:
        """
        :TODO: Merge into one taxonomy method.
        """
        taxonomy = await NCBIClient.__fetch_taxon_docsum(taxon_id)

        return taxonomy["Rank"]

    @staticmethod
    async def fetch_species_taxid(taxid: int) -> int | None:
        """Gets the species taxid for the given lower-rank taxid.

        :TODO: Merge into one taxonomy method.

        :param taxid: NCBI Taxonomy UID
        :return: The NCBI Taxonomy ID of the OTU's species
        """
        taxonomy = await NCBIClient.__fetch_taxon_long(taxid)

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

        :TODO: Use `Enum` for `db`.

        :param name: The OTU name that requires correcting
        :param db: Database to check against. Defaults to 'taxonomy'.
        :return: String containing NCBI-suggested spelling changes
        """
        with Entrez.espell(db=db, term=quote_plus(name)) as f:
            record = Entrez.read(f)

        if record:
            return record["CorrectedQuery"]

        return name

    @staticmethod
    def __get_source_table(raw) -> dict | None:
        for feature in raw[GBSeq.FEATURE_TABLE]:
            if feature["GBFeature_key"] == "source":
                return feature

        raise NCBIParseError(
            keys=NCBIClient.__get_feature_table_keys(raw),
            message="Feature table does not contain source data",
        )

    @staticmethod
    def __feature_table_to_dict(feature: dict) -> dict:
        """Converts the feature table format to a dict"""
        qualifier_dict = {}
        for qualifier in feature["GBFeature_quals"]:
            qual_name = qualifier["GBQualifier_name"]
            qual_value = qualifier["GBQualifier_value"]
            qualifier_dict[qual_name] = qual_value

        return qualifier_dict

    @staticmethod
    def __parse_taxid(source_feature) -> int | None:
        return int(NCBIClient.__split_db_xref(source_feature["db_xref"])[1])

    @staticmethod
    def __split_db_xref(entry: str) -> list[str]:
        xref = entry.split(":")
        if len(xref) != 2:
            raise NCBIParseError

        return xref

    @staticmethod
    def __get_feature_table_keys(raw: dict):
        return [feature["GBFeature_key"] for feature in raw[GBSeq.FEATURE_TABLE]]


class GBSeq(StrEnum):
    ACCESSION = "GBSeq_accession-version"
    DEFINITION = "GBSeq_definition"
    SEQUENCE = "GBSeq_sequence"
    LENGTH = "GBSeq_length"
    COMMENT = "GBSeq_comment"
    FEATURE_TABLE = "GBSeq_feature-table"
