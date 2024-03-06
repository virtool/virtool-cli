import os
from pathlib import Path
from enum import StrEnum
from Bio import Entrez

from structlog import get_logger
from urllib.parse import quote_plus
from urllib.error import HTTPError
from pydantic import ValidationError

from virtool_cli.ncbi.error import NCBIParseError
from virtool_cli.ncbi.model import (
    NCBINuccore,
    NCBISource,
    NCBITaxonomy,
    NCBIDB,
    NCBILineage,
    NCBIRank,
)
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
    def from_repo(cls, repo_path: Path):
        """Initializes the NCBI cache in the default subpath
        under a given repository

        :param repo_path: A path to a reference repository
        :return:
        """
        return NCBIClient(repo_path / ".cache/ncbi")

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

    async def fetch_accessions(
        self,
        accessions: list[str],
        cache_results: bool = True,
        use_cached: bool = True,
    ) -> list[NCBINuccore]:
        """
        Fetch records corresponding to a list of accessions and and validate fetched NCBI Genbank records

        :param accessions: A list of accessions to be fetched
        :param cache_results: If True, caches fetched data as JSON
        :param use_cached: If True, loads data from cache in lieu of fetching if possible
        :return: List of validated records as NCBICNuccore
        """
        if not accessions:
            return []

        logger = base_logger.bind(accessions=accessions)

        if use_cached:
            records = self.cache.load_nuccore_records(accessions)
            if records:
                logger.info("Cached records found")
                return NCBIClient.validate_records(records)

        records = await NCBIClient.fetch_unvalidated_accessions(accessions)
        if records:
            if cache_results:
                self.cache.cache_nuccore_records(records)
                logger.debug("Records cached", n_records=len(records))

            return NCBIClient.validate_records(records)

        return []

    @staticmethod
    async def fetch_unvalidated_accessions(accessions: list[str]) -> list[dict]:
        """
        Take a list of accession numbers, parse the corresponding XML records
        from GenBank using Entrez.Parser and return

        :param accessions: List of accession numbers to fetch from GenBank
        :return: A list of deserialized XML records from NCBI Nucleotide
        """
        logger = base_logger.bind(accessions=accessions)

        try:
            with Entrez.efetch(
                db="nuccore", id=accessions, rettype="gb", retmode="xml"
            ) as f:
                records = Entrez.read(f)

            # Handle cases where not all accessions can be fetched
            if len(records) != len(accessions):
                logger.debug("Partial results fetched, returning results...")

            return records

        except HTTPError as e:
            if e.code == 400:
                logger.error(f"{e}. Bad accessions?")
            else:
                logger.error(e)

        return []

    @staticmethod
    async def fetch_by_taxid(taxid: int) -> list[dict]:
        """Fetch all records linked to a taxonomy record.
        Usable when an OTU does not exist

        :param taxid: A Taxonomy UID
        :return: A list of Entrez-parsed records
        """
        accessions = await NCBIClient.link_accessions(taxid)

        base_logger.debug("Fetching accessions...", taxid=taxid, accessions=accessions)

        return await NCBIClient.fetch_unvalidated_accessions(accessions)

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
                clean_records.append(NCBIClient.validate_record(record))

            except (ValidationError, NCBIParseError) as exc:
                accession = record.get(GBSeq.ACCESSION, "?")
                base_logger.error(f"{exc}", accession=accession)

        return clean_records

    @staticmethod
    def validate_record(raw: dict) -> NCBINuccore:
        """
        Parses an NCBI Genbank record from a Genbank dict to
        a validated NCBINuccore

        :TODO: End function with `return NCBISource(**{...})`. Pydantic can take expanded dict and validate and get rid
               of unwanted fields.

        :param raw: A NCBI Genbank dict record, parsed by Bio.Entrez.Parser
        :return: A validated subset of Genbank record data
        """
        source_dict = NCBIClient._get_source_dict(raw)

        return NCBINuccore(
            accession=raw[GBSeq.ACCESSION],
            definition=raw[GBSeq.DEFINITION],
            sequence=raw[GBSeq.SEQUENCE].upper(),
            comment=raw.get(GBSeq.COMMENT, None),
            source=NCBISource(
                taxid=int(source_dict["db_xref"].split(":")[1]),
                **source_dict,
            ),
        )

    async def fetch_taxonomy(
        self, taxid: int, cache_results: bool = True, use_cached: bool = True
    ) -> NCBITaxonomy | None:
        """
        Fetches and validates a taxonomy record from NCBI Taxonomy.

        If the record rank has an invalid rank (e.g. "no data"),
        does an additional docsum fetch and attempts to extract the rank data from there.

        :param taxid:
        :param cache_results:
        :param use_cached:
        """
        logger = base_logger.bind(taxid=taxid)

        record = None
        if use_cached:
            record = self.cache.load_taxonomy(taxid)
            if record:
                logger.info("Cached record found")

        if record is None:
            with Entrez.efetch(db=NCBIDB.TAXONOMY, id=taxid, rettype="null") as f:
                record = Entrez.read(f)[0]

        if cache_results:
            self.cache.cache_taxonomy(record, taxid)

        try:
            return NCBIClient.validate_taxonomy(record)

        except ValidationError as exc:
            logger.debug("Proper rank not found in record", error=exc.errors())

        logger.debug("Running additional docsum fetch...")
        with Entrez.efetch(
            db=NCBIDB.TAXONOMY, id=taxid, rettype="docsum", retmode="xml"
        ) as f:
            rank_str = Entrez.read(f)[0]["Rank"]
            try:
                rank = NCBIRank(rank_str)
                logger.debug("Valid rank found", rank=rank)
            except ValueError as exc:
                logger.error(
                    "Failed to find a valid rank. Returning empty...", error=exc
                )
                return None
        try:
            return NCBIClient.validate_taxonomy(record, rank)
        except ValidationError as exc:
            logger.error("Failed to find a valid rank. Returning empty...", error=exc)
            return None

    @staticmethod
    def validate_taxonomy(record: dict, rank: NCBIRank | None = None) -> NCBITaxonomy:
        """
        :param record: Taxonomy record in dict form
        :param rank: Optional variable, overrides inline rank data
        :return: A validated NCBI
        """
        species = None
        logger = base_logger.bind(taxid=record["TaxId"])

        # Get lineage list
        lineage = []
        for level_data in record["LineageEx"]:
            level = NCBILineage(
                id=int(level_data["TaxId"]),
                name=level_data["ScientificName"],
                rank=level_data["Rank"],
            )

            if level.rank == NCBIRank.SPECIES:
                species = level

            lineage.append(level)
        logger.debug("Validated lineage", lineage=lineage)

        # Fetch rank if not overwritten
        if rank is None:
            try:
                rank = NCBIRank(record["Rank"])
                if rank is rank.SPECIES:
                    species = NCBILineage(
                        id=int(record["TaxId"]),
                        name=record["ScientificName"],
                        rank=record["Rank"],
                    )
            except ValueError:
                logger.warning("Rank data not included")

        logger.debug("Check rank", rank=rank)

        return NCBITaxonomy(
            id=int(record["TaxId"]),
            species=species,
            lineage=lineage,
            rank=rank,
        )

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
    async def fetch_taxonomy_by_taxid(taxid: int) -> dict:
        """Requests a taxonomy record from NCBI Taxonomy

        :param taxid:
        :return:
        """
        with Entrez.efetch(db=NCBIDB.TAXONOMY, id=taxid, rettype="null") as f:
            record = Entrez.read(f)[0]

        # Handle insufficient rank data
        if record["Rank"] == "no rank":
            with Entrez.efetch(
                db=NCBIDB.TAXONOMY,
                id=taxid,
                rettype="docsum",
                retmode="xml",
            ) as f:
                rank = Entrez.read(f)

            record["Rank"] = rank

        return record

    @staticmethod
    async def get_species_taxid(taxonomy: dict) -> int | None:
        """Gets the species taxid for the given lower-rank taxid.

        :param taxonomy: NCBI Taxonomy UID
        :return: The NCBI Taxonomy ID of the OTU's species
        """
        if taxonomy["Rank"] == "species":
            return int(taxonomy["TaxId"])

        for line in taxonomy["LineageEx"]:
            if line["Rank"] == "species":
                return int(line["TaxId"])

        return None

    @staticmethod
    async def check_spelling(name: str, db: NCBIDB = NCBIDB.TAXONOMY) -> str:
        """Takes the name of an OTU, requests an alternative spelling
        from the Entrez ESpell utility and returns the suggestion

        :param name: The OTU name that requires correcting
        :param db: Database to check against. Defaults to ``taxonomy``.
        :return: String containing NCBI-suggested spelling changes
        """
        with Entrez.espell(db=db, term=quote_plus(name)) as f:
            record = Entrez.read(f)

        if record:
            return record["CorrectedQuery"]

        return name

    @staticmethod
    def _get_source_dict(raw: dict) -> dict:
        """
        :param raw:
        :return:
        """
        source_table = None
        for feature in raw[GBSeq.FEATURE_TABLE]:
            if feature["GBFeature_key"] == "source":
                source_table = feature
                break

        if source_table is not None:
            source_dict = {}

            for qualifier in source_table["GBFeature_quals"]:
                qual_name = qualifier["GBQualifier_name"]
                qual_value = qualifier["GBQualifier_value"]
                source_dict[qual_name] = qual_value

            return source_dict

        raise NCBIParseError(
            keys=[feature["GBFeature_key"] for feature in raw[GBSeq.FEATURE_TABLE]],
            message="Feature table does not contain source data",
        )

    @staticmethod
    def _parse_taxid(source_feature) -> int | None:
        return int(source_feature["db_xref"].split(":")[1])


class GBSeq(StrEnum):
    ACCESSION = "GBSeq_accession-version"
    DEFINITION = "GBSeq_definition"
    SEQUENCE = "GBSeq_sequence"
    LENGTH = "GBSeq_length"
    COMMENT = "GBSeq_comment"
    FEATURE_TABLE = "GBSeq_feature-table"
