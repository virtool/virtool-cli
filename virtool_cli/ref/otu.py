import sys
from collections import defaultdict
from enum import StrEnum
from typing import NamedTuple

import structlog

from virtool_cli.ncbi.client import NCBIClient
from virtool_cli.ncbi.model import NCBIGenbank
from virtool_cli.ref.repo import EventSourcedRepo as Repo
from virtool_cli.ref.resources import EventSourcedRepoOTU as RepoOTU
from virtool_cli.ref.utils import Molecule, IsolateName

base_logger = structlog.get_logger()


class SourceType(StrEnum):
    ISOLATE = "isolate"
    STRAIN = "strain"
    CLONE = "clone"
    REFSEQ = "refseq"


class SourceKey(NamedTuple):
    type: SourceType
    name: str

    def __str__(self):
        return f"{self.type},{self.name}"


def add_otu(repo: Repo, taxid: int) -> RepoOTU:
    logger = base_logger.bind(taxid=taxid)

    ncbi = NCBIClient.from_repo(repo.path, False)

    taxonomy = ncbi.fetch_taxonomy_record(taxid)
    if taxonomy is None:
        logger.fatal(f"Taxonomy ID {taxid} not found")
        sys.exit(1)

    molecule = None

    try:
        otu = repo.create_otu(
            acronym="",
            legacy_id=None,
            name=taxonomy.name,
            molecule=molecule,
            schema=[],
            taxid=taxid,
        )
    except ValueError as e:
        logger.warning(e)
        sys.exit(1)

    logger.debug("Created OTU", id=str(otu.id), name=otu.name, taxid=taxid)

    return otu


class OTUClient:
    """Handles NCBI data for a single OTU."""

    def __init__(self, repo: Repo, otu: RepoOTU, ignore_cache: bool = False):
        self._repo = repo
        self.otu = otu
        self.ignore_cache = ignore_cache
        self.client = NCBIClient.from_repo(self._repo.path, False)

    @classmethod
    def init_from_taxid(
        cls, repo, taxid: int, create_otu: bool = True, ignore_cache: bool = False
    ):
        """
        Initializes a new OTUClient from a Taxonomy ID.

        If create_otu flag is True, fetches the Taxonomy record and adds a new OTU.
        """
        logger = base_logger.bind(taxid=taxid)
        otu_index = repo.index_otus()

        if taxid in otu_index:
            otu = repo.get_otu(otu_index[taxid])

            return OTUClient(repo, otu, ignore_cache)
        else:
            if create_otu:
                logger.info(
                    "This OTU has not been added to the reference yet. Requesting from Taxonomy...",
                )
                otu = add_otu(repo, taxid)

                if otu:
                    return OTUClient(repo, otu, ignore_cache)
            else:
                raise ValueError

        logger.error("OTU could not be found or built.")

        raise ValueError

    @classmethod
    def create_from_taxid(cls, repo, taxid: int, ignore_cache: bool = False):
        """Initializes a new OTU and a new OTUClient from a Taxonomy ID."""
        logger = base_logger.bind(taxid=taxid)
        otu_index = repo.index_otus()
        if taxid in otu_index:
            otu = repo.get_otu(otu_index[taxid])

            logger.error(
                f"Taxonomy ID {taxid} has already been added to this reference.",
                otu_id=str(otu.id),
            )
            raise ValueError

        otu = add_otu(repo, taxid)
        if otu:
            return OTUClient(repo, otu, ignore_cache)

        logger.warning("OTU could not be found or built.")

        raise ValueError

    def update(self):
        """Fetches a full list of Nucleotide accessions associated with the OTU
        and passes the list to the add."""
        ncbi = NCBIClient.from_repo(self._repo.path, self.ignore_cache)

        linked_accessions = ncbi.link_accessions_from_taxid(self.otu.taxid)

        self.add(linked_accessions)

    def add(self, accessions):
        otu_logger = base_logger.bind(
            taxid=self.otu.taxid, otu_id=str(self.otu.id), name=self.otu.name
        )
        fetch_list = list(set(accessions).difference(self.otu.accession_set))
        if not fetch_list:
            otu_logger.info("OTU is up to date.")
            return

        otu_logger.info(
            f"Fetching {len(fetch_list)} accessions...", fetch_list=fetch_list
        )

        records = self.client.fetch_genbank_records(fetch_list)

        record_bins = group_genbank_records_by_isolate(records)

        new_accessions = []

        for isolate_key in record_bins:
            record_bin = record_bins[isolate_key]

            isolate_id = self.otu.get_isolate_id_by_name(
                IsolateName(type=isolate_key.type, value=isolate_key.name)
            )
            if isolate_id is None:
                otu_logger.debug("Creating isolate")
                isolate = self._repo.create_isolate(
                    otu_id=self.otu.id,
                    legacy_id=None,
                    source_name=isolate_key.name,
                    source_type=isolate_key.type,
                )
                isolate_id = isolate.id

            for accession in record_bin:
                record = record_bin[accession]
                sequence = self._repo.create_sequence(
                    otu_id=self.otu.id,
                    isolate_id=isolate_id,
                    accession=record.accession,
                    definition=record.definition,
                    legacy_id=None,
                    segment=record.source.segment,
                    sequence=record.sequence,
                )

                new_accessions.append(sequence.accession)

        if new_accessions:
            otu_logger.info(
                f"Added {len(new_accessions)} sequences to {self.otu.taxid}",
                new_accessions=new_accessions,
            )

        else:
            otu_logger.info(f"No new sequences added to OTU")


def group_genbank_records_by_isolate(records: list[NCBIGenbank]) -> dict:
    """:param records:
    :return:
    """
    isolates = defaultdict(dict)

    for record in records:
        logger = base_logger.bind(
            accession=record.accession,
            definition=record.definition,
            source_data=record.source,
        )

        if record.source.model_fields_set.intersection(
            {SourceType.ISOLATE, SourceType.STRAIN, SourceType.CLONE},
        ):
            for source_type in SourceType:
                if source_type in record.source.model_fields_set:
                    source_key = SourceKey(
                        type=SourceType(source_type),
                        name=record.source.model_dump()[source_type],
                    )

                    isolates[source_key][record.accession] = record

                    break

        elif record.refseq:
            logger.debug(
                "RefSeq record does not contain sufficient source data. "
                + "Edit before inclusion.",
                record=record,
            )

            source_key = SourceKey(
                type=SourceType(SourceType.REFSEQ),
                name=record.accession,
            )

            isolates[source_key][record.accession] = record
        else:
            logger.debug(
                "Record does not contain sufficient source data for inclusion.",
            )

    return isolates


def get_molecule_from_records(records: list[NCBIGenbank]) -> Molecule:
    for record in records:
        if record.refseq:
            return Molecule(
                **{
                    "strandedness": record.strandedness.value,
                    "type": record.moltype.value,
                    "topology": record.topology.value,
                }
            )

    return Molecule(
        **{
            "strandedness": records[0].strandedness.value,
            "type": records[0].moltype.value,
            "topology": records[0].topology.value,
        }
    )
