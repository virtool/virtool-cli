import sys
from collections import defaultdict
from enum import StrEnum
from typing import NamedTuple

import structlog

from virtool_cli.legacy.models import LegacyIsolateSource, LegacySourceType
from virtool_cli.ncbi.client import NCBIClient
from virtool_cli.ncbi.model import NCBIGenbank
from virtool_cli.ref.repo import EventSourcedRepo as Repo
from virtool_cli.ref.resources import EventSourcedRepoOTU as RepoOTU

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
    except ValueError:
        sys.exit(1)

    logger.info("Created OTU", id=str(otu.id), name=otu.name, taxid=taxid)

    return otu


def add_accessions(repo: Repo, otu: RepoOTU, accessions: list[str]):
    otu_logger = base_logger.bind(taxid=otu.taxid, otu_id=str(otu.id))
    fetch_list = list(set(accessions).difference(set(otu.nofetch)))

    ncbi = NCBIClient.from_repo(repo.path, False)

    records = ncbi.fetch_genbank_records(fetch_list)

    record_bins = group_genbank_records_by_isolate(records)

    for isolate_key in record_bins:
        record_bin = record_bins[isolate_key]

        isolate_id = otu.get_isolate_id(type=isolate_key.type, name=isolate_key.name)
        if isolate_id is None:
            otu_logger.debug("Creating isolate")
            isolate = repo.create_isolate(
                otu_id=otu.id,
                legacy_id=None,
                source_name=isolate_key.name,
                source_type=isolate_key.type,
            )
            isolate_id = isolate.id

        for accession in record_bin:
            record = record_bin[accession]
            sequence = repo.create_sequence(
                otu_id=otu.id,
                isolate_id=isolate_id,
                accession=record.accession,
                definition=record.definition,
                legacy_id=None,
                segment=record.source.segment,
                sequence=record.sequence,
            )

        # for accession in record_bin:
        #     otu_logger.debug(record_bin[accession], isolate=str(isolate_key))


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

        logger.debug(
            "Record does not contain sufficient source data for inclusion.",
        )

    return isolates


def extract_isolate_source(
    genbank_records: list[NCBIGenbank],
) -> LegacyIsolateSource:
    """Extract a legacy isolate source from a set of Genbank records associated with the
    isolate.
    """
    for record in genbank_records:
        if record.source.isolate:
            return LegacyIsolateSource(
                name=record.source.isolate,
                type=LegacySourceType.ISOLATE,
            )

        if record.source.strain:
            return LegacyIsolateSource(
                name=record.source.strain,
                type=LegacySourceType.STRAIN,
            )

        if record.source.clone:
            return LegacyIsolateSource(
                name=record.source.clone,
                type=LegacySourceType.CLONE,
            )

    accessions = sorted(
        (record.accession for record in genbank_records if record.accession),
        key=lambda x: int(x.replace("NC_", "").replace(".1", "")),
    )

    return LegacyIsolateSource(
        name=accessions[0].upper(),
        type=LegacySourceType.GENBANK,
    )
