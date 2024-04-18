from collections import defaultdict
from enum import StrEnum
from typing import NamedTuple

import structlog

from virtool_cli.legacy.models import LegacyIsolateSource, LegacySourceType
from virtool_cli.ncbi.model import NCBIGenbank

base_logger = structlog.get_logger()


class SourceType(StrEnum):
    ISOLATE = "isolate"
    STRAIN = "strain"
    CLONE = "clone"
    REFSEQ = "refseq"


class SourceKey(NamedTuple):
    type: SourceType
    name: str


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
