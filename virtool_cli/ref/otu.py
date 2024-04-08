from typing import NamedTuple
from enum import StrEnum

import structlog

from virtool_cli.ncbi.model import NCBIGenbank


base_logger = structlog.get_logger()

UNKNOWN = "unknown"


class SourceType(StrEnum):
    ISOLATE = "isolate"
    STRAIN = "strain"
    CLONE = "clone"
    REFSEQ = "refseq"


class SourceKey(NamedTuple):
    type: SourceType
    name: str


def group_genbank_records_by_isolate(records: list[NCBIGenbank]) -> dict:
    """
    :param records:
    :return:
    """
    isolates = {}

    for record in records:
        logger = base_logger.bind(accession=record.accession)

        if record.source.model_fields_set.intersection(
            {SourceType.ISOLATE, SourceType.STRAIN, SourceType.CLONE}
        ):
            for source_type in SourceType:
                if source_type in record.source.model_fields_set:
                    source_key = SourceKey(
                        type=SourceType(source_type),
                        name=record.source.model_dump()[source_type],
                    )

                    if source_key not in isolates:
                        isolates[source_key] = {}

                    isolates[source_key][record.accession] = record

        else:
            if record.refseq:
                logger.debug(
                    "RefSeq record does not contain sufficient source data. Edit before inclusion.",
                    record=record,
                )

                source_key = SourceKey(
                    type=SourceType(SourceType.REFSEQ), name=record.accession
                )

                if source_key not in isolates:
                    isolates[source_key] = {}

                isolates[source_key][record.accession] = record
            else:
                logger.debug(
                    "Unreviewed record does not contain sufficient source data for inclusion."
                )
                break

    return isolates
