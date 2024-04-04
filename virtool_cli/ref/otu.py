import re
from dataclasses import dataclass
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
                    source_value = record.source.model_dump()[source_type]

                    if source_value not in isolates:
                        isolates[(source_type, source_value)] = {}

                    isolates[(source_type, source_value)][record.accession] = record

        else:
            if record.refseq:
                source_value = record.accession

                if source_value not in isolates:
                    isolates[(SourceType.REFSEQ, source_value)] = {}

                isolates[(SourceType.REFSEQ, source_value)][record.accession] = record
                logger.debug(
                    "RefSeq record does not contain sufficient source data. Edit before inclusion."
                )
            else:
                logger.debug(
                    "Unreviewed record does not contain sufficient source data for inclusion."
                )

    return isolates
