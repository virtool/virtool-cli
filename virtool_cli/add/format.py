from Bio import SeqRecord
from structlog import BoundLogger, get_logger

from virtool_cli.utils.format import format_sequence, get_qualifiers
from virtool_cli.utils.format import check_source_type
from virtool_cli.utils.ncbi import fetch_isolate_metadata
from virtool_cli.add.helpers import find_taxon_id


async def format_record(data: SeqRecord, logger: BoundLogger = get_logger()):
    seq_qualifiers = get_qualifiers(data.features)

    new_sequence = format_sequence(record=data, qualifiers=seq_qualifiers)
    new_sequence["isolate"] = await format_isolate_metadata(qualifiers=seq_qualifiers)

    logger.debug(new_sequence)

    return new_sequence


async def format_isolate_metadata(qualifiers):
    if isolate_type := check_source_type(qualifiers):
        # Isolate metadata contained in qualifiers
        isolate = {
            "source_name": qualifiers.get(isolate_type)[0],
            "source_type": isolate_type,
        }

    else:
        # Extract isolate metadata from NCBI Taxonomy docsum
        isolate = await fetch_isolate_metadata(find_taxon_id(qualifiers["db_xref"]))

    return isolate
