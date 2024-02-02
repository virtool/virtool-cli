
from Bio import SeqRecord

from virtool_cli.add.helpers import find_taxon_id
from virtool_cli.utils.format import check_source_type, format_sequence, get_qualifiers
from virtool_cli.utils.ncbi import fetch_isolate_metadata


async def format_record(record: SeqRecord) -> dict:
    seq_qualifiers = get_qualifiers(record.features)

    return {
        **format_sequence(record=record, qualifiers=seq_qualifiers),
        "isolate": await format_isolate_metadata(qualifiers=seq_qualifiers),
    }


async def format_isolate_metadata(qualifiers) -> dict:
    if isolate_type := check_source_type(qualifiers):
        # Isolate metadata contained in qualifiers
        return {
            "source_name": qualifiers.get(isolate_type)[0],
            "source_type": isolate_type,
        }

    # Extract isolate metadata from NCBI Taxonomy docsum.
    return await fetch_isolate_metadata(find_taxon_id(qualifiers["db_xref"]))
