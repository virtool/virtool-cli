from typing import Tuple

from Bio.SeqIO import SeqRecord
from structlog import BoundLogger, get_logger


async def process_default(
    records: list,
    metadata: dict,
    filter_set: set,
    logger: BoundLogger = get_logger(),
) -> Tuple[list, list]:
    """Format new sequences from NCBI Taxonomy if they do not already exist in the reference.

    :param records: A list of SeqRecords from NCBI Taxonomy
    :param metadata: A deserialized OTU metadata file
    :param filter_set: A set of accessions that should be omitted
    :param logger: Optional entry point for an existing BoundLogger
    :return: A list of processed new sequences/isolates and
        a set of automatically excluded accessions
    """
    auto_excluded = []
    otu_updates = []

    for seq_data in records:
        accession = seq_data.id.split(".")[0]
        seq_qualifier_data = get_qualifiers(seq_data.features)

        if accession in filter_set:
            logger.debug("Accession already exists", accession=seq_data.id)
            continue

        if check_source_type(seq_qualifier_data) is None:
            continue
        isolate = find_isolate_metadata(seq_qualifier_data)

        seq_dict = format_sequence(record=seq_data, qualifiers=seq_qualifier_data)

        if "segment" not in seq_dict:
            schema = metadata.get("schema", [])

            if schema:
                seq_dict["segment"] = schema[0].get("name", "")
            else:
                logger.warning("Missing schema")
                seq_dict["segment"] = ""

        seq_dict["isolate"] = isolate
        otu_updates.append(seq_dict)

    return otu_updates, auto_excluded


def format_sequence(record: SeqRecord, qualifiers: dict) -> dict:
    """Creates a new sequence file for a given isolate

    :param record: Genbank record object for a given accession
    :param qualifiers: Dictionary containing all qualifiers in the source field
        of the features section of a Genbank record
    :return: A dict containing the new sequence data and metadata
    """
    seq_dict = {
        "accession": record.id,
        "definition": record.description,
        "sequence": str(record.seq),
    }

    if (record_host := qualifiers.get("host")) is not None:
        seq_dict["host"] = record_host[0]

    if (record_segment := qualifiers.get("segment")) is not None:
        seq_dict["segment"] = record_segment[0]

    return seq_dict


def format_isolate(source_name: str, source_type: str, isolate_id: str) -> dict:
    """Formats raw isolate data for storage in a reference directory

    :param source_name: Assigned source name for an accession
    :param source_type: Assigned source type for an accession
    :param isolate_id: Unique ID number for this new isolate
    :return: A dict containing the new isolate.json contents
    """
    isolate = {
        "id": isolate_id,
        "source_type": source_type,
        "source_name": source_name,
        "default": False,
    }
    return isolate


def get_qualifiers(seq: list) -> dict:
    """Get relevant qualifiers in a Genbank record

    :param seq: SeqIO features object for a particular accession
    :return: Dictionary containing all qualifiers in the source field of the features section of a Genbank record
    """
    qualifiers = {}

    for feature in [feature for feature in seq if feature.type == "source"]:
        for qual_key in feature.qualifiers:
            qualifiers[qual_key] = feature.qualifiers.get(qual_key)

    return qualifiers


def find_isolate_metadata(qualifiers: dict) -> dict:
    isolate_type = check_source_type(qualifiers)
    isolate_name = qualifiers.get(isolate_type)[0]
    return {"source_name": isolate_name, "source_type": isolate_type}


def check_source_type(qualifiers: dict) -> str | None:
    """Determine the source type in a Genbank record

    :param qualifiers: Dictionary containing qualifiers in a features section of a Genbank record
    :return:
    """
    for qualifier in ["isolate", "strain"]:
        if qualifier in qualifiers:
            return qualifier

    return None
