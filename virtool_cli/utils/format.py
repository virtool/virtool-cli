from Bio.SeqIO import SeqRecord


def format_sequence(record: SeqRecord, qualifiers: dict) -> dict:
    """
    Creates a new sequence file for a given isolate

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

    if (record_host := qualifiers.get('host', None)) is not None:
        seq_dict['host'] = record_host[0]

    if (record_segment := qualifiers.get('segment', None)) is not None:
        seq_dict['segment'] = record_segment[0]

    return seq_dict


def format_isolate(source_name: str, source_type: str, isolate_id: str) -> dict:
    """
    Formats a

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
