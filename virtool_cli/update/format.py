from Bio import SeqIO
from structlog import BoundLogger


def format_sequence(
    record: SeqIO.SeqRecord, qualifiers: dict, 
    logger: BoundLogger
) -> dict:
    """
    Creates a new sequence file for a given isolate

    :param record: Genbank record object for a given accession
    :param qualifiers: Dictionary containing all qualifiers in the source field of the features section of a Genbank record
    :return: A new sequence dictionary if possible, else an empty dict if not
    """
    logger = logger.bind(accession=record.id)
    
    try:
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
    
    except Exception as e:
        logger.exception(e)
        return {}