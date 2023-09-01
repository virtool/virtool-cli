from structlog import BoundLogger


def evaluate_sequence(
    seq_data, 
    seq_qualifier_data, 
    required_parts, 
    logger: BoundLogger
):
    """
    Evaluates the validity of the record using catalogued metadata
    """

    record_segment_list = seq_qualifier_data.get("segment", None)
    logger.debug(record_segment_list)

    segment_name = ''
    if len(required_parts) > 1:
        # This doesn't take into account segment name differences
        if record_segment_list is None:
            logger.debug('No segment name. Discarding record...')
            return False
        
        segment_name = record_segment_list[0]
        
        logger = logger.bind(record_segment=segment_name)
        if segment_name not in required_parts.keys():
            logger.debug('Required segment not found. Moving on...')
            return False

        if required_parts[segment_name] < 0:
            logger.error(
                f'Sequence length parameter for {segment_name} is not in catalog.' + 'Moving on...')
            return False
        
        try:
            listed_length = required_parts[segment_name]
        except KeyError as e:
            logger.exception(e)
    
    else:
        segment_name, listed_length = required_parts.copy().popitem()
        logger.debug(f'{segment_name}: {listed_length}')
    
    logger = logger.bind(segment=segment_name)
    
    if listed_length < 0:
        logger.warning('No valid length available for this segment.')
    
    seq_length = len(seq_data.seq)
    
    if not valid_length(seq_length, listed_length):
        logger.debug('Bad length. Moving on...', 
            seq_length=seq_length, listing_length=listed_length)
        
        return False
    
    return True

def valid_length(seq_length: int, listed_length: int) -> bool:
    """
    Returns true if an integer is within 10% above or below a listed average value.
    Used to check if the length of a new sequence is within acceptable bounds.

    :param seq_length: Length of newly fetched sequence from GenBank
    :param listed_length: Average length of accepted sequences corresponding to this segment
    :return: Dictionary containing all qualifiers in the source field of the features section of a Genbank record
    """
    max_length = listed_length * 1.10
    min_length = listed_length * 0.9
    
    if seq_length > max_length or seq_length < min_length:
        return False
    
    return True

def get_qualifiers(seq: list) -> dict:
    """
    Get relevant qualifiers in a Genbank record

    :param seq: SeqIO features object for a particular accession
    :return: Dictionary containing all qualifiers in the source field of the features section of a Genbank record
    """    
    qualifiers = {}

    for feature in [feature for feature in seq if feature.type == "source"]:
        for qual_key in feature.qualifiers:
            qualifiers[qual_key] = feature.qualifiers.get(qual_key)

    return qualifiers