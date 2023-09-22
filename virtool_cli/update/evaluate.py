from structlog import BoundLogger


def evaluate_sequence(
    seq_data, seq_qualifier_data, required_parts: dict, logger: BoundLogger
):
    """
    Evaluates the validity of the record using catalogued metadata
    """

    record_segment_list = seq_qualifier_data.get("segment", None)
    logger.debug(f"segments={record_segment_list}")

    if len(required_parts) > 1:
        # This doesn't take into account segment name differences
        if record_segment_list is None:
            logger.debug("No segment name. Discarding record...")
            return False

        segment_name = record_segment_list[0]

        logger = logger.bind(record_segment=segment_name)
        if segment_name not in required_parts.keys():
            logger.debug("Required segment not found. Moving on...")
            return False

        if required_parts[segment_name] < 0:
            logger.error(
                f"Sequence length parameter for {segment_name} is not in catalog."
                + "Moving on..."
            )
            return False

        try:
            listed_length = required_parts[segment_name]
        except KeyError as e:
            logger.exception(e)
            listed_length = -1

    else:
        segment_name, listed_length = required_parts.copy().popitem()
        logger.debug(f"{segment_name}: {listed_length}")

    logger = logger.bind(segment=segment_name)

    if listed_length < 0:
        logger.warning("No valid length available for this segment.")

    seq_length = len(seq_data.seq)

    if not valid_length(seq_length, listed_length):
        logger.debug(
            "Bad length. Moving on...",
            seq_length=seq_length,
            listing_length=listed_length,
        )

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


def get_lengthdict_multipartite(schema: list, logger: BoundLogger) -> dict:
    """
    Returns a dict of each required segment and its average length.

    :param schema: Augmented schema from the catalog listing, contains average length data
    :param logger: Optional entry point for an existing BoundLogger
    :return: A dict of required segments and their average lengths
    """
    required_part_lengths = {}

    for part in schema:
        if part["required"]:
            part_length = part.get("length", -1)
            if part_length < 0:
                logger.warning(f"{part['name']} lacks a length listing.")

            required_part_lengths[part["name"]] = part_length

    return required_part_lengths


def get_lengthdict_monopartite(schema: list, logger: BoundLogger) -> dict:
    """
    Returns a dict of the segment name and its average length.
    Given that the OTU is monopartite, a filler segment name can be used if none is found.

    :param schema: Augmented schema from the catalog listing, contains average length data
    :param logger: Optional entry point for an existing BoundLogger
    :return: A dict of the single segment and its average length
    """
    part = schema[0]

    if (part_name := part.get("name", None)) is None:
        part_name = "unlabelled"

    part_length = part.get("length", -1)
    if part_length < 0:
        logger.warning(f"{part['name']} lacks a length listing.")

    return {part_name: part_length}
