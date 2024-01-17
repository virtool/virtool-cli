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


def check_source_type(qualifiers: dict) -> str | None:
    """
    Determine the source type in a Genbank record

    :param qualifiers: Dictionary containing qualifiers in a features section of a Genbank record
    :return: "isolate" or "strain" if found in the record metadata, None if not found
    """
    for qualifier in ["isolate", "strain"]:
        if qualifier in qualifiers:
            return qualifier

    return None


def find_taxon_id(db_xref: list[str]) -> int | None:
    """
    Searches the database cross-reference data for the associated NCBI taxonomy UID.

    :param db_xref: List of NCBI cross-reference information taken from NCBI taxonomy record
    :return: NCBI Taxonomy UID as an integer if found, None if not found
    """
    for xref in db_xref:
        [key, value] = xref.split(":")
        if key == "taxon":
            return int(value)

    return None
