async def is_accession_extant(new_accession: str, accession_list: list) -> bool:
    """
    Check if a new accession already exists in a list of already-assessed accessions.

    :param new_accession: A new accession
    :param accession_list: A list of accessions that should not be added anew
    :return: True if the accession collides with the accession list, False if not
    """
    for existing_accession in accession_list:
        new_accession_stripped = new_accession.split(".")
        extant_accession_stripped = existing_accession.split(".")

        if new_accession_stripped[0] == extant_accession_stripped[0]:
            return True

    return False


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
