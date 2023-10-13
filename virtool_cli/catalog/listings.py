import json
import aiofiles
from pathlib import Path
from structlog import BoundLogger, get_logger


async def parse_listing(path):
    """
    Asynchronously reads listing data from path and returns in dict form

    :param path: Path to the listing file
    :return: Deserialized listing data
    """
    async with aiofiles.open(path, "r") as f:
        contents = await f.read()
        listing = json.loads(contents)

    return listing


def generate_new_listing(otu_id="", name="", taxid=-1) -> dict:
    """
    Generates a new listing from new OTU data regardless of input

    :return: An empty listing dictionary if no input is given OR
        a barebones skeleton for filling out if input is given
    """
    catalog_listing = {
        "_id": otu_id,
        "accessions": {"excluded": [], "included": []},
        "name": name,
        "schema": [],
        "taxid": taxid,
    }
    return catalog_listing


async def update_listing(data: dict, path: Path):
    """
    Overwrites the listing file at path with new data

    :param data: Updated listing data
    :param path: Path to the listing file
    """
    with open(path, "w") as f:
        json.dump(data, f, indent=2, sort_keys=True)


async def write_new_listing(listing: dict, catalog_path: Path) -> Path | None:
    """
    Writes prepared listing data to a listing file under the catalog directory

    :param listing: Deserialized OTU catalog listing
    :param catalog_path: Path to an accession catalog
    """
    listing_path = catalog_path / f"{listing['taxid']}--{listing['_id']}.json"

    await update_listing(data=listing, path=listing_path)

    return listing_path


async def add_new_listing(otu_path: Path, catalog_path: Path, logger=get_logger()):
    """
    Creates a new listing for a newly created OTU with no isolate or accession data

    :param otu_path: Path to an OTU directory
    :param catalog_path: Path to an accession catalog
    :param logger: Optional entry point for an existing BoundLogger
    """
    async with aiofiles.open(otu_path / "otu.json", "r") as f:
        contents = await f.read()
        otu = json.loads(contents)

    listing = generate_new_listing(
        otu_id=otu["_id"], name=otu["name"], taxid=otu["taxid"]
    )

    listing_path = await write_new_listing(listing, catalog_path)
    if listing_path is None:
        logger.error("Listing could not be created under catalog")
        return

    logger.info("Listing written to listing_path", listing_path=str(listing_path))

    return catalog_path


async def generate_listing(
    otu_data: dict,
    accession_list: list,
    sequence_metadata: dict,
    logger: BoundLogger = get_logger(),
) -> dict:
    """
    Takes a fully populated OTU directory and generates a listing based on the contained data

    :param otu_data: OTU data in dict form
    :param accession_list: list of included accesssions
    :param sequence_metadata: dict of sequence metadata, including average length
    :param logger: Optional entry point for an existing BoundLogger
    """
    if (taxid := otu_data.get("taxid", None)) is None:
        taxid = -1

    catalog_listing = generate_new_listing(
        otu_id=otu_data["_id"], name=otu_data["name"], taxid=taxid
    )

    catalog_listing["name"] = otu_data["name"]

    catalog_listing["accessions"]["included"] = accession_list

    schema = otu_data.get("schema", [])
    logger.debug(f"{schema}")

    if schema:
        if len(schema) > 1:
            required_parts = get_required_parts(schema)
            try:
                length_dict = measure_multipartite(sequence_metadata, required_parts)
            except Exception as e:
                logger.exception(e)
                raise e

            for part in schema:
                if part["name"] in length_dict:
                    part["length"] = length_dict[part["name"]]

        else:
            try:
                average_length = measure_monopartite(sequence_metadata)
                logger.debug(f"Average length is {average_length}")
                schema[0]["length"] = average_length
            except Exception as e:
                logger.exception(e)

    catalog_listing["schema"] = schema

    return catalog_listing


def get_required_parts(schema: list) -> list:
    """
    Takes a schema list and returns a list containing the names of required parts only

    :param schema: Schema from the catalog listing, originally taken from the reference data
    :return: List of segment names that are listed as "required"
    """
    required_parts = []
    for part in schema:
        if part["required"]:
            required_parts.append(part["name"])

    return required_parts


def measure_monopartite(sequence_metadata: dict) -> int:
    """
    Takes a dict containing all sequence lengths and returns the average length

    :param sequence_metadata: Dict containing segment data
    :return: Average length of all sequences as an integer
    """
    sequence_lengths = []

    for accession in sequence_metadata:
        metadata = sequence_metadata[accession]
        seq_length = metadata["length"]
        sequence_lengths.append(seq_length)

    average_length = sum(sequence_lengths) / len(sequence_lengths)

    return int(average_length)


def measure_multipartite(sequence_metadata: dict, segment_list: list) -> dict:
    """
    Takes a dict containing all sequence lengths and a list of segments,
    and returns the average length of each constituent segment.

    :param sequence_metadata: Dict containing segment data
    :param segment_list: List of segment names
    :return: A dict keyed by segment name containing the average length of each segment
    """
    part_total_dict = {element: [] for index, element in enumerate(segment_list)}

    for accession in sequence_metadata:
        metadata = sequence_metadata[accession]
        seq_length = metadata["length"]

        segment_name = metadata.get("segment", None)
        if segment_name is None:
            continue
        if segment_name in part_total_dict.keys():
            part_total_dict[segment_name].append(seq_length)

    length_dict = {}
    for segment in part_total_dict:
        if not part_total_dict[segment]:
            length_dict[segment] = 0
            continue
        average_length = sum(part_total_dict[segment]) / len(part_total_dict[segment])
        length_dict[segment] = int(average_length)

    return length_dict
