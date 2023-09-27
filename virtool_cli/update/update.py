import json
from pathlib import Path
import asyncio
from typing import Optional, Tuple
import structlog
from structlog import get_logger, BoundLogger
from urllib.error import HTTPError

from virtool_cli.utils.logging import DEFAULT_LOGGER, DEBUG_LOGGER
from virtool_cli.utils.reference import get_otu_paths, get_isolate_paths
from virtool_cli.utils.id_generator import generate_unique_ids, get_unique_ids
from virtool_cli.utils.ncbi import (
    request_linked_accessions,
    request_from_nucleotide,
    NCBI_REQUEST_INTERVAL,
)
from virtool_cli.update.evaluate import (
    evaluate_sequence,
    get_qualifiers,
    get_lengthdict_monopartite,
    get_lengthdict_multipartite,
)
from virtool_cli.utils.format import format_isolate, format_sequence
from virtool_cli.utils.storage import store_isolate, store_sequence
from virtool_cli.catalog.helpers import search_by_otu_id

DEFAULT_INTERVAL = 0.001


def run(
    otu_path: Path,
    catalog_path: Path,
    auto_evaluate: bool = False,
    debugging: bool = False,
):
    """
    CLI entry point for update.update.run()

    Requests updates for a single OTU directory
    Searches the catalog for a matching catalog listing and requests new updates if it finds one.

    :param otu_path: Path to a OTU directory
    :param catalog_path: Path to a catalog directory
    :param auto_evaluate: Auto-evaluation flag, enables automatic filtering for fetched results
    :param debugging: Enables verbose logs for debugging purposes
    """
    structlog.configure(wrapper_class=DEBUG_LOGGER if debugging else DEFAULT_LOGGER)
    logger = get_logger().bind(otu_path=str(otu_path))

    if auto_evaluate:
        logger.warning(
            "Auto-evaluation is in active development and may produce false negatives."
        )

    [_, otu_id] = otu_path.name.split("--")
    listing_path = search_by_otu_id(otu_id, catalog_path)

    if listing_path:
        asyncio.run(update_otu(otu_path, listing_path, auto_evaluate))

    else:
        logger.critical("Listing not found for this OTU.")
        logger.info("Use virtool acc init to create a new accession catalog.")


async def update_otu(otu_path: Path, listing_path: Path, auto_evaluate: bool = False):
    """
    Requests new records for a single taxon ID
    and writes new data under the corresponding path.

    :param otu_path: Path to a OTU directory
    :param listing_path: Path to a listing in an accession catalog directory
    :param auto_evaluate: Auto-evaluation flag, enables automatic filtering for fetched results
    """
    src_path = otu_path.parent
    listing = json.loads(listing_path.read_text())

    # extract taxon ID and _id hash from listing filename
    [taxid, otu_id] = listing_path.stem.split("--")

    logger = get_logger().bind(taxid=taxid, otu_id=otu_id)

    record_data = await request_new_records(listing)
    if not record_data:
        return

    otu_updates = await process_records(
        records=record_data, listing=listing, auto_evaluate=auto_evaluate, logger=logger
    )
    if not otu_updates:
        return

    # List all isolate and sequence IDs presently in src
    unique_iso, unique_seq = await get_unique_ids(get_otu_paths(src_path))

    await write_data(otu_path, otu_updates, unique_iso, unique_seq, logger=logger)


async def request_new_records(
    listing: dict, logger: BoundLogger = get_logger()
) -> list:
    """
    :param listing: Deserialized OTU catalog listing
    :param logger: Optional entry point for a shared BoundLogger
    """
    try:
        new_accessions = await fetch_upstream_accessions(listing=listing, logger=logger)
        await asyncio.sleep(NCBI_REQUEST_INTERVAL)

    except HTTPError as e:
        logger.error(e)
        return []

    except Exception as e:
        logger.exception(e)
        return []

    try:
        record_data = await request_from_nucleotide(new_accessions)
        await asyncio.sleep(NCBI_REQUEST_INTERVAL)

    except HTTPError as e:
        logger.error(e)
        return []

    except Exception as e:
        logger.exception(e)
        return []

    return record_data


async def process_records(
    records: list,
    listing: dict,
    auto_evaluate: bool = True,
    logger: BoundLogger = get_logger(),
) -> list:
    """
    Takes sequence records and:
        1) Evaluates whether those records should be added to the database,
        2) Formats the records into a smaller dictionary
        3) Returns new formatted dicts in a list

    WARNING: Auto-evaluation is still under active development, especially multipartite filtering

    :param records: SeqRecords retrieved from the NCBI Nucleotide database
    :param listing: Deserialized OTU catalog listing
    :param auto_evaluate: Boolean flag for whether automatic evaluation functions
        should be run
    :param logger: Optional entry point for a shared BoundLogger
    :return: A list of valid sequences formatted for the Virtool reference database
    """
    filter_set = set(listing["accessions"]["included"])
    filter_set.update(listing["accessions"]["excluded"])

    try:
        if auto_evaluate:
            otu_updates, auto_excluded = await process_auto_evaluate(
                records, listing, filter_set, logger
            )
        else:
            otu_updates, auto_excluded = await process_default(
                records, listing, filter_set, logger
            )
    except Exception as e:
        logger.exception(e)
        raise e

    if auto_excluded:
        logger.info(
            "Consider adding these accessions to the catalog exclusion list",
            auto_excluded=auto_excluded,
        )

    if otu_updates:
        return otu_updates

    return []


async def write_data(
    otu_path: Path,
    new_sequences: list,
    unique_iso: set,
    unique_seq: set,
    logger: BoundLogger = get_logger(),
):
    """
    :param otu_path: A path to an OTU directory under a src reference directory
    :param new_sequences: List of new sequences under the OTU
    :param unique_iso: Set of all unique isolate IDs present in the reference
    :param unique_seq: Set of all unique sequence IDs present in the reference
    :param logger: Optional entry point for a shared BoundLogger
    """
    ref_isolates = await label_isolates(otu_path)

    try:
        seq_hashes = generate_unique_ids(
            n=len(new_sequences), excluded=list(unique_seq)
        )
    except Exception as e:
        logger.exception(e)
        return e

    logger.debug(f"Writing {len(new_sequences)} sequences...")

    for seq_data in new_sequences:
        isolate_data = seq_data.pop("isolate")
        isolate_name = isolate_data["source_name"]
        isolate_type = isolate_data["source_type"]

        if isolate_name in ref_isolates:
            iso_id = ref_isolates[isolate_name]["id"]
            logger.debug(
                "Existing isolate name found", iso_name=isolate_name, iso_hash=iso_id
            )

        else:
            try:
                iso_id = generate_unique_ids(n=1, excluded=list(unique_iso)).pop()
            except Exception as e:
                logger.exception(e)
                continue

            logger.debug(
                "Assigning new isolate hash", iso_name=isolate_name, iso_hash=iso_id
            )

            new_isolate = format_isolate(isolate_name, isolate_type, iso_id)

            await store_isolate(new_isolate, iso_id, otu_path)

            unique_iso.add(iso_id)
            ref_isolates[isolate_name] = new_isolate

            logger.info("Created a new isolate directory", path=str(otu_path / iso_id))

        iso_path = otu_path / iso_id

        seq_hash = seq_hashes.pop()
        logger.debug(
            "Assigning new sequence",
            seq_hash=seq_hash,
        )

        await store_sequence(seq_data, seq_hash, iso_path)
        unique_seq.add(seq_hash)

        logger.info(
            f"Wrote new sequence '{seq_hash}'", path=str(iso_path / f"{seq_hash}.json")
        )


async def label_isolates(otu_path: Path) -> dict:
    """
    Return all isolates present in an OTU directory

    :param otu_path: Path to an OTU directory
    :return: A dictionary of isolates indexed by source_name
    """
    if not otu_path.exists():
        raise FileNotFoundError

    isolates = {}

    for iso_path in get_isolate_paths(otu_path):
        with open(iso_path / "isolate.json", "r") as f:
            isolate = json.load(f)
        isolates[isolate.get("source_name")] = isolate

    return isolates


def find_isolate(record_features: dict) -> Optional[str]:
    """
    Determine the source type in a Genbank record

    :param record_features: Dictionary containing qualifiers in a features section of a Genbank record
    :return: "isolate" or "strain" if listed as such in the qualifiers, else None
    """
    for qualifier in ["isolate", "strain"]:
        if qualifier in record_features:
            return qualifier

    return None


async def fetch_upstream_accessions(
    listing: dict, logger: BoundLogger = get_logger()
) -> list:
    """
    Requests a list of all uninspected accessions associated with an OTU's taxon ID

    :param listing: Corresponding catalog listing for this OTU
    :param logger: Optional entry point for an existing BoundLogger
    :return: A list of accessions from NCBI Genbank for the taxon ID,
        sans included and excluded accessions
    """
    taxid = listing.get("taxid")
    included_set = set(listing["accessions"]["included"])
    excluded_set = set(listing["accessions"]["excluded"])

    logger = logger.bind(taxid=taxid)
    logger.debug(
        "Exclude catalogued accessions", included=included_set, excluded=excluded_set
    )

    try:
        upstream_accessions = await request_linked_accessions(taxon_id=taxid)
    except HTTPError:
        logger.error("Could not retrieve new accessions from NCBI.")
        return []

    upstream_set = set(upstream_accessions)

    return list(upstream_set.difference(included_set))


async def process_default(
    records: list, listing: dict, filter_set: set, logger: BoundLogger = get_logger()
) -> Tuple[list, list]:
    """
    Format new sequences from NCBI Taxonomy if they do not already exist in the reference.

    :param records: A list of SeqRecords from NCBI Taxonomy
    :param listing: A deserialized catalog listing for the OTU
    :param filter_set: A set of accessions that should be omitted
    :param logger: Optional entry point for an existing BoundLogger
    :return: A list of processed new sequences/isolates and
        a set of automatically excluded accessions
    """
    auto_excluded = []
    otu_updates = []

    for seq_data in records:
        [accession, _] = seq_data.id.split(".")
        seq_qualifier_data = get_qualifiers(seq_data.features)

        if accession in filter_set:
            logger.debug("Accession already exists", accession=seq_data.id)
            continue

        isolate_type = find_isolate(seq_qualifier_data)
        if isolate_type is None:
            continue
        isolate_name = seq_qualifier_data.get(isolate_type)[0]

        seq_dict = format_sequence(
            record=seq_data, qualifiers=seq_qualifier_data, logger=logger
        )

        if "segment" not in seq_dict:
            seq_dict["segment"] = listing.get("schema")[0]["name"]

        isolate = {"source_name": isolate_name, "source_type": isolate_type}
        seq_dict["isolate"] = isolate
        otu_updates.append(seq_dict)

    return otu_updates, auto_excluded


async def process_auto_evaluate(
    records: list, listing: dict, filter_set: set, logger: BoundLogger = get_logger()
) -> Tuple[list, list]:
    """
    Format new sequences from NCBI Taxonomy if they do not already exist in the reference.
    Automatically evaluate new sequences during formatting.

    :param records: A list of SeqRecords from NCBI Taxonomy
    :param listing: A deserialized catalog listing for the OTU
    :param filter_set: A set of accessions that should be omitted
    :param logger: Optional entry point for an existing BoundLogger
    :return: A list of processed new sequences/isolates and
        a set of automatically excluded accessions
    """
    auto_excluded = []
    otu_updates = []

    if not listing.get("schema", []):
        logger.warning("Missing schema. Moving on...")
        return [], []

    if len(listing["schema"]) < 1:
        required_parts = get_lengthdict_multipartite(listing["schema"], logger)

    else:
        required_parts = get_lengthdict_monopartite(listing["schema"], logger)

    logger = logger.bind(required_parts=required_parts)

    for seq_data in records:
        [accession, _] = seq_data.id.split(".")
        seq_qualifier_data = get_qualifiers(seq_data.features)

        if accession in filter_set:
            logger.debug("Accession already exists", accession=seq_data.id)
            return [], []

        if "_" in seq_data.id:
            logger.warning(
                "This is a RefSeq accession. This data may already exist under a different accession.",
                accession=seq_data.id,
            )

        inclusion_passed = evaluate_sequence(
            seq_data=seq_data,
            seq_qualifier_data=seq_qualifier_data,
            required_parts=required_parts,
            logger=logger,
        )
        if not inclusion_passed:
            logger.debug(
                f"GenBank #{accession} did not pass the curation test...",
                accession=accession,
                including=False,
            )
            auto_excluded.append(accession)
            continue

        isolate_type = find_isolate(seq_qualifier_data)
        if isolate_type is None:
            continue
        isolate_name = seq_qualifier_data.get(isolate_type)[0]

        seq_dict = format_sequence(record=seq_data, qualifiers=seq_qualifier_data)

        if "segment" not in seq_dict:
            seq_dict["segment"] = listing.get("schema")[0]["name"]

        isolate = {"source_name": isolate_name, "source_type": isolate_type}
        seq_dict["isolate"] = isolate
        otu_updates.append(seq_dict)

    return otu_updates, auto_excluded
