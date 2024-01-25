import json
from pathlib import Path
from urllib.error import HTTPError

import click

from virtool_cli.utils.cache import generate_taxid_table
from virtool_cli.utils.id_generator import generate_unique_ids
from virtool_cli.utils.ncbi import fetch_taxonomy_record
from virtool_cli.utils.reference import generate_otu_dirname, get_unique_otu_ids, is_v1

OTU_KEYS = ["_id", "name", "abbreviation", "schema", "taxid"]


async def add_otu(taxid: int, path: Path):
    """Fetch NCBI Taxonomy data for a given taxon ID number and
    write the OTU to the reference directory at src_path

    :param taxid: an NCBI taxonomy id
    :param path: the path to the reference repository
    """
    click.echo("Adding OTU...")

    src_path = path / "src"

    if is_v1(src_path):
        click.echo("Repository is a deprecated v1 reference.", err=True)

    taxid_table = generate_taxid_table(src_path)

    try:
        taxonomy_data = await fetch_taxonomy_record(taxon_id=taxid)
        if not taxonomy_data:
            click.echo(
                "Could not find a record under this taxon ID on NCBI Taxonomy",
                err=True,
            )
            return

    except HTTPError:
        click.echo(
            "Could not find a record under this taxon ID on NCBI Taxonomy",
            err=True,
        )

    if taxid in taxid_table:
        click.echo(
            "An OTU with this taxonomy ID already exists in the reference.",
            err=True,
        )
        return

    new_id = generate_unique_ids(1, excluded=await get_unique_otu_ids(src_path)).pop()

    new_otu = generate_otu(taxonomy_data, new_id)

    if not new_otu:
        click.echo("Received invalid data from NCBI Taxonomy DB", err=True)

    otu_path = await write_otu(new_otu, src_path)

    if otu_path is None:
        click.echo("Could not create OTU directory", err=True)
        return

    click.echo(f"OTU {new_otu['name']} added to reference.")
    click.echo(
        "The OTU directory is empty. Use virtool ref add accessions before building, updating or submitting.",
    )


def generate_otu(taxonomy_data: dict, new_id: str) -> dict:
    """Create a new OTU from a given NCBI Taxonomy unique Id (taxid) and
    a pre-generated Virtool OTU Id

    :param taxonomy_data: Dictionary containing fetched NCBI Taxonomy metadata
    :param new_id: Pre-generated Virtool OTU Id
    :return: The deserialized contents of an OTU metadata file
    """
    for key in ["ScientificName", "Id"]:
        if key not in taxonomy_data:
            return {}

    otu = {
        "_id": new_id,
        "name": taxonomy_data.get("ScientificName", ""),
        "abbreviation": "",
        "schema": [],
        "taxid": int(taxonomy_data.get("Id", "-1")),
    }

    return otu


async def write_otu(otu: dict, src_path: Path) -> Path:
    """Generate a new directory for given OTU metadata, store the metadata under in otu.json
    and return the path to the new directory.

    :param otu: Dict of OTU metadata
    :param src_path: Path to a reference directory
    :return: Path of the new OTU directory under src_path
    """
    dirname = generate_otu_dirname(name=otu["name"], otu_id=otu["_id"])

    otu_path = src_path / dirname
    otu_path.mkdir()

    with open(otu_path / "otu.json", "w") as f:
        json.dump(otu, f, indent=4, sort_keys=True)

    with open(otu_path / "exclusions.json", "w") as f:
        json.dump([], f, indent=4, sort_keys=True)

    return otu_path
