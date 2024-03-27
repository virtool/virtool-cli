import json
from pathlib import Path

import arrow

from virtool_cli.ref.legacy import Repo


def build_json(indent: bool, output_path: Path, path: Path, version: str):
    """Build a Virtool reference JSON file from a data directory.

    :param path: The path to a reference repository
    :param output_path: The path to write the output JSON file to
    :param indent: write indented json when enabled
    :param version: the version string to include in the reference.json file
    """
    repo = Repo(path)

    otus = []

    for otu in repo.iter_otus():
        isolates = []

        for isolate in otu.isolates:
            sequences = [
                {
                    "_id": sequence.id,
                    "accession": sequence.accession,
                    "definition": sequence.definition,
                    "host": sequence.host,
                    "segment": sequence.segment,
                    "sequence": sequence.sequence,
                }
                for sequence in isolate.sequences
            ]

            isolates.append(
                {
                    "id": isolate.id,
                    "default": isolate.default,
                    "sequences": sequences,
                    "source_name": isolate.source_name,
                    "source_type": isolate.source_type,
                },
            )

        isolates.sort(key=lambda x: x["id"])

        otus.append(
            {
                "_id": otu.id,
                "abbreviation": otu.abbreviation,
                "isolates": isolates,
                "name": otu.name,
                "schema": otu.schema,
                "taxid": otu.taxid,
            },
        )

    otus.sort(key=lambda x: x["name"])

    with open(output_path, "w") as f:
        json.dump(
            {
                "created_at": arrow.utcnow().isoformat(),
                "data_type": repo.data_type,
                "name": version,
                "organism": repo.organism,
                "otus": otus,
            },
            f,
            indent=4 if indent else None,
            sort_keys=True,
        )
