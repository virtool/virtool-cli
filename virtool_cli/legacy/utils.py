from dataclasses import dataclass
from pathlib import Path
from typing import Generator

import orjson
from orjson import orjson
from pydantic_core import ErrorDetails

from virtool_cli.ncbi.client import NCBIClient


class HandleErrorContext:
    """A context that is passed to a ``ValidationError`` error handler.

    This context is for one error from ``ValidationError.errors``, not for the
    ``ValidationError`` itself.

    """

    def __init__(
        self,
        error: ErrorDetails,
        fix: bool,
        ncbi_client: NCBIClient,
        otu: dict,
        repaired_otu: dict,
    ):
        self.error = error
        """The details of the Pydantic validation error."""

        self.fix = fix
        """Whether to attempt to fix the error or not."""

        self.ncbi_client = ncbi_client
        """An instance of the NCBI client."""

        self.otu = otu
        """A dictionary representation of the OTU."""

        self.repaired_otu = repaired_otu
        """A copy of the OTU that repairs should be applied to."""

    def update_otu(self, update: dict):
        """Update OTU associated with the error.

        Example:
        -------
        .. code-block:: python

            ctx.update_otu({"name": "new_name"})

        :param update: the update dict

        """
        self.repaired_otu.update(update)

    def update_isolate(self, isolate_index: int, update: dict):
        """Update a legacy isolate on disk given an update dict.

        Example:
        -------
        .. code-block:: python

           ctx.update_isolate(
               isolate_index,
              {"source_name": source.name, "source_type": source.type.lower()},
           )

        :param isolate_index: the index of the isolate in the OTU
        :param update: the update dict

        """
        self.repaired_otu["isolates"][isolate_index].update(update)

    def update_sequence(self, isolate_index: int, sequence_index: int, update: dict):
        """Update a legacy sequence on disk given an update dict.

        Example:
        -------
        .. code-block:: python

            ctx.update_sequence(
                isolate_index,
                sequence_index,
                {"sequence": "ATGCATCGTATGCTYACTAGCTAC"},
            )

        """
        self.repaired_otu["isolates"][isolate_index]["sequences"][
            sequence_index
        ].update(
            update,
        )


@dataclass
class ErrorHandledResult:
    """The return value for a legacy validation ``ValidationError`` handler."""

    message: str
    """The message to display to the user.

    Includes ``rich``-compatible markup.
    """

    fixed: bool = False
    """Whether the error was fixed.

    This value is used to show the user whether the error was fixed or not.
    """


def build_legacy_otu(path: Path) -> dict:
    """Build a dict-based representation of an OTU from an OTU directory.

    :param path: the path to the OTU directory
    :return: a dict-based representation of the OTU

    """
    with open(path / "otu.json") as f:
        otu_data = orjson.loads(f.read())

    otu_data["isolates"] = []

    for isolate_dir_path in sorted(path.iterdir()):
        if isolate_dir_path.name == "otu.json":
            continue

        with open(isolate_dir_path / "isolate.json") as f:
            isolate_data = orjson.loads(f.read())

        isolate_data["sequences"] = []

        for sequence_path in sorted(isolate_dir_path.iterdir()):
            if sequence_path.name == "isolate.json":
                continue

            with open(sequence_path) as f:
                sequence_data = orjson.loads(f.read())

            isolate_data["sequences"].append(sequence_data)

        otu_data["isolates"].append(isolate_data)

    return otu_data


def replace_otu(path: Path, otu: dict):
    """Replace the otu in the path with the data in ``otu``.

    :param path: the path to the otu
    :param otu: the new otu data

    """
    if not path.parent.is_dir():
        raise FileNotFoundError(f"No such OTU parent directory: {path.parent}")

    if not path.is_dir():
        raise FileNotFoundError(f"No such OTU directory: {path}")

    with open(path / "otu.json", "wb") as f:
        f.write(
            orjson.dumps(
                {
                    "_id": otu["_id"],
                    "name": otu["name"],
                    "abbreviation": otu["abbreviation"],
                    "schema": otu["schema"],
                    "taxid": otu["taxid"],
                },
                option=orjson.OPT_INDENT_2,
            ),
        )

    for isolate in sorted(otu["isolates"], key=lambda x: x["id"]):
        isolate_id = isolate["id"]

        with open(path / isolate_id / "isolate.json", "wb") as f:
            f.write(
                orjson.dumps(
                    {
                        "id": isolate["id"],
                        "default": isolate["default"],
                        "source_name": isolate["source_name"],
                        "source_type": isolate["source_type"],
                    },
                    option=orjson.OPT_INDENT_2,
                ),
            )

        for sequence in sorted(isolate["sequences"], key=lambda x: x["_id"]):
            with open(path / isolate_id / f"{sequence['_id']}.json", "wb") as f:
                f.write(
                    orjson.dumps(
                        {
                            "_id": sequence["_id"],
                            "accession": sequence["accession"],
                            "definition": sequence["definition"],
                            "host": sequence["host"],
                            "segment": sequence["segment"],
                            "sequence": sequence["sequence"],
                        },
                        option=orjson.OPT_INDENT_2,
                    ),
                )


def iter_legacy_otus(src_path: Path) -> Generator[dict, None, None]:
    """Iterate over all OTUs in a legacy repository.

    :param src_path: The path to the legacy repository.
    :return: a generator that yields OTU data dictionaries.
    """
    for alpha_dir_path in sorted(src_path.iterdir()):
        if alpha_dir_path.name == "meta.json":
            continue

        for otu_dir_path in sorted(alpha_dir_path.iterdir()):
            yield build_legacy_otu(otu_dir_path)
