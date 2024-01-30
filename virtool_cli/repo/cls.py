"""An access layer for a Virtool reference repository.

This is a work in progress.

TODO: Check if excluded accessions exist in the repo.
TODO: Check that OTUs have only one default isolate.
TODO: Support writing sequences.
TODO: Support writing isolates.

"""

import json
from dataclasses import dataclass
from pathlib import Path

from virtool_cli.utils.id_generator import (
    generate_random_alphanumeric,
)
from virtool_cli.utils.reference import generate_otu_dirname, is_v1


@dataclass
class RepoSequence:
    """Represents a sequence in a Virtool reference repository."""

    id: str
    """The sequence id."""

    accession: str
    """The sequence accession."""

    definition: str
    """The sequence definition."""

    host: str
    """The sequence host."""

    sequence: str
    """The sequence."""

    segment: str
    """The sequence segment."""


class RepoIsolate:
    """Represents an isolate in a Virtool reference repository."""

    def __init__(
        self,
        id: str,
        default: bool,
        otu: "RepoOTU",
        path: Path,
        source_name: str,
        source_type: str,
    ):
        self.id = id
        """The isolate id."""

        self.default = default
        """Indicates whether the isolate is the Virtool default for the OTU."""

        self.otu = otu
        """The isolate's parent OTU."""

        self.path = path
        """The path to the isolate directory."""

        self.source_name = source_name
        """The isolate's source name."""

        self.source_type = source_type
        """The isolate's source type."""

        self.sequences = []

        for sequence_path in (
            p for p in self.path.iterdir() if p.name != "isolate.json"
        ):
            with open(sequence_path) as f:
                sequence_data = json.load(f)

            self.sequences.append(
                RepoSequence(
                    id=sequence_data["_id"],
                    accession=sequence_data["accession"],
                    definition=sequence_data["definition"],
                    host=sequence_data["host"],
                    sequence=sequence_data["sequence"],
                    segment=sequence_data["segment"],
                ),
            )

        self.sequences.sort(key=lambda s: s.accession)

    @property
    def name(self) -> str:
        """The formatted isolate name."""
        return f"{self.source_type} {self.source_name}"


class RepoOTU:
    def __init__(
        self,
        id_: str,
        path: Path,
    ):
        self.id = id_
        """The OTU id."""

        self.isolates: list[RepoIsolate] = []
        """A list of child isolates."""

        self.path = path
        """The path to the OTU directory."""

        self._isolate_name_id_map: dict[tuple[str, str], str] = {}
        """Maps isolate source type-name combos to isolate ids."""

        for isolate_path in [p for p in self.path.iterdir() if p.is_dir()]:
            with open(isolate_path / "isolate.json") as f:
                isolate_data = json.load(f)

            isolate = RepoIsolate(
                id=isolate_data["id"],
                default=isolate_data["default"],
                path=isolate_path,
                otu=self,
                source_name=isolate_data["source_name"],
                source_type=isolate_data["source_type"],
            )

            self.isolates.append(isolate)

            self._isolate_name_id_map[
                (isolate.source_type, isolate.source_name)
            ] = isolate.id

        with open(self.path / "otu.json") as f:
            data = json.load(f)

            if self.id != data["_id"]:
                raise ValueError("otu_id does not match otu.json _id")

            self.abbreviation = data["abbreviation"]
            """The OTU abbreviation (eg. TMV)."""

            self.name = data["name"]
            """The OTU name (eg. Tobacco mosaic virus)."""

            self.schema = data["schema"]

            self.taxid: int = int(data["taxid"])
            """The OTU taxonomy ID."""

        with open(self.path / "exclusions.json") as f:
            self.exclusions = json.load(f)
            """
            A list of accessions that have been blacklisted for this OTU.

            They should not be fetched again or nested in this OTU.
            """

    @property
    def blocked_accessions(self):
        existing_accessions = []

        for isolate in self.isolates:
            for sequence in isolate.sequences:
                existing_accessions.append(sequence.accession)

        return sorted(existing_accessions + self.exclusions)

    def get_isolate_by_id(self, isolate_id: str) -> RepoIsolate:
        """Get an isolate by its id."""
        for isolate in self.isolates:
            if isolate.id == isolate_id:
                return isolate

        raise ValueError(f"isolate_id {isolate_id} not found")

    def get_isolate_by_name(
        self,
        source_type: str,
        source_name: str,
    ) -> RepoIsolate | None:
        """Get an isolate by its source type and name."""
        try:
            isolate_id = self._isolate_name_id_map[(source_type, source_name)]
        except KeyError:
            return None

        return self.get_isolate_by_id(isolate_id)

    def update(
        self,
        abbreviation: str | None = None,
        exclusions: list[str] | None = None,
        name: str | None = None,
    ):
        json_path = self.path / "otu.json"

        with open(json_path) as f:
            otu = json.load(f)

            if abbreviation is not None:
                otu["abbreviation"] = abbreviation

            if name is not None:
                otu["name"] = name

        with open(json_path, "w") as f:
            json.dump(otu, f, indent=4, sort_keys=True)

        with open(self.path / "exclusions.json", "w") as f:
            json.dump(exclusions, f)


class RepoMaps:
    def __init__(self, path: Path):
        if is_v1(path / "src"):
            raise ValueError(
                "This reference database is a deprecated v1 reference. Run 'virtool ref migrate' before trying again.",
            )

        self.otu_id_to_path: dict[str, Path] = {}
        """Maps otu ids to their respective paths."""

        self.sequence_id_to_otu_id: dict[str, str] = {}
        """Maps sequence ids to parent otu ids."""

        self.sequence_id_to_path: dict[str, Path] = {}
        """Maps sequence ids to their respective paths."""

        self.taxid_to_otu_id: dict[int, str] = {}
        """Maps taxids to otu ids."""

        for otu_path in (p for p in (path / "src").iterdir() if p.is_dir()):
            if not otu_path.is_dir():
                continue

            name, otu_id = otu_path.name.split("--")

            self.otu_id_to_path[otu_id] = otu_path

            with open(otu_path / "otu.json") as f:
                taxid = int(json.load(f)["taxid"])

            self.taxid_to_otu_id[taxid] = otu_id

            for isolate_path in (p for p in otu_path.iterdir() if p.is_dir()):
                for sequence_path in [
                    p for p in isolate_path.iterdir() if p.name != "isolate.json"
                ]:
                    sequence_id = sequence_path.stem

                    self.sequence_id_to_otu_id[sequence_id] = otu_id
                    self.sequence_id_to_path[sequence_id] = sequence_path


class Repo:
    def __init__(self, path: Path):
        self.path = path
        """The path to the repo directory."""

        with open(self.path / "src" / "meta.json") as f:
            data = json.load(f)

            self.data_type = data["data_type"]
            """The sequence data type of the repo (eg. genome)."""

            self.organism = data["organism"]
            """The type of organism organized in the reference (eg. virus)."""

        self.maps = RepoMaps(self.path)

    @property
    def src_path(self) -> Path:
        """The path to the repo src directory."""
        return self.path / "src"

    def get_otu_path(self, otu_id: str) -> Path:
        """Get the path to the otu with the given ``otu_id``."""
        try:
            return self.maps.otu_id_to_path[otu_id]
        except KeyError:
            raise ValueError(f"No OTU with ID {otu_id}")

    def get_sequence_path(self, sequence_id: str) -> Path:
        """Get the path to the sequence with the given ``sequence_id``."""
        try:
            return self.maps.sequence_id_to_path[sequence_id]
        except KeyError:
            raise ValueError(f"sequence_id {sequence_id} not found")

    def get_otu_by_id(self, otu_id: str):
        """Get an OTU by its ID."""
        try:
            path = self.get_otu_path(otu_id)
        except ValueError:
            raise ValueError(f"No OTU with ID {otu_id}")

        return RepoOTU(
            otu_id,
            path,
        )

    def get_otu_by_taxid(self, taxid: int) -> RepoOTU:
        try:
            otu_id = self.maps.taxid_to_otu_id[taxid]
        except KeyError:
            raise ValueError(f"No OTU with taxid {taxid}")

        return self.get_otu_by_id(otu_id)

    def create_otu(self, name: str, taxid: int, abbreviation: str = "") -> RepoOTU:
        """Create a new OTU."""
        otu_id = generate_random_alphanumeric(
            length=8,
            excluded=list(self.maps.otu_id_to_path.keys()),
        )

        otu_path = self.src_path / generate_otu_dirname(name, otu_id=otu_id)

        self.maps.otu_id_to_path[otu_id] = otu_path
        self.maps.taxid_to_otu_id[taxid] = otu_id

        try:
            otu_path.mkdir()

            with open(otu_path / "otu.json", "w") as f:
                json.dump(
                    {
                        "_id": otu_id,
                        "abbreviation": abbreviation,
                        "name": name,
                        "schema": [],
                        "taxid": taxid,
                    },
                    f,
                    indent=4,
                    sort_keys=True,
                )

            with open(otu_path / "exclusions.json", "w") as f:
                json.dump([], f, indent=4, sort_keys=True)
        except Exception as e:
            try:
                otu_path.rmdir()
            except FileNotFoundError:
                pass

            try:
                del self.maps.otu_id_to_path[otu_id]
            except KeyError:
                pass

            try:
                del self.maps.taxid_to_otu_id[taxid]
            except KeyError:
                pass

            raise e

        return self.get_otu_by_id(otu_id)

    def iter_otus(self):
        """Iterate over all OTUs."""
        for otu_id in self.maps.otu_id_to_path:
            yield self.get_otu_by_id(otu_id)
