"""Code for the old reference format, not based on event sourcing."""

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import click

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

    isolate: "RepoIsolate"
    """The parent isolate."""

    sequence: str
    """The sequence."""

    segment: str
    """The sequence segment."""


@dataclass
class RepoIsolate:
    """Represents an isolate in a Virtool reference repository."""

    def __init__(
        self,
        id_: str,
        default: bool,
        otu: "RepoOTU",
        path: Path,
        source_name: str,
        source_type: str,
    ):
        self.id = id_
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

        self.sequences: list[RepoSequence] = []
        """A list of child sequences."""

        for sequence_path in (
            p for p in self.path.iterdir() if p.name != "isolate.json"
        ):
            with open(sequence_path) as f:
                sequence_data = json.load(f)
                sequence = RepoSequence(
                    sequence_data["_id"],
                    sequence_data["accession"],
                    sequence_data["definition"],
                    sequence_data["host"],
                    self,
                    sequence_data["sequence"],
                    sequence_data["segment"],
                )
            self.sequences.append(sequence)
        self.sequences.sort(key=lambda s: s.accession)

    @property
    def name(self) -> str:
        """The isolate name."""
        return f"{self.source_type} {self.source_name}"

    def to_dict(self) -> dict[str, Any]:
        sequences = []

        for sequence in self.sequences:
            sequence_dict = asdict(sequence)
            sequence_dict.pop("isolate")

        return {
            "id": self.id,
            "default": self.default,
            "sequences": sequences,
            "source_name": self.source_name,
            "source_type": self.source_type,
        }

    def add_sequence(
        self,
        accession: str,
        definition: str,
        host: str,
        segment: str,
        sequence: str,
    ) -> RepoSequence | None:
        """Add a sequence to the isolate."""
        if accession in self.otu.exclusions:
            click.echo(f"accession {accession} is in exclusions list", err=True)
            return None

        existing_accessions = {
            a.split(".")[0]
            for a in self.otu.repo.maps.sequence_id_to_accession.values()
        }

        if accession.split(".")[0] in existing_accessions:
            click.echo(f"accession {accession} already exists in repo", err=True)
            return None

        sequence_id = generate_random_alphanumeric(
            length=8,
            excluded=set(self.otu.repo.maps.sequence_id_to_path.keys()),
        )

        sequence_path = self.path / f"{sequence_id}.json"

        with open(sequence_path, "w") as f:
            json.dump(
                {
                    "_id": sequence_id,
                    "accession": accession,
                    "definition": definition,
                    "host": host,
                    "segment": segment,
                    "sequence": sequence,
                },
                f,
                indent=4,
                sort_keys=True,
            )

        sequence = RepoSequence(
            sequence_id,
            accession,
            definition,
            host,
            self,
            sequence,
            segment,
        )
        self.otu.repo.maps.sequence_id_to_otu_id[sequence_id] = self.otu.id
        self.otu.repo.maps.sequence_id_to_path[sequence_id] = sequence_path

        self.sequences.append(sequence)

        return sequence


@dataclass
class RepoOTU:
    """Represents an OTU in a Virtool reference repository."""

    def __init__(self, id_: str, path: Path, repo: "Repo"):
        self.id = id_
        """The OTU id."""

        self.path = path
        """The path to the OTU directory."""

        self.repo = repo
        """The OTU's parent repository."""

        self.exclusions = []
        """
        A list of accessions that have been blacklisted for this OTU.

        They should not be fetched again or nested in this OTU.
        """

        self.isolates = []
        """A list of child isolates."""

        self._isolate_name_id_map: dict[tuple[str, str], str] = {}
        """Maps isolate source type-name combos to isolate ids."""

        self.abbreviation = ""
        """The OTU abbreviation (eg. TMV)."""

        self.name = ""
        """The OTU name (eg. Tobacco mosaic virus)."""

        self.schema = []

        self.taxid: int = -1
        """The OTU taxonomy ID."""

        for isolate_path in [p for p in self.path.iterdir() if p.is_dir()]:
            with open(isolate_path / "isolate.json") as f:
                isolate_data = json.load(f)

            isolate = RepoIsolate(
                isolate_data["id"],
                isolate_data["default"],
                self,
                isolate_path,
                isolate_data["source_name"],
                isolate_data["source_type"],
            )

            self.isolates.append(isolate)

            self._isolate_name_id_map[(isolate.source_type, isolate.source_name)] = (
                isolate.id
            )

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

    def to_dict(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "abbreviation": self.abbreviation,
            "name": self.name,
            "schema": self.schema,
            "taxid": self.taxid,
            "isolates": [i.to_dict() for i in self.isolates],
            "exclusions": self.exclusions,
        }

    @property
    def blocked_accessions(self):
        existing_accessions = []

        for isolate in self.isolates:
            for sequence in isolate.sequences:
                existing_accessions.append(sequence.accession)

        return sorted(existing_accessions + self.exclusions)

    def add_isolate(self, source_type: str, source_name: str):
        """Add a new isolate to the OTU."""
        isolate_id = generate_random_alphanumeric(
            length=8,
            excluded=set(self.repo.maps.isolate_id_to_otu_id.keys()),
        )

        isolate_path = self.path / isolate_id
        isolate_path.mkdir()

        with open(isolate_path / "isolate.json", "w") as f:
            json.dump(
                {
                    "id": isolate_id,
                    "default": False,
                    "source_name": source_name,
                    "source_type": source_type,
                },
                f,
                indent=4,
                sort_keys=True,
            )

        isolate = RepoIsolate(
            isolate_id,
            False,
            self,
            isolate_path,
            source_name,
            source_type,
        )

        self.isolates.append(isolate)

        self._isolate_name_id_map[(source_type, source_name)] = isolate_id

        return isolate

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
        self.isolate_id_to_otu_id: dict[str, str] = {}
        """Maps isolate ids to parent otu ids."""

        self.otu_id_to_path: dict[str, Path] = {}
        """Maps otu ids to their respective paths."""

        self.sequence_id_to_accession: dict[str, str] = {}
        """Maps sequence ids to their respective accessions."""

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
                self.isolate_id_to_otu_id[isolate_path.stem] = otu_id

                for sequence_path in [
                    p for p in isolate_path.iterdir() if p.name != "isolate.json"
                ]:
                    sequence_id = sequence_path.stem

                    with open(sequence_path) as f:
                        self.sequence_id_to_accession[sequence_id] = json.load(f)[
                            "accession"
                        ]

                    self.sequence_id_to_otu_id[sequence_id] = otu_id
                    self.sequence_id_to_path[sequence_id] = sequence_path


class Repo:
    def __init__(self, path: Path):
        self.path = path
        """The path to the repo directory."""

        if is_v1(self.path / "src"):
            raise ValueError("Repository is a deprecated v1 reference.")

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
            self,
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
