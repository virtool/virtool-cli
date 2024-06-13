import shutil
from pathlib import Path
from uuid import UUID

from orjson import orjson

from virtool_cli.ref.resources import (
    EventSourcedRepoOTU,
    EventSourcedRepoSequence,
    EventSourcedRepoIsolate,
)


class OTUSnapshot:
    """Manages snapshot data for a single OTU."""

    def __init__(self, path: Path):
        self.path = path
        """The path of this snapshot's directory."""

        self._data_path = self.path / "data"
        """The path of this snapshot's data directory."""

        if not self.path.exists():
            self.path.mkdir()

        if not self._data_path.exists():
            self._data_path.mkdir()

        self._metadata_path = self.path / "metadata.json"
        """The path of this snapshot's metadata file."""

        self.at_event = self._get_at_event()

    @property
    def _otu_path(self):
        return self.path / "otu.json"

    @property
    def _toc_path(self):
        return self.path / "toc.json"

    def _get_at_event(self) -> int | None:
        metadata_dict = self._read_metadata()
        return metadata_dict.get("at_event")

    def _read_metadata(self) -> dict:
        if self._metadata_path.exists():
            with open(self._metadata_path, "rb") as f:
                metadata_dict = orjson.loads(f.read())

            return metadata_dict

        return {}

    def clean(self):
        shutil.rmtree(self.path)
        self.path.mkdir()
        self._data_path.mkdir()

    def cache(
        self, otu: "EventSourcedRepoOTU", at_event: int | None = None, options=None
    ):
        otu_dict = otu.dict(exclude_contents=True)

        with open(self._otu_path, "wb") as f:
            f.write(orjson.dumps(otu_dict, option=options))

        toc = {
            f"{isolate.name}": {
                "id": isolate.id,
                "accessions": {
                    accession: isolate.get_sequence_by_accession(accession).id
                    for accession in sorted(isolate.accessions)
                },
            }
            for isolate in otu.isolates
        }
        with open(self._toc_path, "wb") as f:
            f.write(orjson.dumps(toc, option=options))

        for isolate in otu.isolates:
            with open(self._data_path / f"{isolate.id}.json", "wb") as f:
                f.write(
                    orjson.dumps(isolate.dict(exclude_contents=True), option=options)
                )

            for sequence in isolate.sequences:
                with open(self._data_path / f"{sequence.id}.json", "wb") as f:
                    f.write(orjson.dumps(sequence.dict(), option=options))

        with open(self._metadata_path, "wb") as f:
            f.write(orjson.dumps({"at_event": at_event}, option=options))

    def cache_isolate(self, isolate, options=None):
        with open(self._toc_path, "rb") as f:
            toc = orjson.loads(f.read())
        toc[f"{isolate.name}"] = {
            "id": isolate.id,
            "accessions": {
                accession: isolate.get_sequence_by_accession(accession).id
                for accession in sorted(isolate.accessions)
            },
        }
        with open(self._toc_path, "wb") as f:
            f.write(orjson.dumps(toc, option=options))

        with open(self._data_path / f"{isolate.id}.json", "wb") as f:
            f.write(orjson.dumps(isolate.dict(exclude_contents=True), option=options))

    def cache_sequence(self, sequence, isolate_id: UUID, options=None):
        with open(self._toc_path, "rb") as f:
            toc = orjson.loads(f.read())

        for isolate_name in toc:
            if toc[isolate_name]["id"] == isolate_id:
                toc[isolate_name]["accessions"].append(sequence.accession)
                toc[isolate_name]["accessions"].sort()

        with open(self._toc_path, "wb") as f:
            f.write(orjson.dumps(toc, option=options))

        with open(self._data_path / f"{sequence.id}.json", "wb") as f:
            f.write(orjson.dumps(sequence.dict(), option=options))

    def load(self) -> "EventSourcedRepoOTU":
        """Load an OTU from the snapshot."""
        with open(self._otu_path, "rb") as f:
            otu_dict = orjson.loads(f.read())

        with open(self._toc_path, "rb") as f:
            toc = orjson.loads(f.read())

        isolates = []
        for key in toc:
            isolate_id = UUID(toc[key]["id"])

            with open(self._data_path / f"{isolate_id}.json", "rb") as f:
                isolate_dict = orjson.loads(f.read())

            sequences = []

            for accession in toc[key]["accessions"]:
                sequence_id = UUID(toc[key]["accessions"][accession])

                with open(self._data_path / f"{sequence_id}.json", "rb") as f:
                    sequence_dict = orjson.loads(f.read())

                sequence = EventSourcedRepoSequence(**sequence_dict)

                sequences.append(sequence)

            isolate_dict["uuid"] = UUID(isolate_dict.pop("id"))
            isolate_dict["sequences"] = sequences

            isolate = EventSourcedRepoIsolate(**isolate_dict)

            isolates.append(isolate)

        otu_dict["uuid"] = UUID(otu_dict.pop("id"))
        otu_dict["isolates"] = isolates

        return EventSourcedRepoOTU(**otu_dict)
