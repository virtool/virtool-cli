import shutil
from pathlib import Path
from uuid import UUID
from typing import TYPE_CHECKING

from pydantic import BaseModel, UUID4, Field
from orjson import orjson
from bidict import bidict

if TYPE_CHECKING:
    from virtool_cli.ref.repo import EventSourcedRepo
from virtool_cli.ref.resources import (
    EventSourcedRepoOTU,
    EventSourcedRepoIsolate,
    EventSourcedRepoSequence,
)


class IndexData(BaseModel):
    taxid: int
    otu_id: UUID4 = Field(validation_alias="id")


class OTUSnapshotCache:
    def __init__(self, path: Path):
        self.path = path

        self.path.mkdir(exist_ok=True)

        self._metadata_path = self.path / "otu.json"

        self._toc_path = self.path / "toc.json"

        self._data_path = self.path / "data"

        self._data_path.mkdir(exist_ok=True)

    def clean(self):
        shutil.rmtree(self.path)
        self.path.mkdir()
        self._data_path.mkdir()

    def cache(self, otu: "EventSourcedRepoOTU", options=None):
        otu_dict = otu.dict(exclude_contents=True)

        with open(self._metadata_path, "wb") as f:
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

        # if otu.schema:
        #     with open(self.path / "schema.json", "wb") as f:
        #         f.write(orjson.dumps(otu.schema, option=options))
        #
        # if otu.excluded_accessions:
        #     with open(self.path / "toc_accessions_excluded.json", "wb") as f:
        #         f.write(orjson.dumps(otu.excluded_accessions))

        # write full data
        for isolate in otu.isolates:
            with open(self._data_path / f"{isolate.id}.json", "wb") as f:
                f.write(
                    orjson.dumps(isolate.dict(exclude_contents=True), option=options)
                )

            for sequence in isolate.sequences:
                with open(self._data_path / f"{sequence.id}.json", "wb") as f:
                    f.write(orjson.dumps(sequence.dict(), option=options))

    def load_metadata(self) -> dict:
        with open(self._metadata_path, "rb") as f:
            otu_dict = orjson.loads(f.read())

        return otu_dict

    def load_full(self) -> "EventSourcedRepoOTU":
        with open(self._metadata_path, "rb") as f:
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


class SnapshotCache:
    def __init__(self, repo: "EventSourcedRepo"):
        self._repo = repo

        self.path = self._repo.path / ".cache" / "snapshot"

        self.path.mkdir(exist_ok=True)

        self._meta_path = self.path / "meta.json"
        self._index_path = self.path / "index.json"

        self._excluded_files = [
            self._meta_path.name,
            self._index_path.name,
            "Thumbs.db",
            ".DS_Store",
        ]

    @property
    def index(self) -> bidict:
        """A bidirectional dictionary of OTUs by Taxonomy ID: Virtool UUID.

        To retrieve using Virtool UUID, use index.inverse"""
        if self._index_path.exists():
            return self._load_index()

        return self._build_index()

    def clean(self):
        shutil.rmtree(self.path)
        self.path.mkdir(exist_ok=True)

    def snapshot(self, indent: bool = False):
        """Serializes the current reference data as a snapshot"""
        self.clean()

        options = orjson.OPT_INDENT_2 if indent else None

        reference = {
            "id": self._repo.meta.id,
            "created_at": self._repo.meta.created_at,
            "data_type": self._repo.meta.data_type,
            "name": self._repo.meta.name,
            "organism": self._repo.meta.organism,
        }

        with open(self.path / "meta.json", "wb") as f:
            f.write(orjson.dumps(reference, option=options))

        index = bidict()
        for otu in self._repo.iter_otus():
            self.cache_otu(otu, options=options)
            index[otu.taxid] = otu.id

        self._cache_index(index)

    def cache_otu(self, otu: "EventSourcedRepoOTU", options):
        """Snapshots a single OTU"""
        otu_snap = OTUSnapshotCache(self.path / f"{otu.taxid}")
        otu_snap.cache(otu, options)

    def load_otu(self, taxid: int):
        """Loads an OTU from the most recent repo snapshot"""
        otu_snap = OTUSnapshotCache(self.path / f"{taxid}")

        return otu_snap.load_full()

    def load_otu_metadata(self, taxid: int) -> dict:
        otu_metadata_path = self.path / f"{taxid}" / "otu.json"

        if otu_metadata_path.exists():
            with open(otu_metadata_path, "rb") as f:
                return orjson.loads(f.read())

        raise FileNotFoundError("OTU not found in snapshot")

    def _get_otu_path(self, taxid: int) -> Path:
        return self.path / f"{taxid}"

    def _build_index(self) -> bidict:
        index_by_id = bidict()

        for taxid in [
            int(path.stem)
            for path in self.path.iterdir()
            if path.suffix == ".json" and path.name not in self._excluded_files
        ]:
            index_data = IndexData(**self.load_otu_metadata(taxid))
            index_by_id[index_data.taxid] = index_data.otu_id

        return index_by_id

    def _cache_index(self, index: bidict):
        """Cache the index and taxid bidict as a dictionary with OTU ID as the key"""
        index_by_id = {str(key): index.inverse[key] for key in index.inverse}

        with open(self._index_path, "wb") as f:
            f.write(orjson.dumps(index_by_id))

    def _load_index(self) -> bidict:
        """Serialize the OTU ID, Taxonomy ID dictionary and returns as bidict"""
        with open(self._index_path, "rb") as f:
            index_dict = orjson.loads(f.read())

        return bidict({UUID(key): index_dict[key] for key in index_dict}).inverse
