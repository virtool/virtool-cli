from pathlib import Path
from uuid import UUID

import orjson
import pytest

from virtool_cli.ref.repo import EventSourcedRepo
from virtool_cli.ref.resources import (
    EventSourcedRepoIsolate,
    EventSourcedRepoOTU,
    EventSourcedRepoSequence,
)
from virtool_cli.ref.utils import DataType


@pytest.fixture()
def _repo(tmp_path: Path):
    return EventSourcedRepo.new(
        DataType.GENOME,
        "Generic Viruses",
        tmp_path / "test_repo",
        "virus",
    )


def test_new(_repo: EventSourcedRepo, tmp_path: Path):
    """Test that creating a new ``Repo`` object returns the expected object and creates
    the expected directory structure.
    """
    assert _repo.path == tmp_path / "test_repo"
    assert _repo.last_id == 1

    assert _repo.meta.data_type == DataType.GENOME
    assert _repo.meta.name == "Generic Viruses"
    assert _repo.meta.organism == "virus"


def test_create_otu(_repo: EventSourcedRepo):
    """Test that creating an OTU returns the expected ``RepoOTU`` object and creates
    the expected event file.
    """
    otu = _repo.create_otu("TMV", "Tobacco mosaic virus", [], 12242)

    assert otu == EventSourcedRepoOTU(
        id=otu.id,
        acronym="TMV",
        excluded_accessions=[],
        isolates=[],
        name="Tobacco mosaic virus",
        schema=[],
        taxid=12242,
    )

    with open(_repo.path.joinpath("src", "00000002.json")) as f:
        event = orjson.loads(f.read())

    del event["timestamp"]

    assert event == {
        "data": {
            "id": str(otu.id),
            "acronym": "TMV",
            "excluded_accessions": [],
            "name": "Tobacco mosaic virus",
            "rep_isolate": None,
            "schema": [],
            "taxid": 12242,
        },
        "id": 2,
        "query": {
            "otu_id": str(otu.id),
        },
        "type": "CreateOTU",
    }

    assert _repo.last_id == 2


def test_create_isolate(_repo: EventSourcedRepo):
    """Test that creating an isolate returns the expected ``RepoIsolate`` object and
    creates the expected event file.
    """
    otu = _repo.create_otu("TMV", "Tobacco mosaic virus", [], 12242)

    isolate = _repo.create_isolate(otu.id, "A", "isolate")

    assert isinstance(isolate.id, UUID)
    assert isolate.sequences == []
    assert isolate.source_name == "A"
    assert isolate.source_type == "isolate"

    with open(_repo.path.joinpath("src", "00000003.json")) as f:
        event = orjson.loads(f.read())

    del event["timestamp"]

    assert event == {
        "data": {
            "id": str(isolate.id),
            "source_name": "A",
            "source_type": "isolate",
        },
        "id": 3,
        "query": {
            "otu_id": str(otu.id),
            "isolate_id": str(isolate.id),
        },
        "type": "CreateIsolate",
    }

    assert _repo.last_id == 3


def test_create_sequence(_repo: EventSourcedRepo):
    """Test that creating a sequence returns the expected ``RepoSequence`` object and
    creates the expected event file.
    """
    otu = _repo.create_otu("TMV", "Tobacco mosaic virus", [], 12242)
    isolate = _repo.create_isolate(otu.id, "A", "isolate")

    sequence = _repo.create_sequence(
        otu.id,
        isolate.id,
        "TMVABC.1",
        "TMV",
        "RNA",
        "ACGT",
    )

    assert sequence == EventSourcedRepoSequence(
        id=sequence.id,
        accession="TMVABC.1",
        definition="TMV",
        segment="RNA",
        sequence="ACGT",
    )

    with open(_repo.path.joinpath("src", "00000004.json")) as f:
        event = orjson.loads(f.read())

    del event["timestamp"]

    assert event == {
        "data": {
            "id": str(sequence.id),
            "accession": "TMVABC.1",
            "definition": "TMV",
            "segment": "RNA",
            "sequence": "ACGT",
        },
        "id": 4,
        "query": {
            "otu_id": str(otu.id),
            "isolate_id": str(isolate.id),
            "sequence_id": str(sequence.id),
        },
        "type": "CreateSequence",
    }

    assert _repo.last_id == 4


def test_get_otu(_repo: EventSourcedRepo):
    """Test that getting an OTU returns the expected ``RepoOTU`` object including two
    isolates with one sequence each.
    """
    otu = _repo.create_otu("TMV", "Tobacco mosaic virus", [], 12242)

    isolate_a = _repo.create_isolate(otu.id, "A", "isolate")
    _repo.create_sequence(
        otu.id,
        isolate_a.id,
        "TMVABC.1",
        "TMV",
        "RNA",
        "ACGT",
    )

    isolate_b = _repo.create_isolate(otu.id, "B", "isolate")
    _repo.create_sequence(
        otu.id,
        isolate_b.id,
        "TMVABCB.1",
        "TMV",
        "RNA",
        "ACGTGGAGAGACC",
    )

    otu = _repo.get_otu(otu.id)

    assert otu == EventSourcedRepoOTU(
        id=otu.id,
        acronym="TMV",
        excluded_accessions=[],
        isolates=[
            EventSourcedRepoIsolate(
                id=isolate_a.id,
                source_name="A",
                source_type="isolate",
                sequences=[
                    EventSourcedRepoSequence(
                        id=otu.isolates[0].sequences[0].id,
                        accession="TMVABC.1",
                        definition="TMV",
                        segment="RNA",
                        sequence="ACGT",
                    ),
                ],
            ),
            EventSourcedRepoIsolate(
                id=isolate_b.id,
                source_name="B",
                source_type="isolate",
                sequences=[
                    EventSourcedRepoSequence(
                        id=otu.isolates[1].sequences[0].id,
                        accession="TMVABCB.1",
                        definition="TMV",
                        segment="RNA",
                        sequence="ACGTGGAGAGACC",
                    ),
                ],
            ),
        ],
        name="Tobacco mosaic virus",
        schema=[],
        taxid=12242,
    )

    assert _repo.last_id == 6
