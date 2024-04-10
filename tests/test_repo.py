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
from virtool_cli.ref.utils import DataType, Molecule, Strandedness, MolType, Topology
from virtool_cli.ref.model import IsolateName


@pytest.fixture()
def _repo(tmp_path: Path):
    return EventSourcedRepo.new(
        DataType.GENOME,
        "Generic Viruses",
        tmp_path / "test_repo",
        "virus",
    )


def init_otu(repo: EventSourcedRepo) -> EventSourcedRepoOTU:
    return repo.create_otu(
        "TMV",
        "abcd1234",
        "Tobacco mosaic virus",
        Molecule(Strandedness.SINGLE, MolType.RNA, Topology.LINEAR),
        [],
        12242,
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


class TestCreateOTU:
    def test_ok(self, _repo: EventSourcedRepo):
        """Test that creating an OTU returns the expected ``RepoOTU`` object and creates
        the expected event file.
        """
        otu = _repo.create_otu(
            "TMV",
            "abcd1234",
            "Tobacco mosaic virus",
            Molecule(Strandedness.SINGLE, MolType.RNA, Topology.LINEAR),
            [],
            12242,
        )

        assert otu == EventSourcedRepoOTU(
            id=otu.id,
            acronym="TMV",
            excluded_accessions=[],
            isolates=[],
            legacy_id="abcd1234",
            name="Tobacco mosaic virus",
            molecule=Molecule(Strandedness.SINGLE, MolType.RNA, Topology.LINEAR),
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
                "legacy_id": "abcd1234",
                "name": "Tobacco mosaic virus",
                "molecule": {
                    "strandedness": "single",
                    "type": "RNA",
                    "topology": "linear",
                },
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

    def test_duplicate_name(self, _repo: EventSourcedRepo):
        """Test that creating an OTU with a name that already exists raises a
        ``ValueError``.
        """
        _repo.create_otu(
            "TMV",
            None,
            "Tobacco mosaic virus",
            Molecule(Strandedness.SINGLE, MolType.RNA, Topology.LINEAR),
            [],
            12242,
        )

        with pytest.raises(
            ValueError,
            match="An OTU with the name 'Tobacco mosaic virus' already exists",
        ):
            _repo.create_otu(
                "TMV",
                None,
                "Tobacco mosaic virus",
                Molecule(Strandedness.SINGLE, MolType.RNA, Topology.LINEAR),
                [],
                12242,
            )

    def test_duplicate_legacy_id(self, _repo: EventSourcedRepo):
        """Test that creating an OTU with a legacy ID that already exists raises a
        ``ValueError``.
        """
        _repo.create_otu(
            "TMV",
            "abcd1234",
            "Tobacco mosaic virus",
            Molecule(Strandedness.SINGLE, MolType.RNA, Topology.LINEAR),
            [],
            12242,
        )

        with pytest.raises(
            ValueError,
            match="An OTU with the legacy ID 'abcd1234' already exists",
        ):
            _repo.create_otu(
                "",
                "abcd1234",
                "Abaca bunchy top virus",
                Molecule(Strandedness.SINGLE, MolType.RNA, Topology.LINEAR),
                [],
                438782,
            )


def test_create_isolate(_repo: EventSourcedRepo):
    """Test that creating an isolate returns the expected ``RepoIsolate`` object and
    creates the expected event file.
    """
    otu = init_otu(_repo)

    isolate = _repo.create_isolate(otu.id, None, "A", "isolate")

    assert isinstance(isolate.id, UUID)
    assert isolate.sequences == []
    assert isolate.name.value == "A"
    assert isolate.name.type == "isolate"

    with open(_repo.path.joinpath("src", "00000003.json")) as f:
        event = orjson.loads(f.read())

    del event["timestamp"]

    assert event == {
        "data": {
            "id": str(isolate.id),
            "legacy_id": None,
            "name": {"type": "isolate", "value": "A"},
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
    otu = init_otu(_repo)

    isolate = _repo.create_isolate(otu.id, None, "A", "isolate")

    sequence = _repo.create_sequence(
        otu.id,
        isolate.id,
        "TMVABC.1",
        "TMV",
        None,
        "RNA",
        "ACGT",
    )

    assert sequence == EventSourcedRepoSequence(
        id=sequence.id,
        accession="TMVABC.1",
        definition="TMV",
        legacy_id=None,
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
            "legacy_id": None,
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
    otu = _repo.create_otu(
        "TMV",
        None,
        "Tobacco mosaic virus",
        Molecule(Strandedness.SINGLE, MolType.RNA, Topology.LINEAR),
        [],
        12242,
    )

    isolate_a = _repo.create_isolate(otu.id, None, "A", "isolate")
    _repo.create_sequence(
        otu.id,
        isolate_a.id,
        "TMVABC.1",
        "TMV",
        None,
        "RNA",
        "ACGT",
    )

    isolate_b = _repo.create_isolate(otu.id, None, "B", "isolate")
    _repo.create_sequence(
        otu.id,
        isolate_b.id,
        "TMVABCB.1",
        "TMV",
        None,
        "RNA",
        "ACGTGGAGAGACC",
    )

    otu = _repo.get_otu(otu.id)

    assert otu == EventSourcedRepoOTU(
        id=otu.id,
        acronym="TMV",
        excluded_accessions=[],
        legacy_id=None,
        isolates=[
            EventSourcedRepoIsolate(
                id=isolate_a.id,
                legacy_id=None,
                name=IsolateName(type="isolate", value="A"),
                sequences=[
                    EventSourcedRepoSequence(
                        id=otu.isolates[0].sequences[0].id,
                        accession="TMVABC.1",
                        definition="TMV",
                        legacy_id=None,
                        segment="RNA",
                        sequence="ACGT",
                    ),
                ],
            ),
            EventSourcedRepoIsolate(
                id=isolate_b.id,
                legacy_id=None,
                name=IsolateName(type="isolate", value="B"),
                sequences=[
                    EventSourcedRepoSequence(
                        id=otu.isolates[1].sequences[0].id,
                        accession="TMVABCB.1",
                        definition="TMV",
                        legacy_id=None,
                        segment="RNA",
                        sequence="ACGTGGAGAGACC",
                    ),
                ],
            ),
        ],
        name="Tobacco mosaic virus",
        molecule=Molecule(Strandedness.SINGLE, MolType.RNA, Topology.LINEAR),
        schema=[],
        taxid=12242,
    )

    assert _repo.last_id == 6
