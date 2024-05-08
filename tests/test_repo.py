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
from virtool_cli.ref.utils import (
    DataType,
    IsolateName,
    Molecule,
    MolType,
    Strandedness,
    Topology,
)


@pytest.fixture
def initialized_repo(empty_repo: EventSourcedRepo):
    otu = empty_repo.create_otu(
        "TMV",
        None,
        "Tobacco mosaic virus",
        Molecule(Strandedness.SINGLE, MolType.RNA, Topology.LINEAR),
        [],
        12242,
    )

    isolate_a = empty_repo.create_isolate(otu.id, None, "A", "isolate")
    empty_repo.create_sequence(
        otu.id,
        isolate_a.id,
        "TMVABC",
        "TMV",
        None,
        "RNA",
        "ACGT",
    )

    yield empty_repo


def init_otu(empty_repo: EventSourcedRepo) -> EventSourcedRepoOTU:
    return empty_repo.create_otu(
        "TMV",
        "abcd1234",
        "Tobacco mosaic virus",
        Molecule(Strandedness.SINGLE, MolType.RNA, Topology.LINEAR),
        [],
        12242,
    )


def test_new(empty_repo: EventSourcedRepo, tmp_path: Path):
    """Test that creating a new ``Repo`` object returns the expected object and creates
    the expected directory structure.
    """
    assert empty_repo.path == tmp_path / "test_repo"
    assert empty_repo.last_id == 1

    assert empty_repo.meta.data_type == DataType.GENOME
    assert empty_repo.meta.name == "Generic Viruses"
    assert empty_repo.meta.organism == "virus"


class TestCreateOTU:
    def test_ok(self, empty_repo: EventSourcedRepo):
        """Test that creating an OTU returns the expected ``RepoOTU`` object and creates
        the expected event file.
        """
        otu = empty_repo.create_otu(
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
            excluded_accessions=set(),
            isolates=[],
            legacy_id="abcd1234",
            name="Tobacco mosaic virus",
            molecule=Molecule(Strandedness.SINGLE, MolType.RNA, Topology.LINEAR),
            schema=[],
            taxid=12242,
        )

        with open(empty_repo.path.joinpath("src", "00000002.json")) as f:
            event = orjson.loads(f.read())

        del event["timestamp"]

        assert event == {
            "data": {
                "id": str(otu.id),
                "acronym": "TMV",
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

        assert empty_repo.last_id == 2

    def test_duplicate_name(self, empty_repo: EventSourcedRepo):
        """Test that creating an OTU with a name that already exists raises a
        ``ValueError``.
        """
        empty_repo.create_otu(
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
            empty_repo.create_otu(
                "TMV",
                None,
                "Tobacco mosaic virus",
                Molecule(Strandedness.SINGLE, MolType.RNA, Topology.LINEAR),
                [],
                12242,
            )

    def test_duplicate_legacy_id(self, empty_repo: EventSourcedRepo):
        """Test that creating an OTU with a legacy ID that already exists raises a
        ``ValueError``.
        """
        empty_repo.create_otu(
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
            empty_repo.create_otu(
                "",
                "abcd1234",
                "Abaca bunchy top virus",
                Molecule(Strandedness.SINGLE, MolType.RNA, Topology.LINEAR),
                [],
                438782,
            )


def test_create_isolate(empty_repo: EventSourcedRepo):
    """Test that creating an isolate returns the expected ``RepoIsolate`` object and
    creates the expected event file.
    """
    otu = init_otu(empty_repo)

    isolate = empty_repo.create_isolate(otu.id, None, "A", "isolate")

    assert isinstance(isolate.id, UUID)
    assert isolate.sequences == []
    assert isolate.name.value == "A"
    assert isolate.name.type == "isolate"

    with open(empty_repo.path.joinpath("src", "00000003.json")) as f:
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

    assert empty_repo.last_id == 3


def test_create_sequence(empty_repo: EventSourcedRepo):
    """Test that creating a sequence returns the expected ``RepoSequence`` object and
    creates the expected event file.
    """
    otu = init_otu(empty_repo)

    isolate = empty_repo.create_isolate(otu.id, None, "A", "isolate")

    sequence = empty_repo.create_sequence(
        otu.id,
        isolate.id,
        "TMVABC",
        "TMV",
        None,
        "RNA",
        "ACGT",
    )

    assert sequence == EventSourcedRepoSequence(
        id=sequence.id,
        accession="TMVABC",
        definition="TMV",
        legacy_id=None,
        segment="RNA",
        sequence="ACGT",
    )

    with open(empty_repo.path.joinpath("src", "00000004.json")) as f:
        event = orjson.loads(f.read())

    del event["timestamp"]

    assert event == {
        "data": {
            "id": str(sequence.id),
            "accession": "TMVABC",
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

    assert empty_repo.last_id == 4


class TestRetrieveOTU:
    def test_get_otu(self, empty_repo: EventSourcedRepo):
        """Test that getting an OTU returns the expected ``RepoOTU`` object including two
        isolates with one sequence each.
        """
        otu = empty_repo.create_otu(
            "TMV",
            None,
            "Tobacco mosaic virus",
            Molecule(Strandedness.SINGLE, MolType.RNA, Topology.LINEAR),
            [],
            12242,
        )

        isolate_a = empty_repo.create_isolate(otu.id, None, "A", "isolate")
        empty_repo.create_sequence(
            otu.id,
            isolate_a.id,
            "TMVABC",
            "TMV",
            None,
            "RNA",
            "ACGT",
        )

        isolate_b = empty_repo.create_isolate(otu.id, None, "B", "isolate")
        empty_repo.create_sequence(
            otu.id,
            isolate_b.id,
            "TMVABCB",
            "TMV",
            None,
            "RNA",
            "ACGTGGAGAGACC",
        )

        otu = empty_repo.get_otu(otu.id, ignore_cache=True)

        assert otu == EventSourcedRepoOTU(
            id=otu.id,
            acronym="TMV",
            excluded_accessions=set(),
            legacy_id=None,
            isolates=[
                EventSourcedRepoIsolate(
                    id=isolate_a.id,
                    legacy_id=None,
                    name=IsolateName(type="isolate", value="A"),
                    _sequences_by_accession={
                        "TMVABC": EventSourcedRepoSequence(
                            id=otu.isolates[0].sequences[0].id,
                            accession="TMVABC",
                            definition="TMV",
                            legacy_id=None,
                            segment="RNA",
                            sequence="ACGT",
                        ),
                    },
                ),
                EventSourcedRepoIsolate(
                    id=isolate_b.id,
                    legacy_id=None,
                    name=IsolateName(type="isolate", value="B"),
                    _sequences_by_accession={
                        "TMVABCB": EventSourcedRepoSequence(
                            id=otu.isolates[1].sequences[0].id,
                            accession="TMVABCB",
                            definition="TMV",
                            legacy_id=None,
                            segment="RNA",
                            sequence="ACGTGGAGAGACC",
                        ),
                    },
                ),
            ],
            name="Tobacco mosaic virus",
            molecule=Molecule(Strandedness.SINGLE, MolType.RNA, Topology.LINEAR),
            schema=[],
            taxid=12242,
        )

        assert empty_repo.last_id == 6

    def test_get_accessions(self, initialized_repo: EventSourcedRepo):
        otu = list(initialized_repo.iter_otus())[0]

        assert otu.accession_set == {"TMVABC"}

        isolate_b = initialized_repo.create_isolate(otu.id, None, "B", "isolate")
        initialized_repo.create_sequence(
            otu.id,
            isolate_b.id,
            "TMVABCB",
            "TMV",
            None,
            "RNA",
            "ACGTGGAGAGACC",
        )

        otu = list(initialized_repo.iter_otus())[0]

        assert otu.accession_set == {"TMVABC", "TMVABCB"}

    def test_get_blocked_accessions(self, initialized_repo: EventSourcedRepo):
        otu_id = initialized_repo.index_otus()[12242]
        isolate_b = initialized_repo.create_isolate(otu_id, None, "B", "isolate")
        initialized_repo.create_sequence(
            otu_id,
            isolate_b.id,
            "TMVABCB",
            "TMV",
            None,
            "RNA",
            "ACGTGGAGAGACC",
        )

        initialized_repo.exclude_accession(otu_id, "GROK")
        initialized_repo.exclude_accession(otu_id, "TOK")

        otu = initialized_repo.get_otu(otu_id)

        assert otu.blocked_accession_set == {"TMVABC", "TMVABCB", "GROK", "TOK"}

    def test_get_isolate(self, initialized_repo: EventSourcedRepo):
        otu = list(initialized_repo.iter_otus())[0]

        isolate_ids = {isolate.id for isolate in otu.isolates}

        for isolate_id in isolate_ids:
            assert otu.get_isolate(isolate_id) in otu.isolates

    def test_get_isolate_id_by_name(self, initialized_repo: EventSourcedRepo):
        otu = list(initialized_repo.iter_otus())[0]

        isolate_ids = {isolate.id for isolate in otu.isolates}

        assert (
            otu.get_isolate_id_by_name(IsolateName(type="isolate", value="A"))
            in isolate_ids
        )


def test_exclude_accession(empty_repo: EventSourcedRepo):
    """Test that excluding an accession from an OTU writes the expected event file and
    returns the expected OTU objects.
    """
    otu = empty_repo.create_otu(
        "TMV",
        None,
        "Tobacco mosaic virus",
        Molecule(Strandedness.SINGLE, MolType.RNA, Topology.LINEAR),
        [],
        12242,
    )

    empty_repo.exclude_accession(otu.id, "TMVABC")

    with open(empty_repo.path.joinpath("src", "00000003.json")) as f:
        event = orjson.loads(f.read())

        del event["timestamp"]

        assert event == {
            "data": {
                "accession": "TMVABC",
            },
            "id": 3,
            "query": {
                "otu_id": str(otu.id),
            },
            "type": "ExcludeAccession",
        }

    assert empty_repo.get_otu(otu.id, ignore_cache=True).excluded_accessions == {
        "TMVABC"
    }

    empty_repo.exclude_accession(otu.id, "ABTV")

    assert empty_repo.get_otu(otu.id, ignore_cache=True).excluded_accessions == {
        "TMVABC",
        "ABTV",
    }


class TestEventIndexCache:
    def test_equivalence(self, initialized_repo: EventSourcedRepo):
        assert initialized_repo.last_id == 4

        otu = list(initialized_repo.iter_otus())[0]

        otu_from_cache = initialized_repo.get_otu(otu.id, ignore_cache=False)
        otu_from_scratch = initialized_repo.get_otu(otu.id, ignore_cache=True)
        assert otu_from_cache == otu_from_scratch

        isolate_b = initialized_repo.create_isolate(otu.id, None, "B", "isolate")
        initialized_repo.create_sequence(
            otu.id,
            isolate_b.id,
            "TMVABCB",
            "TMV",
            None,
            "RNA",
            "ACGTGGAGAGACC",
        )

        assert initialized_repo.last_id == 6

        otu_from_cache = initialized_repo.get_otu(otu.id, ignore_cache=False)
        otu_from_scratch = initialized_repo.get_otu(otu.id, ignore_cache=True)
        assert otu_from_cache == otu_from_scratch

    def test_retrieve_nonexistent_otu(self, initialized_repo: EventSourcedRepo):
        assert (
            initialized_repo.get_otu(
                UUID("48b453fb-b60a-4f53-85f1-2941cd3ac5af"), ignore_cache=False
            )
            is None
        )

    def test_get_otu_without_cached_events(self, initialized_repo: EventSourcedRepo):
        otu = list(initialized_repo.iter_otus())[0]

        initialized_repo._event_index_cache.clear_cached_otu_events(otu.id)

        assert initialized_repo.get_otu(otu.id, ignore_cache=False) is not None
