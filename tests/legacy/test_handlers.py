from pathlib import Path

import pytest

from virtool_cli.legacy.utils import ErrorHandledResult
from virtool_cli.legacy.validate import validate_legacy_otu
from virtool_cli.ncbi.client import NCBIClient


def test_ok(legacy_otu: dict, legacy_repo_path: Path, scratch_ncbi_client: NCBIClient):
    """Test that `None` is returned for a valid legacy OTU."""
    assert (
        validate_legacy_otu(
            True,
            scratch_ncbi_client,
            legacy_otu,
        )
        is None
    )


class TestOTU:
    def test_name_too_short(
        self,
        legacy_otu: dict,
        scratch_ncbi_client: NCBIClient,
    ):
        """Test that an empty name returns the expected validation result."""
        legacy_otu["name"] = "Iota"

        otu_result = validate_legacy_otu(
            False,
            scratch_ncbi_client,
            legacy_otu,
        )

        assert otu_result.handler_results == [
            ErrorHandledResult(
                "[bold]name[/bold] must be a string with a minimum length of 5. "
                "it currently has a length of [u]4[/u].",
                False,
            ),
        ]

        assert otu_result.repaired_otu == legacy_otu

    def test_empty_schema(self, legacy_otu: dict, scratch_ncbi_client: NCBIClient):
        """Test that an empty schema raises a ValueError."""
        legacy_otu["schema"] = []

        otu_result = validate_legacy_otu(
            False,
            scratch_ncbi_client,
            legacy_otu,
        )

        assert otu_result.handler_results == [
            ErrorHandledResult(
                "[bold]schema[/bold] must be a List with a minimum length of 1. "
                "it currently has a length of [u]0[/u].",
                False,
            ),
        ]

        assert otu_result.repaired_otu == legacy_otu

    def test_missing_segment(
        self,
        legacy_otu: dict,
        scratch_ncbi_client: NCBIClient,
    ):
        """Test that validation fails if an isolate is missing a required sequence
        segment defined in the OTU schema.
        """
        legacy_otu["isolates"][0]["sequences"].pop()

        otu_result = validate_legacy_otu(
            False,
            scratch_ncbi_client,
            legacy_otu,
        )

        assert otu_result.handler_results == [
            ErrorHandledResult(
                "isolate does not contain all required schema segments",
                False,
            ),
        ]

        assert otu_result.repaired_otu == legacy_otu

    def test_invalid_sequence_segment_name(
        self,
        legacy_otu: dict,
        scratch_ncbi_client: NCBIClient,
    ):
        """Test that validation fails if an isolate has a sequence with a name that is
        not defined in the OTU schema.
        """
        legacy_otu["isolates"][0]["sequences"][0]["segment"] = "DNA Invalid"

        otu_result = validate_legacy_otu(
            False,
            scratch_ncbi_client,
            legacy_otu,
        )

        assert otu_result.handler_results == [
            ErrorHandledResult(
                "sequence contains invalid segment name",
                False,
            ),
        ]

        assert otu_result.repaired_otu == legacy_otu

    def test_duplicate_segment_name(
        self,
        legacy_otu: dict,
        scratch_ncbi_client: NCBIClient,
    ):
        """Test that validation fails when an OTU has a segment with a duplicate name in
        its schema.
        """
        legacy_otu["schema"].append(legacy_otu["schema"][0])

        otu_result = validate_legacy_otu(
            False,
            scratch_ncbi_client,
            legacy_otu,
        )

        assert otu_result.handler_results == [
            ErrorHandledResult(
                "All schema segments must have a unique name",
                False,
            ),
        ]

        assert otu_result.repaired_otu == legacy_otu

    def test_molecule_inconsistency(
        self,
        legacy_otu: dict,
        scratch_ncbi_client: NCBIClient,
    ):
        """Test that validation fails when an OTU has segments with different values
        for their ``molecule`` fields.
        """
        legacy_otu["schema"][0]["molecule"] = "DNA"

        otu_result = validate_legacy_otu(
            False,
            scratch_ncbi_client,
            legacy_otu,
        )

        assert otu_result.handler_results == [
            ErrorHandledResult(
                "All segments in a schema must have the same molecule",
                False,
            ),
        ]

        assert otu_result.repaired_otu == legacy_otu

    def test_too_many_default_isolates(
        self,
        legacy_otu: dict,
        scratch_ncbi_client: NCBIClient,
    ):
        """Test that validation fails when an OTU has more than one default isolate."""
        legacy_otu["isolates"][0]["default"] = True
        legacy_otu["isolates"][1]["default"] = True

        otu_result = validate_legacy_otu(
            False,
            scratch_ncbi_client,
            legacy_otu,
        )

        assert otu_result.handler_results == [
            ErrorHandledResult(
                "Only one isolate can be default",
                False,
            ),
        ]

        assert otu_result.repaired_otu == legacy_otu

    def test_no_default_isolate(
        self,
        legacy_otu: dict,
        scratch_ncbi_client: NCBIClient,
    ):
        """Test that validation fails when an OTU has no default isolate."""
        legacy_otu["isolates"][0]["default"] = False
        legacy_otu["isolates"][1]["default"] = False

        otu_result = validate_legacy_otu(
            False,
            scratch_ncbi_client,
            legacy_otu,
        )

        assert otu_result.handler_results == [
            ErrorHandledResult(
                "At least one isolate must be default",
                False,
            ),
        ]

        assert otu_result.repaired_otu == legacy_otu


class TestIsolate:
    @pytest.mark.parametrize("fix", [True, False], ids=["fix", "no_fix"])
    def test_source_type_invalid(
        self,
        fix: bool,
        legacy_otu: dict,
        scratch_ncbi_client: NCBIClient,
    ):
        """Test that an invalid source type raises a ValueError when auto-fix is
        disabled.
        """
        legacy_otu["isolates"][0]["source_type"] = "unknown"

        otu_result = validate_legacy_otu(
            fix,
            scratch_ncbi_client,
            legacy_otu,
        )

        assert otu_result.handler_results == [
            ErrorHandledResult(
                "[bold]source_type[/bold] must be one of [u]'clone', 'culture', "
                "'isolate', 'genbank', 'genotype', 'serotype', 'strain' or "
                "'variant'[/u], but is [u]'unknown'[/u].",
                fix,
            ),
        ]

        if fix:
            assert otu_result.repaired_otu["isolates"][0]["source_type"] == "isolate"
        else:
            assert otu_result.repaired_otu["isolates"][0]["source_type"] == "unknown"

    @pytest.mark.parametrize("fix", [True, False], ids=["fix", "no_fix"])
    def test_source_name_empty(
        self,
        fix: bool,
        legacy_otu: dict,
        scratch_ncbi_client: NCBIClient,
    ):
        """Test that an isolate with an empty source name raises a ValueError when
        autofix is disabled.
        """
        legacy_otu["isolates"][0]["source_name"] = ""

        otu_result = validate_legacy_otu(
            fix,
            scratch_ncbi_client,
            legacy_otu,
        )

        assert otu_result.handler_results == [
            ErrorHandledResult(
                "[bold]source_name[/bold] must be a string with a minimum length of 1. "
                "it currently has a length of [u]0[/u].",
                fix,
            ),
        ]

        if fix:
            assert otu_result.repaired_otu["isolates"][0]["source_name"] == "Q767"
        else:
            assert otu_result.repaired_otu == legacy_otu

    def test_source_name_empty_failed_fix(
        self,
        mocker,
        legacy_otu: dict,
    ):
        """Test that an isolate with an empty source name raises a ValueError."""
        legacy_otu["isolates"][0]["source_name"] = ""

        bad_ncbi_client = mocker.Mock(spec=NCBIClient)
        bad_ncbi_client.fetch_genbank_records.return_value = []

        otu_result = validate_legacy_otu(
            True,
            bad_ncbi_client,
            legacy_otu,
        )

        assert otu_result.handler_results == [
            ErrorHandledResult(
                "[bold]source_name[/bold] must be a string with a minimum length of 1. "
                "it currently has a length of [u]0[/u]. we could not find a source "
                "name on genbank.",
                False,
            ),
        ]

        assert otu_result.repaired_otu["isolates"][0]["source_name"] == ""

    @pytest.mark.parametrize("fix", [True, False], ids=["fix", "no_fix"])
    def test_too_few_sequences(
        self,
        fix: bool,
        legacy_otu: dict,
        scratch_ncbi_client: NCBIClient,
    ):
        """Test that an isolate with no sequences raises a ValueError."""
        legacy_otu["isolates"][0]["sequences"] = []

        otu_result = validate_legacy_otu(
            fix,
            scratch_ncbi_client,
            legacy_otu,
        )

        assert otu_result.handler_results == [
            ErrorHandledResult(
                "[bold]sequences[/bold] must be a List with a minimum length of 1. "
                "it currently has a length of [u]0[/u].",
                False,
            ),
        ]


class TestSequence:
    @pytest.mark.parametrize("fix", [True, False], ids=["fix", "no_fix"])
    def test_invalid_accession(
        self,
        fix: bool,
        legacy_otu: dict,
        scratch_ncbi_client: NCBIClient,
    ):
        """Test that validation fails when an isolate has a sequence containing invalid
        characters.
        """
        legacy_otu["isolates"][0]["sequences"][0]["accession"] = "NC_010317"

        otu_result = validate_legacy_otu(
            fix,
            scratch_ncbi_client,
            legacy_otu,
        )

        assert otu_result.handler_results == [
            ErrorHandledResult(
                "[bold]accession[/bold] value [u]NC_010317[/u] does not contain a "
                "[u].[/u] and needs a valid version suffix.",
                fix,
            ),
        ]

        if fix:
            assert (
                otu_result.repaired_otu["isolates"][0]["sequences"][0]["accession"]
                == "NC_010317.1"
            )
        else:
            assert otu_result.repaired_otu == legacy_otu

    @pytest.mark.parametrize(
        "field,input_",
        [
            ("accession", "A.2"),
            ("definition", "Too short"),
            ("sequence", "ATGCAT"),
        ],
    )
    def test_field_too_short(
        self,
        field: str,
        input_: str,
        legacy_otu: dict,
        scratch_ncbi_client: NCBIClient,
    ):
        """Test that a sequence with a definition that is too short raises a ValueError."""
        legacy_otu["isolates"][0]["sequences"][0][field] = input_

        otu_result = validate_legacy_otu(
            False,
            scratch_ncbi_client,
            legacy_otu,
        )

        assert len(otu_result.handler_results) == 1

        message = otu_result.handler_results[0].message

        assert (
            f"[bold]{field}[/bold] must be a string with a minimum length of "
            in message
        )
        assert f"it currently has a length of [u]{len(input_)}[/u]." in message

        assert otu_result.repaired_otu == legacy_otu
