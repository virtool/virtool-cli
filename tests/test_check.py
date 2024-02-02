import json
import shutil
from pathlib import Path

import pytest

from virtool_cli.check.otu import check_otu
from virtool_cli.utils.reference import get_isolate_paths


class TestCheckOTU:
    @pytest.mark.parametrize(
        "otu_dirname",
        [
            "abaca_bunchy_top_virus--c93ec9a9",
            "babaco_mosaic_virus--xcl20vqt",
            "cabbage_leaf_curl_jamaica_virus--d226290f",
            "faba_bean_necrotic_stunt_alphasatellite_1--6444acf3",
        ],
    )
    def test_ok(self, log, otu_dirname: str, src_scratch_path: Path):
        """Test that the check passes when an OTU is valid."""
        check_otu(src_scratch_path / otu_dirname)
        assert log.has("OTU is valid")

    def test_missing_otu_json(self, log, src_scratch_path: Path):
        """Test that the check fails when an OTU is missing its metadata file."""
        otu_path = src_scratch_path / "abaca_bunchy_top_virus--c93ec9a9"
        (otu_path / "otu.json").unlink()

        assert check_otu(otu_path) is False

        assert log.has("Could not find otu.json")

    def test_missing_isolate_json(self, log, src_scratch_path: Path):
        """Test that the check fails when an isolate is missing its metadata file."""
        otu_path = src_scratch_path / "babaco_mosaic_virus--xcl20vqt"

        (get_isolate_paths(otu_path)[0] / "isolate.json").unlink()

        assert check_otu(otu_path) is False

        assert log.has("Could not find isolate.json")

    @pytest.mark.parametrize(
        "field",
        [
            "id",
            "default",
            "source_name",
            "source_type",
        ],
    )
    def test_missing_isolate_json_value(
        self,
        field: str,
        log,
        src_scratch_path: Path,
    ):
        """Test that the check fails when an isolate.json file is missing a critical
        value.
        """
        otu_path = (
            src_scratch_path / "faba_bean_necrotic_stunt_alphasatellite_1--6444acf3"
        )

        for isolate_path in get_isolate_paths(otu_path):
            isolate = json.loads((isolate_path / "isolate.json").read_text())
            isolate[field] = ""

            with open(isolate_path / "isolate.json", "w") as f:
                json.dump(isolate, f, indent=4)

            assert check_otu(otu_path) is False

            assert log.has(f"Invalid or missing {field} field in isolate.json")

    def test_no_isolates(self, log, src_scratch_path: Path):
        """Test that the check fails when an OTU has no isolates."""
        otu_path = src_scratch_path / "abaca_bunchy_top_virus--c93ec9a9"

        for isolate_path in get_isolate_paths(otu_path):
            shutil.rmtree(isolate_path)

        assert check_otu(otu_path) is False

        assert log.has("No accession data in OTU")

    def test_no_sequences(self, log, src_scratch_path: Path):
        """Test that the check fails when an isolate has no sequences."""
        otu_path = src_scratch_path / "cabbage_leaf_curl_jamaica_virus--d226290f"

        isolate_path = get_isolate_paths(otu_path)[0]

        for sequence_path in isolate_path.iterdir():
            sequence_path.unlink()

        assert check_otu(otu_path) is False

        assert log.has("No sequences in isolate")
