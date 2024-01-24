from pathlib import Path

import pytest

from virtool_cli.utils.storage import get_otu_accessions


@pytest.mark.parametrize(
    "otu_dir_name,result",
    [
        (
            "abaca_bunchy_top_virus--c93ec9a9",
            [
                "EF546802.1",
                "EF546803.1",
                "EF546804.1",
                "EF546805.1",
                "EF546806.1",
                "EF546807.1",
                "NC_010314",
                "NC_010315",
                "NC_010316",
                "NC_010317",
                "NC_010318",
                "NC_010319",
            ],
        ),
        ("babaco_mosaic_virus--xcl20vqt", ["NC_036587"]),
        (
            "cabbage_leaf_curl_jamaica_virus--d226290f",
            ["DQ178610", "DQ178611", "DQ178613", "DQ178614"],
        ),
        ("faba_bean_necrotic_stunt_alphasatellite_1--6444acf3", ["NC_023881"]),
    ],
)
async def test_utils_get_otu_accessions(
    otu_dir_name: str, result: list[str], src_test_path: Path
):
    assert get_otu_accessions(src_test_path / otu_dir_name) == result
