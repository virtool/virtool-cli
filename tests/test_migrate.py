from pathlib import Path

import json
import shutil
import subprocess


def convert_to_dict(otu_list: list):
    otu_dict = {}
    for otu in otu_list:
        otu_id = otu["_id"]
        isolates = otu.pop("isolates")

        iso_dict = {}
        for isolate in isolates:
            iso_id = isolate["id"]

            sequences = isolate.pop("sequences")
            seq_dict = {}
            for sequence in sequences:
                seq_id = sequence["_id"]
                seq_dict[seq_id] = sequence

            isolate["sequences"] = seq_dict
            iso_dict[iso_id] = isolate

        otu["isolates"] = iso_dict
        otu_dict[otu_id] = otu

    return otu_dict


def test_migrate(test_files_path: Path, tmp_path):
    """
    Test that the generated catalog creates all the same filenames as
    the control catalog
    """
    migrate_path = tmp_path / "migration"

    shutil.copytree(test_files_path / "src_v1", migrate_path)

    assert len([migrate_path.glob("[a-z]")]) > 0

    reference_v1 = json.loads(open(test_files_path / "reference_v1.json").read())

    v1_otu_set = convert_to_dict(reference_v1["otus"])

    ref_v2_path = tmp_path / "reference_v2.json"

    subprocess.call(["virtool", "ref", "migrate", "-src", str(migrate_path)])

    subprocess.call(
        ["virtool", "ref", "build", "-src", str(migrate_path), "-o", str(ref_v2_path)]
    )

    reference_v2 = json.loads(ref_v2_path.read_text())

    v2_otu_set = convert_to_dict(reference_v2["otus"])

    assert v1_otu_set == v2_otu_set
