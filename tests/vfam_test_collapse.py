import filecmp
import pytest


from virtool_cli.vfam_collapse import *
from virtool_cli.vfam_curate import *
from pathlib import Path


TEST_PATH = "vfam_input"
LARGE_INPUT_PATH = "vfam_large_input"
CLUSTER_FILE_PATH = "vfam_intermediates/generic_input.clstr"
CLUSTER_FASTA_PATH = "vfam_intermediates/clustered_file.faa"
BLAST_FILE_PATH = "vfam_intermediates/generic_input_all_by_all_blast.br"


@pytest.fixture()
def output(tmp_path):
    output_path = tmp_path / "dir"
    output_path.mkdir()
    return output_path


@pytest.fixture()
def cluster_results(output):
    input_paths = get_input_paths(Path(TEST_PATH))
    no_phages = remove_phages(input_paths)
    no_dupes = remove_dupes(no_phages, output)

    return generate_clusters(no_dupes, output, None, 1.0)


@pytest.fixture()
def blast_results(output):
    input_paths = get_input_paths(Path(TEST_PATH))
    no_phages = remove_phages(input_paths)
    no_dupes = remove_dupes(no_phages, output)

    # doesn't do filter by name step
    cdhit_result = generate_clusters(no_dupes, output, None, 1.0)
    return all_by_all_blast(cdhit_result, output, 8)


def test_cluster_sequences(output, cluster_results):
    """Test that cluster output file from cluster sequences matches desired output from original vfam"""
    assert filecmp.cmp(cluster_results, Path(CLUSTER_FASTA_PATH))
    cluster_file = str(cluster_results) + ".clstr"
    assert filecmp.cmp(Path(cluster_file), Path(CLUSTER_FILE_PATH), shallow=True)


def test_polyprotein_list(output, cluster_results):
    no_polyproteins = []
    filtered_by_name = polyprotein_name_check(cluster_results, output)
    recs_in_filtered = 0
    with Path(filtered_by_name) as handle:
        for record in SeqIO.parse(filtered_by_name, "fasta"):
            recs_in_filtered += 1
            no_polyproteins.append(record.id)

    polyproteins = []
    recs_in_unfiltered = 0
    with Path(cluster_results) as handle:
        for record in SeqIO.parse(cluster_results, "fasta"):
            recs_in_unfiltered += 1
            if "polyprotein" in record.description:
                polyproteins.append(record.id)

    assert recs_in_filtered == recs_in_unfiltered - len(polyproteins)
    for record in no_polyproteins:
        assert record not in polyproteins


def test_all_by_all_blast(output, blast_results):
    """Test that cluster output file from cluster sequences matches desired output from original vfam"""
    assert filecmp.cmp(blast_results, Path(BLAST_FILE_PATH))


if __name__ == "__main__":
    pytest.main([__file__, "-k", "test", "-v", "-s"])
