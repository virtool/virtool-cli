import pytest
import filecmp
import os


from virtool_cli.vfam_curate import get_input_paths, group_input_paths, remove_dupes
from virtool_cli.vfam_collapse import generate_clusters, all_by_all_blast
from virtool_cli.vfam_polyprotein import find_polyproteins
from virtool_cli.vfam_markov import write_abc, blast_to_mcl, mcl_to_fasta
from pathlib import Path


DUPES_INPUT = "vfam_input/Dupes"
GENERIC_INPUT = "vfam_input/Generic"
LARGE_INPUT = "vfam_input/Large"

BLAST_DUPES = "vfam_og_intermediates/Dupes/blast"
BLAST_GENERIC = "vfam_og_intermediates/Generic/blast"
BLAST_LARGE = "vfam_og_intermediates/Large/blast"

DUPES_CLUSTERS = "vfam_og_intermediates/Dupes/fasta_files"
GENERIC_CLUSTERS = "vfam_og_intermediates/Generic/fasta_files"
LARGE_CLUSTERS = "vfam_og_intermediates/Large/fasta_files"


@pytest.fixture()
def output(tmp_path):
    output_path = tmp_path / "dir"
    output_path.mkdir()
    return output_path


@pytest.mark.parametrize("input_dir, blast_path", [(DUPES_INPUT, BLAST_DUPES), (GENERIC_INPUT, BLAST_GENERIC),
                                                   (LARGE_INPUT, BLAST_LARGE)])
def test_abc(input_dir, blast_path, output):
    input_paths = get_input_paths(Path(input_dir))
    no_phages = group_input_paths(input_paths)
    no_dupes = remove_dupes(no_phages, output, 1)

    clustered_file = generate_clusters(no_dupes, None, 1.0)
    blast_file = all_by_all_blast(clustered_file, 8)
    polyproteins = find_polyproteins(blast_file)

    result = write_abc(blast_file, polyproteins)
    expected = Path(str(blast_path) + '.abc')

    assert filecmp.cmp(result, expected)


@pytest.mark.parametrize("input_dir, blast_path", [(DUPES_INPUT, BLAST_DUPES), (GENERIC_INPUT, BLAST_GENERIC),
                                                   (LARGE_INPUT, BLAST_LARGE)])
def test_blast_to_mcl(input_dir, blast_path, output):
    input_paths = get_input_paths(Path(input_dir))
    no_phages = group_input_paths(input_paths)
    no_dupes = remove_dupes(no_phages, output, 1)

    clustered_file = generate_clusters(no_dupes, None, 1.0)
    blast_file = all_by_all_blast(clustered_file, 8)
    polyproteins = find_polyproteins(blast_file)

    result = blast_to_mcl(blast_file, polyproteins, None)
    expected = str()


@pytest.mark.parametrize("input_dir, cluster_dir", [(DUPES_INPUT, DUPES_CLUSTERS), (GENERIC_INPUT, GENERIC_CLUSTERS),
                                                   (LARGE_INPUT, LARGE_CLUSTERS)])
def test_mcl_to_fasta(input_dir, cluster_dir, output):

    input_paths = get_input_paths(Path(input_dir))
    no_phages = group_input_paths(input_paths)
    no_dupes = remove_dupes(no_phages, output, 1)

    clustered_file = generate_clusters(no_dupes, None, 1.0)
    blast_file = all_by_all_blast(clustered_file, 8)
    polyproteins = find_polyproteins(blast_file)

    mcl_file = blast_to_mcl(blast_file, polyproteins, None)

    result_files = mcl_to_fasta(mcl_file, clustered_file)
    expected_files = os.listdir(cluster_dir)

    for x in range(len(result_files)):
        assert filecmp.cmp(result_files[x], expected_files[x], shallow=True)


if __name__ == "__main__":
    pytest.main([__file__, "-k", "test", "-v", "-s"])


