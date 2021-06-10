import pytest
import filecmp
import os


from pathlib import Path
from virtool_cli.vfam_curate import get_input_paths, group_input_paths, remove_dupes
from virtool_cli.vfam_collapse import generate_clusters, all_by_all_blast
from virtool_cli.vfam_polyprotein import find_polyproteins
from virtool_cli.vfam_markov import write_abc, blast_to_mcl, mcl_to_fasta


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


@pytest.fixture()
def clustered_file(input_dir, output):
    input_paths = get_input_paths(Path(input_dir))
    no_phages = group_input_paths(input_paths)
    no_dupes = remove_dupes(no_phages, output, None, 1)

    return generate_clusters(no_dupes, None, None, 1.0)


@pytest.fixture()
def blast_results(output, clustered_file):
    return all_by_all_blast(clustered_file, None, 8)


@pytest.fixture()
def polyproteins(output, blast_results):
    return find_polyproteins(blast_results)


@pytest.fixture()
def mcl_results(output, blast_results, polyproteins):
    return blast_to_mcl(blast_results, polyproteins, None, None)


@pytest.mark.parametrize("input_dir, blast_path", [(DUPES_INPUT, BLAST_DUPES), (GENERIC_INPUT, BLAST_GENERIC),
                                                   (LARGE_INPUT, BLAST_LARGE)])
def test_abc(input_dir, blast_path, output, blast_results, polyproteins):
    """test that write_abc produces an abc files to match original vfam output"""
    result = write_abc(blast_results, polyproteins, None)
    expected = Path(str(blast_path) + '.abc')

    assert filecmp.cmp(result, expected)


@pytest.mark.parametrize("input_dir, blast_path", [(DUPES_INPUT, BLAST_DUPES), (GENERIC_INPUT, BLAST_GENERIC),
                                                   (LARGE_INPUT, BLAST_LARGE)])
def test_blast_to_mcl(input_dir, blast_path, output, blast_results, polyproteins):
    """test that blast_to_mcl produces .mcl files that match original vfam output"""
    result = blast_to_mcl(blast_results, polyproteins, None, None)
    expected = str(blast_path) + ".mcl"
    assert filecmp.cmp(result, expected, shallow=True)


@pytest.mark.parametrize("input_dir, blast_path", [(DUPES_INPUT, BLAST_DUPES), (GENERIC_INPUT, BLAST_GENERIC),
                                                   (LARGE_INPUT, BLAST_LARGE)])
def test_mci(input_dir, blast_path, output, blast_results, polyproteins):
    """test that blast_to_mcl produces .mci files that match original vfam output"""
    result = blast_to_mcl(blast_results, polyproteins, None, None).parent / "blast.mci"
    expected = str(blast_path) + ".mci"
    assert filecmp.cmp(result, expected, shallow=True)


@pytest.mark.parametrize("input_dir, cluster_dir", [(DUPES_INPUT, DUPES_CLUSTERS), (GENERIC_INPUT, GENERIC_CLUSTERS)])
def test_mcl_to_fasta(input_dir, cluster_dir, output, clustered_file, mcl_results):
    """test that mcl_to_fasta produces fasta files that match original vfam output"""
    result_files = mcl_to_fasta(mcl_results, clustered_file, None)
    result_files.sort()

    expected_files = []
    for file in os.listdir(cluster_dir):
        expected_files.append(cluster_dir / Path(file))
    expected_files.sort()

    for x in range(len(expected_files)):
        assert filecmp.cmp(expected_files[x], result_files[x])





