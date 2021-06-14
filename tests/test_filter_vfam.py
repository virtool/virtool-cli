import os
import pytest
import filecmp


from Bio import SeqIO
from pathlib import Path
from virtool_cli.vfam_curate import get_input_paths, group_input_paths, remove_dupes
from virtool_cli.vfam_collapse import generate_clusters, all_by_all_blast
from virtool_cli.vfam_polyprotein import find_polyproteins
from virtool_cli.vfam_markov import blast_to_mcl, mcl_to_fasta
from virtool_cli.vfam_filter import filter_on_coverage, filter_on_number


DUPES_INPUT = Path(__file__).parent.parent / "tests" / "vfam_input" / "Dupes"
GENERIC_INPUT = Path(__file__).parent.parent / "tests" / "vfam_input" / "Generic"
LARGE_INPUT = Path(__file__).parent.parent / "tests" / "vfam_input" / "Large"

DUPES_COVERAGE = Path(__file__).parent.parent / "tests" / "vfam_og_intermediates" / "Dupes" / "filtered_on_coverage"
GENERIC_COVERAGE = Path(__file__).parent.parent / "tests" / "vfam_og_intermediates" / "Generic" / "filtered_on_coverage"
LARGE_COVERAGE = Path(__file__).parent.parent / "tests" / "vfam_og_intermediates" / "Large" / "filtered_on_coverage"

DUPES_REMOVED_CVG = Path(__file__).parent.parent / "tests" / "vfam_og_intermediates" / "Dupes" / "removed_on_coverage"
GENERIC_REMOVED_CVG = Path(__file__).parent.parent / "tests" / "vfam_og_intermediates" / "Generic" / \
                      "removed_on_coverage"
LARGE_REMOVED_CVG = Path(__file__).parent.parent / "tests" / "vfam_og_intermediates" / "Large" / "removed_on_coverage"

DUPES_FILTERED_BY_NUM = Path(__file__).parent.parent / "tests" / "vfam_og_intermediates" / "Dupes" / \
                        "filtered_by_number"
GENERIC_FILTERED_BY_NUM = Path(__file__).parent.parent / "tests" / "vfam_og_intermediates" / "Generic" / \
                          "filtered_by_number"


@pytest.fixture()
def output(tmp_path):
    output_path = tmp_path / "dir"
    output_path.mkdir()
    return output_path


@pytest.fixture()
def fasta_files(input_dir, output):
    input_paths = get_input_paths(Path(input_dir))
    no_phages = group_input_paths(input_paths)
    no_dupes = remove_dupes(no_phages, output, None, 1)

    clustered_file = generate_clusters(no_dupes, None, None, 1.0)
    blast_results = all_by_all_blast(clustered_file, None, 8)
    polyproteins = find_polyproteins(blast_results)
    mcl_file = blast_to_mcl(blast_results, polyproteins, None, None)

    return mcl_to_fasta(mcl_file, clustered_file, None)


@pytest.mark.parametrize("input_dir, filtered_files", [(DUPES_INPUT, DUPES_COVERAGE), (GENERIC_INPUT, GENERIC_COVERAGE),
                                                       (LARGE_INPUT, LARGE_COVERAGE)])
def test_on_cvg_files(fasta_files, input_dir, filtered_files):
    """Test that filter_on_coverage filters the same files as original vfam."""
    with Path(filtered_files).open('r') as handle:
        expected = [line.split()[0] for line in handle]
    expected.sort()
    result = filter_on_coverage(fasta_files)
    assert len(expected) == len(result)


@pytest.mark.parametrize("input_dir, removed_recs", [(DUPES_INPUT, DUPES_REMOVED_CVG),
                                                     (GENERIC_INPUT, GENERIC_REMOVED_CVG),
                                                     (LARGE_INPUT, LARGE_REMOVED_CVG)])
def test_on_cvg_recs(input_dir, removed_recs, fasta_files):
    """Test that filter_on_coverage removes the same records as original vfam."""
    with Path(removed_recs).open('r') as handle:
        removed_records = [line.split()[0] for line in handle]

    result = filter_on_coverage(fasta_files)
    for file in result:
        for record in SeqIO.parse(file.open('r'), "fasta"):
            assert record.id not in removed_records


@pytest.mark.parametrize("input_dir, filtered_files", [(DUPES_INPUT, DUPES_FILTERED_BY_NUM),
                                                       (GENERIC_INPUT, GENERIC_FILTERED_BY_NUM)])
def test_on_number(input_dir, filtered_files, fasta_files):
    """Test that filter_on_number filters the same files as original vfam."""
    filtered_on_coverage = filter_on_coverage(fasta_files)
    result = filter_on_number(filtered_on_coverage, 2)
    result.sort()

    filtered_names = os.listdir(filtered_files)
    expected = []
    for name in filtered_names:
        expected.append(filtered_files / Path(name))
    expected.sort()

    assert len(expected) == len(result)
    for r_file, e_file in zip(result, expected):
        assert filecmp.cmp(r_file, e_file)
