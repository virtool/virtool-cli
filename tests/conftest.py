import pytest

from pathlib import Path
from virtool_cli.vfam_curate import get_input_paths, group_input_paths, remove_dupes, remove_phages
from virtool_cli.vfam_collapse import generate_clusters, all_by_all_blast
from virtool_cli.vfam_filter import filter_on_coverage, filter_on_number
from virtool_cli.vfam_markov import blast_to_mcl, write_abc, mcl_to_fasta
from virtool_cli.vfam_msa import batch_muscle_call, batch_hmm_call
from virtool_cli.vfam_polyprotein import find_polyproteins


@pytest.fixture()
def output(tmp_path):
    output_path = tmp_path / "dir"
    output_path.mkdir()
    return output_path


@pytest.fixture()
def input_paths(input_dir):
    return get_input_paths(Path(input_dir))


@pytest.fixture()
def no_phages(input_paths):
    return remove_phages(input_paths)


@pytest.fixture()
def group_records(input_paths):
    return group_input_paths(input_paths)


@pytest.fixture()
def no_dupes(group_records, output):
    return remove_dupes(group_records, output, None, 1)


@pytest.fixture()
def clustered_result(output, no_dupes):
    return generate_clusters(no_dupes, None, None, 1.0)


@pytest.fixture()
def blast_result(clustered_result):
    return all_by_all_blast(clustered_result, None, 8)


@pytest.fixture()
def polyproteins(output, blast_result):
    return find_polyproteins(blast_result)


@pytest.fixture()
def mcl_results(output, blast_result, polyproteins):
    return blast_to_mcl(blast_result, polyproteins, None, None)


@pytest.fixture()
def abc_file(blast_result, polyproteins):
    return write_abc(blast_result, polyproteins, None)


@pytest.fixture()
def fasta_files(output, mcl_results, clustered_result):
    return mcl_to_fasta(mcl_results, clustered_result, None)


@pytest.fixture()
def filter_on_cvg(fasta_files):
    return filter_on_coverage(fasta_files)


@pytest.fixture()
def filter_on_num(filter_on_cvg):
    return filter_on_number(filter_on_cvg, 2)


@pytest.fixture()
def get_msa(output, fasta_files):
    fasta_files = filter_on_coverage(fasta_files)
    fasta_files = filter_on_number(fasta_files, 2)
    return batch_muscle_call(fasta_files)


@pytest.fixture()
def get_hmm(get_msa):
    return batch_hmm_call(get_msa)
