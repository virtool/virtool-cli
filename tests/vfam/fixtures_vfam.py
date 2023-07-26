from pathlib import Path

import pytest
from Bio import SeqIO

from virtool_cli.vfam.collapse import generate_clusters, blast_all_by_all
from virtool_cli.vfam.curate import get_input_paths, group_input_paths
from virtool_cli.vfam.filter import filter_on_coverage, filter_on_number
from virtool_cli.vfam.markov import blast_to_mcl, write_abc, mcl_to_fasta
from virtool_cli.vfam.msa import batch_muscle_call, batch_hmm_call
from virtool_cli.vfam.polyprotein import find_polyproteins


@pytest.fixture()
def output(tmp_path):
    output_path = tmp_path / "dir"
    output_path.mkdir()
    return output_path


@pytest.fixture()
def input_paths(input_dir):
    return get_input_paths(Path(input_dir))


@pytest.fixture()
def group_records(input_paths):
    return group_input_paths(input_paths, False, 2)


@pytest.fixture()
def curated_recs(group_records, output):
    output_dir = output / Path("intermediate_files")
    if not output_dir.exists():
        output_dir.mkdir()

    output_name = "curated_records.faa"
    output_path = output_dir / Path(output_name)

    SeqIO.write(group_records, Path(output_path), "fasta")

    return output_path


@pytest.fixture()
def clustered_result(output, curated_recs):
    return generate_clusters(curated_recs, 1.0)


@pytest.fixture()
def blast_result(clustered_result):
    return blast_all_by_all(clustered_result, 8)


@pytest.fixture()
def polyproteins(output, blast_result):
    return find_polyproteins(blast_result)


@pytest.fixture()
def mcl_results(output, blast_result, polyproteins):
    return blast_to_mcl(blast_result, polyproteins)


@pytest.fixture()
def abc_file(blast_result, polyproteins):
    return write_abc(blast_result, polyproteins)


@pytest.fixture()
def fasta_files(output, mcl_results, clustered_result):
    return mcl_to_fasta(mcl_results, clustered_result)


@pytest.fixture()
def filter_on_cvg(fasta_files):
    return filter_on_coverage(fasta_files)


@pytest.fixture()
def filter_on_num(filter_on_cvg):
    return filter_on_number(filter_on_cvg, 2)


@pytest.fixture()
def filtered_msa(output, fasta_files):
    fasta_files = filter_on_coverage(fasta_files)
    fasta_files = filter_on_number(fasta_files, 2)
    return batch_muscle_call(fasta_files)


@pytest.fixture()
def filtered_hmm(filtered_msa):
    return batch_hmm_call(filtered_msa)
