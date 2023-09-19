import json
import re
from pathlib import Path
from typing import Optional


def is_v1(src_path: Path) -> bool:
    """
    Returns True if the given reference directory is formatted under version 1 guidelines.
    Determines if virtool ref migrate should be run before using this version of virtool-cli
    
    :param src_path: Path to a src database directory
    :return: Boolean value depending on whether alphabetized bins are found.
    """
    alpha_bins = list(src_path.glob('[a-z]'))
    
    if alpha_bins:
        return True
    
    return False

def get_subdir(path: Path) -> list:
    """
    Generates a list of all subdirectories under a path

    :param path: Path to a src database directory
    :return: List of paths to all OTU in a src directory
    """
    return [subdir for subdir in path.iterdir() if subdir.is_dir()]

def get_otu_paths(src_path: Path) -> list:
    """
    Generates a list of paths to all OTUs in a src directory.

    :param src_path: Path to a src database directory
    :return: List of paths to all OTU in a src directory
    """
    return [otu for otu in src_path.iterdir() if otu.is_dir()]

def get_isolate_paths(otu_path: Path) -> list:
    """
    Generates a list of paths to all isolates in an OTU directory.

    :param otu_path: Path to an OTU directory under a reference directory
    :return: List of paths to all OTU in a src directory
    """
    return [iso_path for iso_path in otu_path.iterdir() if iso_path.is_dir()]

def get_otu_accessions(otu_path: Path) -> list:
    """
    Gets all accessions from an OTU directory and returns a list

    :param otu_path: Path to an OTU directory under a reference directory
    :return: A list of all accessions under an OTU
    """
    accessions = []
    
    for isolate_path in get_isolate_paths(otu_path):
        for sequence_path in get_sequence_paths(isolate_path):
            sequence = json.loads(sequence_path.read_text())
            accessions.append(sequence['accession'])

    return accessions

def get_sequence_paths(isolate_path: Path) -> list:
    """
    Generates a list of paths to all sequences in an isolate directory.

    :param isolate_path: Path to an isolate directory under an OTU directory
    :return: A list of paths to all sequence files in an isolate directory
    """
    sequence_ids = [
        i for i in isolate_path.glob('*.json') if i.stem != "isolate"
    ]

    return sequence_ids

def search_otu_by_id(src_path: Path, id: str) -> Optional[Path]:
    """
    Searches filenames in the database by unique id and
    returns the path if it finds a match

    :param src_path: Path to a reference database directory
    :param otu_id: ID hash of OTU
    """
    for path in src_path.glob(f'*--{id}'):
        if path.is_dir():
            return path
    
    return None

def read_otu(path: Path) -> dict:
    """
    Returns a json file in dict form

    :param path: Path to an OTU directory under a reference source
    :return: Deserialized OTU data in dict form
    """
    with open(path / "otu.json", "r") as f:
        otu = json.load(f)

    return otu

def generate_otu_dirname(name: str, id: str = '') -> str:
    """
    Takes in a human-readable string, replaces whitespace and symbols
    and adds the Virtool hash id as a suffix
    
    :param name: Human-readable, searchable name of the OTU
    :param id: ID hash of OTU 
    :return: A directory name in the form of 'converted_otu_name--taxid'
    """
    no_plus = name.replace('+', 'plus ')
    no_symbols = re.split(r'[():/-]+', no_plus)
    joined = ' '.join(no_symbols)
    no_whitespace = re.sub(r'[\s]+', "_", joined)

    dirname = no_whitespace.lower()
    dirname += '--' + id

    return dirname

def get_sequence_metadata(sequence_path: Path) -> dict:
    """
    Gets the accession length and segment name from a sequence file
    and returns it in a dict

    :param sequence_path: Path to a sequence file
    :return: A dict containing the sequence accession, sequence length and segment name if present
    """
    sequence = json.loads(sequence_path.read_text())

    sequence_metadata = {
        'accession': sequence['accession'],
        'length': len(sequence['sequence'])
    }

    segment = sequence.get('segment', None)
    if segment is not None:
        sequence_metadata['segment'] = segment

    return sequence_metadata

def get_otu_accessions_metadata(otu_path) -> dict:
    """
    Returns sequence metadata for all sequences present under an OTU

    :param otu_path: Path to an OTU directory under a reference directory
    :return: An accession-keyed dict containing all constituent sequence metadata
    """
    # get length and segment metadata from sequences
    all_metadata = {}

    for isolate_path in get_isolate_paths(otu_path):
        for sequence_path in get_sequence_paths(isolate_path):
            
            sequence_metadata = get_sequence_metadata(sequence_path)
            accession = sequence_metadata['accession']
            all_metadata[accession] = sequence_metadata
    
    return all_metadata