import json
import shutil
from pathlib import Path


class NCBICache:
    """Manages caching functionality for NCBI data"""

    def __init__(self, path: Path):
        """
        :param path: A directory that will store cached data
        """
        self.path = path

        self.nuccore = self.path / "nuccore"
        self.taxonomy = self.path / "taxonomy"

        self.path.mkdir(exist_ok=True)
        self.nuccore.mkdir(exist_ok=True)
        self.taxonomy.mkdir(exist_ok=True)

    def clear(self):
        """Clear the cache."""
        shutil.rmtree(self.path)
        self.path.mkdir()

    def cache_records(self, records: list[dict], filestem: str):
        """Add a list of NCBI Nucleotide records to the cache."""
        with open(self._get_nuccore_path(filestem), "w") as f:
            json.dump(records, f)

    def load_records(self, filestem: str) -> list[dict] | None:
        """
        Retrieve a list of NCBI Nucleotide records from the cache.
        Returns None if the records are not found in the cache.
        """
        try:
            with open(self._get_nuccore_path(filestem), "r") as f:
                return json.load(f)
        except FileNotFoundError:
            return None

    def cache_taxonomy(self, taxonomy: dict, taxon_id: int):
        """Add a NCBI Taxonomy record to the cache"""
        with open(self._get_taxonomy_path(taxon_id), "w") as f:
            json.dump(taxonomy, f)

    def load_taxonomy(self, taxon_id: int) -> dict | None:
        """Load data from a cached record fetch"""
        try:
            with open(self._get_taxonomy_path(taxon_id), "r") as f:
                return json.load(f)
        except FileNotFoundError:
            return None

    def _get_nuccore_path(self, otu_id: str) -> Path:
        """Returns a standardized path for a set of cached NCBI Nucleotide records"""
        return self.nuccore / f"{otu_id}.json"

    def _get_taxonomy_path(self, taxid: int) -> Path:
        """Returns a standardized path for a cached NCBI Taxonomy record"""
        return self.taxonomy / f"{taxid}.json"
