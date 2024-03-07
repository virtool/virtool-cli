import json
import shutil
from pathlib import Path

from structlog import get_logger


base_logger = get_logger()


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
        """Clear and reset the cache."""
        shutil.rmtree(self.path)
        self.path.mkdir()

        self.nuccore.mkdir()
        self.taxonomy.mkdir()

    def cache_nuccore_record(
        self, data: dict, accession: str, no_overwrite: bool = False
    ):
        """
        Add a Genbank record from NCBI Nucleotide to the cache.

        :param data: A data from a Genbank record corresponding
        :param accession: The NCBI accession of the record
        :param no_overwrite: If True, raise a FileExistsError
        """
        cached_record_path = self._get_nuccore_path(f"{accession}")
        if no_overwrite and cached_record_path.exists():
            raise FileExistsError

        with open(cached_record_path, "w") as f:
            json.dump(data, f)

            if not cached_record_path.exists():
                raise FileNotFoundError

    def load_nuccore_record(self, accession: str) -> dict | None:
        """
        Retrieve a list of NCBI Nucleotide records from the cache.
        Returns None if the records are not found in the cache.

        :param accession: The NCBI accession of the record
        :return: Deserialized Genbank data if file is found in cache, else None
        """

        try:
            with open(self._get_nuccore_path(accession), "r") as f:
                return json.load(f)

        except FileNotFoundError:
            return None

    def cache_taxonomy(self, data: dict, taxid: int, no_overwrite: bool = False):
        """
        Add a NCBI Taxonomy record to the cache

        :param data: NCBI Taxonomy record data
        :param taxid: The UID of a NCBI Taxonomy record
        :param no_overwrite: If True, raise a FileExistsError
        """
        cached_taxonomy_path = self._get_taxonomy_path(taxid)

        if no_overwrite and cached_taxonomy_path.exists():
            raise FileExistsError

        with open(cached_taxonomy_path, "w") as f:
            json.dump(data, f)

        if not cached_taxonomy_path.exists():
            raise FileNotFoundError

    def load_taxonomy(self, taxid: int) -> dict | None:
        """Load data from a cached record fetch

        :param taxid: The UID of a NCBI Taxonomy record
        :return: Deserialized Taxonomy data if file is found in cache, else None
        """
        try:
            with open(self._get_taxonomy_path(taxid), "r") as f:
                return json.load(f)
        except FileNotFoundError:
            return None

    def _get_nuccore_path(self, accession: str) -> Path:
        """Returns a standardized path for a set of cached NCBI Nucleotide records

        :param accession: The NCBI accession of a Genbank record
        :return: A properly-formatted path to a cached record
        """
        return self.nuccore / f"{accession}.json"

    def _get_taxonomy_path(self, taxid: int) -> Path:
        """Returns a standardized path for a cached NCBI Taxonomy record

        :param taxid: The UID of a NCBI Taxonomy record
        :return: A properly-formatted path to a cached record
        """
        return self.taxonomy / f"{taxid}.json"
