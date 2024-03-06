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
        """Clear and reset the cache."""
        shutil.rmtree(self.path)
        self.path.mkdir()

        self.nuccore.mkdir()
        self.taxonomy.mkdir()

    def cache_nuccore_records(
        self, accessions: list[dict], overwrite_enabled: bool = True
    ):
        """Add a list of NCBI Nucleotide records to the cache."""
        for record in accessions:
            self.cache_nuccore_record(
                record, record["GBSeq_primary-accession"], overwrite_enabled
            )

    def cache_nuccore_record(
        self, record: dict, accession: str, overwrite_enabled: bool = True
    ):
        """
        :param record:
        :param accession:
        """
        cached_record_path = self._get_nuccore_path(f"{accession}")
        if not overwrite_enabled and cached_record_path.exists():
            raise FileExistsError

        with open(cached_record_path, "w") as f:
            json.dump(record, f)

            if not cached_record_path.exists():
                raise FileNotFoundError

    def load_nuccore_records(self, accessions: list[str]) -> list[dict] | None:
        """
        Retrieve a list of NCBI Nucleotide records from the cache.
        Returns None if the records are not found in the cache.

        :param accessions:
        :return:
        """
        records = []
        for accession in accessions:
            record = self.load_nuccore_record(accession)
            if record is not None:
                records.append(record)

        return records

    def load_nuccore_record(self, accession: str) -> dict | None:
        """
        Retrieve a list of NCBI Nucleotide records from the cache.
        Returns None if the records are not found in the cache.

        :param accession:
        :return:
        """

        try:
            with open(self._get_nuccore_path(accession), "r") as f:
                return json.load(f)

        except FileNotFoundError:
            return None

    def cache_taxonomy(
        self, taxonomy: dict, taxon_id: int, overwrite_enabled: bool = True
    ):
        """Add a NCBI Taxonomy record to the cache"""
        cached_taxonomy_path = self._get_taxonomy_path(taxon_id)
        if overwrite_enabled and cached_taxonomy_path.exists():
            raise FileExistsError

        with open(cached_taxonomy_path, "w") as f:
            json.dump(taxonomy, f)

        if not cached_taxonomy_path.exists():
            raise FileNotFoundError

    def load_taxonomy(self, taxon_id: int) -> dict | None:
        """Load data from a cached record fetch"""
        try:
            with open(self._get_taxonomy_path(taxon_id), "r") as f:
                return json.load(f)
        except FileNotFoundError:
            return None

    def _get_nuccore_path(self, accession: str) -> Path:
        """Returns a standardized path for a set of cached NCBI Nucleotide records"""
        return self.nuccore / f"{accession}.json"

    def _get_taxonomy_path(self, taxid: int) -> Path:
        """Returns a standardized path for a cached NCBI Taxonomy record"""
        return self.taxonomy / f"{taxid}.json"
