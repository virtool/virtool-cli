import json
import shutil
from pathlib import Path


class NCBICache:
    def __init__(self, path: Path):
        path.mkdir(exist_ok=True)
        self.path = path

        self.subpaths = {
            "nuccore": self.path / "nuccore",
            "taxonomy": self.path / "taxonomy",
        }

        for key in self.subpaths:
            self.subpaths[key].mkdir(exist_ok=True)

    def clear(self):
        """Clear the cache."""
        shutil.rmtree(self.path)
        self.path.mkdir()

    async def cache_records(self, records: list[dict], filestem: str, indent=False):
        """Add a list of NCBI records to the cache"""
        with open(self._get_record_path(filestem), "w") as f:
            json.dump(records, f, indent=2 if indent else None)

    async def decache_records(self, filestem: str) -> list[dict] | None:
        """Retrieve a list of NCBI records from the cache"""
        try:
            with open(self._get_record_path(filestem), "r") as f:
                return json.load(f)
        except FileNotFoundError:
            return None

    async def cache_taxonomy(self, taxonomy: dict, taxon_id: int, indent=False):
        with open(self._get_taxonomy_path(taxon_id), "w") as f:
            json.dump(taxonomy, f, indent=2 if indent else None)

    async def decache_taxonomy(self, taxon_id: int) -> dict | None:
        try:
            with open(self._get_taxonomy_path(taxon_id), "r") as f:
                return json.load(f)
        except FileNotFoundError:
            return None

    def _get_record_path(self, otu_id):
        return self.subpaths["nuccore"] / f"{otu_id}.json"

    def _get_taxonomy_path(self, taxid):
        return self.subpaths["taxonomy"] / f"{taxid}.json"
