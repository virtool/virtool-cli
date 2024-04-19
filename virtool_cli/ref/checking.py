from typing import TYPE_CHECKING

import orjson

if TYPE_CHECKING:
    from virtool_cli.ref.repo import EventSourcedRepo


class Checker:
    def __init__(self, repo: "EventSourcedRepo"):
        self._repo = repo

        self._cache_path = self._repo.path / ".cache" / "repo"
        self._cache_path.mkdir(exist_ok=True)

    def _cache_by_last_id(self, func, key: str):
        last_id = self._repo.last_id

        try:
            with open(self._cache_path / f"{key}_{last_id}.json") as f:
                return orjson.loads(f.read())
        except FileNotFoundError:
            data = func()

            with open(self._cache_path / f"{key}_{last_id}.json", "wb") as f:
                f.write(orjson.dumps(data))

            return data

    def check_legacy_id_exists(self, legacy_id: str | None):
        if legacy_id is None:
            return

        legacy_otu_ids = self._cache_by_last_id(
            lambda: [
                otu.legacy_id
                for otu in self._repo.iter_otus()
                if otu.legacy_id is not None
            ],
            "legacy_otu_ids",
        )

        if legacy_id in legacy_otu_ids:
            raise ValueError(f"An OTU with the legacy ID '{legacy_id}' already exists")

    def check_otu_name_exists(self, name: str):
        lowered_otu_names = self._cache_by_last_id(
            lambda: [otu.name.lower() for otu in self._repo.iter_otus()],
            "lowered_otu_names",
        )

        if name.lower() in set(lowered_otu_names):
            raise ValueError(f"An OTU with the name '{name}' already exists")
