import json
from pathlib import Path


class RepoMaps:
    def __init__(self, path: Path):
        self.isolate_id_to_otu_id: dict[str, str] = {}
        """Maps isolate ids to parent otu ids."""

        self.otu_id_to_path: dict[str, Path] = {}
        """Maps otu ids to their respective paths."""

        self.sequence_id_to_accession: dict[str, str] = {}
        """Maps sequence ids to their respective accessions."""

        self.sequence_id_to_otu_id: dict[str, str] = {}
        """Maps sequence ids to parent otu ids."""

        self.sequence_id_to_path: dict[str, Path] = {}
        """Maps sequence ids to their respective paths."""

        self.taxid_to_otu_id: dict[int, str] = {}
        """Maps taxids to otu ids."""

        for otu_path in (p for p in (path / "src").iterdir() if p.is_dir()):
            if not otu_path.is_dir():
                continue

            name, otu_id = otu_path.name.split("--")

            self.otu_id_to_path[otu_id] = otu_path

            with open(otu_path / "otu.json") as f:
                taxid = int(json.load(f)["taxid"])

            self.taxid_to_otu_id[taxid] = otu_id

            for isolate_path in (p for p in otu_path.iterdir() if p.is_dir()):
                self.isolate_id_to_otu_id[isolate_path.stem] = otu_id

                for sequence_path in [
                    p for p in isolate_path.iterdir() if p.name != "isolate.json"
                ]:
                    sequence_id = sequence_path.stem

                    with open(sequence_path) as f:
                        self.sequence_id_to_accession[sequence_id] = json.load(f)[
                            "accession"
                        ]

                    self.sequence_id_to_otu_id[sequence_id] = otu_id
                    self.sequence_id_to_path[sequence_id] = sequence_path
