from enum import StrEnum

from pydantic import BaseModel


class NCBIRank(StrEnum):
    SPECIES = "species"
    ISOLATE = "isolate"


class NCBISourceType(StrEnum):
    ISOLATE = "isolate"
    STRAIN = "strain"
    CLONE = "clone"
    GENOTYPE = "genotype"


class NCBISource(BaseModel):
    taxid: int
    host: str | None = None
    segment: str | None = None
    isolate: str | None = None
    strain: str | None = None
    clone: str | None = None


class NCBINuccore(BaseModel):
    accession: str
    definition: str
    sequence: str
    source: NCBISource
    comment: str = ""


class NCBITaxonomy(BaseModel):
    id: int
    accessions: list[str]
    lineage: list[str]
    rank: NCBIRank


class NCBIDB(StrEnum):
    NUCCORE = "nuccore"
    TAXONOMY = "taxonomy"
