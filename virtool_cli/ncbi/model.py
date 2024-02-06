from enum import Enum

from pydantic import BaseModel
from dataclasses import dataclass


class NCBIRank(str, Enum):
    SPECIES = "species"
    ISOLATE = "isolate"


class NCBISourceType(str, Enum):
    CLONE = "clone"
    GENOTYPE = "genotype"
    ISOLATE = "isolate"
    STRAIN = "strain"


class NCBISource(BaseModel):
    type: NCBISourceType
    name: str
    host: str
    segment: str | None
    taxid: int


class NCBIAccession(BaseModel):
    accession: str
    definition: str
    sequence: str
    comment: str


class NCBITaxonomy(BaseModel):
    id: int
    accessions: list[str]
    lineage: list[str]
    rank: NCBIRank
