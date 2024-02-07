from enum import Enum

from pydantic import BaseModel
from dataclasses import dataclass


class NCBIRank(str, Enum):
    SPECIES = "species"
    ISOLATE = "isolate"


class NCBISourceType(str, Enum):
    ISOLATE = "isolate"
    STRAIN = "strain"
    CLONE = "clone"
    GENOTYPE = "genotype"


class NCBISource(BaseModel):
    type: NCBISourceType
    name: str
    taxid: int
    host: str = ""
    segment: str = ""


class NCBIAccession(BaseModel):
    accession: str
    definition: str
    sequence: str
    comment: str = ""


class NCBITaxonomy(BaseModel):
    id: int
    accessions: list[str]
    lineage: list[str]
    rank: NCBIRank
