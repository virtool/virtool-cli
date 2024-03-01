from enum import StrEnum
from typing import Optional

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
    type: NCBISourceType
    name: str
    taxid: int
    host: str = ""
    segment: str = ""


class NCBIAccession(BaseModel):
    accession: str
    definition: str
    sequence: str
    comment: Optional[str] = ""


class NCBINuccore(BaseModel):
    accession: str
    definition: str
    sequence: str
    taxid: int
    host: str = ""
    segment: str = ""
    isolate: str = ""
    strain: str = ""
    clone: str = ""
    comment: str = ""


class NCBITaxonomy(BaseModel):
    id: int
    accessions: list[str]
    lineage: list[str]
    rank: NCBIRank
