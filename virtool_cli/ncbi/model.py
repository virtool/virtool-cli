from enum import StrEnum

from pydantic import BaseModel, Field, field_validator


class NCBIRank(StrEnum):
    SPECIES = "species"
    ISOLATE = "isolate"


class NCBISourceType(StrEnum):
    ISOLATE = "isolate"
    STRAIN = "strain"
    CLONE = "clone"
    GENOTYPE = "genotype"


class NCBISource(BaseModel):
    taxid: int = Field(validation_alias="db_xref")
    host: str = ""
    segment: str = ""
    isolate: str = ""
    strain: str = ""
    clone: str = ""

    @field_validator("taxid", mode="before")
    @classmethod
    def db_xref_to_taxid(cls, raw: str) -> int:
        return int(raw.split(":")[1])


class NCBINuccore(BaseModel):
    accession: str
    definition: str
    sequence: str
    source: NCBISource
    comment: str = ""

    @field_validator("sequence", mode="after")
    @classmethod
    def sequence_to_upper(cls, raw: str) -> str:
        return raw.upper()


class NCBILineage(BaseModel):
    id: int = Field(validation_alias="TaxId")
    name: str = Field(validation_alias="ScientificName")
    rank: str = Field(validation_alias="Rank")


class NCBITaxonomy(BaseModel):
    id: int
    rank: NCBIRank
    lineage: list[NCBILineage]
    species: NCBILineage


class NCBIDB(StrEnum):
    NUCCORE = "nuccore"
    TAXONOMY = "taxonomy"
