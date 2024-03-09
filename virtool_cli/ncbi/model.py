from enum import StrEnum

from pydantic import BaseModel, Field, AliasChoices, field_validator, model_validator


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

    # @model_validator(mode="after")
    # def check_source_type(self):
    #     for source_type in ("isolate", "strain", "clone"):
    #         if getattr(self, source_type) != "":
    #             return self
    #
    #     raise ValueError("No source type included in data")


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
    id: int = Field(validation_alias=AliasChoices("TaxId", "id"))
    name: str = Field(validation_alias=AliasChoices("ScientificName", "name"))
    rank: str = Field(validation_alias=AliasChoices("Rank", "rank"))


class NCBITaxonomy(BaseModel):
    id: int
    rank: NCBIRank
    lineage: list[NCBILineage]
    species: NCBILineage


class NCBIDB(StrEnum):
    NUCCORE = "nuccore"
    TAXONOMY = "taxonomy"
