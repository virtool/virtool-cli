from enum import StrEnum

from pydantic import BaseModel, Field, AliasChoices, field_validator, model_validator
from pydantic import ValidationError


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
    accession: str = Field(validation_alias="GBSeq_primary-accession")
    definition: str = Field(validation_alias="GBSeq_definition")
    sequence: str = Field(validation_alias="GBSeq_sequence")
    comment: str = Field("", validation_alias="GBSeq_comment")
    source: NCBISource = Field(validation_alias="GBSeq_feature-table")

    @field_validator("sequence", mode="after")
    @classmethod
    def to_uppercase(cl, raw: str) -> str:
        return raw.upper()

    @field_validator("source", mode="before")
    @classmethod
    def create_source(cls, raw: list) -> NCBISource:
        source_table = None
        for feature in raw:
            if feature["GBFeature_key"] == "source":
                source_table = feature
                break
        if source_table is not None:
            source_dict = {}

            for qualifier in source_table["GBFeature_quals"]:
                qual_name = qualifier["GBQualifier_name"]
                qual_value = qualifier["GBQualifier_value"]
                source_dict[qual_name] = qual_value

            return NCBISource(**source_dict)


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
