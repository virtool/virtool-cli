from enum import StrEnum
from typing import Annotated

from pydantic import (
    AliasChoices,
    BaseModel,
    ConfigDict,
    Field,
    computed_field,
    field_validator,
    model_validator,
)


class NCBIDatabase(StrEnum):
    """NCBI Databases used by NCBIClient"""

    NUCCORE = "nuccore"
    TAXONOMY = "taxonomy"


class NCBIRank(StrEnum):
    """Relevant OTU rank types"""

    SPECIES = "species"
    ISOLATE = "isolate"


class NCBISourceMolType(StrEnum):
    """The in vivo molecule type of a sequence
    Based on the INSDC controlled vocabulary list for the /mol_type qualifier

    Reference:
    https://www.insdc.org/submitting-standards/controlled-vocabulary-moltype-qualifier/
    """

    GENOMIC_DNA = "genomic DNA"
    OTHER_DNA = "other DNA"
    UNASSIGNED_DNA = "unassigned DNA"

    GENOMIC_RNA = "genomic RNA"
    MRNA = "mRNA"
    TRNA = "tRNA"
    TRANSCRIBED_RNA = "transcribed RNA"
    VIRAL_CRNA = "viral cRNA"
    OTHER_RNA = "other RNA"


class NCBISource(BaseModel):
    taxid: Annotated[int, Field(validation_alias="db_xref")]
    organism: str
    mol_type: NCBISourceMolType
    isolate: str = ""
    host: str = ""
    segment: str = ""
    strain: str = ""
    clone: str = ""

    @field_validator("taxid", mode="before")
    @classmethod
    def db_xref_to_taxid(cls, raw: str) -> int:
        return int(raw.split(":")[1])


class NCBIMolType(StrEnum):
    """The in vivo molecule type of a sequence, corresponds to Genbank's moltype field"""

    DNA = "DNA"
    RNA = "RNA"
    TRNA = "tRNA"
    MRNA = "mRNA"
    CRNA = "cRNA"


class NCBIStrandedness(StrEnum):
    SINGLE = "single"
    DOUBLE = "double"


class NCBITopology(StrEnum):
    LINEAR = "linear"
    CIRCULAR = "circular"


class NCBIGenbank(BaseModel):
    accession: Annotated[str, Field(validation_alias="GBSeq_primary-accession")]
    accession_version: Annotated[str, Field(validation_alias="GBSeq_accession-version")]
    strandedness: Annotated[
        NCBIStrandedness,
        Field(validation_alias="GBSeq_strandedness"),
    ]
    moltype: Annotated[NCBIMolType, Field(validation_alias="GBSeq_moltype")]
    topology: Annotated[NCBITopology, Field(validation_alias="GBSeq_topology")]
    definition: Annotated[str, Field(validation_alias="GBSeq_definition")]
    organism: Annotated[str, Field(validation_alias="GBSeq_organism")]
    sequence: Annotated[
        str,
        Field(validation_alias="GBSeq_sequence", pattern=r"^[gatc]+$"),
    ]
    source: Annotated[NCBISource, Field(validation_alias="GBSeq_feature-table")]
    comment: Annotated[str, Field("", validation_alias="GBSeq_comment")]

    @computed_field()
    def refseq(self) -> bool:
        return self.accession.startswith("NC_")

    @field_validator("sequence", mode="after")
    @classmethod
    def to_uppercase(cls, raw: str) -> str:
        return raw.upper()

    @field_validator("source", mode="before")
    @classmethod
    def create_source(cls, raw: list) -> NCBISource:
        for feature in raw:
            if feature["GBFeature_key"] == "source":
                source_dict = {
                    qual["GBQualifier_name"]: qual.get("GBQualifier_value", "")
                    for qual in feature["GBFeature_quals"]
                }
                return NCBISource(**source_dict)

        raise ValueError("Feature table contains no ``source`` table.")

    @model_validator(mode="after")
    def check_source(self):
        if self.source.organism != self.organism:
            raise ValueError("Non-matching organism fields on record and source")

        return self


class NCBILineage(BaseModel):
    model_config = ConfigDict(populate_by_name=True)

    id: Annotated[int, Field(validation_alias="TaxId")]
    name: Annotated[str, Field(validation_alias="ScientificName")]
    rank: Annotated[str, Field(validation_alias="Rank")]


class NCBITaxonomyOtherNames(BaseModel):
    acronym: Annotated[list[str], Field([], validation_alias="Acronym")]
    genbank_acronym: Annotated[list[str], Field([], validation_alias="GenbankAcronym")]
    equivalent_name: Annotated[list[str], Field([], validation_alias="EquivalentName")]
    synonym: Annotated[list[str], Field([], validation_alias="Synonym")]
    includes: Annotated[list[str], Field([], validation_alias="Includes")]


class NCBITaxonomy(BaseModel):
    model_config = ConfigDict(populate_by_name=True)

    id: Annotated[int, Field(validation_alias="TaxId")]
    name: Annotated[str, Field(validation_alias="ScientificName")]
    other_names: Annotated[
        NCBITaxonomyOtherNames,
        Field(NCBITaxonomyOtherNames(), validation_alias="OtherNames"),
    ]
    lineage: Annotated[list[NCBILineage], Field(validation_alias="LineageEx")]
    rank: Annotated[NCBIRank, Field(validation_alias=AliasChoices("rank", "Rank"))]

    @computed_field
    def species(self) -> NCBILineage:
        if self.rank is NCBIRank.SPECIES:
            return NCBILineage(id=self.id, name=self.name, rank=self.rank)

        for item in self.lineage:
            if item.rank == "species":
                return item

        raise ValueError("No species level taxon found in lineage")
