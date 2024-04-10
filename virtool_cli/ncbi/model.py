import re
from enum import StrEnum

from pydantic import BaseModel, Field, ConfigDict, AliasChoices
from pydantic import field_validator, model_validator, computed_field

GATC = re.compile("[gatc]+")


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
    taxid: int = Field(validation_alias="db_xref")
    organism: str
    mol_type: NCBISourceMolType
    isolate: str = ""
    host: str = ""
    segment: str = ""
    strain: str = ""
    clone: str = ""

    @field_validator("mol_type", mode="before")
    @classmethod
    def validate_moltype(cls, raw: str) -> NCBISourceMolType:
        return NCBISourceMolType(raw)

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
    accession: str = Field(validation_alias="GBSeq_primary-accession")
    accession_version: str = Field(validation_alias="GBSeq_accession-version")
    strandedness: NCBIStrandedness = Field(validation_alias="GBSeq_strandedness")
    moltype: NCBIMolType = Field(validation_alias="GBSeq_moltype")
    topology: NCBITopology = Field(validation_alias="GBSeq_topology")
    definition: str = Field(validation_alias="GBSeq_definition")
    organism: str = Field(validation_alias="GBSeq_organism")
    sequence: str = Field(validation_alias="GBSeq_sequence")
    source: NCBISource = Field(validation_alias="GBSeq_feature-table")
    comment: str = Field("", validation_alias="GBSeq_comment")

    @computed_field()
    def refseq(self) -> bool:
        return self.accession.startswith("NC_")

    @field_validator("moltype", mode="before")
    @classmethod
    def validate_moltype(cls, raw: str) -> NCBIMolType:
        return NCBIMolType(raw)

    @field_validator("topology", mode="before")
    @classmethod
    def validate_topology(cls, raw: str) -> NCBITopology:
        return NCBITopology(raw)

    @field_validator("sequence", mode="before")
    @classmethod
    def validate_sequence(cls, raw: str) -> str:
        if GATC.match(raw):
            return raw

        raise ValueError("Invalid sequence code")

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

    id: int = Field(validation_alias="TaxId")
    name: str = Field(validation_alias="ScientificName")
    rank: str = Field(validation_alias="Rank")


class NCBITaxonomyOtherNames(BaseModel):
    acronym: list[str] = Field([], validation_alias="Acronym")
    genbank_acronym: list[str] = Field([], validation_alias="GenbankAcronym")
    equivalent_name: list[str] = Field([], validation_alias="EquivalentName")
    synonym: list[str] = Field([], validation_alias="Synonym")
    includes: list[str] = Field([], validation_alias="Includes")


class NCBITaxonomy(BaseModel):
    model_config = ConfigDict(populate_by_name=True)

    id: int = Field(validation_alias="TaxId")
    name: str = Field(validation_alias="ScientificName")
    other_names: NCBITaxonomyOtherNames = Field(
        NCBITaxonomyOtherNames(), validation_alias="OtherNames"
    )

    lineage: list[NCBILineage] = Field(validation_alias="LineageEx")
    rank: NCBIRank = Field(validation_alias=AliasChoices("rank", "Rank"))

    @computed_field
    def species(self) -> NCBILineage:
        if self.rank is NCBIRank.SPECIES:
            return NCBILineage(id=self.id, name=self.name, rank=self.rank)

        for item in self.lineage:
            if item.rank == "species":
                return item

        raise ValueError("No species level taxon found in lineage")

    @field_validator("lineage", mode="before")
    @classmethod
    def create_lineage_list(cls, raw: list[dict]):
        return [NCBILineage(**level_data) for level_data in raw]

    @field_validator("rank", mode="before")
    @classmethod
    def validate_rank(cls, raw: str):
        return NCBIRank(raw)
