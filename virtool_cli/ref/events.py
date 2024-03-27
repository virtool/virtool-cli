import datetime

from pydantic import UUID4, BaseModel, Field, computed_field


class EventQuery(BaseModel):
    """A base class for representing the query targeting an event at a specific
    resource.
    """


class RepoQuery(EventQuery):
    """Represents the query for targeting an event at a repository."""

    repository_id: UUID4


class OTUQuery(EventQuery):
    """Represents the query for targeting an event at an OTU."""

    otu_id: UUID4


class IsolateQuery(OTUQuery):
    """Represents the query for targeting an event at an isolate."""

    isolate_id: UUID4


class SequenceQuery(IsolateQuery):
    """Represents the query for targeting an event at a sequence."""

    sequence_id: UUID4


class EventData(BaseModel):
    """Represents the data for an event.

    Different event data classes are used
    """


class Event(BaseModel):
    id: int
    data: EventData
    query: EventQuery
    timestamp: datetime.datetime

    @computed_field
    @property
    def type(self) -> str:
        return self.__class__.__name__


class CreateRepoData(EventData):
    """Represents the data for the creation of a new repository."""

    id: UUID4
    data_type: str
    name: str
    organism: str


class CreateRepo(Event):
    """Represents the creation of a new repository."""

    data: CreateRepoData
    query: RepoQuery


class CreateOTUData(EventData):
    """Represents the data for the creation of a new OTU."""

    id: UUID4
    abbreviation: str
    excluded_accessions: list[str]
    rep_isolate: UUID4 | None
    name: str
    otu_schema: list = Field(alias="schema")
    taxid: int


class CreateOTU(Event):
    """Represents the creation of a new OTU."""

    data: CreateOTUData
    query: OTUQuery


class CreateIsolateData(EventData):
    """Represents the data for the creation of a new isolate."""

    id: UUID4
    source_name: str
    source_type: str


class CreateIsolate(Event):
    """Represents the creation of a new isolate."""

    data: CreateIsolateData
    query: IsolateQuery


class CreateSequenceData(EventData):
    """Represents the data for the creation of a new sequence."""

    id: UUID4
    accession: str
    definition: str
    segment: str
    sequence: str


class CreateSequence(Event):
    """Represents the creation of a new sequence."""

    data: CreateSequenceData
    query: SequenceQuery
