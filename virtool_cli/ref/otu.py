import sys
from collections import defaultdict

import structlog

from virtool_cli.ncbi.client import NCBIClient
from virtool_cli.ncbi.model import NCBIGenbank
from virtool_cli.ref.repo import EventSourcedRepo
from virtool_cli.ref.resources import EventSourcedRepoOTU
from virtool_cli.ref.utils import Molecule, IsolateName, IsolateNameType

base_logger = structlog.get_logger()


def add_otu(
    repo: EventSourcedRepo, taxid: int, ignore_cache: bool = False
) -> EventSourcedRepoOTU:
    """Fetch a Taxonomy record and add the OTU to the given repo.

    If the OTU cannot be added, terminate the command."""
    logger = base_logger.bind(taxid=taxid)

    ncbi = NCBIClient.from_repo(repo.path, ignore_cache)

    taxonomy = ncbi.fetch_taxonomy_record(taxid)
    if taxonomy is None:
        logger.fatal(f"Taxonomy ID {taxid} not found")
        sys.exit(1)

    molecule = None

    try:
        otu = repo.create_otu(
            acronym="",
            legacy_id=None,
            name=taxonomy.name,
            molecule=molecule,
            schema=[],
            taxid=taxid,
        )
    except ValueError as e:
        logger.warning(e)
        sys.exit(1)

    logger.debug("Created OTU", id=str(otu.id), name=otu.name, taxid=taxid)

    return otu


def get_otu_from_taxid(repo, taxid: int) -> EventSourcedRepoOTU:
    """
    Initialize a new OTU from a Taxonomy ID.

    If create_otu flag is True, fetches the Taxonomy record and adds a new OTU.
    """
    otu_index = repo.index_otus()

    if taxid in otu_index:
        otu = repo.get_otu(otu_index[taxid])
        return otu

    raise ValueError(f"No OTU found for Taxonomy id {taxid} in this repo.")


def create_otu_from_taxid(
    repo, taxid: int, ignore_cache: bool = False
) -> EventSourcedRepoOTU:
    """Initialize a new OTU from a Taxonomy ID."""
    logger = base_logger.bind(taxid=taxid)

    otu_index = repo.index_otus()

    if taxid in otu_index:
        otu = repo.get_otu(otu_index[taxid])

        logger.error(
            f"Taxonomy ID {taxid} has already been added to this reference.",
            otu_id=str(otu.id),
        )
        raise ValueError(
            f"Taxonomy ID {taxid} has already been added to this reference under OTU Id {str(otu.id)}."
        )

    otu = add_otu(repo, taxid, ignore_cache)
    if otu:
        return otu

    raise ValueError(
        f"Could not add Taxonomy id {taxid} to this repo due to record formatting issues."
    )


def update_otu(
    repo: EventSourcedRepo, otu: EventSourcedRepoOTU, ignore_cache: bool = False
):
    """Fetch a full list of Nucleotide accessions associated with the OTU
    and pass the list to the add method."""
    ncbi = NCBIClient.from_repo(repo.path, ignore_cache)

    linked_accessions = ncbi.link_accessions_from_taxid(otu.taxid)

    add_sequences(repo, otu, linked_accessions)


def group_genbank_records_by_isolate(records: list[NCBIGenbank]) -> dict:
    """Indexes Genbank records by isolate name"""
    isolates = defaultdict(dict)

    for record in records:
        logger = base_logger.bind(
            accession=record.accession,
            definition=record.definition,
            source_data=record.source,
        )

        if record.source.model_fields_set.intersection(
            {IsolateNameType.ISOLATE, IsolateNameType.STRAIN, IsolateNameType.CLONE},
        ):
            for source_type in IsolateNameType:
                if source_type in record.source.model_fields_set:
                    isolate_name = IsolateName(
                        **{
                            "type": IsolateNameType(source_type),
                            "value": record.source.model_dump()[source_type],
                        },
                    )

                    isolates[isolate_name.frozen][record.accession] = record

                    break

        elif record.refseq:
            logger.debug(
                "RefSeq record does not contain sufficient source data. "
                + "Edit before inclusion.",
                record=record,
            )

            isolate_name = IsolateName(
                **{
                    "type": IsolateNameType(IsolateNameType.REFSEQ),
                    "value": record.accession,
                },
            )

            isolates[isolate_name.frozen][record.accession] = record

        logger.debug(
            "Record does not contain sufficient source data for inclusion.",
        )

    return isolates


def add_sequences(repo: EventSourcedRepo, otu: EventSourcedRepoOTU, accessions: list):
    """Take a list of accessions, filter for eligible accessions and
    add new sequences to the OTU"""
    client = NCBIClient.from_repo(repo.path, False)

    otu_logger = base_logger.bind(taxid=otu.taxid, otu_id=str(otu.id), name=otu.name)
    fetch_list = list(set(accessions).difference(otu.blocked_accessions))
    if not fetch_list:
        otu_logger.info("OTU is up to date.")
        return

    otu_logger.info(f"Fetching {len(fetch_list)} accessions...", fetch_list=fetch_list)

    records = client.fetch_genbank_records(fetch_list)

    if records and not otu.molecule:
        # TODO: Upcoming UpdateMolecule event?
        molecule = get_molecule_from_records(records)
        otu_logger.debug("Retrieved new molecule data", molecule=molecule)

    record_bins = group_genbank_records_by_isolate(records)

    new_accessions = []

    for isolate_key in record_bins:
        record_bin = record_bins[isolate_key]

        isolate_id = otu.get_isolate_id_by_name(isolate_key)
        if isolate_id is None:
            otu_logger.debug(
                f"Creating isolate for {isolate_key.type}, {isolate_key.value}"
            )
            isolate = repo.create_isolate(
                otu_id=otu.id,
                legacy_id=None,
                source_name=isolate_key.value,
                source_type=isolate_key.type,
            )
            isolate_id = isolate.id

        for accession in record_bin:
            record = record_bin[accession]
            sequence = repo.create_sequence(
                otu_id=otu.id,
                isolate_id=isolate_id,
                accession=record.accession,
                definition=record.definition,
                legacy_id=None,
                segment=record.source.segment,
                sequence=record.sequence,
            )

            new_accessions.append(sequence.accession)

    if new_accessions:
        otu_logger.info(
            f"Added {len(new_accessions)} sequences to {otu.taxid}",
            new_accessions=new_accessions,
        )

    else:
        otu_logger.info(f"No new sequences added to OTU")


def get_molecule_from_records(records: list[NCBIGenbank]) -> Molecule:
    for record in records:
        if record.refseq:
            return Molecule(
                **{
                    "strandedness": record.strandedness.value,
                    "type": record.moltype.value,
                    "topology": record.topology.value,
                }
            )

    return Molecule(
        **{
            "strandedness": records[0].strandedness.value,
            "type": records[0].moltype.value,
            "topology": records[0].topology.value,
        }
    )
