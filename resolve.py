"""Glue between grounding and UniProt.

* ``from_pubtator`` / ``from_gilda`` convert each path's raw output into the
  shared ``ProteinMention`` type (with a confidence level).
* ``resolve`` turns a mention's grounding into a concrete UniProt entry,
  deterministically, since the input is an identifier rather than a fuzzy name.
"""
from __future__ import annotations

import models
import uniprot
from grounding import GroundedMention
from models import ProteinMention
from pubtator import GeneMention


# --------------------------------------------------------------------------- #
# Raw path output -> ProteinMention
# --------------------------------------------------------------------------- #
def from_pubtator(mentions: list[GeneMention]) -> list[ProteinMention]:
    # PubTator3 entities are research-grade normalized -> high confidence.
    return [
        ProteinMention(
            label=m.best_label(),
            count=m.count,
            source="PubTator3",
            db="ncbigene",
            ident=m.ncbi_gene_id,
            confidence=models.HIGH,
        )
        for m in mentions
    ]


def from_gilda(mentions: list[GroundedMention]) -> list[ProteinMention]:
    return [
        ProteinMention(
            label=m.best_label(),
            count=m.count,
            source="Gilda",
            db=m.db,                       # "hgnc" or "uniprot"
            ident=m.ident,
            confidence=m.confidence(),
            symbol_hint=m.name,            # gilda's approved symbol for HGNC
            score=round(m.top_score, 3),
        )
        for m in mentions
    ]


# --------------------------------------------------------------------------- #
# Grounding -> UniProt entry
# --------------------------------------------------------------------------- #
def resolve(mention: ProteinMention) -> uniprot.ProteinHit | None:
    """Resolve a mention's grounding to a UniProt entry."""
    db = mention.db.lower()
    if db in ("ncbigene", "ncbi_gene", "geneid", "entrez"):
        return uniprot.by_gene_id(mention.ident)
    if db == "hgnc":
        # gilda gives the approved symbol; gene_exact + human is deterministic.
        hit = uniprot.by_gene_symbol(mention.symbol_hint or mention.label, taxon_id=9606)
        return hit or uniprot.by_gene_symbol(mention.symbol_hint or mention.label, taxon_id=None)
    if db == "uniprot":
        return uniprot.by_accession(mention.ident)
    return None
