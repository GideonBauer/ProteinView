"""Shared data model.

``ProteinMention`` is the common currency between the two extraction paths
(PubTator3 and gilda) and the resolution/UI layers. It carries a *grounding* —
a namespace + identifier — rather than a raw string, plus a confidence level so
the UI can be calibrated: clean and silent when confident, flagged when not.
"""
from __future__ import annotations

from dataclasses import dataclass

# Confidence levels, ordered.
HIGH = "high"
MEDIUM = "medium"
LOW = "low"

_BADGE = {HIGH: "●●●", MEDIUM: "●●○", LOW: "●○○"}


@dataclass
class ProteinMention:
    label: str          # best surface form for display (as it appeared in text)
    count: int          # number of mentions in the document
    source: str         # "PubTator3" | "Gilda"
    db: str             # grounding namespace: "ncbigene" | "hgnc" | "uniprot"
    ident: str          # identifier within that namespace
    confidence: str     # HIGH | MEDIUM | LOW
    symbol_hint: str = ""      # approved symbol, when known (gilda HGNC groundings)
    score: float | None = None  # gilda match score, when applicable

    def badge(self) -> str:
        return _BADGE.get(self.confidence, "")
