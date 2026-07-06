"""gilda-based extraction — the fallback path for papers not in PMC.

gilda does dictionary-based NER *and* grounding with machine-learned
disambiguation, so it replaces both the recognizer and the fuzzy lookup — no
regex, no scispaCy. It grounds human genes/proteins to HGNC and non-human ones
to UniProt, each with a score. We keep only gene/protein groundings and derive a
confidence level from the score and how ambiguous the match was.

gilda: https://github.com/gyorilab/gilda  (pip install gilda)
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field

import models

# gilda namespaces that correspond to a single protein we can fetch a structure
# for. (Family/complex groundings like FamPlex are intentionally excluded —
# they don't map to one structure.)
_GENE_NAMESPACES = {"hgnc", "uniprot"}


@dataclass
class GroundedMention:
    db: str
    ident: str
    name: str
    surface_forms: Counter = field(default_factory=Counter)
    count: int = 0
    top_score: float = 0.0
    ambiguous: bool = False

    def best_label(self) -> str:
        return self.surface_forms.most_common(1)[0][0] if self.surface_forms else self.name

    def confidence(self) -> str:
        if self.top_score >= 0.7 and not self.ambiguous:
            return models.HIGH
        if self.top_score >= 0.5 and not self.ambiguous:
            return models.MEDIUM
        return models.LOW


def extract_gene_mentions(text: str, min_score: float = 0.5) -> list[GroundedMention]:
    """Recognize + ground gene/protein mentions in free text."""
    import gilda  # imported lazily; only needed on the fallback path

    grounded: dict[str, GroundedMention] = {}
    for ann in gilda.annotate(text):
        matches = getattr(ann, "matches", None) or []
        if not matches:
            continue
        top = matches[0]
        db = top.term.db.lower()
        if db not in _GENE_NAMESPACES or top.score < min_score:
            continue

        key = f"{db}:{top.term.id}"
        gm = grounded.get(key)
        if gm is None:
            gm = GroundedMention(db, top.term.id, top.term.entry_name)
            grounded[key] = gm
        gm.count += 1
        gm.surface_forms[getattr(ann, "text", "")] += 1
        gm.top_score = max(gm.top_score, top.score)
        # A close-scoring runner-up means the surface form is ambiguous.
        if len(matches) > 1 and (top.score - matches[1].score) < 0.1:
            gm.ambiguous = True

    return sorted(grounded.values(), key=lambda g: g.count, reverse=True)
