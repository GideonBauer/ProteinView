"""Protein-mention extraction from free text.

Two backends are available:

* ``scispacy`` — biomedical NER (recommended). Uses the ``en_ner_bionlp13cg_md``
  model and keeps ``GENE_OR_GENE_PRODUCT`` entities, which covers proteins.
* ``regex`` — a dependency-free fallback that over-generates gene-symbol-shaped
  tokens (e.g. ``TP53``, ``BRCA1``, ``Cas9``). It is intentionally noisy; the
  UniProt resolution step in the app filters out tokens that don't map to a
  real protein, so the two stages together stay usable without scispaCy.

``backend="auto"`` uses scispaCy if it's installed, otherwise the regex fallback.
"""
from __future__ import annotations

import importlib.util
import re
from collections import Counter
from functools import lru_cache

# scispaCy model + the entity label that corresponds to proteins.
# Swap for "en_ner_jnlpba_md" / {"protein"} if you prefer that model's labels.
_SCISPACY_MODEL = "en_ner_bionlp13cg_md"
_PROTEIN_LABELS = {"GENE_OR_GENE_PRODUCT"}


# --------------------------------------------------------------------------- #
# Public API
# --------------------------------------------------------------------------- #
def extract_protein_candidates(
    text: str, backend: str = "auto", top_k: int = 40
) -> tuple[list[tuple[str, int]], str]:
    """Return ``(candidates, backend_used)``.

    ``candidates`` is a list of ``(name, mention_count)`` tuples, most frequent
    first, capped at ``top_k``.
    """
    if backend == "auto":
        backend = "scispacy" if scispacy_available() else "regex"

    if backend == "scispacy":
        counts = _extract_scispacy(text)
    else:
        counts = _extract_regex(text)

    return counts.most_common(top_k), backend


def scispacy_available() -> bool:
    """True if scispaCy and the NER model are importable."""
    try:
        import scispacy  # noqa: F401
        import spacy  # noqa: F401
    except Exception:
        return False
    return importlib.util.find_spec(_SCISPACY_MODEL) is not None


# --------------------------------------------------------------------------- #
# scispaCy backend
# --------------------------------------------------------------------------- #
@lru_cache(maxsize=1)
def _load_scispacy():
    import spacy  # imported lazily so the app runs even without scispaCy

    return spacy.load(_SCISPACY_MODEL)


def _extract_scispacy(text: str) -> Counter:
    nlp = _load_scispacy()
    counts: Counter = Counter()
    # spaCy has a per-doc length cap; process the document in chunks.
    for chunk in _chunks(text, 90_000):
        doc = nlp(chunk)
        for ent in doc.ents:
            if ent.label_ in _PROTEIN_LABELS:
                name = ent.text.strip()
                if _looks_valid(name):
                    counts[name] += 1
    return counts


# --------------------------------------------------------------------------- #
# Regex fallback
# --------------------------------------------------------------------------- #
_CANDIDATE_RE = re.compile(r"\b[A-Za-z][A-Za-z0-9]{1,9}\b")

# Common all-caps tokens that look like gene symbols but aren't proteins.
# Extend this as you notice noise in your own documents.
_STOPWORDS = {
    "DNA", "RNA", "MRNA", "CDNA", "PCR", "QPCR", "ELISA", "PAGE", "SDS",
    "PBS", "ATP", "GTP", "ADP", "NADH", "NADPH", "UV", "USA", "PDF", "FIG",
    "AND", "THE", "FOR", "WITH", "III", "II", "IV", "PH", "OD", "RT",
}


def _extract_regex(text: str) -> Counter:
    counts: Counter = Counter()
    for tok in _CANDIDATE_RE.findall(text):
        if tok.upper() in _STOPWORDS:
            continue
        # Keep gene-symbol-shaped tokens: all-caps (EGFR) or containing a
        # digit (TP53, Cas9). Plain lowercase/capitalized words are dropped.
        if tok.isupper() or any(c.isdigit() for c in tok):
            counts[tok] += 1
    return counts


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def _looks_valid(name: str) -> bool:
    return len(name) >= 2 and not name.isdigit()


def _chunks(text: str, size: int):
    """Yield ``text`` in <= ``size`` char slices, breaking on whitespace so
    entities aren't split across chunk boundaries."""
    start, n = 0, len(text)
    while start < n:
        end = min(start + size, n)
        if end < n:
            space = text.rfind(" ", start, end)
            if space > start:
                end = space
        yield text[start:end]
        start = end
