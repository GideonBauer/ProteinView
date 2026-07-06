"""UniProt REST API helpers.

Resolves free-text protein names to UniProt accessions and pulls the metadata
that powers the app's current and planned features:

* ``search_uniprot``     — name -> candidate entries (used for disambiguation)
* ``get_entry``          — full UniProtKB JSON for one accession (cached)
* ``extract_features``   — functional regions (domains, sites) -> residue ranges
* ``extract_links``      — useful external links for an entry

Docs: https://www.uniprot.org/help/api_queries
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import requests

SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
ENTRY_URL = "https://rest.uniprot.org/uniprotkb/{acc}"

# Feature types worth surfacing/highlighting. Full list:
# https://www.uniprot.org/help/sequence_annotation
DEFAULT_FEATURE_TYPES = {
    "Domain", "Region", "Active site", "Binding site", "Motif", "Repeat",
}


@dataclass
class ProteinHit:
    accession: str
    entry_name: str
    protein_name: str
    gene_names: str
    organism: str
    reviewed: bool
    taxon_id: int | None = None

    def label(self) -> str:
        """Human-readable one-liner for a selectbox."""
        star = "★" if self.reviewed else "○"
        gene = f" ({self.gene_names})" if self.gene_names else ""
        return f"{star} {self.accession} — {self.protein_name}{gene} · {self.organism}"


# --------------------------------------------------------------------------- #
# Search + resolution
# --------------------------------------------------------------------------- #
def resolve_protein(
    name: str, prefer_taxon_id: int | None = 9606, max_alternatives: int = 8
) -> tuple[ProteinHit | None, list[ProteinHit]]:
    """Resolve a free-text protein name to a *single* best UniProt entry, plus a
    ranked list of alternatives (orthologs / name collisions) for optional
    override.

    The ranking is the important part: a plain free-text search for e.g. "TIA1"
    also returns entries that merely *mention* TIA1 in their annotation (like
    FASTK), and UniProt's default relevance can float those above the real one.
    We therefore pull candidates via an exact gene-symbol match *and* free text,
    then score so that an exact gene match in the preferred organism wins.

    ``prefer_taxon_id`` defaults to human (9606); set to ``None`` for no
    organism preference.
    """
    # Gather reviewed candidates from two angles, de-duplicated by accession.
    seen: dict[str, ProteinHit] = {}
    for query in (
        f'gene_exact:"{name}" AND reviewed:true',
        f'"{name}" AND reviewed:true',
    ):
        for hit in _search(query, size=10):
            seen.setdefault(hit.accession, hit)
    candidates = list(seen.values())

    # Last resort: nothing reviewed matched, so allow unreviewed (TrEMBL).
    if not candidates:
        candidates = _search(f'"{name}"', size=10)
    if not candidates:
        return None, []

    name_u = name.upper()

    def score(h: ProteinHit) -> tuple:
        gene_u = h.gene_names.upper()
        if gene_u == name_u:
            gene_score = 2  # exact gene-symbol match
        elif name_u in gene_u.replace(";", " ").split():
            gene_score = 1  # matches one of several gene names
        else:
            gene_score = 0  # only an annotation/name mention
        organism_score = 1 if (
            prefer_taxon_id is not None and h.taxon_id == prefer_taxon_id
        ) else 0
        reviewed_score = 1 if h.reviewed else 0
        return (gene_score, organism_score, reviewed_score)

    candidates.sort(key=score, reverse=True)
    return candidates[0], candidates[:max_alternatives]


def search_uniprot(
    name: str, size: int = 5, reviewed_only: bool = True
) -> list[ProteinHit]:
    """Plain search of UniProtKB for ``name`` (no ranking). Prefers reviewed
    entries, falling back to unreviewed if a reviewed search returns nothing.
    Kept for direct use; ``resolve_protein`` is preferred for picking one."""
    query = f'"{name}" AND reviewed:true' if reviewed_only else f'"{name}"'
    hits = _search(query, size=size)
    if not hits and reviewed_only:
        return _search(f'"{name}"', size=size)
    return hits


def _search(query: str, size: int = 10) -> list[ProteinHit]:
    params = {
        "query": query,
        "format": "json",
        "size": size,
        "fields": "accession,id,protein_name,gene_names,organism_name,reviewed",
    }
    resp = requests.get(SEARCH_URL, params=params, timeout=30)
    resp.raise_for_status()
    return [_parse_hit(r) for r in resp.json().get("results", [])]


def _parse_hit(r: dict) -> ProteinHit:
    accession = r.get("primaryAccession", "")
    entry_name = r.get("uniProtkbId", "")

    desc = r.get("proteinDescription", {})
    protein_name = (
        desc.get("recommendedName", {}).get("fullName", {}).get("value", "")
    )
    if not protein_name:
        subs = desc.get("submissionNames", [])
        if subs:
            protein_name = subs[0].get("fullName", {}).get("value", "")
    protein_name = protein_name or "(unnamed)"

    genes = r.get("genes", [])
    gene_names = genes[0].get("geneName", {}).get("value", "") if genes else ""

    organism_obj = r.get("organism", {})
    organism = organism_obj.get("scientificName", "")
    taxon_id = organism_obj.get("taxonId")
    entry_type = r.get("entryType", "").lower()
    reviewed = "reviewed" in entry_type and "unreviewed" not in entry_type

    return ProteinHit(
        accession, entry_name, protein_name, gene_names, organism, reviewed,
        taxon_id,
    )


# --------------------------------------------------------------------------- #
# Single-entry fetch + derived data
# --------------------------------------------------------------------------- #
@lru_cache(maxsize=128)
def get_entry(accession: str) -> dict:
    """Full UniProtKB JSON for one accession. Cached to avoid refetching when
    both features and links are requested for the same protein."""
    resp = requests.get(
        ENTRY_URL.format(acc=accession), params={"format": "json"}, timeout=30
    )
    resp.raise_for_status()
    return resp.json()


def extract_features(entry: dict, feature_types: set[str] | None = None) -> list[dict]:
    """Return functional features with residue ranges:
    ``[{"type", "description", "start", "end"}, ...]``."""
    wanted = feature_types or DEFAULT_FEATURE_TYPES
    out: list[dict] = []
    for f in entry.get("features", []):
        ftype = f.get("type", "")
        if ftype not in wanted:
            continue
        loc = f.get("location", {})
        start = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value")
        if start is None or end is None:
            continue
        out.append(
            {
                "type": ftype,
                "description": f.get("description", "") or ftype,
                "start": int(start),
                "end": int(end),
            }
        )
    return out


def extract_links(entry: dict) -> dict[str, str]:
    """Return a small set of useful external links for an entry."""
    acc = entry.get("primaryAccession", "")
    links = {
        "UniProt": f"https://www.uniprot.org/uniprotkb/{acc}/entry",
        "AlphaFold": f"https://alphafold.ebi.ac.uk/entry/{acc}",
    }
    # Add the first experimental PDB structure, if any are cross-referenced.
    for xref in entry.get("uniProtKBCrossReferences", []):
        if xref.get("database") == "PDB":
            pdb_id = xref.get("id", "")
            if pdb_id:
                links["RCSB PDB"] = f"https://www.rcsb.org/structure/{pdb_id}"
                break
    return links