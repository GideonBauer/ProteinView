"""UniProt REST API helpers.

Turns a *grounded identifier* (an NCBI Gene ID, an HGNC symbol, or a UniProt
accession) into a concrete UniProt entry. Because the input is already a
resolved identifier rather than a fuzzy name, these lookups are deterministic —
the FASTK-style false positives from free-text search can't occur here.

Docs: https://www.uniprot.org/help/query-fields
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import requests

SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
ENTRY_URL = "https://rest.uniprot.org/uniprotkb/{acc}"

DEFAULT_FEATURE_TYPES = {
    "Domain", "Region", "Active site", "Binding site", "Motif", "Repeat",
}

_FIELDS = "accession,id,protein_name,gene_names,organism_name,reviewed"


@dataclass
class ProteinHit:
    accession: str
    entry_name: str
    protein_name: str
    gene_names: str
    organism: str
    reviewed: bool
    taxon_id: int | None = None


# --------------------------------------------------------------------------- #
# Identifier -> entry (deterministic)
# --------------------------------------------------------------------------- #
def by_gene_id(ncbi_gene_id: str) -> ProteinHit | None:
    """Map an NCBI (Entrez) Gene ID to its UniProt entry. The Gene ID is
    species-specific, so this is unambiguous; prefer the reviewed entry."""
    hits = _search(f"(xref:geneid-{ncbi_gene_id}) AND reviewed:true", size=5)
    if not hits:
        hits = _search(f"(xref:geneid-{ncbi_gene_id})", size=5)
    return _best_reviewed(hits)


def by_gene_symbol(symbol: str, taxon_id: int | None = 9606) -> ProteinHit | None:
    """Map an *approved* gene symbol to its reviewed UniProt entry. Uses
    ``gene_exact`` so only entries whose gene name is exactly ``symbol`` match."""
    query = f'gene_exact:"{symbol}" AND reviewed:true'
    if taxon_id is not None:
        query += f" AND organism_id:{taxon_id}"
    hits = _search(query, size=5)
    if not hits and taxon_id is not None:  # widen to any organism
        hits = _search(f'gene_exact:"{symbol}" AND reviewed:true', size=5)
    return _best_reviewed(hits)


def by_accession(accession: str) -> ProteinHit | None:
    """Fetch a UniProt entry directly by accession."""
    hits = _search(f"accession:{accession}", size=1)
    return hits[0] if hits else None


# --------------------------------------------------------------------------- #
# Single-entry fetch + derived data (features, links)
# --------------------------------------------------------------------------- #
@lru_cache(maxsize=128)
def get_entry(accession: str) -> dict:
    """Full UniProtKB JSON for one accession. Cached so features and links for
    the same protein don't refetch."""
    resp = requests.get(
        ENTRY_URL.format(acc=accession), params={"format": "json"}, timeout=30
    )
    resp.raise_for_status()
    return resp.json()


def extract_features(entry: dict, feature_types: set[str] | None = None) -> list[dict]:
    """Functional features with residue ranges: ``[{type, description, start, end}]``."""
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
    """A few useful external links for an entry."""
    acc = entry.get("primaryAccession", "")
    links = {
        "UniProt": f"https://www.uniprot.org/uniprotkb/{acc}/entry",
        "AlphaFold": f"https://alphafold.ebi.ac.uk/entry/{acc}",
    }
    for xref in entry.get("uniProtKBCrossReferences", []):
        if xref.get("database") == "PDB":
            pdb_id = xref.get("id", "")
            if pdb_id:
                links["RCSB PDB"] = f"https://www.rcsb.org/structure/{pdb_id}"
                break
    return links


# --------------------------------------------------------------------------- #
# Internals
# --------------------------------------------------------------------------- #
def _search(query: str, size: int = 5) -> list[ProteinHit]:
    params = {"query": query, "format": "json", "size": size, "fields": _FIELDS}
    resp = requests.get(SEARCH_URL, params=params, timeout=30)
    resp.raise_for_status()
    return [_parse_hit(r) for r in resp.json().get("results", [])]


def _best_reviewed(hits: list[ProteinHit]) -> ProteinHit | None:
    if not hits:
        return None
    hits.sort(key=lambda h: (h.reviewed,), reverse=True)  # reviewed first
    return hits[0]


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
