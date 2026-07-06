"""Identify the paper: DOI -> PubMed / PMC IDs.

Using regex *here* is appropriate — a DOI has a strict, well-defined syntax
(unlike a protein name). We take the first DOI in the front matter (the
article's own, not one from the reference list) and resolve it through NCBI's
ID Converter, which only returns IDs for articles in PubMed/PMC. That "in PMC?"
signal is exactly the branch we want: in PMC -> use PubTator3; otherwise -> gilda.

ID Converter: https://www.ncbi.nlm.nih.gov/pmc/tools/id-converter-api/
"""
from __future__ import annotations

import re

import requests

# Crossref's recommended DOI pattern.
_DOI_RE = re.compile(r"10\.\d{4,9}/[-._;()/:A-Za-z0-9]+", re.IGNORECASE)
_IDCONV_URL = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"


def extract_doi(text: str) -> str | None:
    """Return the article DOI, searching the front matter first."""
    # The article's own DOI appears on the first page; references (which also
    # contain DOIs) come later, so bias toward the start of the document.
    head = text[:6000]
    match = _DOI_RE.search(head) or _DOI_RE.search(text)
    if not match:
        return None
    return match.group(0).rstrip(".,;)")


def doi_to_ids(
    doi: str,
    tool: str = "protein-explorer",
    email: str = "protein-explorer@example.com",
) -> dict:
    """Resolve a DOI to ``{"pmid", "pmcid", "doi"}``. Returns ``{}`` if the
    article isn't in PubMed/PMC (the caller then falls back to gilda)."""
    params = {"ids": doi, "format": "json", "tool": tool, "email": email}
    resp = requests.get(_IDCONV_URL, params=params, timeout=30)
    resp.raise_for_status()
    records = resp.json().get("records", [])
    if not records:
        return {}
    rec = records[0]
    if rec.get("status") == "error" or rec.get("errmsg"):
        return {}
    return {
        "pmid": rec.get("pmid"),
        "pmcid": rec.get("pmcid"),
        "doi": rec.get("doi"),
    }
