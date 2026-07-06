"""PubTator3 client — the preferred, high-trust extraction path.

For any paper in the PMC open-access subset, NCBI has already run research-grade
gene/protein recognition *and* normalization (AIONER + GNorm2) and tagged every
mention with a stable NCBI Gene ID. We just read those annotations. This is far
more trustworthy than anything computed locally, and the Gene IDs are
species-specific, so orthologs separate cleanly instead of colliding.

API: https://www.ncbi.nlm.nih.gov/research/pubtator3/api  (limit: 3 req/sec)
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field

import requests

_BASE = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api"


@dataclass
class GeneMention:
    ncbi_gene_id: str
    surface_forms: Counter = field(default_factory=Counter)
    count: int = 0

    def best_label(self) -> str:
        return self.surface_forms.most_common(1)[0][0] if self.surface_forms else self.ncbi_gene_id


def get_gene_mentions(pmid: str | None = None, pmcid: str | None = None) -> list[GeneMention]:
    """Return gene mentions (normalized to NCBI Gene IDs) for a paper, most
    frequent first. Prefers full text (PMC) when available."""
    data = _export(pmid=pmid, pmcid=pmcid)
    genes: dict[str, GeneMention] = {}
    for passage in _iter_passages(data):
        for ann in passage.get("annotations", []):
            infons = ann.get("infons", {})
            if infons.get("type") != "Gene":
                continue
            raw_id = str(infons.get("identifier") or infons.get("NCBI Gene") or "")
            text = (ann.get("text") or "").strip()
            for gid in _split_ids(raw_id):
                gm = genes.setdefault(gid, GeneMention(gid))
                gm.count += 1
                if text:
                    gm.surface_forms[text] += 1
    return sorted(genes.values(), key=lambda g: g.count, reverse=True)


# --------------------------------------------------------------------------- #
# Internals
# --------------------------------------------------------------------------- #
def _export(pmid: str | None = None, pmcid: str | None = None) -> object:
    if pmcid:
        url = f"{_BASE}/publications/pmc_export/biocjson"
        params = {"pmcids": pmcid, "full": "true"}
    elif pmid:
        url = f"{_BASE}/publications/export/biocjson"
        params = {"pmids": pmid, "full": "true"}
    else:
        raise ValueError("provide a pmid or pmcid")
    resp = requests.get(url, params=params, timeout=60)
    resp.raise_for_status()
    return resp.json()


def _iter_passages(data: object):
    """BioC-JSON may arrive as a single document, a collection with a
    ``documents`` list, or a bare list of documents. Normalize and yield
    passages."""
    docs: list = []
    if isinstance(data, list):
        docs = data
    elif isinstance(data, dict):
        if isinstance(data.get("documents"), list):
            docs = data["documents"]
        elif isinstance(data.get("PubTator3"), list):
            docs = data["PubTator3"]
        elif "passages" in data:
            docs = [data]
    for doc in docs:
        if isinstance(doc, dict):
            for passage in doc.get("passages", []):
                yield passage


def _split_ids(raw_id: str):
    """A single annotation can list several gene IDs (';' or ',' separated)."""
    for part in raw_id.replace(",", ";").split(";"):
        part = part.strip()
        if part and part.lower() != "none":
            yield part
