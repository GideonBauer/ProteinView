"""Protein-protein interactions via the STRING database.

This module is a working starting point for your planned "connections to other
proteins" feature. It is intentionally NOT wired into the UI yet — call
``get_interaction_partners`` from ``app.py`` when you're ready, and render the
result as a table or a network graph (e.g. with ``streamlit-agraph`` or an
``st.graphviz_chart``).

STRING API docs: https://string-db.org/help/api/
"""
from __future__ import annotations

import requests

_STRING_PARTNERS_URL = "https://string-db.org/api/json/interaction_partners"

# STRING species (NCBI taxonomy) IDs. Add whatever organisms you work with.
HUMAN = 9606
MOUSE = 10090


def get_interaction_partners(
    identifier: str, species: int = HUMAN, limit: int = 10
) -> list[dict]:
    """Return the top interaction partners for a protein.

    ``identifier`` can be a gene symbol (e.g. ``"TP53"``) or a UniProt accession.
    Each item: ``{"partner", "score"}`` where ``score`` is STRING's combined
    confidence in [0, 1].
    """
    params = {
        "identifiers": identifier,
        "species": species,
        "limit": limit,
        "caller_identity": "protein_explorer",
    }
    resp = requests.get(_STRING_PARTNERS_URL, params=params, timeout=30)
    resp.raise_for_status()

    partners: list[dict] = []
    for row in resp.json():
        partners.append(
            {
                "partner": row.get("preferredName_B", ""),
                "score": float(row.get("score", 0.0)),
            }
        )
    return partners
