"""3D protein structure: fetch predicted models and render them.

Structures come from the AlphaFold Protein Structure Database, which has a
predicted model for essentially every UniProt entry. We resolve the file via the
AlphaFold *API* (``/api/prediction/{accession}``) rather than guessing the file
URL, because the DB is versioned (v6 as of 2025) and hardcoded
``.../model_v4.pdb`` links now 404 for many entries. The API returns the current
``pdbUrl`` regardless of version. Rendering uses py3Dmol, whose ``_make_html()``
output embeds cleanly in Streamlit via ``st.components.v1.html``.

AlphaFold DB: https://alphafold.ebi.ac.uk/  ·  API docs: https://alphafold.ebi.ac.uk/api-docs
"""
from __future__ import annotations

from functools import lru_cache

import py3Dmol
import requests

# The API is keyed on UniProt accession and returns the current file URLs, so we
# never hardcode a model version. (No API key needed for read access.)
_ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"


@lru_cache(maxsize=64)
def fetch_alphafold_pdb(accession: str) -> str | None:
    """Return PDB-format text for a UniProt accession, or ``None`` if AlphaFold
    has no model for it.

    Two hops: ask the API for the prediction metadata, then download the
    ``pdbUrl`` it reports. This is robust across AlphaFold DB version bumps.
    """
    api_resp = requests.get(_ALPHAFOLD_API.format(acc=accession), timeout=30)
    if not api_resp.ok:
        return None

    predictions = api_resp.json()
    if isinstance(predictions, dict):  # defensive: some errors return an object
        predictions = [predictions]
    if not predictions:
        return None

    entry = predictions[0]
    pdb_url = entry.get("pdbUrl")
    if not pdb_url:
        # Defensive fallback if the field is ever renamed: take any URL-valued
        # field that points at a .pdb file.
        for key, val in entry.items():
            if (
                key.lower().endswith("url")
                and isinstance(val, str)
                and val.endswith(".pdb")
            ):
                pdb_url = val
                break
    if not pdb_url:
        return None

    pdb_resp = requests.get(pdb_url, timeout=30)
    return pdb_resp.text if pdb_resp.ok else None


def render_structure_html(
    pdb_data: str,
    color_mode: str = "spectrum",
    highlights: list[tuple[int, int, str]] | None = None,
    width: int = 720,
    height: int = 520,
) -> str:
    """Render a PDB structure to standalone HTML.

    ``color_mode``:
      * ``"spectrum"``   — rainbow N-terminus -> C-terminus.
      * ``"confidence"`` — color by AlphaFold pLDDT (stored in the B-factor
        column); low confidence -> high confidence. Approximates AlphaFold's
        official orange->blue scheme.

    ``highlights``: optional ``[(start_resi, end_resi, color), ...]`` overlaid on
    top of the base coloring (used for domains / functional regions).
    """
    view = py3Dmol.view(width=width, height=height)
    view.addModel(pdb_data, "pdb")

    if color_mode == "confidence":
        view.setStyle(
            {},
            {
                "cartoon": {
                    "colorscheme": {
                        "prop": "b",
                        "gradient": "roygb",
                        "min": 50,
                        "max": 90,
                    }
                }
            },
        )
    else:
        view.setStyle({}, {"cartoon": {"color": "spectrum"}})

    for start, end, color in highlights or []:
        view.setStyle(
            {"resi": f"{start}-{end}"}, {"cartoon": {"color": color}}
        )

    view.zoomTo()
    view.setBackgroundColor("white")
    return view._make_html()