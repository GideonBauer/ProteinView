"""Protein Explorer — Streamlit app.

Flow:
    upload PDF -> extract text -> find protein mentions -> resolve to a
    UniProt entry -> fetch the AlphaFold structure -> render it in 3D.

Run with:  streamlit run app.py

The heavy/network steps are wrapped in Streamlit caches so they don't re-run on
every widget interaction (Streamlit reruns the whole script each time).
"""
from __future__ import annotations

import streamlit as st

import ner
import structure
import uniprot
from pdf_reader import extract_text

# A small palette for highlighting functional regions on the structure.
_REGION_COLORS = [
    "red", "orange", "green", "magenta", "cyan", "yellow", "purple", "pink",
]


# --------------------------------------------------------------------------- #
# Cached wrappers (keep the modules themselves framework-agnostic)
# --------------------------------------------------------------------------- #
@st.cache_resource(show_spinner=False)
def _warm_ner_model():
    # Loads the scispaCy model once per session (no-op for the regex backend).
    if ner.scispacy_available():
        ner._load_scispacy()  # noqa: SLF001
    return True


@st.cache_data(show_spinner=False)
def _cached_text(data: bytes) -> str:
    return extract_text(data)


@st.cache_data(show_spinner=False)
def _cached_candidates(text: str):
    return ner.extract_protein_candidates(text)


@st.cache_data(show_spinner=False)
def _cached_resolve(name: str):
    return uniprot.resolve_protein(name)


@st.cache_data(show_spinner=False)
def _cached_entry(accession: str) -> dict:
    return uniprot.get_entry(accession)


@st.cache_data(show_spinner=False)
def _cached_pdb(accession: str):
    return structure.fetch_alphafold_pdb(accession)


# --------------------------------------------------------------------------- #
# UI
# --------------------------------------------------------------------------- #
st.set_page_config(page_title="Protein Explorer", page_icon="🧬", layout="wide")
st.title("🧬 Protein Explorer")
st.caption(
    "Upload a PDF, extract the proteins it mentions, and explore their "
    "3D structures."
)

uploaded = st.file_uploader("Upload a PDF", type=["pdf"])
if uploaded is None:
    st.info("Upload a PDF to get started.")
    st.stop()

# 1. Text ------------------------------------------------------------------- #
with st.spinner("Reading PDF…"):
    text = _cached_text(uploaded.getvalue())

if not text.strip():
    st.warning(
        "No extractable text found. This looks like a scanned/image-only PDF — "
        "it would need OCR (e.g. pytesseract) before proteins can be extracted."
    )
    st.stop()

# 2. Protein mentions ------------------------------------------------------- #
_warm_ner_model()
with st.spinner("Finding protein mentions…"):
    candidates, backend = _cached_candidates(text)

if backend == "regex":
    st.info(
        "Using the built-in regex extractor. Install scispaCy (see README) for "
        "higher-quality biomedical NER."
    )

if not candidates:
    st.warning("No protein mentions detected in this document.")
    st.stop()

st.subheader("Detected protein mentions")
options = {
    f"{name}  ·  {count} mention{'s' if count > 1 else ''}": name
    for name, count in candidates
}
chosen_label = st.selectbox("Pick a protein to explore", list(options))
chosen_name = options[chosen_label]

# 3. Resolve to a single UniProt entry ------------------------------------- #
with st.spinner(f"Looking up “{chosen_name}” in UniProt…"):
    hit, alternatives = _cached_resolve(chosen_name)

if hit is None:
    st.warning(
        f"“{chosen_name}” didn’t match a UniProt entry. Try another mention — "
        "the extractor casts a wide net and not every hit is a real protein."
    )
    st.stop()

# One protein, one structure by default. The chooser is only shown (collapsed)
# when there are genuine alternatives — e.g. the same protein in another species.
others = [h for h in alternatives if h.accession != hit.accession]
if others:
    with st.expander(f"Not {hit.gene_names or hit.protein_name}? "
                     f"Pick a different match ({len(others)} other"
                     f"{'s' if len(others) > 1 else ''})"):
        labels = {h.label(): h for h in alternatives}
        default_idx = next(
            i for i, h in enumerate(alternatives) if h.accession == hit.accession
        )
        picked = st.selectbox(
            "UniProt entry", list(labels), index=default_idx,
            label_visibility="collapsed",
        )
        hit = labels[picked]

# 4. Structure + rendering options ----------------------------------------- #
left, right = st.columns([3, 2], gap="large")

with right:
    st.markdown(f"### {hit.protein_name}")
    st.markdown(
        f"**Accession:** `{hit.accession}`  \n"
        f"**Gene:** {hit.gene_names or '—'}  \n"
        f"**Organism:** {hit.organism}"
    )

    color_mode = st.radio(
        "Color by",
        ["spectrum", "confidence"],
        format_func=lambda m: {
            "spectrum": "Rainbow (N→C)",
            "confidence": "AlphaFold confidence (pLDDT)",
        }[m],
        horizontal=True,
    )
    show_regions = st.checkbox("Highlight functional regions")

    # Links (foundation for your planned "important links" feature).
    entry = _cached_entry(hit.accession)
    st.markdown("**Links**")
    st.markdown(
        "  ·  ".join(
            f"[{name}]({url})"
            for name, url in uniprot.extract_links(entry).items()
        )
    )

# Build region highlights + a legend if requested.
highlights: list[tuple[int, int, str]] = []
features: list[dict] = []
if show_regions:
    features = uniprot.extract_features(entry)
    for i, feat in enumerate(features):
        color = _REGION_COLORS[i % len(_REGION_COLORS)]
        highlights.append((feat["start"], feat["end"], color))
        feat["color"] = color

with left:
    with st.spinner("Fetching AlphaFold structure…"):
        pdb = _cached_pdb(hit.accession)

    if pdb is None:
        st.warning(
            "AlphaFold has no predicted structure for this accession. "
            "You could fall back to an experimental PDB structure here."
        )
    else:
        html = structure.render_structure_html(
            pdb, color_mode=color_mode, highlights=highlights
        )
        st.components.v1.html(html, height=540)

if show_regions:
    if features:
        st.markdown("**Highlighted regions**")
        st.dataframe(
            [
                {
                    "Color": f["color"],
                    "Type": f["type"],
                    "Residues": f"{f['start']}–{f['end']}",
                    "Description": f["description"],
                }
                for f in features
            ],
            use_container_width=True,
            hide_index=True,
        )
    else:
        st.caption("No annotated functional regions found for this entry.")