"""Protein Explorer — Streamlit app.

Pipeline:
    upload PDF -> extract text + DOI
    -> if the paper is in PMC: read PubTator3's normalized gene annotations
       else: recognize + ground mentions locally with gilda
    -> the reader picks a protein (each carries a confidence badge)
    -> resolve its grounding to a UniProt entry (deterministic)
    -> fetch the AlphaFold structure and render it in 3D.

The key idea: we anchor to grounded identifiers, not fuzzy name search, so the
one protein shown is correct by construction (or clearly flagged when it isn't).

Run with:  streamlit run app.py
"""
from __future__ import annotations

import models
import streamlit as st

import grounding
import paper
import pubtator
import resolve
import structure
import uniprot
from pdf_reader import extract_text

_REGION_COLORS = [
    "red", "orange", "green", "magenta", "cyan", "yellow", "purple", "pink",
]

_CONF_HELP = {
    models.HIGH: "High confidence — grounded to a specific gene/protein identifier.",
    models.MEDIUM: "Medium confidence — grounded, but the name is somewhat ambiguous.",
    models.LOW: "Low confidence — the name is ambiguous; verify this is the intended protein.",
}


# --------------------------------------------------------------------------- #
# Cached analysis (runs once per uploaded file)
# --------------------------------------------------------------------------- #
@st.cache_data(show_spinner=False)
def _analyze(pdf_bytes: bytes) -> dict:
    text = extract_text(pdf_bytes)
    if not text.strip():
        return {"text_empty": True}

    doi = paper.extract_doi(text)
    ids = {}
    if doi:
        try:
            ids = paper.doi_to_ids(doi)
        except Exception:
            ids = {}

    mentions, source = [], None
    if ids.get("pmcid") or ids.get("pmid"):
        try:
            gm = pubtator.get_gene_mentions(pmid=ids.get("pmid"), pmcid=ids.get("pmcid"))
            if gm:
                mentions, source = resolve.from_pubtator(gm), "PubTator3"
        except Exception:
            pass

    if not mentions:  # not in PMC, or PubTator unavailable -> gilda fallback
        gm = grounding.extract_gene_mentions(text)
        mentions, source = resolve.from_gilda(gm), "Gilda"

    return {"doi": doi, "ids": ids, "mentions": mentions, "source": source}


@st.cache_data(show_spinner=False)
def _cached_resolve(db: str, ident: str, symbol_hint: str):
    return resolve.resolve(
        models.ProteinMention("", 0, "", db, ident, models.HIGH, symbol_hint)
    )


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
    "Upload a paper, and explore the 3D structures of the proteins it discusses."
)

uploaded = st.file_uploader("Upload a PDF", type=["pdf"])
if uploaded is None:
    st.info("Upload a PDF to get started.")
    st.stop()

with st.spinner("Reading the paper and identifying proteins…"):
    result = _analyze(uploaded.getvalue())

if result.get("text_empty"):
    st.warning(
        "No extractable text found — this looks like a scanned/image-only PDF, "
        "which would need OCR before proteins can be extracted."
    )
    st.stop()

mentions: list[models.ProteinMention] = result.get("mentions", [])
source = result.get("source")

# Provenance banner — tell the reader where the grounding came from.
if source == "PubTator3":
    st.success(
        "Proteins identified from **PubTator3** — NCBI's research-grade "
        "literature annotations for this paper. Each is normalized to a specific "
        "gene, so what you see is grounded, not guessed."
    )
elif source == "Gilda":
    st.info(
        "This paper isn't in PubMed Central, so proteins were recognized and "
        "grounded locally with **gilda**. Most are reliable; anything ambiguous "
        "is flagged below."
    )

if not mentions:
    st.warning("No protein mentions could be grounded in this document.")
    st.stop()

# --- Protein picker --------------------------------------------------------- #
st.subheader("Proteins in this paper")
options = {
    f"{m.badge()}  {m.label}  ·  {m.count} mention{'s' if m.count > 1 else ''}": m
    for m in mentions
}
chosen_label = st.selectbox(
    "Pick a protein to explore  (● = grounding confidence)", list(options)
)
mention = options[chosen_label]

if mention.confidence != models.HIGH:
    st.caption(f"{mention.badge()}  {_CONF_HELP[mention.confidence]}")

# --- Resolve to a UniProt entry -------------------------------------------- #
with st.spinner(f"Resolving “{mention.label}” to a UniProt entry…"):
    hit = _cached_resolve(mention.db, mention.ident, mention.symbol_hint)

if hit is None:
    st.warning(
        f"“{mention.label}” is grounded ({mention.db}:{mention.ident}) but has no "
        "matching UniProt entry, so there's no structure to show."
    )
    st.stop()

# --- Structure + details ---------------------------------------------------- #
left, right = st.columns([3, 2], gap="large")

with right:
    st.markdown(f"### {hit.protein_name}")
    st.markdown(
        f"**Accession:** `{hit.accession}`  \n"
        f"**Gene:** {hit.gene_names or '—'}  \n"
        f"**Organism:** {hit.organism}"
    )
    provenance = (
        f"NCBI Gene {mention.ident}" if source == "PubTator3"
        else f"{mention.db}:{mention.ident}"
        + (f" · score {mention.score}" if mention.score is not None else "")
    )
    st.caption(f"{mention.badge()} Grounded via {source} · {provenance}")

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

    entry = _cached_entry(hit.accession)
    st.markdown("**Links**")
    st.markdown(
        "  ·  ".join(
            f"[{name}]({url})" for name, url in uniprot.extract_links(entry).items()
        )
    )

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
            "An experimental PDB structure could be shown here as a fallback."
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
