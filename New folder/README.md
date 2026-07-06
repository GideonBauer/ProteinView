# 🧬 Protein Explorer

Upload a PDF → extract the proteins it mentions → explore their 3D structures.

## What it does

1. **Reads** a PDF's text (`pdf_reader.py`, via `pdfplumber`).
2. **Finds protein mentions** (`ner.py`) using biomedical NER, with a
   dependency-free fallback.
3. **Resolves** each mention to a **UniProt** entry (`uniprot.py`).
4. **Fetches** the **AlphaFold** predicted structure and **renders it in 3D**
   (`structure.py`, via `py3Dmol`).

It also already includes the groundwork for your planned features:
external **links**, functional-**region highlighting**, and a stub for
protein–protein **interactions** (`interactions.py`, STRING API).

## Quick start

```bash
pip install -r requirements.txt
streamlit run app.py
```

That's enough to run the app immediately — it uses the built-in regex extractor
for protein mentions. For better extraction, add scispaCy below.

## Optional: high-quality NER with scispaCy (recommended)

The regex extractor casts a wide net (the UniProt step filters false positives).
For cleaner biomedical NER, install scispaCy and its NER model. The app detects
it automatically and switches backends.

```bash
pip install "scispacy>=0.6.2" "spacy>=3.7,<3.9" "numpy<2.0"
pip install https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.4/en_ner_bionlp13cg_md-0.5.4.tar.gz
```

Notes:
- On Python 3.10+ scispaCy uses `nmslib-metabrainz`, which ships prebuilt
  wheels — no C compiler needed. (Older `nmslib` on 3.9 sometimes required one.)
- `en_ner_bionlp13cg_md` tags `GENE_OR_GENE_PRODUCT` entities. To use the
  alternative `en_ner_jnlpba_md` model (explicit `protein` label), install that
  tarball instead and update `_SCISPACY_MODEL` / `_PROTEIN_LABELS` in `ner.py`.

## Project layout

| File               | Responsibility                                             |
|--------------------|------------------------------------------------------------|
| `app.py`           | Streamlit UI + caching; wires everything together          |
| `pdf_reader.py`    | PDF bytes → text                                           |
| `ner.py`           | Text → protein-name candidates (scispaCy or regex)         |
| `uniprot.py`       | Name → UniProt entry; features & links                     |
| `structure.py`     | Accession → AlphaFold PDB; render with py3Dmol             |
| `interactions.py`  | STRING interaction partners (stub, not yet in the UI)      |

Each module avoids importing Streamlit, so they stay reusable and testable;
caching lives in `app.py`.

## Extending it (your roadmap)

- **Important links** — `uniprot.extract_links()` already returns UniProt /
  AlphaFold / RCSB links. UniProt entries carry many more cross-references under
  `entry["uniProtKBCrossReferences"]` (PDB, Pfam, InterPro, Ensembl, …) to add.
- **Highlight important regions** — wired up via the "Highlight functional
  regions" checkbox. `uniprot.extract_features()` returns residue ranges that
  `structure.render_structure_html(highlights=...)` colors on the model. Widen
  `DEFAULT_FEATURE_TYPES` to surface more annotation types.
- **Connections to other proteins** — `interactions.get_interaction_partners()`
  calls the STRING API. Call it from `app.py` and render the result as a table
  or a network graph (e.g. `streamlit-agraph` or `st.graphviz_chart`).

## Known limitations

- **Scanned PDFs** have no text layer; the app detects this and tells you OCR
  (e.g. `pytesseract`) would be needed first.
- **Ambiguity** — a name like `TP53` can match several organisms/isoforms; the
  app lets you pick the exact UniProt entry (★ marks reviewed Swiss-Prot ones).
- **Network access** — UniProt, AlphaFold, and STRING are live public APIs, so
  the app needs internet access at runtime.
