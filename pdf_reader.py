"""PDF text extraction.

Keeps a single responsibility: turn PDF bytes into plain text. Deliberately
framework-agnostic (no Streamlit imports) so it can be reused or unit-tested.
"""
from __future__ import annotations

import io

import pdfplumber


def extract_text(file_bytes: bytes) -> str:
    """Extract all text from a PDF given as raw bytes.

    Returns an empty string for scanned/image-only PDFs (no embedded text
    layer). Callers should treat empty output as "needs OCR".
    """
    parts: list[str] = []
    with pdfplumber.open(io.BytesIO(file_bytes)) as pdf:
        for page in pdf.pages:
            parts.append(page.extract_text() or "")
    return "\n".join(parts)
