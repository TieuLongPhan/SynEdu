"""Checks for reusable notebook presentation helpers."""

from synedu.Utils.rxn_vis import (
    render_code_html,
    render_html_heading,
    render_mapping_agreement,
)
from synedu.Utils.vis import render_smiles_annotation, tokenize_smiles


def test_smiles_annotation_escapes_input_and_classifies_tokens() -> None:
    tokens = tokenize_smiles("C[C@H](O)Cl")
    assert any(token["label"] == "Bracketed atom / group" for token in tokens)
    assert "&lt;unsafe&gt;" in render_smiles_annotation("C<unsafe>")


def test_mapping_presentation_helpers_escape_and_report_status() -> None:
    assert "&lt;heading&gt;" in render_html_heading("<heading>")
    assert "&lt;code&gt;" in render_code_html("<code>")
    assert "All three mappers agree" in render_mapping_agreement(True)
    assert "Mappers disagree" in render_mapping_agreement(False)
