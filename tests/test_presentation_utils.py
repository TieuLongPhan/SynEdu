"""Checks for reusable notebook presentation helpers."""

import math

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

from synedu.Utils.rxn_vis import (
    _draw_reaction_graph_svg,
    _ensure_2d,
    _svg_has_non_finite_coordinates,
    render_code_html,
    render_html_heading,
    render_mapping_agreement,
    render_reaction_gallery,
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


def test_reaction_gallery_scales_every_svg_to_its_card() -> None:
    gallery = render_reaction_gallery(
        "CCBr>>CCO",
        ["CCBr>>CCN", "CCBr>>CCCl"],
        title="Candidates",
        description="Showing {shown} of {total}; {hits} hits.",
    )

    assert gallery.count('class="synedu-rxn-svg"') == 3
    assert gallery.count('preserveAspectRatio="xMidYMid meet"') == 3
    assert gallery.count("width:100%;height:auto;max-width:100%") == 3
    assert "overflow-x:auto" not in gallery


def test_reaction_renderer_replaces_non_finite_coordinates() -> None:
    molecule = Chem.MolFromSmiles("CCO")
    AllChem.Compute2DCoords(molecule)
    molecule.GetConformer().SetAtomPosition(0, Point3D(float("nan"), 0.0, 0.0))

    _ensure_2d(molecule)

    conformer = molecule.GetConformer()
    assert all(
        math.isfinite(value)
        for atom_index in range(molecule.GetNumAtoms())
        for value in (
            conformer.GetAtomPosition(atom_index).x,
            conformer.GetAtomPosition(atom_index).y,
            conformer.GetAtomPosition(atom_index).z,
        )
    )


def test_reaction_renderer_rejects_non_finite_svg_paths() -> None:
    assert _svg_has_non_finite_coordinates("<path d='M nan,nan L 2,3'/>")
    assert not _svg_has_non_finite_coordinates("<path d='M 1,1 L 2,3'/>")


def test_graph_fallback_produces_finite_svg() -> None:
    svg = _draw_reaction_graph_svg("[CH3:1][Br:2]>>[CH3:1][OH:2]", "Fallback")

    assert "<svg" in svg
    assert "Fallback" in svg
    assert not _svg_has_non_finite_coordinates(svg)
