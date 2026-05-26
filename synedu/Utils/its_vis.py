from __future__ import annotations

from typing import Any, Dict, List, Optional, Set, Tuple

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx
from rdkit import Chem
from rdkit.Chem import AllChem

from .conversion import graph_to_mol

# ── CPK element palette (fill, border) ───────────────────────────────────
_ELEMENT_PALETTE: Dict[str, Tuple[str, str]] = {
    "C": ("#636363", "#3d3d3d"),
    "O": ("#E8524A", "#b83830"),
    "N": ("#5B8DD9", "#3a65b0"),
    "S": ("#E8A838", "#c07a10"),
    "Cl": ("#3DBE6C", "#1e8a46"),
    "F": ("#5BC8AF", "#2a9178"),
    "Br": ("#A0522D", "#6b3118"),
    "I": ("#8C54C8", "#5e2fa0"),
    "P": ("#E878C8", "#b84898"),
    "H": ("#C8C8C8", "#909090"),
    "Na": ("#AB5CF2", "#7b34c8"),
    "Mg": ("#8AFF00", "#58b000"),
    "Si": ("#F0C8A0", "#b88860"),
}
_DEFAULT_FILL = "#A0A0A0"
_DEFAULT_BORDER = "#606060"

# ── edge-change type colors ───────────────────────────────────────────────
_EDGE_STYLES = {
    "broken": {
        "color": "#D62728",
        "lw_mult": 1.4,
        "style": "-",
        "alpha": 0.95,
        "label": "Broken  (b_r>b_p)",
    },
    "formed": {
        "color": "#2CA02C",
        "lw_mult": 1.4,
        "style": "-",
        "alpha": 0.95,
        "label": "Formed  (b_r<b_p)",
    },
    "changed": {
        "color": "#FF7F0E",
        "lw_mult": 1.4,
        "style": "-",
        "alpha": 0.95,
        "label": "Changed (b_r≠b_p)",
    },
    "unchanged": {
        "color": "#888888",
        "lw_mult": 0.7,
        "style": "--",
        "alpha": 0.45,
        "label": "Unchanged",
    },
}


def _luminance(hex_color: str) -> float:
    h = hex_color.lstrip("#")
    r, g, b = (int(h[i : i + 2], 16) / 255.0 for i in (0, 2, 4))
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def _fill(el: str) -> str:
    return _ELEMENT_PALETTE.get(el, (_DEFAULT_FILL, _DEFAULT_BORDER))[0]


def _border(el: str) -> str:
    return _ELEMENT_PALETTE.get(el, (_DEFAULT_FILL, _DEFAULT_BORDER))[1]


def _pos_from_mol(
    mol: Chem.Mol, G: nx.Graph
) -> Optional[Dict[Any, Tuple[float, float]]]:
    """Return {node_id: (x, y)} from RDKit 2D conformer matched via atom_map attribute."""
    if mol.GetNumConformers() == 0:
        AllChem.Compute2DCoords(mol)
    conf = mol.GetConformer(0)
    mol_pos: Dict[int, Tuple[float, float]] = {
        a.GetAtomMapNum(): (
            conf.GetAtomPosition(a.GetIdx()).x,
            conf.GetAtomPosition(a.GetIdx()).y,
        )
        for a in mol.GetAtoms()
        if a.GetAtomMapNum() > 0
    }
    out: Dict[Any, Tuple[float, float]] = {}
    for n in G.nodes():
        am = G.nodes[n].get("atom_map", n)
        if am in mol_pos:
            out[n] = mol_pos[am]
        elif n in mol_pos:
            out[n] = mol_pos[n]
    return out if len(out) == G.number_of_nodes() else None


def _side_index(coordinate: str) -> Optional[int]:
    key = coordinate.lower().strip()
    if key in {"its", "auto", "none"}:
        return None
    if key in {"reactant", "reactants", "r", "left", "lhs"}:
        return 0
    if key in {"product", "products", "p", "right", "rhs"}:
        return 1
    raise ValueError(
        f"coordinate must be 'its', 'reactant', or 'product' (got {coordinate!r})"
    )


def _side_value(value: Any, idx: int) -> Any:
    if isinstance(value, (tuple, list)) and len(value) > idx:
        return value[idx]
    return value


def _as_order(value: Any, idx: int) -> float:
    value = _side_value(value, idx)
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def its_to_side_graph(its: nx.Graph, coordinate: str = "reactant") -> nx.Graph:
    """
    Project an ITS graph onto the reactant or product molecular graph.

    ITS edges store bond orders as ``(b_reactant, b_product)``.  This helper
    keeps only bonds present on the requested side and converts tuple-valued
    node attributes, such as formal charge or hydrogen count, to that side.
    """
    idx = _side_index(coordinate)
    if idx is None:
        raise ValueError("its_to_side_graph requires 'reactant' or 'product'")

    G = nx.Graph()
    for node, data in its.nodes(data=True):
        attrs = {key: _side_value(value, idx) for key, value in data.items()}
        attrs.setdefault("atom_map", node)
        if not attrs.get("atom_map"):
            attrs["atom_map"] = node
        G.add_node(node, **attrs)

    for u, v, data in its.edges(data=True):
        order = _as_order(data.get("order", (0.0, 0.0)), idx)
        if order <= 1e-6:
            continue
        attrs = {key: _side_value(value, idx) for key, value in data.items()}
        attrs["order"] = order
        attrs["aromatic"] = (
            bool(attrs.get("aromatic", False)) or abs(order - 1.5) < 1e-6
        )
        G.add_edge(u, v, **attrs)

    return G


def _set_atom_maps_from_graph(mol: Chem.Mol, G: nx.Graph) -> Chem.Mol:
    for atom, node in zip(mol.GetAtoms(), G.nodes()):
        atom_map = G.nodes[node].get("atom_map", node)
        try:
            atom_map = int(atom_map)
        except (TypeError, ValueError):
            continue
        if atom_map > 0:
            atom.SetAtomMapNum(atom_map)
    return mol


def its_to_side_mol(
    its: nx.Graph, coordinate: str = "reactant", sanitize: bool = True
) -> Chem.Mol:
    """
    Convert an ITS graph to the reactant-side or product-side RDKit molecule.

    This is the molecular counterpart to :func:`its_to_side_graph`.
    """
    side_graph = its_to_side_graph(its, coordinate)
    mol = graph_to_mol(side_graph, sanitize=sanitize)
    return _set_atom_maps_from_graph(mol, side_graph)


def _pos_from_its_coordinate(
    its: nx.Graph, coordinate: str
) -> Optional[Dict[Any, Tuple[float, float]]]:
    """Return RDKit 2D coordinates from a reactant/product ITS projection."""
    try:
        side_graph = its_to_side_graph(its, coordinate)
    except Exception:
        return None

    try:
        mol = its_to_side_mol(its, coordinate, sanitize=True)
    except Exception:
        try:
            mol = its_to_side_mol(its, coordinate, sanitize=False)
        except Exception:
            return None

    if mol.GetNumConformers() == 0:
        AllChem.Compute2DCoords(mol)
    conf = mol.GetConformer(0)
    nodes = list(side_graph.nodes())
    if len(nodes) != mol.GetNumAtoms():
        return None
    return {
        node: (
            conf.GetAtomPosition(atom_idx).x,
            conf.GetAtomPosition(atom_idx).y,
        )
        for atom_idx, node in enumerate(nodes)
    }


def its_coordinate_layout(
    its: nx.Graph, coordinate: str = "reactant"
) -> Dict[Any, Tuple[float, float]]:
    """
    Return RDKit 2D coordinates for an ITS graph on one reaction side.

    Use ``coordinate="reactant"`` to lay atoms out from reactant-side bonds, or
    ``coordinate="product"`` to lay them out from product-side bonds.
    """
    pos = _pos_from_its_coordinate(its, coordinate)
    if pos is None:
        raise ValueError(f"Could not compute {coordinate!r} coordinates for ITS graph")
    return pos


def _classify_edge(br, bp) -> str:
    br, bp = float(br), float(bp)
    if abs(br - bp) < 1e-6:
        return "unchanged"
    if br > 0 and bp < 1e-6:
        return "broken"
    if br < 1e-6 and bp > 0:
        return "formed"
    return "changed"


def visualize_its(  # noqa: C901
    its: nx.Graph,
    *,
    mol: Optional[Chem.Mol] = None,
    ax: Optional[plt.Axes] = None,
    title: Optional[str] = None,
    pos: Optional[Dict] = None,
    coordinate: str = "its",
    layout: str = "kamada_kawai",  # "spring" | "kamada_kawai" | "circular"
    node_size: int = 700,
    font_size: int = 9,
    edge_width: float = 2.4,
    show_edge_labels: bool = True,
    show_unchanged_edge_labels: bool = False,
    show_node_labels: bool = True,
    show_atom_map: bool = False,  # show :mapnum on label
    show_legend: bool = True,
    rc_ring_color: str = "#FFD700",  # gold ring on reaction-center nodes
) -> plt.Axes:
    """
    Visualize an ITS graph with CPK node colors, reaction-center rings,
    and colored edge types (broken/formed/changed/unchanged).

    Edge 'order' attribute expected as (b_r, b_p) tuple.

    :param coordinate: ``"its"`` uses a graph layout, ``"reactant"`` lays the
        ITS over the reactant-side molecular graph, and ``"product"`` lays it
        over the product-side molecular graph.
    """
    created_fig = False
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 5), facecolor="white")
        created_fig = True

    ax.set_facecolor("white")
    ax.set_axis_off()
    if title:
        ax.set_title(
            title, fontsize=font_size + 2, fontweight="bold", pad=8, color="#1a1a1a"
        )

    # ── layout ───────────────────────────────────────────────────────────
    if pos is None:
        coord_idx = _side_index(coordinate)
        if coord_idx is not None:
            pos = _pos_from_its_coordinate(its, coordinate)
        if pos is None and mol is not None:
            pos = _pos_from_mol(mol, its)
        if pos is None:
            if layout == "spring":
                pos = nx.spring_layout(its, seed=0, k=1.0)
            elif layout == "circular":
                pos = nx.circular_layout(its)
            else:
                try:
                    pos = nx.kamada_kawai_layout(its)
                except Exception:
                    pos = nx.spring_layout(its, seed=0)

    # ── classify edges ───────────────────────────────────────────────────
    rc_nodes: Set = set()
    edge_groups: Dict[str, List[Tuple]] = {k: [] for k in _EDGE_STYLES}
    edge_labels: Dict[str, Dict[Tuple, str]] = {k: {} for k in _EDGE_STYLES}

    for u, v, d in its.edges(data=True):
        br, bp = d.get("order", (0.0, 0.0))
        etype = _classify_edge(br, bp)
        edge_groups[etype].append((u, v))
        lbl = f"({float(br):g},{float(bp):g})"
        edge_labels[etype][(u, v)] = lbl
        if etype != "unchanged":
            rc_nodes.update([u, v])

    # ── node styling ─────────────────────────────────────────────────────
    nodelist = list(its.nodes())
    node_colors = []
    node_borders = []
    label_colors = []

    for n in nodelist:
        el = its.nodes[n].get("element", "?")
        fill = _fill(el)
        bord = _border(el)
        node_colors.append(fill)
        node_borders.append(bord)
        label_colors.append("white" if _luminance(fill) < 0.50 else "#1a1a1a")

    # ── draw reaction-center node glow ───────────────────────────────────
    rc_list = [n for n in nodelist if n in rc_nodes]
    if rc_list:
        nc = nx.draw_networkx_nodes(
            its,
            pos,
            nodelist=rc_list,
            node_size=int(node_size * 1.9),
            node_color=rc_ring_color,
            edgecolors="none",
            linewidths=0,
            ax=ax,
            alpha=0.20,
        )
        nc.set_zorder(1)

    # ── draw edges in order (unchanged first, then colored) ──────────────
    for etype in ("unchanged", "broken", "formed", "changed"):
        eg = edge_groups[etype]
        if not eg:
            continue
        st = _EDGE_STYLES[etype]
        nx.draw_networkx_edges(
            its,
            pos,
            ax=ax,
            edgelist=eg,
            width=edge_width * st["lw_mult"],
            edge_color=st["color"],
            alpha=st["alpha"],
            style=st["style"],
            arrows=False,
        )

    # ── draw nodes ───────────────────────────────────────────────────────
    nc = nx.draw_networkx_nodes(
        its,
        pos,
        nodelist=nodelist,
        node_size=node_size,
        node_color=node_colors,
        edgecolors=node_borders,
        linewidths=max(1.2, node_size**0.5 * 0.055),
        ax=ax,
    )
    nc.set_zorder(3)

    # RC ring outline
    if rc_list:
        nc = nx.draw_networkx_nodes(
            its,
            pos,
            nodelist=rc_list,
            node_size=int(node_size * 1.32),
            node_color="none",
            edgecolors=rc_ring_color,
            linewidths=max(2.0, node_size**0.5 * 0.10),
            ax=ax,
            alpha=0.9,
        )
        nc.set_zorder(4)

    # ── node labels ──────────────────────────────────────────────────────
    if show_node_labels:
        for i, n in enumerate(nodelist):
            el = its.nodes[n].get("element", "?")
            am = its.nodes[n].get("atom_map", 0)
            if show_atom_map and am:
                lbl = f"{el}:{am}"
            else:
                lbl = el
            x, y = pos[n]
            ax.text(
                x,
                y,
                lbl,
                ha="center",
                va="center",
                fontsize=font_size,
                fontweight="bold",
                color=label_colors[i],
                zorder=9,
                path_effects=[pe.withStroke(linewidth=1.5, foreground="none")],
            )

    # ── edge labels ──────────────────────────────────────────────────────
    if show_edge_labels:
        for etype in ("broken", "formed", "changed"):
            lbls = edge_labels[etype]
            if lbls:
                nx.draw_networkx_edge_labels(
                    its,
                    pos,
                    ax=ax,
                    edge_labels=lbls,
                    font_size=font_size - 1,
                    font_color=_EDGE_STYLES[etype]["color"],
                    bbox=dict(
                        boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85
                    ),
                )
        if show_unchanged_edge_labels and edge_labels["unchanged"]:
            nx.draw_networkx_edge_labels(
                its,
                pos,
                ax=ax,
                edge_labels=edge_labels["unchanged"],
                font_size=font_size - 2,
                font_color="#888888",
            )

    # ── legend ───────────────────────────────────────────────────────────
    if show_legend:
        legend_handles = []
        for etype in ("broken", "formed", "changed", "unchanged"):
            if not edge_groups[etype]:
                continue
            st = _EDGE_STYLES[etype]
            legend_handles.append(
                Line2D(
                    [0],
                    [0],
                    color=st["color"],
                    lw=2.2 * st["lw_mult"],
                    linestyle=st["style"],
                    alpha=min(st["alpha"] + 0.1, 1.0),
                    label=st["label"],
                )
            )
        if legend_handles:
            ax.legend(
                handles=legend_handles,
                loc="upper right",
                frameon=True,
                framealpha=0.92,
                fontsize=font_size - 1,
                edgecolor="#cccccc",
                handlelength=2.0,
                labelspacing=0.4,
            )

    if created_fig:
        plt.tight_layout()
        plt.show()

    return ax
