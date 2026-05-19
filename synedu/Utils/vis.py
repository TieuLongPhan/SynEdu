from __future__ import annotations

from typing import Optional, Tuple, Dict, Set
import math

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import networkx as nx
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

_ELEMENT_COLORS = {
    "C": "#4D4D4D",
    "O": "#D62728",
    "N": "#1F77B4",
    "S": "#FF7F0E",
    "Cl": "#2CA02C",
    "F": "#2CA02C",
    "Br": "#8C564B",
    "I": "#9467BD",
    "P": "#E377C2",
    "H": "#BBBBBB",
}


def _ensure_2d(mol: Chem.Mol) -> None:
    if mol.GetNumConformers() == 0:
        AllChem.Compute2DCoords(mol)


def _luminance(hex_color: str) -> float:
    h = hex_color.lstrip("#")
    r, g, b = (int(h[i : i + 2], 16) / 255.0 for i in (0, 2, 4))
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def _avg_edge_length(pos: Dict[int, Tuple[float, float]], G: nx.Graph) -> float:
    if G.number_of_edges() == 0:
        return 1.0
    ds = []
    for u, v in G.edges():
        x1, y1 = pos[int(u)]
        x2, y2 = pos[int(v)]
        ds.append(math.hypot(x2 - x1, y2 - y1))
    return sum(ds) / len(ds)


def _perp_offset(
    p1: Tuple[float, float], p2: Tuple[float, float], offset: float
) -> Tuple[float, float]:
    dx = p2[0] - p1[0]
    dy = p2[1] - p1[1]
    L = math.hypot(dx, dy)
    if L == 0:
        return (0.0, 0.0)
    ux, uy = dx / L, dy / L
    return (-uy * offset, ux * offset)


def _index_offset_vec(
    n: int,
    G: nx.Graph,
    pos: Dict[int, Tuple[float, float]],
    *,
    base: float,
) -> Tuple[float, float]:
    """
    Small local offset for index label near node n.

    Direction: away from the centroid of neighbors (keeps label near atom
    while avoiding bond lines). Fallback: slight upward offset.
    """
    x, y = pos[n]
    nbrs = [int(m) for m in G.neighbors(n)]
    if not nbrs:
        return (0.0, base)

    cx = sum(pos[m][0] for m in nbrs) / len(nbrs)
    cy = sum(pos[m][1] for m in nbrs) / len(nbrs)

    dx = x - cx
    dy = y - cy
    L = math.hypot(dx, dy)
    if L == 0:
        return (0.0, base)

    ux, uy = dx / L, dy / L
    return (ux * base, uy * base)


def _draw_bond_lines(
    ax: plt.Axes,
    p1: Tuple[float, float],
    p2: Tuple[float, float],
    *,
    order: int,
    aromatic: bool,
    aromatic_style: str,
    offset: float,
    lw: float,
) -> None:
    if aromatic and aromatic_style == "dashed":
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], linestyle="--", color="k", linewidth=lw)
        return

    # aromatic "circle" style: draw as single bond (circle drawn separately)
    if aromatic:
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], linestyle="-", color="k", linewidth=lw)
        return

    if order <= 1:
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], linestyle="-", color="k", linewidth=lw)
        return

    # double/triple as parallel lines
    dx, dy = _perp_offset(p1, p2, offset)
    if order == 2:
        ax.plot(
            [p1[0] + dx, p2[0] + dx], [p1[1] + dy, p2[1] + dy], color="k", linewidth=lw
        )
        ax.plot(
            [p1[0] - dx, p2[0] - dx], [p1[1] - dy, p2[1] - dy], color="k", linewidth=lw
        )
    elif order == 3:
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color="k", linewidth=lw * 0.95)
        ax.plot(
            [p1[0] + dx, p2[0] + dx],
            [p1[1] + dy, p2[1] + dy],
            color="k",
            linewidth=lw * 0.95,
        )
        ax.plot(
            [p1[0] - dx, p2[0] - dx],
            [p1[1] - dy, p2[1] - dy],
            color="k",
            linewidth=lw * 0.95,
        )
    else:
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], linestyle="-", color="k", linewidth=lw)


def _draw_aromatic_circles(
    ax: plt.Axes, G: nx.Graph, pos: Dict[int, Tuple[float, float]], scale: float
) -> None:
    # draw a circle for each aromatic cycle found by cycle_basis
    for cyc in nx.cycle_basis(G):
        if len(cyc) < 5:  # avoid tiny cycles
            continue
        if not all(bool(G.nodes[int(n)].get("aromatic", False)) for n in cyc):
            continue

        # check edges aromatic along the cycle
        ok = True
        for i in range(len(cyc)):
            u = int(cyc[i])
            v = int(cyc[(i + 1) % len(cyc)])
            if not bool(G.edges[u, v].get("aromatic", False)):
                ok = False
                break
        if not ok:
            continue

        xs = [pos[int(n)][0] for n in cyc]
        ys = [pos[int(n)][1] for n in cyc]
        cx = sum(xs) / len(xs)
        cy = sum(ys) / len(ys)
        rs = [math.hypot(x - cx, y - cy) for x, y in zip(xs, ys)]
        r = (sum(rs) / len(rs)) * scale
        ax.add_patch(mpatches.Circle((cx, cy), r, fill=False, linewidth=1.1, color="k"))


def draw_molecular_graph(
    G: nx.Graph,
    mol: Optional[Chem.Mol] = None,
    *,
    ax: Optional[plt.Axes] = None,
    include_mol: bool = False,
    label_mode: str = "hetero",  # "all" | "hetero" | "none"
    show_indices: bool = False,
    indices_for_carbons: bool = True,
    show_bond_labels: bool = False,
    aromatic_style: str = "circle",  # "circle" | "dashed"
    node_size: int = 260,
    bond_lw: float = 1.6,
    figsize: Tuple[float, float] = (10, 5),
    # --- highlighting (e.g. MCS) ---
    highlight_nodes: Optional[Set[int]] = None,
    highlight_edges: Optional[Set[Tuple[int, int]]] = None,
    highlight_color: str = "#FFB000",
    highlight_alpha: float = 0.9,
    # --- typography ---
    element_fontsize: int = 10,
    index_fontsize: int = 12,  # bigger by default
) -> plt.Axes | Tuple[plt.Figure, Tuple[plt.Axes, plt.Axes]]:
    """
    Visualize a molecular NetworkX graph (and optionally an RDKit 2D depiction).

    - Stable node ordering ensures colors/labels align.
    - Indices are placed close to atoms using a *local* offset away from neighbors.
    """
    aromatic_style = aromatic_style.lower()
    if aromatic_style not in {"circle", "dashed"}:
        raise ValueError("aromatic_style must be 'circle' or 'dashed'")

    label_mode = label_mode.lower()
    if label_mode not in {"all", "hetero", "none"}:
        raise ValueError("label_mode must be 'all', 'hetero', or 'none'")

    # normalize highlight edges (undirected)
    hl_edges_norm: Set[Tuple[int, int]] = set()
    if highlight_edges:
        for u, v in highlight_edges:
            uu, vv = int(u), int(v)
            hl_edges_norm.add((min(uu, vv), max(uu, vv)))

    created_fig = False
    if include_mol:
        fig, (ax_mol, ax_g) = plt.subplots(1, 2, figsize=figsize)
        created_fig = True
    else:
        if ax is None:
            fig, ax_g = plt.subplots(figsize=(figsize[1], figsize[1]))  # square-ish
            created_fig = True
        else:
            ax_g = ax
            fig = ax_g.figure

    # ---------------------------
    # stable node order (critical)
    # ---------------------------
    nodelist = sorted(int(n) for n in G.nodes())

    # ---------------------------
    # positions: prefer RDKit 2D
    # ---------------------------
    pos: Dict[int, Tuple[float, float]] = {}
    if mol is not None:
        _ensure_2d(mol)
        conf = mol.GetConformer(0)
        for n in nodelist:
            p = conf.GetAtomPosition(int(n))
            pos[int(n)] = (p.x, p.y)
    else:
        pos = {
            int(k): (float(v[0]), float(v[1]))
            for k, v in nx.spring_layout(G, seed=42).items()
        }

    avg_len = _avg_edge_length(pos, G)
    bond_offset = avg_len * 0.08
    idx_offset = avg_len * 0.10  # keeps indices close to nodes

    # ---------------------------
    # nodes + element labels
    # ---------------------------
    node_colors: list[str] = []
    element_labels: Dict[int, str] = {}
    label_color: Dict[int, str] = {}

    for n in nodelist:
        data = G.nodes[n]
        el = str(data.get("element", "C"))
        c = _ELEMENT_COLORS.get(el, "#7F7F7F")
        node_colors.append(c)

        if label_mode == "none":
            txt = ""
        elif label_mode == "hetero" and el == "C":
            txt = ""
        else:
            txt = el
        element_labels[n] = txt

        label_color[n] = "white" if _luminance(c) < 0.55 else "black"

    nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=nodelist,
        node_size=node_size,
        node_color=node_colors,
        linewidths=0.0,
        ax=ax_g,
    )

    # highlight nodes as an outline overlay
    if highlight_nodes:
        hl = [int(n) for n in highlight_nodes if int(n) in G]
        if hl:
            nx.draw_networkx_nodes(
                G,
                pos,
                nodelist=hl,
                node_size=int(node_size * 1.18),
                node_color="none",
                edgecolors=highlight_color,
                linewidths=2.6,
                ax=ax_g,
                alpha=highlight_alpha,
            )

    # ---------------------------
    # edges
    # ---------------------------
    for u, v, data in G.edges(data=True):
        u, v = int(u), int(v)
        p1 = pos[u]
        p2 = pos[v]

        aromatic = bool(data.get("aromatic", False))
        if aromatic:
            order = 1
        else:
            try:
                order = int(round(abs(float(data.get("order", 1.0)))))
            except Exception:
                order = 1

        # highlight edge as thick underlay (so black bond stays visible)
        if hl_edges_norm and ((min(u, v), max(u, v)) in hl_edges_norm):
            ax_g.plot(
                [p1[0], p2[0]],
                [p1[1], p2[1]],
                color=highlight_color,
                linewidth=bond_lw * 3.0,
                alpha=highlight_alpha,
                solid_capstyle="round",
                zorder=0,
            )

        _draw_bond_lines(
            ax_g,
            p1,
            p2,
            order=order,
            aromatic=aromatic,
            aromatic_style=aromatic_style,
            offset=bond_offset,
            lw=bond_lw,
        )

        if show_bond_labels and (not aromatic):
            mx, my = 0.5 * (p1[0] + p2[0]), 0.5 * (p1[1] + p2[1])
            ax_g.text(
                mx,
                my,
                str(order),
                fontsize=8,
                ha="center",
                va="center",
                color="k",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.9),
                zorder=8,
            )

    if aromatic_style == "circle":
        _draw_aromatic_circles(ax_g, G, pos, scale=0.55)

    # element labels (centered)
    for n in nodelist:
        txt = element_labels.get(n, "")
        if not txt:
            continue
        x, y = pos[n]
        ax_g.text(
            x,
            y,
            txt,
            ha="center",
            va="center",
            fontsize=element_fontsize,
            color=label_color[n],
            zorder=9,
        )

    # ---------------------------
    # index labels (near nodes, bigger font)
    # ---------------------------
    if show_indices:
        for n in nodelist:
            el = str(G.nodes[n].get("element", "C"))
            if (label_mode == "hetero") and (el == "C") and (not indices_for_carbons):
                continue

            x, y = pos[n]
            dx, dy = _index_offset_vec(n, G, pos, base=idx_offset)
            ix, iy = x + dx, y + dy

            ax_g.text(
                ix,
                iy,
                str(n),
                fontsize=index_fontsize,
                ha="center",
                va="center",
                color="k",
                bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.95),
                zorder=10,
            )

    ax_g.set_axis_off()
    ax_g.set_aspect("equal")

    # ---------------------------
    # optional RDKit panel
    # ---------------------------
    if include_mol:
        ax_mol.set_axis_off()
        if mol is not None:
            _ensure_2d(mol)
            # RDKit panel can show atom indices too
            try:
                dopt = Draw.MolDrawOptions()
                dopt.addAtomIndices = bool(show_indices)
                img = Draw.MolToImage(
                    mol, size=(500, 500), kekulize=False, options=dopt
                )
            except Exception:
                img = Draw.MolToImage(mol, size=(500, 500), kekulize=False)
            ax_mol.imshow(img)

        if created_fig:
            fig.tight_layout()
        return fig, (ax_mol, ax_g)

    if created_fig:
        fig.tight_layout()
    return ax_g
