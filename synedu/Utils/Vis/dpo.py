from __future__ import annotations

import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from dataclasses import dataclass
from typing import Dict, Tuple, Optional, Any


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


# ============================================================
# ITS visualizer (your original; unchanged)
# ============================================================
def visualize_its(
    its: nx.Graph,
    *,
    ax=None,
    title: str | None = None,
    pos: dict | None = None,
    layout: str = "kamada_kawai",  # "spring" | "kamada_kawai" | "circular"
    node_size: int = 900,
    font_size: int = 10,
    edge_width: float = 2.8,
    show_edge_labels: bool = True,
    show_unchanged_edge_labels: bool = False,
    show_node_labels: bool = True,
    show_legends: bool = False,
):
    """
    Visualize an ITS graph.

    Nodes:
      - colored by element (its.nodes[n]['element'])
      - thicker outline if in reaction center

    Edges (its[u][v]['order'] == (br, bp)):
      - br > bp : broken  (red)
      - br < bp : formed  (green)
      - br = bp : unchanged (black, thinner)
    """
    created_fig = False
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))
        created_fig = True

    ax.set_axis_off()
    if title:
        ax.set_title(title, fontsize=font_size + 1)

    if pos is None:
        if layout == "spring":
            pos = nx.spring_layout(its, seed=0, k=0.9)
        elif layout == "circular":
            pos = nx.circular_layout(its)
        else:
            pos = nx.kamada_kawai_layout(its)

    broken, formed, unchanged = [], [], []
    lbl_b, lbl_f, lbl_u = {}, {}, {}
    rc_nodes = set()

    for u, v, d in its.edges(data=True):
        br, bp = d.get("order", (0.0, 0.0))
        if br > bp:
            broken.append((u, v))
            lbl_b[(u, v)] = f"({br:g},{bp:g})"
            rc_nodes.update([u, v])
        elif br < bp:
            formed.append((u, v))
            lbl_f[(u, v)] = f"({br:g},{bp:g})"
            rc_nodes.update([u, v])
        else:
            unchanged.append((u, v))
            lbl_u[(u, v)] = f"({br:g},{bp:g})"

    node_colors, node_edgecolors, node_linewidths = [], [], []
    elems_present = []
    for n, d in its.nodes(data=True):
        elem = d.get("element", "?")
        if elem not in elems_present:
            elems_present.append(elem)

        node_colors.append(_ELEMENT_COLORS.get(elem, "#CCCCCC"))
        if n in rc_nodes:
            node_edgecolors.append("#000000")
            node_linewidths.append(2.2)
        else:
            node_edgecolors.append("#666666")
            node_linewidths.append(1.2)

    labels = None
    if show_node_labels:
        labels = {n: f"{n}:{its.nodes[n].get('element','?')}" for n in its.nodes()}

    nx.draw_networkx_nodes(
        its,
        pos,
        ax=ax,
        node_size=node_size,
        node_color=node_colors,
        edgecolors=node_edgecolors,
        linewidths=node_linewidths,
    )

    nx.draw_networkx_edges(
        its,
        pos,
        ax=ax,
        edgelist=unchanged,
        width=max(1.1, edge_width * 0.55),
        edge_color="black",
        alpha=0.55,
    )
    nx.draw_networkx_edges(
        its,
        pos,
        ax=ax,
        edgelist=broken,
        width=edge_width * 1.25,
        edge_color="red",
        alpha=0.95,
    )
    nx.draw_networkx_edges(
        its,
        pos,
        ax=ax,
        edgelist=formed,
        width=edge_width * 1.25,
        edge_color="green",
        alpha=0.95,
    )

    if labels is not None:
        nx.draw_networkx_labels(its, pos, ax=ax, labels=labels, font_size=font_size)

    if show_edge_labels:
        if lbl_b:
            nx.draw_networkx_edge_labels(
                its, pos, ax=ax, edge_labels=lbl_b, font_size=font_size - 1
            )
        if lbl_f:
            nx.draw_networkx_edge_labels(
                its, pos, ax=ax, edge_labels=lbl_f, font_size=font_size - 1
            )
        if show_unchanged_edge_labels and lbl_u:
            nx.draw_networkx_edge_labels(
                its, pos, ax=ax, edge_labels=lbl_u, font_size=font_size - 3
            )

    if show_legends:
        edge_legend = [
            Line2D([0], [0], color="red", lw=3, label="broken (br>bp)"),
            Line2D([0], [0], color="green", lw=3, label="formed (br<bp)"),
            Line2D([0], [0], color="black", lw=2, alpha=0.6, label="unchanged (br=bp)"),
        ]
        elem_legend = [
            Patch(
                facecolor=_ELEMENT_COLORS.get(e, "#CCCCCC"),
                edgecolor="#666666",
                label=e,
            )
            for e in elems_present
        ]

        leg1 = ax.legend(handles=edge_legend, loc="upper left", frameon=False)
        if elem_legend:
            leg2 = ax.legend(
                handles=elem_legend,
                loc="lower left",
                frameon=False,
                ncol=min(6, len(elem_legend)),
            )
            ax.add_artist(leg1)

    if created_fig:
        plt.tight_layout()
        plt.show()

    return ax


# ============================================================
# DPO utilities
# ============================================================
@dataclass(frozen=True)
class DPODecomp:
    K_nodes: set
    K_edges: set  # preserved edges (same order)
    L_only_edges: set  # removed/changed (delete in L)
    R_only_edges: set  # added/changed (add in R)
    L_orders: Dict[Tuple[int, int], float]
    R_orders: Dict[Tuple[int, int], float]


def _ekey(u: int, v: int) -> Tuple[int, int]:
    return (u, v) if u < v else (v, u)


def _edge_orders(G: nx.Graph) -> Dict[Tuple[int, int], float]:
    out: Dict[Tuple[int, int], float] = {}
    for u, v, d in G.edges(data=True):
        out[_ekey(u, v)] = float(d.get("order", 1.0))
    return out


def dpo_decompose_atom_conserving(
    L: nx.Graph,
    R: nx.Graph,
    *,
    order_tol: float = 1e-9,
) -> DPODecomp:
    """
    Atom-conserving DPO-like decomposition:
      - K_nodes: nodes present on both sides (by node id)
      - K_edges: edges present on both sides with same 'order'
      - L_only_edges: edges deleted or order-changed (treated as delete)
      - R_only_edges: edges created or order-changed (treated as add)
    """
    K_nodes = set(L.nodes()) & set(R.nodes())
    L_orders = _edge_orders(L)
    R_orders = _edge_orders(R)

    common = set(L_orders) & set(R_orders)
    K_edges = {e for e in common if abs(L_orders[e] - R_orders[e]) <= order_tol}

    L_only = set(L_orders) - K_edges
    R_only = set(R_orders) - K_edges

    return DPODecomp(
        K_nodes=K_nodes,
        K_edges=K_edges,
        L_only_edges=L_only,
        R_only_edges=R_only,
        L_orders=L_orders,
        R_orders=R_orders,
    )


def build_its_from_LR(
    L: nx.Graph, R: nx.Graph, *, dec: Optional[DPODecomp] = None
) -> nx.Graph:
    """Build ITS with edge attribute order=(br,bp) from L and R."""
    if dec is None:
        dec = dpo_decompose_atom_conserving(L, R)

    its = nx.Graph()

    for n in set(L.nodes()) | set(R.nodes()):
        if n in L.nodes:
            its.add_node(n, **L.nodes[n])
        else:
            its.add_node(n, **R.nodes[n])

    all_edges = set(dec.L_orders) | set(dec.R_orders)
    for u, v in all_edges:
        br = dec.L_orders.get((u, v), 0.0)
        bp = dec.R_orders.get((u, v), 0.0)
        its.add_edge(u, v, order=(br, bp))

    return its


def _layout_pos(G: nx.Graph, layout: str, seed: int = 0) -> Dict[Any, Any]:
    if layout == "spring":
        return nx.spring_layout(G, seed=seed, k=0.9)
    if layout == "circular":
        return nx.circular_layout(G)
    if layout == "kamada_kawai":
        return nx.kamada_kawai_layout(G)
    raise ValueError("layout must be one of: 'spring', 'kamada_kawai', 'circular'")


# ============================================================
# DPO visualizer (NO change-hint texts anymore)
# ============================================================
def visualize_dpo_rule(
    L: nx.Graph,
    R: nx.Graph,
    *,
    use_its: bool = False,
    ax=None,  # None -> create a new 1x3 figure, otherwise pass axes array-like of length 3
    title: str | None = None,
    layout: str = "kamada_kawai",
    pos: dict | None = None,
    seed: int = 0,
    node_size: int = 900,
    font_size: int = 10,
    edge_width: float = 2.8,
    show_edge_labels: bool = True,
    show_node_labels: bool = True,
    show_legends: bool = True,
    show_unchanged_edge_labels: bool = False,  # only for ITS
    order_tol: float = 1e-9,
):
    """
    Visualize a DPO rule as L | K(or ITS) | R, sharing the same layout.

    - use_its=False: middle panel is K (context), L/R highlight changes:
        L: removed/changed edges = dashed red
        R: added/changed edges   = dashed green
    - use_its=True: middle panel is ITS (br,bp) with your visualize_its semantics.

    IMPORTANT: No "dashed = ..." hint texts are drawn anywhere.
    """
    created_fig = False
    if ax is None:
        fig, axes = plt.subplots(1, 3, figsize=(13, 4))
        created_fig = True
    else:
        axes = ax

    dec = dpo_decompose_atom_conserving(L, R, order_tol=order_tol)

    # build K (context) graph: preserved nodes+edges only
    K = nx.Graph()
    for n in dec.K_nodes:
        # prefer L attrs, fallback to R
        if n in L.nodes:
            K.add_node(n, **L.nodes[n])
        else:
            K.add_node(n, **R.nodes[n])

    for u, v in dec.K_edges:
        # order is same in L/R; keep from L if possible
        if L.has_edge(u, v):
            K.add_edge(u, v, **L.get_edge_data(u, v))
        else:
            K.add_edge(u, v, **R.get_edge_data(u, v))

    its = build_its_from_LR(L, R, dec=dec) if use_its else None

    # shared layout on union for stability
    if pos is None:
        union = nx.compose(L, R)
        pos = _layout_pos(union, layout=layout, seed=seed)

    # reaction center nodes = nodes incident to any changed edge
    rc_nodes = set()
    for e in dec.L_only_edges | dec.R_only_edges:
        rc_nodes.update(e)

    def _node_color(G: nx.Graph, n: Any) -> str:
        elem = G.nodes[n].get("element", "?")
        return _ELEMENT_COLORS.get(elem, "#CCCCCC")

    def _node_label(G: nx.Graph, n: Any) -> str:
        elem = G.nodes[n].get("element", "?")
        amap = G.nodes[n].get("atom_map", n)
        return f"{amap}:{elem}"

    def _draw_panel(
        G: nx.Graph,
        axp,
        *,
        panel_title: str,
        dashed_edges: set[Tuple[int, int]] = set(),
        dashed_color: str = "red",
    ):
        axp.set_axis_off()
        axp.set_title(panel_title, fontsize=font_size + 1)

        nodes = list(G.nodes())
        node_colors = [_node_color(G, n) for n in nodes]

        edgecolors, linewidths = [], []
        for n in nodes:
            if n in rc_nodes:
                edgecolors.append("#000000")
                linewidths.append(2.6)
            else:
                edgecolors.append("#666666")
                linewidths.append(1.2)

        nx.draw_networkx_nodes(
            G,
            pos,
            ax=axp,
            node_size=node_size,
            node_color=node_colors,
            edgecolors=edgecolors,
            linewidths=linewidths,
        )

        if show_node_labels:
            labels = {n: _node_label(G, n) for n in nodes}
            nx.draw_networkx_labels(G, pos, ax=axp, labels=labels, font_size=font_size)

        # split edges into preserved (solid) vs changed (dashed)
        solid = []
        dashed = []
        for u, v in G.edges():
            if _ekey(u, v) in dashed_edges:
                dashed.append((u, v))
            else:
                solid.append((u, v))

        if solid:
            nx.draw_networkx_edges(
                G,
                pos,
                ax=axp,
                edgelist=solid,
                width=max(1.1, edge_width * 0.65),
                edge_color="black",
                alpha=0.60,
            )
        if dashed:
            nx.draw_networkx_edges(
                G,
                pos,
                ax=axp,
                edgelist=dashed,
                width=edge_width * 1.35,
                edge_color=dashed_color,
                style="dashed",
                alpha=0.95,
            )

        if show_edge_labels:
            elabs = {}
            for u, v, d in G.edges(data=True):
                o = float(d.get("order", 1.0))
                s = str(int(o)) if o.is_integer() else f"{o:g}"
                elabs[(u, v)] = s
            if elabs:
                nx.draw_networkx_edge_labels(
                    G, pos, ax=axp, edge_labels=elabs, font_size=font_size - 1
                )

    if title and created_fig:
        fig.suptitle(title, fontsize=font_size + 2)

    # L panel
    _draw_panel(
        L,
        axes[0],
        panel_title="L pattern",
        dashed_edges=dec.L_only_edges,
        dashed_color="red",
    )

    # middle panel: K or ITS
    if use_its:
        visualize_its(
            its,  # type: ignore[arg-type]
            ax=axes[1],
            title="ITS graph",
            pos=pos,
            layout=layout,
            node_size=node_size,
            font_size=font_size,
            edge_width=edge_width,
            show_edge_labels=show_edge_labels,
            show_unchanged_edge_labels=show_unchanged_edge_labels,
            show_node_labels=show_node_labels,
            show_legends=False,  # we control legend globally below
        )
    else:
        _draw_panel(K, axes[1], panel_title="K (context)")

    # R panel
    _draw_panel(
        R,
        axes[2],
        panel_title="R pattern",
        dashed_edges=dec.R_only_edges,
        dashed_color="green",
    )

    # legends (optional)
    if show_legends:
        if use_its:
            edge_handles = [
                Line2D([0], [0], color="red", lw=3, label="broken (br>bp)"),
                Line2D([0], [0], color="green", lw=3, label="formed (br<bp)"),
                Line2D(
                    [0], [0], color="black", lw=2, alpha=0.6, label="unchanged (br=bp)"
                ),
            ]
        else:
            edge_handles = [
                Line2D([0], [0], color="black", lw=2, alpha=0.6, label="preserved"),
                Line2D(
                    [0],
                    [0],
                    color="red",
                    lw=3,
                    linestyle="--",
                    label="removed/changed (in L)",
                ),
                Line2D(
                    [0],
                    [0],
                    color="green",
                    lw=3,
                    linestyle="--",
                    label="added/changed (in R)",
                ),
            ]

        elems = []
        for _, d in nx.compose(L, R).nodes(data=True):
            e = d.get("element", "?")
            if e not in elems:
                elems.append(e)

        elem_handles = [
            Patch(
                facecolor=_ELEMENT_COLORS.get(e, "#CCCCCC"),
                edgecolor="#666666",
                label=e,
            )
            for e in elems
        ]

        # put both legends on the middle axis
        leg1 = axes[1].legend(handles=edge_handles, loc="upper left", frameon=False)
        if elem_handles:
            leg2 = axes[1].legend(
                handles=elem_handles,
                loc="lower left",
                frameon=False,
                ncol=min(6, len(elem_handles)),
            )
            axes[1].add_artist(leg1)

    if created_fig:
        plt.tight_layout()
        plt.show()

    return axes, dec, pos, its
