import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

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
    Visualize an ITS graph on a given Matplotlib axis.

    Nodes:
      - colored by element (ITS.nodes[n]['element'])
      - thicker black outline if in reaction center

    Edges (ITS[u][v]['order'] == (br, bp)):
      - br > bp : broken  (red)
      - br < bp : formed  (green)
      - br = bp : unchanged (black, thinner)
    """
    # axis handling (supports subplots)
    created_fig = False
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))
        created_fig = True

    ax.set_axis_off()
    if title:
        ax.set_title(title, fontsize=font_size + 1)

    # layout
    if pos is None:
        if layout == "spring":
            pos = nx.spring_layout(its, seed=0, k=0.9)
        elif layout == "circular":
            pos = nx.circular_layout(its)
        else:
            pos = nx.kamada_kawai_layout(its)

    # classify edges + reaction center nodes
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

    # node styling
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

    # node labels
    labels = None
    if show_node_labels:
        labels = {n: f"{n}:{its.nodes[n].get('element', '?')}" for n in its.nodes()}

    # draw nodes
    nx.draw_networkx_nodes(
        its,
        pos,
        ax=ax,
        node_size=node_size,
        node_color=node_colors,
        edgecolors=node_edgecolors,
        linewidths=node_linewidths,
    )

    # draw edges (unchanged behind)
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

    # labels
    if labels is not None:
        nx.draw_networkx_labels(its, pos, ax=ax, labels=labels, font_size=font_size)

    # edge labels
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

    # legends
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

        # two legends: add one, then re-add as artist
        leg1 = ax.legend(handles=edge_legend, loc="upper left", frameon=False)
        if elem_legend:
            ax.legend(
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
