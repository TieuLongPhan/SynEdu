"""Graph clustering and hashing utilities for SynEdu.

Provides WL graph hashing (fast isomorphism pre-filter) and exact
graph-isomorphism clustering of ITS / reaction-center graphs.
"""

from __future__ import annotations

import copy
import random as _random
from collections import defaultdict
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import networkx as nx
from networkx.algorithms import isomorphism as iso

AttrSpec = Union[None, str, Sequence[str]]


# ---------------------------------------------------------------------------
# WL hashing (from S07)
# ---------------------------------------------------------------------------


def _ensure_defaults_nodes(G: nx.Graph, attr: str) -> None:
    for _, d in G.nodes(data=True):
        d.setdefault(attr, "")


def _ensure_defaults_edges(G: nx.Graph, attr: str) -> None:
    for _, _, d in G.edges(data=True):
        d.setdefault(attr, "")


def _combine_node_attrs(G: nx.Graph, attrs: Sequence[str], name: str) -> None:
    for _, d in G.nodes(data=True):
        d[name] = "|".join(str(d.get(a, "")) for a in attrs)


def _combine_edge_attrs(G: nx.Graph, attrs: Sequence[str], name: str) -> None:
    for _, _, d in G.edges(data=True):
        d[name] = "|".join(str(d.get(a, "")) for a in attrs)


def _prepare_graph(
    graph: nx.Graph,
    node_attrs: AttrSpec = None,
    edge_attrs: AttrSpec = None,
) -> Tuple[nx.Graph, Optional[str], Optional[str]]:
    G = copy.deepcopy(graph)

    if isinstance(node_attrs, str):
        n_name = node_attrs
        _ensure_defaults_nodes(G, n_name)
    elif node_attrs:
        nlist = list(node_attrs)
        n_name = nlist[0] if len(nlist) == 1 else "_wl_node"
        if len(nlist) == 1:
            _ensure_defaults_nodes(G, n_name)
        else:
            _combine_node_attrs(G, nlist, n_name)
    else:
        n_name = None

    if isinstance(edge_attrs, str):
        e_name = edge_attrs
        _ensure_defaults_edges(G, e_name)
    elif edge_attrs:
        elist = list(edge_attrs)
        e_name = elist[0] if len(elist) == 1 else "_wl_edge"
        if len(elist) == 1:
            _ensure_defaults_edges(G, e_name)
        else:
            _combine_edge_attrs(G, elist, e_name)
    else:
        e_name = None

    return G, n_name, e_name


def wl_hash(
    graph: nx.Graph,
    node_attrs: AttrSpec = ("element", "charge"),
    edge_attrs: AttrSpec = ("order",),
    iterations: int = 2,
    digest_size: int = 16,
) -> str:
    """
    Compute a Weisfeiler-Lehman graph hash for fast isomorphism pre-filtering.

    :param graph: NetworkX graph to hash.
    :param node_attrs: Node attribute(s) to include in the hash.
    :param edge_attrs: Edge attribute(s) to include in the hash.
    :param iterations: Number of WL iterations.
    :param digest_size: Hash digest size in bytes.
    :returns: Hexadecimal hash string.
    :rtype: str
    """
    G, n_attr, e_attr = _prepare_graph(graph, node_attrs, edge_attrs)
    return nx.weisfeiler_lehman_graph_hash(
        G,
        node_attr=n_attr,
        edge_attr=e_attr,
        iterations=iterations,
        digest_size=digest_size,
    )


# ---------------------------------------------------------------------------
# Exact isomorphism clustering
# ---------------------------------------------------------------------------


def cluster_its_graphs(
    graphs: List[nx.Graph],
    node_attrs: Sequence[str] = ("element",),
    edge_attrs: Sequence[str] = ("order",),
) -> List[int]:
    """
    Cluster ITS / reaction-center graphs by exact graph isomorphism.

    Uses a pairwise isomorphism check with WL hashing as a pre-filter.
    Returns a list of integer class labels (one per input graph).

    :param graphs: List of NetworkX graphs to cluster.
    :param node_attrs: Node attributes used for node matching.
    :param edge_attrs: Edge attributes used for edge matching.
    :returns: List of class IDs aligned with the input list.
    :rtype: List[int]
    """
    node_attrs = tuple(node_attrs)
    edge_attrs = tuple(edge_attrs)

    def _nm(n1: Dict, n2: Dict) -> bool:
        return all(n1.get(a) == n2.get(a) for a in node_attrs)

    def _em(e1: Dict, e2: Dict) -> bool:
        return all(e1.get(a) == e2.get(a) for a in edge_attrs)

    hashes = [wl_hash(G, node_attrs=node_attrs, edge_attrs=edge_attrs) for G in graphs]
    labels = [-1] * len(graphs)
    next_cls = 0

    for i, G in enumerate(graphs):
        if labels[i] != -1:
            continue
        labels[i] = next_cls
        for j in range(i + 1, len(graphs)):
            if labels[j] != -1 or hashes[i] != hashes[j]:
                continue
            gm = iso.GraphMatcher(G, graphs[j], node_match=_nm, edge_match=_em)
            if gm.is_isomorphic():
                labels[j] = next_cls
        next_cls += 1

    return labels


# ---------------------------------------------------------------------------
# Stratified sampling
# ---------------------------------------------------------------------------


def stratified_sample(
    records: List[Any],
    key: str,
    n_per_class: int = 1,
    seed: Optional[int] = None,
) -> List[Any]:
    """
    Return up to *n_per_class* records per unique value of *key*.

    :param records: List of dicts (or objects with ``__getitem__`` for *key*).
    :param key: Dict key whose value defines the class/group.
    :param n_per_class: Maximum records to sample per class.
    :param seed: Optional random seed for reproducibility.
    :returns: Sampled subset.
    :rtype: list
    """
    groups: Dict[Any, List[Any]] = defaultdict(list)
    for rec in records:
        groups[rec.get(key) if hasattr(rec, "get") else rec[key]].append(rec)

    rng = _random.Random(seed)
    result: List[Any] = []
    for group in groups.values():
        result.extend(rng.sample(group, min(n_per_class, len(group))))
    return result
