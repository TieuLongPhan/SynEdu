from __future__ import annotations

from collections import defaultdict
import math
from typing import Any, Dict, Iterable, List, Mapping, Optional, Tuple, Union
import warnings

from networkx.algorithms import isomorphism as iso
import networkx as nx


def node_match(n1, n2):
    """Return True if two nodes match on element, formal charge, and aromaticity."""
    return (
        n1.get("element") == n2.get("element")
        and int(n1.get("formal_charge", 0)) == int(n2.get("formal_charge", 0))
        and bool(n1.get("aromatic", False)) == bool(n2.get("aromatic", False))
    )


def edge_match(e1, e2):
    """Return True if two edges match on bond order and aromaticity."""

    def normalized_order(value):
        if isinstance(value, (tuple, list)):
            return tuple(float(item) for item in value)
        return float(value)

    return normalized_order(e1.get("order", 1.0)) == normalized_order(
        e2.get("order", 1.0)
    ) and bool(e1.get("aromatic", False)) == bool(e2.get("aromatic", False))


def enumerate_automorphisms(G: nx.Graph):
    """Return all automorphisms of G as a list of node-mapping dicts."""
    GM_self = iso.GraphMatcher(G, G, node_match=node_match, edge_match=edge_match)
    return list(GM_self.isomorphisms_iter())


def compute_orbits_from_automorphisms(G: nx.Graph, automorphisms=None):
    """Group nodes into symmetry orbits using union-find over automorphisms."""
    if automorphisms is None:
        automorphisms = enumerate_automorphisms(G)

    parent = {v: v for v in G.nodes()}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for auto in automorphisms:
        for v, fv in auto.items():
            union(v, fv)

    orbits = {}
    for v in G.nodes():
        r = find(v)
        orbits.setdefault(r, set()).add(v)

    return sorted(orbits.values(), key=lambda s: min(s))


def nx_subgraph_matches(
    host_G: nx.Graph,
    pattern_G: nx.Graph,
    *,
    invert: bool = True,
    node_match_fn=node_match,
    edge_match_fn=edge_match,
) -> List[Dict[Any, Any]]:
    """
    Enumerate labeled subgraph monomorphisms of ``pattern_G`` inside ``host_G``.

    NetworkX yields mappings in the ``host_node -> pattern_node`` direction.
    With ``invert=True`` (default), this helper returns the teaching-friendly
    ``pattern_node -> host_node`` mapping used throughout the talktorials. Host
    edges between matched nodes that are absent from the pattern are allowed.
    """
    matcher = iso.GraphMatcher(
        host_G,
        pattern_G,
        node_match=node_match_fn,
        edge_match=edge_match_fn,
    )
    out: List[Dict[Any, Any]] = []
    for host_to_pattern in matcher.subgraph_monomorphisms_iter():
        if invert:
            out.append({pattern: host for host, pattern in host_to_pattern.items()})
        else:
            out.append(dict(host_to_pattern))
    return out


def dedup_by_host_image_with_orbits(
    matches: Iterable[Dict[Any, Any]],
    *,
    orbits: Union[Iterable[Iterable[Any]], Mapping[Any, Any], None] = None,
    host_autos: Optional[Iterable[Mapping[Any, Any]]] = None,
    pattern_node_order: Optional[Iterable[Any]] = None,
    return_groups: bool = False,
) -> Union[List[Dict[Any, Any]], Dict[Tuple[Any, ...], List[Dict[Any, Any]]]]:
    """
    Deduplicate ``pattern -> host`` mappings by their canonical host image.

    ``host_autos`` may be supplied to identify whole host images related by one
    host automorphism. Individual vertex orbits are insufficient for this:
    orbits of vertices do not determine orbits of vertex sets.
    """
    matches = list(matches)
    if not matches:
        return {} if return_groups else []

    if orbits is not None and host_autos is None:
        warnings.warn(
            "The 'orbits' argument alone cannot canonicalize multi-node host "
            "images and is ignored; pass complete host automorphisms via "
            "'host_autos'.",
            DeprecationWarning,
            stacklevel=2,
        )
    autos = [dict(auto) for auto in host_autos] if host_autos is not None else []
    if pattern_node_order is None:
        pattern_node_order = tuple(sorted(matches[0].keys()))
    else:
        pattern_node_order = tuple(pattern_node_order)

    def canonical_host_image(mapping: Mapping[Any, Any]) -> Tuple[Any, ...]:
        image = tuple(sorted(mapping.values()))
        if not autos:
            return image
        transformed_images = []
        for index, auto in enumerate(autos):
            missing = [node for node in image if node not in auto]
            if missing:
                raise ValueError(
                    f"Host automorphism {index} does not cover image node(s) "
                    f"{missing!r}"
                )
            transformed_images.append(tuple(sorted(auto[node] for node in image)))
        return min(transformed_images)

    buckets: Dict[Tuple[Any, ...], List[Dict[Any, Any]]] = defaultdict(list)
    for mapping in matches:
        buckets[canonical_host_image(mapping)].append(mapping)

    if return_groups:
        return {key: buckets[key] for key in sorted(buckets.keys())}

    representatives: List[Dict[Any, Any]] = []
    for key in sorted(buckets.keys()):
        group = buckets[key]
        representatives.append(
            min(
                group, key=lambda mapping: tuple(mapping[p] for p in pattern_node_order)
            )
        )
    return representatives


def dedup_by_pattern_image_with_orbits(
    matches: Iterable[Dict[Any, Any]],
    *,
    pattern_autos: Iterable[Dict[Any, Any]],
    pattern_node_order: Optional[Iterable[Any]] = None,
    return_groups: bool = False,
) -> Union[List[Dict[Any, Any]], Dict[Tuple[Any, ...], List[Dict[Any, Any]]]]:
    """Deduplicate ``pattern -> host`` mappings modulo pattern automorphisms."""
    matches = list(matches)
    if not matches:
        return {} if return_groups else []

    if pattern_node_order is None:
        pattern_node_order = tuple(sorted(matches[0].keys()))
    else:
        pattern_node_order = tuple(pattern_node_order)

    autos = list(pattern_autos)

    def canonical_tuple(mapping: Dict[Any, Any]) -> Tuple[Any, ...]:
        best = None
        for auto in autos:
            value = tuple(mapping[auto[p]] for p in pattern_node_order)
            if best is None or value < best:
                best = value
        return (
            best if best is not None else tuple(mapping[p] for p in pattern_node_order)
        )

    buckets: Dict[Tuple[Any, ...], List[Dict[Any, Any]]] = defaultdict(list)
    for mapping in matches:
        buckets[canonical_tuple(mapping)].append(mapping)

    if return_groups:
        return {key: buckets[key] for key in sorted(buckets.keys())}

    representatives: List[Dict[Any, Any]] = []
    for key in sorted(buckets.keys()):
        group = buckets[key]
        representatives.append(
            min(
                group, key=lambda mapping: tuple(mapping[p] for p in pattern_node_order)
            )
        )
    return representatives


def rdkit_subgraph_matches(host_mol, pattern_mol, *, uniquify: bool = False):
    """Enumerate RDKit substructure matches as ``pattern_atom -> host_atom`` dicts."""
    matches = host_mol.GetSubstructMatches(pattern_mol, uniquify=uniquify)
    return [
        {pattern_idx: host_idx for pattern_idx, host_idx in enumerate(match)}
        for match in matches
    ]


def extract_induced_subgraph(
    G: nx.Graph, nodes: Iterable[Any], *, copy: bool = True
) -> nx.Graph:
    """Return the induced subgraph on ``nodes`` after checking that all exist."""
    nodes = set(nodes)
    missing = nodes - set(G.nodes)
    if missing:
        raise KeyError(f"Nodes not found in graph: {sorted(missing)}")
    H = G.subgraph(nodes)
    return H.copy() if copy else H


def print_graph_attributes(G: nx.Graph) -> None:
    """
    Print a human-readable summary of all node and edge attributes.

    :param G: NetworkX graph to inspect.
    :type G: networkx.Graph
    """
    print(f"Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    print("Nodes:")
    for node, data in G.nodes(data=True):
        attrs = ", ".join(f"{k}={v!r}" for k, v in data.items())
        print(f"  {node}: {attrs}")
    print("Edges:")
    for u, v, data in G.edges(data=True):
        attrs = ", ".join(f"{k}={v!r}" for k, v in data.items())
        print(f"  ({u}, {v}): {attrs}")


_ALLOWED_VALENCE = {
    "C": {4},
    "N": {3, 4},
    "O": {2},
    "S": {2, 4, 6},
    "P": {3, 5},
}

_CHARGED_ALLOWED_VALENCE = {
    ("C", -1): {3},
    ("C", 1): {3},
    ("N", -1): {2, 3},
    ("N", 1): {4},
    ("O", -1): {1},
    ("O", 1): {3},
    ("S", -1): {1, 3, 5},
    ("S", 1): {3, 5},
    ("P", 1): {4, 6},
}

_VALENCE_ELECTRONS = {
    "H": 1,
    "C": 4,
    "N": 5,
    "O": 6,
    "F": 7,
    "P": 5,
    "S": 6,
    "Cl": 7,
    "Br": 7,
    "I": 7,
    "Si": 4,
    "B": 3,
}


def build_be_matrix(G: nx.Graph):
    """
    Build the bond-electron matrix used in S01.

    Off-diagonal entries are bond orders. Diagonal entries approximate leftover
    valence electrons from element valence, formal charge, bonds, and hydrogens.
    """
    import numpy as np

    ordered = sorted(G.nodes())
    idx = {node: i for i, node in enumerate(ordered)}
    matrix = np.zeros((len(ordered), len(ordered)))

    for u, v, data in G.edges(data=True):
        order = float(data.get("order", 1.0))
        matrix[idx[u], idx[v]] = order
        matrix[idx[v], idx[u]] = order

    for node, data in G.nodes(data=True):
        element = data.get("element", "C")
        charge = int(data.get("formal_charge", 0))
        hcount = int(data.get("hcount", 0))
        if element not in _VALENCE_ELECTRONS:
            raise ValueError(
                f"No valence-electron model is defined for element {element!r}"
            )
        valence_electrons = _VALENCE_ELECTRONS[element]
        bond_electrons = (
            sum(
                float(G.edges[node, nbr].get("order", 1.0)) for nbr in G.neighbors(node)
            )
            + hcount
        )
        matrix[idx[node], idx[node]] = max(
            valence_electrons - charge - bond_electrons, 0
        )

    return matrix, ordered


def add_wildcard(G: nx.Graph, indices: Iterable[Any]) -> nx.Graph:
    """Add wildcard atoms using the smallest chemically allowed open valence."""
    H = G.copy()
    next_idx = max(H.nodes, default=0) + 1
    for atom in indices:
        if atom not in H:
            continue
        data = H.nodes[atom]
        element = data.get("element")
        if element not in _ALLOWED_VALENCE:
            continue
        formal_charge = int(data.get("formal_charge", 0))
        hcount = int(data.get("hcount", 0))
        bond_valence = sum(
            float(H.edges[atom, nbr].get("order", 1.0)) for nbr in H.neighbors(atom)
        )
        current_valence = hcount + bond_valence
        allowed = _CHARGED_ALLOWED_VALENCE.get(
            (element, formal_charge), _ALLOWED_VALENCE[element]
        )
        targets = sorted(
            valence for valence in allowed if valence + 1e-9 >= current_valence
        )
        if not targets:
            continue
        deficit = float(targets[0] - current_valence)
        if deficit <= 1e-9:
            continue
        bond_order = min(3, math.floor(deficit + 1e-9))
        if bond_order < 1:
            continue
        wildcard = next_idx
        next_idx += 1
        H.add_node(
            wildcard,
            element="*",
            aromatic=False,
            hcount=0,
            formal_charge=0,
            wildcard=True,
        )
        H.add_edge(atom, wildcard, order=float(bond_order), aromatic=False)
    return H


def find_wildcards(G: nx.Graph) -> List[Any]:
    """Return graph nodes labelled as wildcard atoms."""
    return [n for n, data in G.nodes(data=True) if data.get("element") == "*"]


def _get_unique_wildcard(G: nx.Graph, label: str) -> Any:
    wildcards = find_wildcards(G)
    if len(wildcards) != 1:
        raise ValueError(f"{label} must contain exactly one wildcard")
    return wildcards[0]


def merge_graphs_on_wildcard(
    G1: nx.Graph,
    G2: nx.Graph,
    *,
    wc1: Any = None,
    wc2: Any = None,
) -> nx.Graph:
    """Glue two molecular graphs by removing one wildcard atom from each graph."""
    wc1 = wc1 if wc1 is not None else _get_unique_wildcard(G1, "G1")
    wc2 = wc2 if wc2 is not None else _get_unique_wildcard(G2, "G2")

    nbrs1 = list(G1.neighbors(wc1))
    nbrs2 = list(G2.neighbors(wc2))
    if len(nbrs1) != 1 or len(nbrs2) != 1:
        raise ValueError("Only 1-1 wildcard gluing is supported")

    u1, u2 = nbrs1[0], nbrs2[0]
    merged = nx.Graph()

    for node, data in G1.nodes(data=True):
        if node != wc1:
            merged.add_node(node, **data)
    for u, v, data in G1.edges(data=True):
        if wc1 not in (u, v):
            merged.add_edge(u, v, **data)

    offset = max(merged.nodes, default=-1) + 1
    node_map = {}
    for node, data in G2.nodes(data=True):
        if node == wc2:
            continue
        node_map[node] = node + offset if isinstance(node, int) else f"{node}_{offset}"
        merged.add_node(node_map[node], **data)
    for u, v, data in G2.edges(data=True):
        if wc2 not in (u, v):
            merged.add_edge(node_map[u], node_map[v], **data)

    bond_order = min(
        float(G1.edges[wc1, u1].get("order", 1.0)),
        float(G2.edges[wc2, u2].get("order", 1.0)),
    )
    merged.add_edge(u1, node_map[u2], order=bond_order)
    return merged


def h_to_explicit(G: nx.Graph) -> nx.Graph:
    """
    Expand implicit hydrogen counts into explicit H nodes.

    Reads the ``hcount`` node attribute and adds one H node (connected by a
    single bond) per implicit hydrogen.  The ``hcount`` attribute is reset to 0
    on the heavy atom after expansion.

    :param G: Molecular graph with ``hcount`` node attributes.
    :type G: networkx.Graph
    :returns: New graph with explicit H nodes appended.
    :rtype: networkx.Graph
    """
    H = G.copy()
    next_id = max(H.nodes(), default=-1) + 1
    for node, data in list(G.nodes(data=True)):
        n_h = int(data.get("hcount", 0))
        for _ in range(n_h):
            H.add_node(
                next_id,
                element="H",
                formal_charge=0,
                aromatic=False,
                hcount=0,
                atom_map=0,
            )
            H.add_edge(node, next_id, order=1.0, aromatic=False)
            next_id += 1
        H.nodes[node]["hcount"] = 0
    return H


def h_to_implicit(G: nx.Graph) -> nx.Graph:
    """
    Collapse explicit H nodes back into ``hcount`` attributes on heavy atoms.

    An H node is absorbed only when it has exactly one neighbour and that
    neighbour is a heavy (non-H) atom.

    :param G: Molecular graph that may contain explicit H nodes.
    :type G: networkx.Graph
    :returns: New graph with H nodes removed and ``hcount`` updated.
    :rtype: networkx.Graph
    """
    H = G.copy()
    h_nodes = {n for n, d in G.nodes(data=True) if d.get("element") == "H"}
    for h in h_nodes:
        neighbors = list(G.neighbors(h))
        if len(neighbors) == 1:
            heavy = neighbors[0]
            if heavy not in h_nodes:
                H.nodes[heavy]["hcount"] = int(H.nodes[heavy].get("hcount", 0)) + 1
                H.remove_node(h)
    return H
