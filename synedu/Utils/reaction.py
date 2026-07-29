"""Reaction-level graph utilities for SynEdu.

Converts between reaction SMILES and graph representations, and builds
Imaginary Transition State (ITS) graphs for reaction analysis.
"""

from __future__ import annotations

from collections import Counter
from typing import Any, Dict, List, Literal, Optional, Tuple

import networkx as nx
from rdkit import Chem

from .conversion import smiles_to_graph


def rsmi_to_graph(
    rsmi: str,
    drop_non_aam: bool = False,
    use_index_as_atom_map: bool = False,
) -> Tuple[nx.Graph, nx.Graph]:
    """
    Convert a reaction SMILES string to a pair of reactant / product graphs.

    Each side of the reaction is parsed into a single (possibly disconnected)
    NetworkX graph where multiple components appear as separate connected
    components.  Node IDs default to atom-map numbers when present, falling
    back to globally-unique atom indices.

    :param rsmi: Reaction SMILES string (``reactants>>products`` or
        ``reactants>agents>products``).
    :param drop_non_aam: If True, skip atoms without atom-map numbers.
    :param use_index_as_atom_map: If True, use atom index as node ID and
        ``atom_map`` attribute; otherwise prefer the SMILES atom-map number.
    :returns: Tuple ``(reactant_graph, product_graph)``.
    :rtype: Tuple[networkx.Graph, networkx.Graph]
    """
    parts = rsmi.split(">")
    reactant_smi = parts[0]
    product_smi = parts[-1]
    kwargs = dict(
        drop_non_aam=drop_non_aam, use_index_as_atom_map=use_index_as_atom_map
    )
    r = smiles_to_graph(reactant_smi, **kwargs)
    p = smiles_to_graph(product_smi, **kwargs)
    return r, p


def graph_to_rsmi(reactants: nx.Graph, products: nx.Graph) -> str:
    """
    Convert reactant and product graphs to an atom-mapped reaction SMILES.

    Node ``atom_map`` attributes are used to set RDKit atom-map numbers so
    the resulting SMILES carries full AAM information.  Disconnected graphs
    produce multi-component SMILES joined with ``"."``.

    :param reactants: NetworkX graph of reactants.
    :param products: NetworkX graph of products.
    :returns: Reaction SMILES string ``"r_smi>>p_smi"``.
    :rtype: str
    """
    from .conversion import graph_to_smi

    r_smi = graph_to_smi(reactants)
    p_smi = graph_to_smi(products)
    return f"{r_smi}>>{p_smi}"


AUX_LIBRARY: Dict[str, Counter[str]] = {
    "O": Counter({"H": 2, "O": 1}),
    "[H+]": Counter({"H": 1, "__charge__": 1}),
    "[OH-]": Counter({"H": 1, "O": 1, "__charge__": -1}),
}


def parse_reaction_smiles(rsmi: str) -> Tuple[List[Chem.Mol], List[Chem.Mol]]:
    """Parse a reaction SMILES string into reactant and product molecule lists."""
    try:
        lhs, rhs = rsmi.split(">>")
    except ValueError as exc:
        raise ValueError("Reaction SMILES must contain '>>'") from exc

    reactants = [Chem.MolFromSmiles(s) for s in lhs.split(".") if s]
    products = [Chem.MolFromSmiles(s) for s in rhs.split(".") if s]
    if any(m is None for m in reactants + products):
        raise ValueError("Invalid molecule in reaction SMILES")
    return reactants, products


def molecular_formula(mol: Chem.Mol) -> Counter[str]:
    """Return element counts and net formal charge, including hydrogens."""
    mol_h = Chem.AddHs(mol)
    formula: Counter[str] = Counter()
    for atom in mol_h.GetAtoms():
        formula[atom.GetSymbol()] += 1
        formula["__charge__"] += int(atom.GetFormalCharge())
    if formula["__charge__"] == 0:
        del formula["__charge__"]
    return formula


def total_formula(mols: List[Chem.Mol]) -> Counter[str]:
    """Return the summed molecular formula for a list of molecules."""
    total: Counter[str] = Counter()
    for mol in mols:
        for key, value in molecular_formula(mol).items():
            total[key] += value
    if total["__charge__"] == 0:
        del total["__charge__"]
    return total


def signed_difference(a: Counter[str], b: Counter[str]) -> Counter[str]:
    """Compute ``a - b`` while preserving negative element counts."""
    delta: Counter[str] = Counter()
    for elem in set(a) | set(b):
        value = a.get(elem, 0) - b.get(elem, 0)
        if value:
            delta[elem] = value
    return delta


def reaction_imbalance(
    rsmi: str,
) -> Tuple[
    Counter[str],
    Literal[
        "balanced",
        "missing_in_products",
        "missing_in_reactants",
        "missing_in_both_sides",
    ],
]:
    """
    Compute atom-count and net-charge imbalance as reactants minus products.

    Positive counts indicate atoms missing from the products; negative counts
    indicate atoms missing from the reactants.
    """
    reactants, products = parse_reaction_smiles(rsmi)
    delta = signed_difference(total_formula(reactants), total_formula(products))
    if not delta:
        return delta, "balanced"

    atom_delta = {
        element: count for element, count in delta.items() if element != "__charge__"
    }
    direction_values = atom_delta.values() if atom_delta else delta.values()
    has_pos = any(value > 0 for value in direction_values)
    has_neg = any(value < 0 for value in direction_values)
    if has_pos and not has_neg:
        return delta, "missing_in_products"
    if has_neg and not has_pos:
        return delta, "missing_in_reactants"
    return delta, "missing_in_both_sides"


def can_impute_with_aux(
    delta: Counter[str],
    aux_formula: Counter[str],
) -> Optional[int]:
    """Return stoichiometric count if ``delta`` equals ``c * aux_formula``."""
    if not delta:
        return 0

    ratios = []
    for elem, count in aux_formula.items():
        if elem not in delta or count == 0:
            return None
        if delta[elem] % count != 0:
            return None
        ratios.append(delta[elem] // count)

    if not ratios or len(set(ratios)) != 1:
        return None

    coefficient = ratios[0]
    if coefficient <= 0:
        return None

    remainder = signed_difference(
        delta, Counter({k: coefficient * v for k, v in aux_formula.items()})
    )
    if remainder:
        return None
    return coefficient


def impute_rsmi_rule(
    rsmi: str,
    delta: Counter[str],
    status: str,
    library: Dict[str, Counter[str]] = AUX_LIBRARY,
) -> Optional[str]:
    """Impute a small auxiliary species on the side indicated by imbalance."""
    if status == "balanced":
        return rsmi
    if status == "missing_in_both_sides":
        return None

    lhs, rhs = rsmi.split(">>")
    if status == "missing_in_products":
        target_delta = delta
        side = "products"
    elif status == "missing_in_reactants":
        target_delta = -delta
        side = "reactants"
    else:
        return None

    for aux_smiles, formula in library.items():
        coefficient = can_impute_with_aux(target_delta, formula)
        if coefficient is None:
            continue
        aux_block = ".".join([aux_smiles] * coefficient)
        if side == "products":
            rhs = f"{rhs}.{aux_block}" if rhs else aux_block
        else:
            lhs = f"{lhs}.{aux_block}" if lhs else aux_block
        return f"{lhs}>>{rhs}"
    return None


def impute_smi_mcs(
    rsmi: str,
    status: str,
    missing_smi: Optional[str] = None,
) -> Optional[str]:
    """Add an MCS-inferred missing species to the side indicated by imbalance."""
    if status == "balanced" or missing_smi is None:
        return rsmi
    if status == "missing_in_both_sides":
        return None

    lhs, rhs = rsmi.split(">>")
    if status == "missing_in_products":
        rhs = f"{rhs}.{missing_smi}" if rhs else missing_smi
    elif status == "missing_in_reactants":
        lhs = f"{lhs}.{missing_smi}" if lhs else missing_smi
    else:
        return None
    return f"{lhs}>>{rhs}"


# ---------------------------------------------------------------------------
# ITS helpers (adapted from talktorial S04)
# ---------------------------------------------------------------------------


def _get_edge_order(G: nx.Graph, u: Any, v: Any) -> float:
    ed = G.get_edge_data(u, v)
    if not ed:
        return 0.0
    val = ed.get("order", 0.0)
    if isinstance(val, (tuple, list)):
        return float(val[0]) if val else 0.0
    try:
        return float(val)
    except (TypeError, ValueError):
        return 0.0


def build_its(
    G_r: nx.Graph,
    G_p: nx.Graph,
    *,
    atom_map_key: str = "atom_map",
    node_attrs: Optional[List[str]] = None,
    include_zero_edges: bool = False,
    require_element_match: bool = True,
) -> nx.Graph:
    """
    Build an Imaginary Transition State (ITS) graph from reactant/product graphs.

    Nodes are keyed by atom-map number.  Each edge carries an ``order``
    attribute that is a ``(b_r, b_p)`` tuple of bond orders in reactant and
    product respectively.

    :param G_r: Reactant graph (nodes must have ``atom_map`` attribute).
    :param G_p: Product graph.
    :param atom_map_key: Node attribute used for alignment.
    :param node_attrs: Extra node attributes to propagate as ``(val_r, val_p)``
        tuples.  Defaults to all attributes except element and atom_map.
    :param include_zero_edges: If True, include (0,0) edges between all pairs.
    :param require_element_match: Raise if the same atom-map maps to different
        elements in reactant vs product.
    :returns: ITS graph.
    :rtype: networkx.Graph
    """
    map_r = {
        data[atom_map_key]: n
        for n, data in G_r.nodes(data=True)
        if atom_map_key in data and data[atom_map_key]
    }
    map_p = {
        data[atom_map_key]: n
        for n, data in G_p.nodes(data=True)
        if atom_map_key in data and data[atom_map_key]
    }
    all_maps = sorted(set(map_r.keys()) | set(map_p.keys()))

    if node_attrs is None:
        keys: set = set()
        for _, d in G_r.nodes(data=True):
            keys.update(d.keys())
        for _, d in G_p.nodes(data=True):
            keys.update(d.keys())
        keys.discard(atom_map_key)
        keys.discard("element")
        node_attrs = sorted(keys)

    ITS = nx.Graph()
    ITS.add_nodes_from(all_maps)

    for m in all_maps:
        r_node = map_r.get(m)
        p_node = map_p.get(m)
        r_data = G_r.nodes[r_node] if r_node is not None else {}
        p_data = G_p.nodes[p_node] if p_node is not None else {}

        elem_r = r_data.get("element")
        elem_p = p_data.get("element")
        if elem_r is not None and elem_p is not None and elem_r != elem_p:
            if require_element_match:
                raise ValueError(
                    f"Element mismatch for atom_map {m}: {elem_r!r} vs {elem_p!r}"
                )
            ITS.nodes[m]["element"] = (elem_r, elem_p)
        else:
            ITS.nodes[m]["element"] = elem_r if elem_r is not None else elem_p

        for key in node_attrs:
            ITS.nodes[m][key] = (r_data.get(key), p_data.get(key))

    for idx_i, i in enumerate(all_maps):
        for j in all_maps[idx_i + 1 :]:
            r_u, r_v = map_r.get(i), map_r.get(j)
            br = (
                _get_edge_order(G_r, r_u, r_v)
                if (r_u is not None and r_v is not None)
                else 0.0
            )
            p_u, p_v = map_p.get(i), map_p.get(j)
            bp = (
                _get_edge_order(G_p, p_u, p_v)
                if (p_u is not None and p_v is not None)
                else 0.0
            )
            if include_zero_edges or (br != 0.0 or bp != 0.0):
                ITS.add_edge(i, j, order=(br, bp))

    return ITS


def _get_reaction_center(its: nx.Graph) -> nx.Graph:
    """Extract the reaction-center subgraph (edges where bond order changes)."""
    rc_nodes: set = set()
    rc_edges = []
    for u, v, d in its.edges(data=True):
        br, bp = d.get("order", (0.0, 0.0))
        if float(br) != float(bp):
            rc_edges.append((u, v))
            rc_nodes.update([u, v])
    rc = nx.Graph()
    rc.add_nodes_from((n, its.nodes[n]) for n in rc_nodes)
    rc.add_edges_from((u, v, its.edges[u, v]) for u, v in rc_edges)
    return rc


def rsmi_to_its(rsmi: str, core: bool = False) -> nx.Graph:
    """
    Convert a reaction SMILES to an ITS (or reaction-center) graph.

    :param rsmi: Atom-mapped reaction SMILES.
    :param core: If True, return only the reaction-center subgraph (atoms and
        bonds where bond order changes).
    :returns: ITS or reaction-center NetworkX graph.
    :rtype: networkx.Graph
    """
    r, p = rsmi_to_graph(rsmi)
    its = build_its(r, p)
    if core:
        return _get_reaction_center(its)
    return its
