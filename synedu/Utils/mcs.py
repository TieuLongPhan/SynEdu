"""Maximum-common-substructure helpers used in S03 and later talktorials."""

from __future__ import annotations

from itertools import combinations
from typing import Any, Dict, Iterable, Set, Tuple

import networkx as nx
from networkx.algorithms import isomorphism as iso
from rdkit import Chem
from rdkit.Chem import rdFMCS


def rdkit_mcs(
    mols: Iterable[Chem.Mol | str],
    timeout: int = 10,
    ringMatchesRingOnly: bool = True,
    atomCompare=rdFMCS.AtomCompare.CompareElements,
    bondCompare=rdFMCS.BondCompare.CompareOrder,
):
    """Return the MCS SMARTS string and the RDKit ``FindMCS`` result object."""
    mol_objs = [Chem.MolFromSmiles(m) if isinstance(m, str) else m for m in mols]
    result = rdFMCS.FindMCS(
        mol_objs,
        atomCompare=atomCompare,
        bondCompare=bondCompare,
        ringMatchesRingOnly=ringMatchesRingOnly,
        timeout=timeout,
    )
    return result.smartsString, result


def mcs_size(mol_a: Chem.Mol, mol_b: Chem.Mol, timeout: int = 3) -> int:
    """Return the atom count of the RDKit MCS between two molecules."""
    try:
        smarts, result = rdkit_mcs([mol_a, mol_b], timeout=timeout)
        if result.canceled or not smarts:
            return 0
        query = Chem.MolFromSmarts(smarts)
        return query.GetNumAtoms() if query is not None else 0
    except Exception:
        return 0


def mcs_class(mcs_atoms: int, ref_atoms: int) -> str:
    """Classify an MCS size relative to a reference molecule atom count."""
    if mcs_atoms == ref_atoms:
        return "full reference MCS"
    if mcs_atoms == ref_atoms - 1:
        return "large partial MCS"
    if mcs_atoms == 1:
        return "single-atom MCS"
    return "weak/no MCS"


def nx_mcs(
    G: nx.Graph,
    H: nx.Graph,
    *,
    node_match,
    edge_match,
) -> tuple[nx.Graph, Dict[Any, Any]]:
    """
    Return one maximum common induced subgraph of ``G`` and a mapping to ``H``.

    This brute-force helper is intentionally small and transparent for teaching.
    It is exponential and should be used only for small graphs.
    """
    nodes = list(G.nodes())
    for size in range(len(nodes), 0, -1):
        for subset in combinations(nodes, size):
            subgraph = G.subgraph(subset).copy()
            matcher = iso.GraphMatcher(
                H, subgraph, node_match=node_match, edge_match=edge_match
            )
            if matcher.subgraph_is_isomorphic():
                h_to_sg = next(matcher.subgraph_isomorphisms_iter())
                return subgraph, {sg: h for h, sg in h_to_sg.items()}
    return nx.Graph(), {}


def map_edges(
    subgraph: nx.Graph, sg_to_h: Dict[Any, Any], H: nx.Graph
) -> Set[Tuple[Any, Any]]:
    """Map edges from a subgraph into a host graph using a node mapping."""
    out: Set[Tuple[Any, Any]] = set()
    for u, v in subgraph.edges():
        hu, hv = sg_to_h[u], sg_to_h[v]
        if H.has_edge(hu, hv):
            out.add((hu, hv))
    return out
