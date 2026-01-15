from rdkit import Chem
import networkx as nx
from typing import Dict, Any


def mol_to_graph(mol: Chem.Mol) -> nx.Graph:
    """
    Convert an RDKit Mol to a lightweight heavy-atom NetworkX graph.

    :param mol: Sanitized RDKit molecule.
    :type mol: rdkit.Chem.Mol
    :returns: Graph keyed by atom index.
    :rtype: networkx.Graph
    """
    G = nx.Graph()
    for atom in mol.GetAtoms():
        node_id = atom.GetIdx()
        attrs: Dict[str, Any] = {
            "element": atom.GetSymbol(),
            "formal_charge": int(atom.GetFormalCharge()),
            "aromatic": bool(atom.GetIsAromatic()),
            "hcount": int(atom.GetTotalNumHs()),
        }
        G.add_node(node_id, **attrs)

    for bond in mol.GetBonds():
        u = bond.GetBeginAtomIdx()
        v = bond.GetEndAtomIdx()
        G.add_edge(
            u,
            v,
            order=float(bond.GetBondTypeAsDouble()),
            aromatic=bool(bond.GetIsAromatic()),
        )
    return G


def graph_to_mol(
    G: nx.Graph,
    sanitize: bool = True,
    use_h_count: bool = False,
) -> Chem.Mol:
    """
    Reconstruct an RDKit Mol from a lightweight NetworkX graph.

    :param G: Graph keyed by atom index.
    :type G: networkx.Graph
    :param sanitize: If True, sanitize molecule.
    :type sanitize: bool
    :param use_h_count: If True, add explicit H atoms according to node ``hcount`` (default: False).
    :type use_h_count: bool
    :returns: Reconstructed RDKit molecule.
    :rtype: rdkit.Chem.Mol
    """
    rw = Chem.RWMol()
    node_to_idx: Dict[Any, int] = {}

    # 1) add heavy atoms (defer aromatic perception)
    for node, data in G.nodes(data=True):
        element = data.get("element", "C")
        charge = int(data.get("formal_charge", 0))
        atom = Chem.Atom(element)
        atom.SetFormalCharge(charge)
        atom.SetIsAromatic(False)
        idx = rw.AddAtom(atom)
        node_to_idx[node] = idx

    # 2) add heavy-heavy bonds
    for u, v, data in G.edges(data=True):
        i = node_to_idx[u]
        j = node_to_idx[v]
        if bool(data.get("aromatic", False)):
            btype = Chem.BondType.AROMATIC
        else:
            try:
                order = int(round(abs(float(data.get("order", 1.0)))))
            except Exception:
                order = 1
            btype = {
                1: Chem.BondType.SINGLE,
                2: Chem.BondType.DOUBLE,
                3: Chem.BondType.TRIPLE,
            }.get(order, Chem.BondType.SINGLE)
        rw.AddBond(i, j, btype)

    # 3) optionally add explicit H atoms
    if use_h_count:
        for node, data in list(G.nodes(data=True)):
            try:
                n_h = int(data.get("hcount", 0))
            except Exception:
                continue
            if n_h <= 0:
                continue
            heavy_idx = node_to_idx[node]
            heavy_atom = rw.GetAtomWithIdx(heavy_idx)
            heavy_atom.SetNoImplicit(True)
            heavy_atom.SetNumExplicitHs(int(n_h))

    mol = rw.GetMol()

    if sanitize:
        Chem.SanitizeMol(
            mol,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
            ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
            ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE,
        )
        Chem.SetAromaticity(mol)

    return mol
