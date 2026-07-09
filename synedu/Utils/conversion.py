from rdkit import Chem
import networkx as nx
from typing import Any, Dict, Optional


def mol_to_graph(
    mol: Chem.Mol,
    drop_non_aam: bool = False,
    use_index_as_atom_map: bool = False,
) -> nx.Graph:
    """
    Convert an RDKit Mol to a lightweight heavy-atom NetworkX graph.

    :param mol: RDKit molecule.
    :type mol: rdkit.Chem.Mol
    :param drop_non_aam: If True, skip atoms that carry no atom-map number.
    :type drop_non_aam: bool
    :param use_index_as_atom_map: If True, use atom index as node ID and set
        ``atom_map`` attribute to the index. If False, use atom-map number
        when present, falling back to a unique atom index.
    :type use_index_as_atom_map: bool
    :returns: Graph with element/formal_charge/aromatic/hcount/atom_map node
        attributes and order/aromatic edge attributes.
    :rtype: networkx.Graph
    """
    G: nx.Graph = nx.Graph()
    local: Dict[int, Any] = {}
    seen_aam: Dict[int, int] = {}

    # Existing atom-map numbers, used to avoid collisions with fallback IDs
    used_aam = {
        atom.GetAtomMapNum() for atom in mol.GetAtoms() if atom.GetAtomMapNum() > 0
    }

    next_fallback_id = 1

    def get_unique_fallback_id() -> int:
        nonlocal next_fallback_id
        while next_fallback_id in used_aam:
            next_fallback_id += 1
        out = next_fallback_id
        used_aam.add(out)
        next_fallback_id += 1
        return out

    for atom in mol.GetAtoms():
        atom_idx = atom.GetIdx()
        am = int(atom.GetAtomMapNum())

        if drop_non_aam and am == 0:
            continue

        if use_index_as_atom_map:
            node_id = atom_idx + 1
            node_am = atom_idx + 1
        elif am > 0:
            if am in seen_aam:
                raise ValueError(
                    "Duplicate atom-map number "
                    f"{am} on atoms {seen_aam[am]} and {atom_idx}"
                )
            seen_aam[am] = atom_idx
            node_id = am
            node_am = am
        else:
            node_id = get_unique_fallback_id()
            node_am = 0

        local[atom_idx] = node_id

        G.add_node(
            node_id,
            element=atom.GetSymbol(),
            formal_charge=int(atom.GetFormalCharge()),
            aromatic=bool(atom.GetIsAromatic()),
            hcount=int(atom.GetTotalNumHs()),
            atom_map=node_am,
            rdkit_idx=atom_idx,
        )

    for bond in mol.GetBonds():
        u = local.get(bond.GetBeginAtomIdx())
        v = local.get(bond.GetEndAtomIdx())

        if u is None or v is None:
            continue

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
    use_atom_map: bool = True,
    preserve_aromatic: bool = True,
) -> Chem.Mol:
    """
    Reconstruct an RDKit Mol from a NetworkX graph.

    Supports:
    - element
    - formal_charge
    - atom_map
    - bond order
    - aromatic bonds
    - optional explicit H count
    """
    rw = Chem.RWMol()
    node_to_idx: Dict[Any, int] = {}

    for node, data in G.nodes(data=True):
        element = data.get("element", "C")
        charge = int(data.get("formal_charge", 0))

        atom = Chem.Atom(element)
        atom.SetFormalCharge(charge)

        if use_atom_map:
            am = int(data.get("atom_map", 0) or 0)
            if am:
                atom.SetAtomMapNum(am)

        atom.SetIsAromatic(False)

        idx = rw.AddAtom(atom)
        node_to_idx[node] = idx

    for u, v, data in G.edges(data=True):
        i = node_to_idx[u]
        j = node_to_idx[v]

        if preserve_aromatic and bool(data.get("aromatic", False)):
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

    if use_h_count:
        for node, data in G.nodes(data=True):
            try:
                n_h = int(data.get("hcount", 0))
            except Exception:
                continue

            if n_h <= 0:
                continue

            atom = rw.GetAtomWithIdx(node_to_idx[node])
            atom.SetNoImplicit(True)
            atom.SetNumExplicitHs(n_h)

    mol = rw.GetMol()

    if sanitize:
        Chem.SanitizeMol(
            mol,
            #     sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
            #     ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
            #     ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE,
        )
        Chem.SetAromaticity(mol)

    return mol


def smiles_to_graph(
    smiles: str,
    drop_non_aam: bool = False,
    use_index_as_atom_map: bool = False,
    sanitize: bool = True,
) -> nx.Graph:
    """
    Convert a SMILES string, including multi-component SMILES, to a NetworkX graph.

    This function delegates graph construction to :func:`mol_to_graph`.

    :param smiles: SMILES string. Atoms may carry atom-map numbers, e.g. ``[CH3:1]``.
        Multi-component SMILES such as ``CCCl.[OH-]`` are supported.
    :type smiles: str
    :param drop_non_aam: If True, skip atoms that carry no atom-map number.
    :type drop_non_aam: bool
    :param use_index_as_atom_map: If True, ignore existing atom-map numbers and use
        RDKit atom index + 1 as both node ID and ``atom_map``.
    :type use_index_as_atom_map: bool
    :param sanitize: If True, sanitize the RDKit molecule before conversion.
    :type sanitize: bool
    :returns: NetworkX graph with node attributes ``element``, ``formal_charge``,
        ``aromatic``, ``hcount``, and ``atom_map``; and edge attributes ``order``
        and ``aromatic``.
    :rtype: networkx.Graph
    :raises ValueError: If RDKit cannot parse the SMILES string.
    """
    mol = Chem.MolFromSmiles(smiles, sanitize=False)

    if mol is None:
        raise ValueError(f"Could not parse SMILES: {smiles!r}")

    if sanitize:
        Chem.SanitizeMol(mol)

    return mol_to_graph(
        mol,
        drop_non_aam=drop_non_aam,
        use_index_as_atom_map=use_index_as_atom_map,
    )


def _graph_component_to_mol(subgraph: nx.Graph) -> Optional[Chem.Mol]:
    """Convert a connected NetworkX subgraph to an RDKit Mol (with atom maps)."""
    rw = Chem.RWMol()
    node_to_idx: Dict[Any, int] = {}
    for node, data in subgraph.nodes(data=True):
        atom = Chem.Atom(data.get("element", "C"))
        atom.SetFormalCharge(int(data.get("formal_charge", 0)))
        am = int(data.get("atom_map", 0))
        if am:
            atom.SetAtomMapNum(am)
        idx = rw.AddAtom(atom)
        node_to_idx[node] = idx
    for u, v, data in subgraph.edges(data=True):
        order = int(round(abs(float(data.get("order", 1.0)))))
        btype = {
            1: Chem.BondType.SINGLE,
            2: Chem.BondType.DOUBLE,
            3: Chem.BondType.TRIPLE,
        }.get(order, Chem.BondType.SINGLE)
        rw.AddBond(node_to_idx[u], node_to_idx[v], btype)
    mol = rw.GetMol()
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        pass
    return mol


def graph_to_smi(G: nx.Graph) -> str:
    """
    Convert a NetworkX molecular graph to a canonical SMILES string.

    Handles disconnected graphs (multi-component molecules) by joining
    components with ``"."``.

    :param G: Molecular graph with ``element`` node attribute.
    :type G: networkx.Graph
    :returns: Canonical SMILES string.
    :rtype: str
    """
    mol = graph_to_mol(
        G, sanitize=True, use_atom_map=True, preserve_aromatic=True, use_h_count=True
    )
    if mol is None:
        return ""
    return Chem.MolToSmiles(mol)
