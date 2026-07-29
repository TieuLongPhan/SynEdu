"""Focused regression tests for graph and chemistry invariants."""

import networkx as nx
import pytest
from rdkit import Chem
from rdkit.Chem import rdFMCS

from synedu.Utils import (
    add_wildcard,
    build_be_matrix,
    cluster_its_graphs,
    dedup_by_host_image_with_orbits,
    edge_match,
    enumerate_automorphisms,
    graph_to_mol,
    graph_to_smi,
    h_to_implicit,
    impute_rsmi_rule,
    mcs_class,
    nx_subgraph_matches,
    reaction_imbalance,
    smiles_to_graph,
    wl_hash,
)
from synedu.Utils.mcs import rdkit_mcs


def _plain_graph(edges):
    graph = nx.Graph()
    graph.add_edges_from(edges)
    nx.set_node_attributes(graph, "C", "element")
    nx.set_edge_attributes(graph, 1.0, "order")
    return graph


def test_networkx_subgraph_matching_is_monomorphic() -> None:
    """A path embeds in a triangle even though the host has an extra edge."""
    pattern = _plain_graph([(0, 1), (1, 2)])
    host = _plain_graph([(0, 1), (1, 2), (2, 0)])

    assert len(nx_subgraph_matches(host, pattern)) == 6


def test_host_image_dedup_uses_whole_automorphisms() -> None:
    """A vertex-transitive prism still has two distinct edge orbits."""
    host = nx.circular_ladder_graph(3)
    pattern = _plain_graph([(0, 1)])
    nx.set_node_attributes(host, "C", "element")
    nx.set_edge_attributes(host, 1.0, "order")
    matches = nx_subgraph_matches(host, pattern)

    representatives = dedup_by_host_image_with_orbits(
        matches,
        host_autos=enumerate_automorphisms(host),
    )

    assert len(representatives) == 2


def test_host_image_dedup_degrades_safely_for_legacy_orbits() -> None:
    matches = [{0: 2, 1: 3}, {0: 3, 1: 2}]

    with pytest.warns(DeprecationWarning, match="orbits.*ignored"):
        representatives = dedup_by_host_image_with_orbits(
            matches,
            orbits=[{2, 3}],
        )

    assert len(representatives) == 1


def test_host_image_dedup_rejects_incomplete_automorphisms_clearly() -> None:
    with pytest.raises(ValueError, match="does not cover image node"):
        dedup_by_host_image_with_orbits(
            [{0: 2, 1: 3}],
            host_autos=[{2: 3}],
        )


def test_edge_match_distinguishes_aromatic_and_single_orders() -> None:
    assert not edge_match(
        {"order": 1.5, "aromatic": False},
        {"order": 1.0, "aromatic": False},
    )


def test_wildcards_use_smallest_satisfied_valence() -> None:
    methylamine = smiles_to_graph("CN")
    dimethyl_sulfide = smiles_to_graph("CSC")

    assert all(
        data.get("element") != "*"
        for _, data in add_wildcard(methylamine, list(methylamine.nodes)).nodes(
            data=True
        )
    )
    assert all(
        data.get("element") != "*"
        for _, data in add_wildcard(
            dimethyl_sulfide, list(dimethyl_sulfide.nodes)
        ).nodes(data=True)
    )


def test_aromatic_fragment_wildcards_use_supported_bond_orders() -> None:
    benzene = smiles_to_graph("c1ccccc1")
    fragment = benzene.subgraph([0, 1, 2]).copy()
    wildcard_fragment = add_wildcard(fragment, list(fragment.nodes))

    wildcard_orders = [
        data["order"]
        for u, v, data in wildcard_fragment.edges(data=True)
        if wildcard_fragment.nodes[u].get("element") == "*"
        or wildcard_fragment.nodes[v].get("element") == "*"
    ]
    assert wildcard_orders
    assert set(wildcard_orders) <= {1.0, 2.0, 3.0}
    assert graph_to_smi(wildcard_fragment)


def test_wl_hash_and_clustering_preserve_formal_charge() -> None:
    neutral = smiles_to_graph("CCN")
    charged = smiles_to_graph("CC[NH3+]")

    assert wl_hash(neutral) != wl_hash(charged)
    assert cluster_its_graphs([neutral, charged]) == [0, 1]


def test_unabsorbable_hydrogen_is_preserved() -> None:
    hydrogen = smiles_to_graph("[H][H]")

    collapsed = h_to_implicit(hydrogen)

    assert collapsed.number_of_nodes() == 2
    assert collapsed.number_of_edges() == 1


def test_graph_to_mol_respects_zero_hydrogens_and_zero_order_edges() -> None:
    graph = nx.Graph()
    graph.add_node(1, element="C", formal_charge=0, hcount=0)
    graph.add_node(2, element="C", formal_charge=0, hcount=3)
    graph.add_node(3, element="C", formal_charge=0, hcount=3)
    graph.add_edge(1, 2, order=1.0)
    graph.add_edge(1, 3, order=1.0)

    molecule = graph_to_mol(graph, sanitize=False, use_h_count=True)
    assert molecule.GetAtomWithIdx(0).GetNoImplicit()
    assert molecule.GetAtomWithIdx(0).GetNumExplicitHs() == 0

    graph.edges[1, 2]["order"] = 0.0
    molecule = graph_to_mol(graph, sanitize=False, use_h_count=True)
    assert molecule.GetNumBonds() == 1


def test_graph_to_mol_rejects_unsupported_bond_orders() -> None:
    graph = _plain_graph([(0, 1)])
    graph.edges[0, 1]["order"] = 4.0

    with pytest.raises(ValueError, match="Unsupported bond order"):
        graph_to_mol(graph, sanitize=False)


def test_be_matrix_rejects_unknown_element_models() -> None:
    graph = nx.Graph()
    graph.add_node(0, element="Fe", formal_charge=2, hcount=0)

    with pytest.raises(ValueError, match="No valence-electron model"):
        build_be_matrix(graph)


def test_mcs_input_validation_and_small_reference_classes() -> None:
    with pytest.raises(ValueError, match="Could not parse"):
        rdkit_mcs(["not-a-smiles", "CC"])

    assert mcs_class(0, 1) == "weak/no MCS"
    assert mcs_class(1, 2) == "single-atom MCS"
    assert mcs_class(2, 3) == "large partial MCS"

    assert Chem.MolFromSmiles("CC") is not None


def test_mcs_preserves_historical_positional_comparison_arguments() -> None:
    smarts, result = rdkit_mcs(
        ["CC", "CCC"],
        3,
        True,
        rdFMCS.AtomCompare.CompareElements,
        rdFMCS.BondCompare.CompareOrder,
    )

    assert smarts
    assert not result.canceled


def test_reaction_balance_includes_net_charge() -> None:
    delta, status = reaction_imbalance("[Na+]>>[Na-]")

    assert delta == {"__charge__": 2}
    assert status == "missing_in_products"


def test_anionic_auxiliary_imputation_uses_atom_direction_and_charge() -> None:
    rsmi = "CO>>[CH3+]"
    delta, status = reaction_imbalance(rsmi)

    assert status == "missing_in_products"
    assert delta == {"H": 1, "O": 1, "__charge__": -1}
    assert impute_rsmi_rule(rsmi, delta, status) == "CO>>[CH3+].[OH-]"
