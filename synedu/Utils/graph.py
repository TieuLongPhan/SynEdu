from networkx.algorithms import isomorphism as iso
import networkx as nx


def node_match(n1, n2):
    return (
        n1.get("element") == n2.get("element")
        and int(n1.get("formal_charge", 0)) == int(n2.get("formal_charge", 0))
        and bool(n1.get("aromatic", False)) == bool(n2.get("aromatic", False))
    )


def edge_match(e1, e2):
    return int(e1.get("order", 1)) == int(e2.get("order", 1)) and bool(
        e1.get("aromatic", False)
    ) == bool(e2.get("aromatic", False))


def enumerate_automorphisms(G: nx.Graph):
    GM_self = iso.GraphMatcher(G, G, node_match=node_match, edge_match=edge_match)
    return list(GM_self.isomorphisms_iter())


def compute_orbits_from_automorphisms(G: nx.Graph, automorphisms=None):
    # Copied from S01 but kept small here for self-containment.
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
