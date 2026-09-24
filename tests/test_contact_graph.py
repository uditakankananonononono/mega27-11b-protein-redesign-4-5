"""Hermetic contact-graph tests on synthetic coordinates (parser-test style
fixtures; not experimental data)."""
import numpy as np

from src.contact_graph import (AA20, build_graph, mutated_features, one_hot,
                               residue_index)
from src.pdb_fetch import Atom


def _ca(serial, res, chain, seq, x, y, z):
    return Atom(serial=serial, name="CA", res_name=res, chain=chain,
                res_seq=seq, x=x, y=y, z=z, element="C")


# 4 residues on a line, 4 A apart -> cutoff 5 gives a path graph 0-1-2-3.
ATOMS = [
    _ca(1, "GLY", "A", 1, 0.0, 0.0, 0.0),
    _ca(2, "ALA", "A", 2, 4.0, 0.0, 0.0),
    _ca(3, "VAL", "A", 3, 8.0, 0.0, 0.0),
    _ca(4, "LEU", "A", 4, 12.0, 0.0, 0.0),
]


def test_adjacency_path_graph():
    g = build_graph(ATOMS, chain="A", cutoff=5.0)
    expected = np.array([
        [0, 1, 0, 0],
        [1, 0, 1, 0],
        [0, 1, 0, 1],
        [0, 0, 1, 0],
    ], dtype=np.float64)
    assert np.array_equal(g["adj"], expected)
    assert g["edge_index"].shape[1] == 6  # 3 undirected edges, both directions
    assert np.allclose(g["edge_dist"], 4.0)


def test_cutoff_excludes_distant_pairs():
    g = build_graph(ATOMS, chain="A", cutoff=3.9)
    assert g["adj"].sum() == 0.0


def test_normalized_adjacency_math():
    g = build_graph(ATOMS, chain="A", cutoff=5.0)
    a_hat = g["adj"] + np.eye(4)
    deg = a_hat.sum(axis=1)  # [2, 3, 3, 2]
    d = np.power(deg, -0.5)
    expected = d[:, None] * a_hat * d[None, :]
    assert np.allclose(g["norm_adj"], expected)
    # spot-check one entry by hand: (0,1) = 1/sqrt(2*3)
    assert abs(g["norm_adj"][0, 1] - 1 / np.sqrt(6)) < 1e-12


def test_feature_shapes_and_one_hot():
    g = build_graph(ATOMS, chain="A", cutoff=5.0)
    assert g["features"].shape == (4, len(AA20) + 1)
    assert one_hot("GLY").sum() == 1.0
    assert one_hot("UNK").sum() == 0.0
    assert g["features"][0, AA20.index("GLY")] == 1.0


def test_mutation_replaces_only_target_row():
    g = build_graph(ATOMS, chain="A", cutoff=5.0)
    idx = residue_index(g, "A", 2)
    mf = mutated_features(g, "A", 2, "TRP")
    assert mf[idx, AA20.index("TRP")] == 1.0
    assert mf[idx, AA20.index("ALA")] == 0.0
    keep = [i for i in range(4) if i != idx]
    assert np.array_equal(mf[keep], g["features"][keep])
    assert np.array_equal(g["norm_adj"], g["norm_adj"])  # graph untouched


def test_local_subgraph_selection_and_determinism():
    from src.contact_graph import local_subgraph
    # 40 residues on a line, 5 A apart
    atoms = [
        Atom(serial=i + 1, name="CA", res_name="ALA", chain="A",
             res_seq=i + 1, x=5.0 * i, y=0.0, z=0.0, element="C")
        for i in range(40)
    ]
    g = local_subgraph(atoms, "A", 20, min_nodes=8, max_nodes=12,
                       start_radius=6.0, radius_step=5.0, max_radius=60.0)
    # center always present, cap respected, centered on residue 20
    seqs = [r[2] for r in g["residues"]]
    assert 20 in seqs and len(seqs) <= 12
    assert abs(seqs.index(20) - (len(seqs) - 1) / 2) <= 1  # roughly centered
    assert g["residues"][g["center_local_idx"]][2] == 20
    assert g["parent_chain_size"] == 40
    # determinism
    g2 = local_subgraph(atoms, "A", 20, min_nodes=8, max_nodes=12,
                        start_radius=6.0, radius_step=5.0, max_radius=60.0)
    assert [r[2] for r in g2["residues"]] == seqs
    assert np.array_equal(g["norm_adj"], g2["norm_adj"])
    # residue_index still locates the center in the subgraph
    assert residue_index(g, "A", 20) == g["center_local_idx"]


def test_local_subgraph_missing_center_raises():
    from src.contact_graph import local_subgraph
    atoms = [Atom(serial=1, name="CA", res_name="ALA", chain="A",
                  res_seq=1, x=0.0, y=0.0, z=0.0, element="C")]
    try:
        local_subgraph(atoms, "A", 99)
        assert False, "should have raised"
    except KeyError:
        pass
