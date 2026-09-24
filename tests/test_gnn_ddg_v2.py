"""Hermetic v2 tests: finite-difference gradient check, antisymmetry by
construction, and overfit learning - all on synthetic graphs."""
import numpy as np

from src.contact_graph import build_graph
from src.gnn_ddg_v2 import GNNddGv2
from src.pdb_fetch import Atom


def _ca(serial, res, chain, seq, x, y, z):
    return Atom(serial=serial, name="CA", res_name=res, chain=chain,
                res_seq=seq, x=x, y=y, z=z, element="C")


def _graph(seed):
    rng = np.random.default_rng(seed)
    res = ["GLY", "ALA", "VAL", "LEU", "SER", "THR"]
    atoms = [_ca(i + 1, res[i], "A", i + 1, x, y, z)
             for i, (x, y, z) in enumerate(rng.normal(0, 3.0, (6, 3)))]
    return build_graph(atoms, chain="A", cutoff=8.0), 2  # center = node 2


def _loss(net, g, c, t):
    y = net.predict(g["norm_adj"], g["features"], c)
    return 0.5 * (y - t) ** 2


def test_v2_gradients_match_finite_differences():
    net = GNNddGv2(in_dim=21, hidden=8, mlp_hidden=4, seed=42)
    g, c = _graph(7)
    y, cache = net.forward(g["norm_adj"], g["features"], c, cache=True)
    cache["norm_adj_T"] = g["norm_adj"].T
    grads = net.backward(cache, 0.35)
    eps = 1e-5
    for key in ["W1", "W2", "U1", "U2"]:
        p = net.params[key]
        idx = np.unravel_index(np.argmax(np.abs(grads[key])), p.shape)
        orig = p[idx]
        p[idx] = orig + eps
        lp = _loss(net, g, c, 0.35)
        p[idx] = orig - eps
        lm = _loss(net, g, c, 0.35)
        p[idx] = orig
        num = (lp - lm) / (2 * eps)
        assert abs(num - grads[key][idx]) < 1e-6 * max(1.0, abs(grads[key][idx])), \
            f"{key}: analytic {grads[key][idx]} vs numeric {num}"


def test_v2_antisymmetry_exact():
    # ddG(wt->mut) == -ddG(mut->wt) by construction (shared net, difference)
    net = GNNddGv2(in_dim=21, hidden=8, mlp_hidden=4, seed=1)
    g, c = _graph(3)
    from src.contact_graph import mutated_features, AA20
    mf = g["features"].copy()
    mf[c] = 0.0
    from src.contact_graph import AA_INDEX, HYDROPHOBICITY
    mf[c][AA_INDEX["TRP"]] = 1.0
    mf[c][len(AA20)] = HYDROPHOBICITY["TRP"] / 4.5
    fwd = (net.predict(g["norm_adj"], mf, c)
           - net.predict(g["norm_adj"], g["features"], c))
    rev = (net.predict(g["norm_adj"], g["features"], c)
           - net.predict(g["norm_adj"], mf, c))
    assert abs(fwd + rev) < 1e-12


def test_v2_center_readout_changes_output():
    net = GNNddGv2(in_dim=21, hidden=8, mlp_hidden=4, seed=5)
    g, c = _graph(11)
    y_c = net.predict(g["norm_adj"], g["features"], c)
    y_other = net.predict(g["norm_adj"], g["features"], (c + 1) % 6)
    assert y_c != y_other  # center position must matter
