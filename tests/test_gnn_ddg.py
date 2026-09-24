"""Hermetic GNN tests: analytic gradients checked against finite
differences, and end-to-end learning verified by overfitting a tiny
synthetic dataset. All data here is synthetic test data."""
import numpy as np

from src.contact_graph import build_graph
from src.gnn_ddg import GNNddG
from src.pdb_fetch import Atom


def _ca(serial, res, chain, seq, x, y, z):
    return Atom(serial=serial, name="CA", res_name=res, chain=chain,
                res_seq=seq, x=x, y=y, z=z, element="C")


def _tiny_graph(seed):
    rng = np.random.default_rng(seed)
    atoms = []
    pos = rng.normal(0, 3.0, (6, 3))
    res = ["GLY", "ALA", "VAL", "LEU", "SER", "THR"]
    for i, (x, y, z) in enumerate(pos):
        atoms.append(_ca(i + 1, res[i], "A", i + 1, x, y, z))
    return build_graph(atoms, chain="A", cutoff=8.0)


def _example(seed, target):
    g = _tiny_graph(seed)
    return {"norm_adj": g["norm_adj"], "features": g["features"],
            "target": target}


def test_analytic_gradients_match_finite_differences():
    net = GNNddG(in_dim=21, hidden=8, mlp_hidden=4, seed=42)
    ex = _example(seed=7, target=0.35)
    y, cache = net.forward(ex["norm_adj"], ex["features"], cache=True)
    cache["norm_adj_T"] = ex["norm_adj"].T
    grads = net.backward(cache, ex["target"])

    eps = 1e-5
    for key in ["W1", "W2", "U1", "U2"]:
        p = net.params[key]
        flat_idx = np.unravel_index(np.argmax(np.abs(grads[key])), p.shape)
        orig = p[flat_idx]
        p[flat_idx] = orig + eps
        yp = net.forward(ex["norm_adj"], ex["features"])
        lp = 0.5 * (yp - ex["target"]) ** 2
        p[flat_idx] = orig - eps
        ym = net.forward(ex["norm_adj"], ex["features"])
        lm = 0.5 * (ym - ex["target"]) ** 2
        p[flat_idx] = orig
        num = (lp - lm) / (2 * eps)
        ana = grads[key][flat_idx]
        assert abs(num - ana) < 1e-6 * max(1.0, abs(ana)), \
            f"{key}: analytic {ana} vs numeric {num}"


def test_overfits_tiny_synthetic_batch():
    # 4 synthetic graphs with arbitrary targets; a working optimizer must
    # drive the loss far below its initial value.
    batch = [_example(seed=s, target=t)
             for s, t in [(1, -1.0), (2, 0.5), (3, 1.5), (4, -0.25)]]
    net = GNNddG(in_dim=21, hidden=16, mlp_hidden=8, seed=0)
    initial = net.train_step(batch, lr=5e-3)
    loss = initial
    for _ in range(400):
        loss = net.train_step(batch, lr=5e-3)
    assert loss < 0.05 * initial, \
        f"loss did not collapse: {initial:.4f} -> {loss:.4f}"
    for ex in batch:
        pred = net.predict(ex["norm_adj"], ex["features"])
        assert abs(pred - ex["target"]) < 0.3


def test_forward_deterministic_and_scalar():
    net = GNNddG(in_dim=21, hidden=8, mlp_hidden=4, seed=3)
    ex = _example(seed=5, target=0.0)
    y1 = net.predict(ex["norm_adj"], ex["features"])
    y2 = net.predict(ex["norm_adj"], ex["features"])
    assert isinstance(y1, float) and y1 == y2
