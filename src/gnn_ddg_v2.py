"""GNN ddG scorer v2 - center-aware readout.

v1 weakness (preserved in results/model_v1_bias_analysis.md): mean-pooling
over ~256 nodes dilutes the single-node mutation signal, letting the head
shortcut on residue identity (aromatic bias). v2 readout concatenates the
mean-pooled embedding with the mutated site's own node vector:

    Z1 = A_hat X W1 + b1          H1 = relu(Z1)
    Z2 = A_hat H1 W2 + b2         H2 = relu(Z2)
    g  = [ mean_i H2[i] , H2[c] ]               (c = mutation site)
    h3 = relu(g U1 + d1)
    y  = h3 U2 + d2

ddG = f(mutant graph) - f(wild-type graph) through the SAME network, so the
antisymmetry property of v1 is preserved exactly: swapping mutant and
wild-type negates the prediction by construction.

Backpropagation adds only the readout split vs v1:
    dH2 = broadcast(dg_mean / N) + one_hot(c) outer dg_ctr
"""
from __future__ import annotations

import numpy as np

from src.gnn_ddg import _relu, _relu_grad


class GNNddGv2:
    """Two-layer GCN with [mean-pool | center-node] readout + MLP head."""

    def __init__(self, in_dim: int, hidden: int = 64, mlp_hidden: int = 32,
                 seed: int = 0):
        rng = np.random.default_rng(seed)
        self.hidden = hidden
        self.params = {
            "W1": rng.normal(0.0, np.sqrt(2.0 / in_dim), (in_dim, hidden)),
            "b1": np.zeros(hidden),
            "W2": rng.normal(0.0, np.sqrt(2.0 / hidden), (hidden, hidden)),
            "b2": np.zeros(hidden),
            "U1": rng.normal(0.0, np.sqrt(2.0 / (2 * hidden)),
                             (2 * hidden, mlp_hidden)),
            "d1": np.zeros(mlp_hidden),
            "U2": rng.normal(0.0, np.sqrt(2.0 / mlp_hidden), (mlp_hidden, 1)),
            "d2": np.zeros(1),
        }
        self._adam_m = {k: np.zeros_like(v) for k, v in self.params.items()}
        self._adam_v = {k: np.zeros_like(v) for k, v in self.params.items()}
        self._adam_t = 0

    def forward(self, norm_adj, feats, center: int, cache: bool = False):
        a_x = norm_adj @ feats
        Z1 = a_x @ self.params["W1"] + self.params["b1"]
        H1 = _relu(Z1)
        a_h1 = norm_adj @ H1
        Z2 = a_h1 @ self.params["W2"] + self.params["b2"]
        H2 = _relu(Z2)
        g = np.concatenate([H2.mean(axis=0), H2[center]])
        h3_pre = g @ self.params["U1"] + self.params["d1"]
        h3 = _relu(h3_pre)
        y = float(h3 @ self.params["U2"][:, 0] + self.params["d2"][0])
        if not cache:
            return y
        return y, {"a_x": a_x, "Z1": Z1, "H1": H1, "a_h1": a_h1, "Z2": Z2,
                   "H2": H2, "g": g, "h3_pre": h3_pre, "h3": h3, "y": y,
                   "n": feats.shape[0], "center": center}

    def backward(self, c: dict, target: float) -> dict:
        """Analytic gradients of 0.5*(y - t)^2 (dL/dy = y - t)."""
        dy = c["y"] - target
        grads = {"U2": np.outer(c["h3"], [dy]), "d2": np.array([dy])}
        dh3 = _relu_grad(c["h3_pre"]) * (self.params["U2"][:, 0] * dy)
        grads["U1"] = np.outer(c["g"], dh3)
        grads["d1"] = dh3
        dg = dh3 @ self.params["U1"].T
        dg_mean, dg_ctr = dg[:self.hidden], dg[self.hidden:]
        dH2 = np.repeat(dg_mean[None, :] / c["n"], c["n"], axis=0)
        dH2[c["center"]] += dg_ctr
        dZ2 = _relu_grad(c["Z2"]) * dH2
        grads["W2"] = c["a_h1"].T @ dZ2
        grads["b2"] = dZ2.sum(axis=0)
        dH1 = c["norm_adj_T"] @ dZ2 @ self.params["W2"].T
        dZ1 = _relu_grad(c["Z1"]) * dH1
        grads["W1"] = c["a_x"].T @ dZ1
        grads["b1"] = dZ1.sum(axis=0)
        return grads

    def predict(self, norm_adj, feats, center: int) -> float:
        return self.forward(norm_adj, feats, center, cache=False)
