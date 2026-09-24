"""Pure-numpy graph neural network for protein stability (ddG) scoring.

Architecture
------------
Node features X (N x F) and symmetric normalized adjacency A_hat with
self-loops, A_hat = D^-1/2 (A + I) D^-1/2:

    Z1 = A_hat X W1 + b1          H1 = relu(Z1)
    Z2 = A_hat H1 W2 + b2         H2 = relu(Z2)
    g  = (1/N) sum_i H2[i]                    (mean-pool readout)
    h3 = relu(g U1 + d1)
    y  = h3 U2 + d2                           (scalar ddG prediction)

Loss: mean squared error over the batch, L = mean((y - t)^2).

Backpropagation (per example; A_hat symmetric so A_hat^T = A_hat):

    dL/dy = 2(y - t)
    dU2 = h3^T dL/dy                      dd2 = dL/dy
    dh3 = relu'(h3) * (dL/dy U2^T)
    dU1 = g^T dh3                         dd1 = dh3
    dg  = dh3 U1^T
    dH2 = dg / N  (broadcast to every node)
    dZ2 = relu'(Z2) * dH2
    dW2 = (A_hat H1)^T dZ2                db2 = sum_rows dZ2
    dH1 = A_hat (dZ2 W2^T)   (A_hat @ dZ2 @ W2^T)
    dZ1 = relu'(Z1) * dH1
    dW1 = (A_hat X)^T dZ1                 db1 = sum_rows dZ1

Training uses Adam. Gradient correctness is verified hermetically against
finite differences in tests/test_gnn_ddg.py.
"""
from __future__ import annotations

import numpy as np


def _relu(x: np.ndarray) -> np.ndarray:
    return np.maximum(x, 0.0)


def _relu_grad(pre: np.ndarray) -> np.ndarray:
    return (pre > 0.0).astype(np.float64)


class GNNddG:
    """Two-layer GCN with mean-pool readout and a 2-layer MLP head."""

    def __init__(self, in_dim: int, hidden: int = 32, mlp_hidden: int = 16,
                 seed: int = 0):
        rng = np.random.default_rng(seed)
        s1 = np.sqrt(2.0 / in_dim)
        s2 = np.sqrt(2.0 / hidden)
        s3 = np.sqrt(2.0 / hidden)
        s4 = np.sqrt(2.0 / mlp_hidden)
        self.params = {
            "W1": rng.normal(0.0, s1, (in_dim, hidden)),
            "b1": np.zeros(hidden),
            "W2": rng.normal(0.0, s2, (hidden, hidden)),
            "b2": np.zeros(hidden),
            "U1": rng.normal(0.0, s3, (hidden, mlp_hidden)),
            "d1": np.zeros(mlp_hidden),
            "U2": rng.normal(0.0, s4, (mlp_hidden, 1)),
            "d2": np.zeros(1),
        }
        self._adam_m = {k: np.zeros_like(v) for k, v in self.params.items()}
        self._adam_v = {k: np.zeros_like(v) for k, v in self.params.items()}
        self._adam_t = 0

    def forward(self, norm_adj: np.ndarray, feats: np.ndarray,
                cache: bool = False):
        a_x = norm_adj @ feats                    # (N, F)
        Z1 = a_x @ self.params["W1"] + self.params["b1"]
        H1 = _relu(Z1)
        a_h1 = norm_adj @ H1                      # (N, hidden)
        Z2 = a_h1 @ self.params["W2"] + self.params["b2"]
        H2 = _relu(Z2)
        g = H2.mean(axis=0)                       # (hidden,)
        h3_pre = g @ self.params["U1"] + self.params["d1"]
        h3 = _relu(h3_pre)
        y = float(h3 @ self.params["U2"][:, 0] + self.params["d2"][0])
        if not cache:
            return y
        return y, {
            "a_x": a_x, "Z1": Z1, "H1": H1, "a_h1": a_h1, "Z2": Z2,
            "H2": H2, "g": g, "h3_pre": h3_pre, "h3": h3, "y": y,
            "n": feats.shape[0],
        }

    def backward(self, c: dict, target: float) -> dict:
        """Analytic gradients of 0.5*(y - t)^2 w.r.t. all parameters.

        (0.5 factor so dL/dy = y - t; keeps the gradient-check math clean.)
        """
        dy = c["y"] - target
        grads = {}
        grads["U2"] = np.outer(c["h3"], [dy])
        grads["d2"] = np.array([dy])
        dh3 = _relu_grad(c["h3_pre"]) * (self.params["U2"][:, 0] * dy)
        grads["U1"] = np.outer(c["g"], dh3)
        grads["d1"] = dh3
        dg = dh3 @ self.params["U1"].T
        dH2 = np.repeat(dg[None, :] / c["n"], c["n"], axis=0)
        dZ2 = _relu_grad(c["Z2"]) * dH2
        grads["W2"] = c["a_h1"].T @ dZ2
        grads["b2"] = dZ2.sum(axis=0)
        dH1 = c["norm_adj_T"] @ dZ2 @ self.params["W2"].T
        dZ1 = _relu_grad(c["Z1"]) * dH1
        grads["W1"] = c["a_x"].T @ dZ1
        grads["b1"] = dZ1.sum(axis=0)
        return grads

    def predict(self, norm_adj: np.ndarray, feats: np.ndarray) -> float:
        return self.forward(norm_adj, feats, cache=False)

    def train_step(self, batch: list[dict], lr: float = 1e-3,
                   beta1: float = 0.9, beta2: float = 0.999,
                   eps: float = 1e-8) -> float:
        """One Adam step over a batch of {norm_adj, features, target}."""
        accum = {k: np.zeros_like(v) for k, v in self.params.items()}
        loss = 0.0
        for ex in batch:
            y, c = self.forward(ex["norm_adj"], ex["features"], cache=True)
            c["norm_adj_T"] = ex["norm_adj"].T
            loss += 0.5 * (y - ex["target"]) ** 2
            g = self.backward(c, ex["target"])
            for k in accum:
                accum[k] += g[k]
        m = max(len(batch), 1)
        for k in accum:
            accum[k] /= m
        loss /= m
        self._adam_t += 1
        for k, p in self.params.items():
            self._adam_m[k] = beta1 * self._adam_m[k] + (1 - beta1) * accum[k]
            self._adam_v[k] = beta2 * self._adam_v[k] + (1 - beta2) * accum[k] ** 2
            m_hat = self._adam_m[k] / (1 - beta1 ** self._adam_t)
            v_hat = self._adam_v[k] / (1 - beta2 ** self._adam_t)
            p -= lr * m_hat / (np.sqrt(v_hat) + eps)
        return loss
