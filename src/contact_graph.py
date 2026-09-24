"""Residue contact-graph construction from real PDB structures.

Nodes are residues (C-alpha representatives), edges connect residues whose
C-alpha atoms are within a cutoff distance. Produces the normalized adjacency
D^-1/2 (A + I) D^-1/2 used by the GNN, plus node feature matrices.
"""
from __future__ import annotations

import numpy as np

from src.pdb_fetch import Atom, ca_trace

AA20 = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
]
AA_INDEX = {aa: i for i, aa in enumerate(AA20)}
# Kyte-Doolittle hydrophobicity, standard scale.
HYDROPHOBICITY = {
    "ALA": 1.8, "ARG": -4.5, "ASN": -3.5, "ASP": -3.5, "CYS": 2.5,
    "GLN": -3.5, "GLU": -3.5, "GLY": -0.4, "HIS": -3.2, "ILE": 4.5,
    "LEU": 3.8, "LYS": -3.9, "MET": 1.9, "PHE": 2.8, "PRO": -1.6,
    "SER": -0.8, "THR": -0.7, "TRP": -0.9, "TYR": -1.3, "VAL": 4.2,
}


def one_hot(res_name: str) -> np.ndarray:
    """20-dim one-hot for standard amino acids; unknown -> all zeros."""
    v = np.zeros(len(AA20), dtype=np.float64)
    idx = AA_INDEX.get(res_name)
    if idx is not None:
        v[idx] = 1.0
    return v


def build_graph(atoms: list[Atom], chain: str | None = None,
                cutoff: float = 10.0) -> dict:
    """Build the residue contact graph for one chain (or all chains).

    Returns dict with keys:
      residues  list of (res_name, chain, res_seq) per node
      coords    (N, 3) C-alpha coordinates
      edge_index (2, E) directed edges, both directions stored
      edge_dist (E,)   C-alpha distances per edge
      adj       (N, N) binary adjacency without self-loops
      norm_adj  (N, N) symmetric normalized adjacency with self-loops
      features  (N, F) node features: one-hot(20) + hydrophobicity
    """
    trace = ca_trace(atoms, chain=chain)
    n = len(trace)
    if n == 0:
        raise ValueError("no C-alpha atoms found (empty chain?)")
    coords = np.array([[a.x, a.y, a.z] for a in trace], dtype=np.float64)
    residues = [(a.res_name, a.chain, a.res_seq) for a in trace]

    return assemble_graph(residues, coords, cutoff)


def assemble_graph(residues, coords, cutoff: float = 10.0) -> dict:
    """Dense contact graph from a residue list and their coordinates."""
    n = len(residues)
    diff = coords[:, None, :] - coords[None, :, :]
    dist = np.sqrt((diff ** 2).sum(axis=-1))
    adj = ((dist < cutoff) & (dist > 0.0)).astype(np.float64)

    ei = np.array(np.nonzero(adj), dtype=np.int64)
    ed = dist[adj > 0]

    a_hat = adj + np.eye(n)
    deg = a_hat.sum(axis=1)
    d_inv_sqrt = np.power(deg, -0.5)
    norm_adj = d_inv_sqrt[:, None] * a_hat * d_inv_sqrt[None, :]

    hp = np.array([HYDROPHOBICITY.get(r[0], 0.0) for r in residues],
                  dtype=np.float64).reshape(-1, 1) / 4.5  # scale to ~[-1, 1]
    feats = np.concatenate([np.stack([one_hot(r[0]) for r in residues]), hp],
                           axis=1)
    return {
        "residues": residues,
        "coords": coords,
        "edge_index": ei,
        "edge_dist": ed,
        "adj": adj,
        "norm_adj": norm_adj,
        "features": feats,
    }


def local_subgraph(atoms: list[Atom], chain: str, center_seq: int,
                   min_nodes: int = 64, max_nodes: int = 256,
                   start_radius: float = 8.0, radius_step: float = 2.0,
                   max_radius: float = 30.0, cutoff: float = 10.0) -> dict:
    """Contact graph of the local neighborhood around one residue.

    ddG effects are local: residues are selected by C-alpha distance to the
    mutation site - the radius grows from `start_radius` until at least
    `min_nodes` are inside (or `max_radius`), and the set is capped at the
    `max_nodes` nearest residues. Deterministic: the same inputs always
    select the same nodes, and the center residue is always included.
    This keeps memory bounded for very large chains (a dense N x N graph on
    a 5888-residue chain needs hundreds of MB; a 256-node cap needs ~1 MB).
    """
    trace = ca_trace(atoms, chain=chain)
    coords = np.array([[a.x, a.y, a.z] for a in trace], dtype=np.float64)
    center = None
    for i, a in enumerate(trace):
        if a.res_seq == center_seq:
            center = i
            break
    if center is None:
        raise KeyError(f"residue {chain}:{center_seq} not in chain")
    d = np.sqrt(((coords - coords[center]) ** 2).sum(axis=1))
    radius = start_radius
    while (d <= radius).sum() < min_nodes and radius < max_radius:
        radius += radius_step
    sel = np.nonzero(d <= radius)[0]
    if len(sel) > max_nodes:
        sel = np.argsort(d, kind="stable")[:max_nodes]
    sel = np.sort(sel)  # canonical order: sequence order
    residues = [(a.res_name, a.chain, a.res_seq)
                for a in (trace[i] for i in sel)]
    g = assemble_graph(residues, coords[sel], cutoff)
    g["center_local_idx"] = int(np.nonzero(sel == center)[0][0])
    g["subgraph_radius"] = float(radius)
    g["parent_chain_size"] = len(trace)
    return g


def residue_index(graph: dict, chain: str, res_seq: int) -> int:
    """Node index of a residue by chain and sequence number."""
    for i, (_, ch, sq) in enumerate(graph["residues"]):
        if ch == chain and sq == res_seq:
            return i
    raise KeyError(f"residue {chain}:{res_seq} not in graph")


def mutated_features(graph: dict, chain: str, res_seq: int,
                     new_res: str) -> np.ndarray:
    """Feature matrix with one residue's identity features replaced.

    WT coordinates are kept (standard fixed-backbone approximation for
    ddG scoring); only the identity-dependent node features change.
    """
    idx = residue_index(graph, chain, res_seq)
    feats = graph["features"].copy()
    row = np.zeros_like(feats[idx])
    row[:len(AA20)] = one_hot(new_res)
    row[len(AA20)] = HYDROPHOBICITY.get(new_res, 0.0) / 4.5
    feats[idx] = row
    return feats
