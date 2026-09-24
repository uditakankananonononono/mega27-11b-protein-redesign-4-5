"""Hermetic smoke test for tools/gnn_ddg.py (no network, no real ckpt)."""
import json
import os
import pickle
import subprocess
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from src.gnn_ddg_v2 import GNNddGv2


def _write_pdb(path, n=12):
    """Synthetic 12-residue ALA chain, gently curved CA trace."""
    lines = []
    for i in range(n):
        x, y, z = 3.8 * i * 0.4, 3.8 * np.sin(i * 0.6), 3.8 * np.cos(i * 0.6)
        lines.append(
            f"ATOM  {i+1:>5}  CA  ALA A{i+1:>4}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C  ")
    lines.append("TER\nEND")
    path.write_text("\n".join(lines) + "\n")


def test_cli_scan_smoke(tmp_path):
    pdb_dir = tmp_path / "pdb"
    pdb_dir.mkdir()
    _write_pdb(pdb_dir / "TEST.pdb")
    net = GNNddGv2(in_dim=23, hidden=64, mlp_hidden=32, seed=1)
    ckpt = tmp_path / "ckpt.pkl"
    with open(ckpt, "wb") as fh:
        pickle.dump({"params": net.params, "epoch": 0}, fh)
    out = tmp_path / "scan.json"
    r = subprocess.run(
        [sys.executable, os.path.join(ROOT, "tools", "gnn_ddg.py"),
         "scan", "--pdb", "TEST", "--chain", "A", "--top", "5",
         "--ckpt", str(ckpt), "--bias", str(tmp_path / "no-bias.json"),
         "--pdb-dir", str(pdb_dir), "--out", str(out)],
        capture_output=True, text=True, cwd=ROOT)
    assert r.returncode == 0, r.stderr
    d = json.load(open(out))
    assert d["n_positions"] == 12
    assert d["n_scanned"] == 12 * 19
    assert len(d["candidates"]) == 5
    assert all(not r_["ood_wt_pro"] for r_ in d["candidates"])


def test_cli_predict_smoke(tmp_path):
    pdb_dir = tmp_path / "pdb"
    pdb_dir.mkdir()
    _write_pdb(pdb_dir / "TEST.pdb")
    net = GNNddGv2(in_dim=23, hidden=64, mlp_hidden=32, seed=2)
    ckpt = tmp_path / "ckpt.pkl"
    with open(ckpt, "wb") as fh:
        pickle.dump({"params": net.params, "epoch": 0}, fh)
    r = subprocess.run(
        [sys.executable, os.path.join(ROOT, "tools", "gnn_ddg.py"),
         "predict", "--pdb", "TEST", "--chain", "A", "--mutation", "A5V",
         "--ckpt", str(ckpt), "--bias", str(tmp_path / "no-bias.json"),
         "--pdb-dir", str(pdb_dir)],
        capture_output=True, text=True, cwd=ROOT)
    assert r.returncode == 0, r.stderr
    d = json.loads(r.stdout)
    assert d["mutation"] == "A5V"
    assert d["wt"] == "ALA" and d["mut"] == "VAL"
