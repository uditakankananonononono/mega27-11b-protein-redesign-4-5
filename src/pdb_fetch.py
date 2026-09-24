"""PDB download and parsing utilities for the MEGA27-11b redesign pipeline.

All reported results use experimentally determined structures downloaded from
RCSB (https://www.rcsb.org) and cached under data/pdb/. Parser unit tests use
small synthetic fixtures constructed in the test module; those fixtures are
never used as experimental data.
"""
from __future__ import annotations

import math
import os
import time
import urllib.request
from dataclasses import dataclass

RCSB_DOWNLOAD = "https://files.rcsb.org/download/{pdb_id}.pdb"


def fetch_pdb(pdb_id: str, cache_dir: str, retries: int = 3, timeout: int = 30) -> str:
    """Download a PDB file from RCSB with local caching; return the file path."""
    pdb_id = pdb_id.upper()
    os.makedirs(cache_dir, exist_ok=True)
    path = os.path.join(cache_dir, f"{pdb_id}.pdb")
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return path
    url = RCSB_DOWNLOAD.format(pdb_id=pdb_id)
    last_err: Exception | None = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as resp:
                data = resp.read()
            if not data[:6].split()[0] in (b"HEADER", b"REMARK", b"ATOM", b"CRYST1", b"TITLE"):
                raise ValueError(f"unexpected payload for {pdb_id}")
            with open(path, "wb") as fh:
                fh.write(data)
            return path
        except Exception as err:  # re-raised after retries
            last_err = err
            time.sleep(2 ** attempt)
    raise RuntimeError(f"failed to fetch {pdb_id} from RCSB: {last_err}")


@dataclass
class Atom:
    serial: int
    name: str
    res_name: str
    chain: str
    res_seq: int
    x: float
    y: float
    z: float
    element: str

    def distance(self, other: "Atom") -> float:
        return math.sqrt(
            (self.x - other.x) ** 2
            + (self.y - other.y) ** 2
            + (self.z - other.z) ** 2
        )


def parse_pdb(text: str) -> list[Atom]:
    """Parse ATOM records from PDB text (fixed columns, whitespace fallback)."""
    atoms: list[Atom] = []
    for line in text.splitlines():
        if not line.startswith("ATOM"):
            continue
        try:
            atoms.append(Atom(
                serial=int(line[6:11]),
                name=line[12:16].strip(),
                res_name=line[17:20].strip(),
                chain=line[21].strip() or "_",
                res_seq=int(line[22:26]),
                x=float(line[30:38]),
                y=float(line[38:46]),
                z=float(line[46:54]),
                element=line[76:78].strip() or line[12:16].strip()[0],
            ))
        except (ValueError, IndexError):
            parts = line.split()
            if len(parts) < 9:
                continue
            atoms.append(Atom(
                serial=int(parts[1]), name=parts[2], res_name=parts[3],
                chain=parts[4], res_seq=int(parts[5]),
                x=float(parts[6]), y=float(parts[7]), z=float(parts[8]),
                element=parts[-1],
            ))
    return atoms


def load_structure(pdb_id: str, cache_dir: str) -> list[Atom]:
    """Fetch (cached) and parse a structure by PDB id."""
    path = fetch_pdb(pdb_id, cache_dir)
    with open(path) as fh:
        return parse_pdb(fh.read())


def ca_trace(atoms: list[Atom], chain: str | None = None) -> list[Atom]:
    """Return C-alpha atoms ordered by (chain, residue number)."""
    out = [a for a in atoms if a.name == "CA" and (chain is None or a.chain == chain)]
    out.sort(key=lambda a: (a.chain, a.res_seq))
    return out


def distance_matrix(coords: list[tuple[float, float, float]]) -> list[list[float]]:
    """Pairwise Euclidean distance matrix (list of lists), symmetric."""
    n = len(coords)
    dm = [[0.0] * n for _ in range(n)]
    for i in range(n):
        xi, yi, zi = coords[i]
        for j in range(i + 1, n):
            xj, yj, zj = coords[j]
            d = math.sqrt((xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2)
            dm[i][j] = dm[j][i] = d
    return dm
