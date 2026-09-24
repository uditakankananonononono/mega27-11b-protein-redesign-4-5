"""Download every PDB structure required by the benchmark datasets.
Usage: python3 scripts/download_pdbs.py
Idempotent: skips structures already cached in data/pdb/."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.dataset import load_bench_tsv, unique_structures
from src.pdb_fetch import fetch_pdb

BENCH = "data/bench/protddg-bench"
FILES = [
    f"{BENCH}/SSYM/train-s2648-test-ssym.tsv",
    f"{BENCH}/SSYM/ssym-5fold.tsv",
    f"{BENCH}/P53/p53.tsv",
]


def main():
    wanted = set()
    for path in FILES:
        wanted |= set(unique_structures(load_bench_tsv(path)))
    # study structures (WT references + Y220C mutants + SOD1)
    wanted |= {("2XWR", "A"), ("1TSR", "B"), ("2J1X", "A"), ("2VUK", "A"),
               ("2OCJ", "A"), ("1SPD", "A"), ("2C9V", "A")}
    ids = sorted({pdb for pdb, _ in wanted})
    ok, failed = 0, []
    for i, pdb in enumerate(ids, 1):
        try:
            fetch_pdb(pdb, "data/pdb", retries=2, timeout=20)
            ok += 1
            if i % 20 == 0:
                print(f"[{i}/{len(ids)}] ok so far", flush=True)
        except Exception as err:  # noqa: BLE001
            failed.append((pdb, str(err)))
            print(f"FAILED {pdb}: {err}", flush=True)
    print(f"DONE: {ok}/{len(ids)} structures cached; {len(failed)} failed")
    for pdb, err in failed:
        print(f"  {pdb}: {err}")


if __name__ == "__main__":
    main()
