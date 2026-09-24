#!/bin/sh
# Re-fetch all external data (benchmark sets + PDB structures).
# Provenance: ProtDDG-Bench (github.com/protddg-bench/protddg-bench,
# CC BY-NC-SA); structures from RCSB (https://www.rcsb.org).
set -e
cd "$(dirname "$0")/.."
if [ ! -d data/bench/protddg-bench ]; then
  git clone --depth 1 https://github.com/protddg-bench/protddg-bench.git data/bench/protddg-bench
  rm -rf data/bench/protddg-bench/.git
fi
python3 scripts/download_pdbs.py
echo "data ready"
