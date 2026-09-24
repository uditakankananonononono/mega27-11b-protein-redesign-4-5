"""Benchmark dataset loading (ProtDDG-Bench: S2648, Ssym, P53 sets).

Data provenance (cited in the paper):
  S2648 - Dehouck et al. 2011 (PMID 21569468), 2648 variants / 132 structures
  Ssym  - Usmanova et al. 2018 (PMID 29718106), 342 mutations + reverses
  P53   - 42 variants on 2OCJ chain A (PMID 24281696)
Files come from github.com/protddg-bench/protddg-bench (CC BY-NC-SA).
"""
from __future__ import annotations

import csv
from dataclasses import dataclass

AA1_TO_3 = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
    "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
    "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}
AA3_TO_1 = {v: k for k, v in AA1_TO_3.items()}


@dataclass
class MutationRecord:
    pdb_id: str       # e.g. "1AMQ"
    chain: str        # e.g. "A"
    wt_aa3: str       # wild-type residue, 3-letter
    position: int     # residue sequence number
    mut_aa3: str      # mutant residue, 3-letter
    ddg: float        # experimental ddG, kcal/mol (positive = stabilizing per benchmark sign convention)
    direction: str    # "DIR" or "INV" (Ssym), "" otherwise


def split_pdb_chain(code: str) -> tuple[str, str]:
    """'1AMQA' -> ('1AMQ', 'A'). Handles 4-letter PDB + 1-letter chain."""
    code = code.strip()
    if len(code) < 5:
        raise ValueError(f"cannot split PDB+chain code: {code!r}")
    return code[:4].upper(), code[4]


def parse_mutation(mut: str) -> tuple[str, int, str]:
    """'C191Y' -> ('CYS', 191, 'TYR')."""
    mut = mut.strip()
    wt, pos, mt = mut[0], mut[1:-1], mut[-1]
    if wt not in AA1_TO_3 or mt not in AA1_TO_3:
        raise ValueError(f"unknown amino acid in mutation {mut!r}")
    return AA1_TO_3[wt], int(pos), AA1_TO_3[mt]


def load_bench_tsv(path: str, skipped: list | None = None) -> list[MutationRecord]:
    """Load a ProtDDG-Bench TSV (comment header starts with #).

    Rows whose mutation string carries a PDB insertion code (e.g. L27CN on
    1LVE) cannot be mapped by the fixed-column parser, which discards
    iCodes. They are skipped and recorded in `skipped` (if a list is given)
    as (row, reason) pairs - documented data cleaning, never silent.
    """
    records = []
    with open(path) as fh:
        reader = csv.reader((ln for ln in fh if not ln.startswith("#")),
                            delimiter="\t")
        for row in reader:
            if len(row) < 5:
                continue
            row = [c.strip() for c in row]
            try:
                pdb_id, chain = split_pdb_chain(row[2])
                wt3, pos, mt3 = parse_mutation(row[3])
                ddg = float(row[4])
            except ValueError as err:
                if skipped is not None:
                    skipped.append((row, str(err)))
                continue
            direction = row[8] if len(row) > 8 else ""
            records.append(MutationRecord(
                pdb_id=pdb_id, chain=chain, wt_aa3=wt3, position=pos,
                mut_aa3=mt3, ddg=ddg, direction=direction,
            ))
    return records


def unique_structures(records: list[MutationRecord]) -> list[tuple[str, str]]:
    """Sorted unique (pdb_id, chain) pairs needed by a record set."""
    return sorted({(r.pdb_id, r.chain) for r in records})
