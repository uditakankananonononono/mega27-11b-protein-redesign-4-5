"""Group-held-out validation split: no PDB leakage, deterministic."""
from src.dataset import MutationRecord, split_by_group


def _recs():
    recs = []
    for p in range(6):
        for i in range(5):
            recs.append(MutationRecord(f"PDB{p}", "A", "ALA", i + 1,
                                       "CYS", -1.0, ""))
    return recs


def test_no_pdb_overlap_between_sides():
    recs = _recs()
    tr, va = split_by_group(recs, 0.34, seed=0)
    tr_pdbs = {recs[i].pdb_id for i in tr}
    va_pdbs = {recs[i].pdb_id for i in va}
    assert tr_pdbs.isdisjoint(va_pdbs)
    assert len(tr) + len(va) == len(recs)
    assert 0 < len(va) < len(recs)


def test_deterministic_same_seed():
    recs = _recs()
    assert split_by_group(recs, 0.34, seed=7) == split_by_group(recs, 0.34,
                                                                seed=7)


def test_zero_frac_gives_empty_val():
    recs = _recs()
    tr, va = split_by_group(recs, 0.0, seed=0)
    assert va == [] and len(tr) == len(recs)
