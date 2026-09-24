

def test_ca_trace_dedupes_altloc():
    """Altloc duplicate CA records collapse to one node per residue."""
    from src.pdb_fetch import load_structure
    from src.contact_graph import ca_trace
    import os
    pdb = "data/pdb/2VUK.pdb"
    if not os.path.exists(pdb):
        import pytest
        pytest.skip("2VUK not cached")
    tr = ca_trace(load_structure("2VUK", "data/pdb"), chain="A")
    seqs = [a.res_seq for a in tr]
    assert len(seqs) == len(set(seqs)) == 195
