"""Hermetic parser tests. Fixtures are synthetic minimal ATOM records used
only to exercise the parser; they are not experimental data and are never
used in any reported result."""
import math

from src.pdb_fetch import ca_trace, distance_matrix, parse_pdb


def _atom(serial, name, res, chain, seq, x, y, z, elem):
    # Fixed-column PDB ATOM record (cols: 1-6 rec, 7-11 serial, 13-16 name,
    # 18-20 resName, 22 chain, 23-26 resSeq, 31-38/39-46/47-54 xyz,
    # 55-60 occ, 61-66 temp, 77-78 element).
    return (
        f"ATOM  {serial:5d} {name:>4s} {res:>3s} {chain}{seq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{20.00:6.2f}          {elem:>2s}"
    )


FIXTURE = "\n".join([
    "HEADER    SYNTHETIC PARSER FIXTURE - NOT EXPERIMENTAL DATA",
    _atom(1, "N", "GLY", "A", 1, -1.2, 0.4, 0.0, "N"),
    _atom(2, "CA", "GLY", "A", 1, 0.0, 0.0, 0.0, "C"),
    _atom(3, "C", "GLY", "A", 1, 0.5, -1.4, 0.3, "C"),
    _atom(4, "CA", "ALA", "A", 2, 3.8, 0.0, 0.0, "C"),
    _atom(5, "CA", "VAL", "B", 1, 0.0, 5.1, 0.0, "C"),
    "END",
])


def test_parse_counts_and_fields():
    atoms = parse_pdb(FIXTURE)
    assert len(atoms) == 5
    ca = atoms[1]
    assert ca.name == "CA" and ca.res_name == "GLY"
    assert ca.chain == "A" and ca.res_seq == 1
    assert ca.element == "C"
    assert abs(ca.x - 0.0) < 1e-9 and abs(ca.z - 0.0) < 1e-9


def test_ca_trace_chain_filter_and_order():
    atoms = parse_pdb(FIXTURE)
    trace_a = ca_trace(atoms, chain="A")
    assert [a.res_seq for a in trace_a] == [1, 2]
    assert all(a.chain == "A" for a in trace_a)
    assert len(ca_trace(atoms)) == 3  # both chains


def test_distance_matrix_symmetric_and_exact():
    coords = [(0.0, 0.0, 0.0), (3.8, 0.0, 0.0), (0.0, 5.1, 0.0)]
    dm = distance_matrix(coords)
    assert abs(dm[0][1] - 3.8) < 1e-9
    assert abs(dm[0][2] - 5.1) < 1e-9
    assert abs(dm[1][2] - math.sqrt(3.8**2 + 5.1**2)) < 1e-9
    for i in range(3):
        for j in range(3):
            assert dm[i][j] == dm[j][i] and dm[i][i] == 0.0


def test_whitespace_fallback_parsing():
    mangled = "ATOM 7 CA ALA A 3 1.0 2.0 3.0 C"
    atoms = parse_pdb(mangled)
    assert len(atoms) == 1
    assert atoms[0].res_seq == 3 and abs(atoms[0].y - 2.0) < 1e-9
