"""Hermetic dataset-parsing tests with synthetic TSV content (no real
benchmark rows needed)."""
import textwrap

from src.dataset import (load_bench_tsv, parse_mutation, split_pdb_chain,
                         unique_structures)


def test_split_pdb_chain():
    assert split_pdb_chain("1AMQA") == ("1AMQ", "A")
    assert split_pdb_chain("2ocjB") == ("2OCJ", "B")


def test_parse_mutation():
    assert parse_mutation("C191Y") == ("CYS", 191, "TYR")
    assert parse_mutation("A4V") == ("ALA", 4, "VAL")
    try:
        parse_mutation("X4V")
        assert False, "should have raised"
    except ValueError:
        pass


def test_load_bench_tsv(tmp_path):
    content = textwrap.dedent("""\
        #SET\tCLID\tPDB\tMUT\tDDG\tpH\tT\tID\tDIR/INV
        SET_0\t1AMQA\t1AMQA\tC191Y\t-2.30\t7.50\t25.00\t1\tDIR
        SET_0\t1AMQA\t1QIRA\tY191C\t2.30\t7.50\t25.00\t1\tINV
    """)
    p = tmp_path / "toy.tsv"
    p.write_text(content)
    recs = load_bench_tsv(str(p))
    assert len(recs) == 2
    r = recs[0]
    assert (r.pdb_id, r.chain, r.wt_aa3, r.position, r.mut_aa3) == \
        ("1AMQ", "A", "CYS", 191, "TYR")
    assert abs(r.ddg + 2.30) < 1e-9 and r.direction == "DIR"
    assert unique_structures(recs) == [("1AMQ", "A"), ("1QIR", "A")]


def test_unparseable_rows_are_logged_and_skipped(tmp_path):
    content = "#H\nSET_0\t1LVEA\t1LVEA\tL27CN\t-1.0\t7.0\t25.0\t9\tDIR\nSET_0\t1AMQA\t1AMQA\tC191Y\t-2.30\t7.5\t25.0\t1\tDIR\n"
    p = tmp_path / "toy2.tsv"
    p.write_text(content)
    skipped = []
    recs = load_bench_tsv(str(p), skipped=skipped)
    assert len(recs) == 1 and recs[0].position == 191
    assert len(skipped) == 1 and "27C" in skipped[0][1]
