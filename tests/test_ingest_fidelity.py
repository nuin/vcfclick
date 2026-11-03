"""Row-building fidelity: what one VCF record becomes in the tables.

The unit tests here read small plain-text VCFs with cyvcf2 in-process
(no storage session) and call the row builders in ingest.vcf_rows
directly. The end-to-end tests go through the CLI on both backends.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from cyvcf2 import VCF

from ingest._arrow import (
    GENOTYPES_ARROW_SCHEMA,
    GENOTYPES_COLUMNS,
    VARIANTS_ARROW_SCHEMA,
    VARIANTS_COLUMNS,
    write_parquet,
)
from ingest.vcf_rows import build_record_rows
from tests.conftest import bgzip_vcf, run_cli

HEADER = """\
##fileformat=VCFv4.3
##contig=<ID=chr1,length=1000000>
##INFO=<ID=SB,Number=4,Type=Integer,Description="Strand bias">
##INFO=<ID=AC,Number=.,Type=Integer,Description="Allele count">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihoods">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Min depth">
##FORMAT=<ID=SBX,Number=4,Type=Integer,Description="Strand counts">
##FORMAT=<ID=FX,Number=1,Type=Float,Description="A float">
##FORMAT=<ID=FS,Number=1,Type=String,Description="A string">
"""


def _write_vcf(path: Path, samples: list[str], records: list[str]) -> Path:
    cols = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    body = "\n".join([cols + "\t" + "\t".join(samples), *records]) + "\n"
    path.write_text(HEADER + body)
    return path


def _records(path: Path):
    vcf = VCF(str(path))
    return list(vcf.samples), list(vcf)


def _rows(path: Path, keep_reference: bool = False):
    samples, variants = _records(path)
    from ingest.routing import classify_header

    extra = classify_header(VCF(str(path)))["extra_format"]
    out = []
    for v in variants:
        out.append(build_record_rows(v, samples, extra, "i1", keep_reference))
    return samples, out


def _gt_dict(genotype_rows: list[list]) -> dict[str, dict]:
    idx = {c: i for i, c in enumerate(GENOTYPES_COLUMNS)}
    return {r[idx["sample_id"]]: dict(zip(GENOTYPES_COLUMNS, r)) for r in genotype_rows}


def _variant_dict(variant_row: list) -> dict:
    return dict(zip(VARIANTS_COLUMNS, variant_row))


# ─── finding 6: numeric FORMAT overflow is stored as the numpy repr ────


def test_numeric_format_extra_is_comma_joined(tmp_path):
    """MIN_DP=5 is '5', SBX=1,2,3,4 is '1,2,3,4', a missing value is '.',
    and a shorter vector drops cyvcf2's end-of-vector padding."""
    vcf = _write_vcf(
        tmp_path / "x.vcf",
        ["A", "B", "C"],
        [
            "chr1\t100\t.\tA\tG\t50\tPASS\tDP=9\tGT:MIN_DP:SBX:FX:FS"
            "\t0/1:5:1,2,3,4:1.5:foo\t1/1:.:.:.:.\t0/1:7:1,2:2.5:bar"
        ],
    )
    _, rows = _rows(vcf)
    gts = _gt_dict(rows[0][1])
    assert gts["A"]["format_extra"] == {
        "MIN_DP": "5",
        "SBX": "1,2,3,4",
        "FX": "1.5",
        "FS": "foo",
    }
    assert gts["B"]["format_extra"] == {"MIN_DP": ".", "SBX": ".", "FX": ".", "FS": "."}
    assert gts["C"]["format_extra"]["SBX"] == "1,2"
    assert gts["C"]["format_extra"]["FX"] == "2.5"


@pytest.mark.parametrize("backend", ["duckdb", "chdb"])
def test_numeric_format_extra_queryable_after_ingest(vcfclick_home, tmp_path, backend):
    """CAST(format_extra['MIN_DP'] AS INT) works on a stored value."""
    src = _write_vcf(
        tmp_path / "fx.vcf",
        ["A", "B"],
        ["chr1\t100\t.\tA\tG\t50\tPASS\tDP=9\tGT:MIN_DP:SBX\t0/1:5:1,2,3,4\t1/1:.:."],
    )
    gz = bgzip_vcf(src, tmp_path / "fx.vcf.gz")
    run_cli(vcfclick_home, backend, "db", "create", "fx")
    run_cli(vcfclick_home, backend, "db", "ingest", "fx", str(gz), "--serial")
    r = run_cli(
        vcfclick_home,
        backend,
        "db",
        "query",
        "fx",
        "SELECT sample_id, format_extra['MIN_DP'], format_extra['SBX'] "
        "FROM genotypes ORDER BY sample_id",
        "--format",
        "JSONCompact",
    )
    assert json.loads(r.stdout)["data"] == [["A", "5", "1,2,3,4"], ["B", ".", "."]]


# ─── finding 36: a multi-value INFO field aimed at a scalar column ──────


def test_multi_value_info_scalar_goes_to_info_extra(tmp_path):
    """INFO SB declared Number=4 (GATK) cannot fill the float info_SB
    column: the tuple lands in info_extra and the row still writes."""
    vcf = _write_vcf(
        tmp_path / "sb.vcf",
        ["A"],
        ["chr1\t100\t.\tA\tG\t50\tPASS\tSB=1,2,3,4;AC=1,0;DP=9\tGT\t0/1"],
    )
    _, rows = _rows(vcf)
    variant = _variant_dict(rows[0][0])
    assert variant["info_SB"] is None
    assert variant["info_AC"] is None
    assert variant["info_DP"] == 9
    assert variant["info_extra"] == {"SB": "1,2,3,4", "AC": "1,0"}
    write_parquet([rows[0][0]], VARIANTS_ARROW_SCHEMA, tmp_path / "v.parquet")


# ─── finding 37: FORMAT integers above the uint16 column range ──────────


def test_wide_format_integers_clamp_to_uint16(tmp_path):
    """DP=70000 and AD=80000,5 are clamped to 65535, so the batch writes
    instead of aborting the ingest with an ArrowInvalid."""
    vcf = _write_vcf(
        tmp_path / "dp.vcf",
        ["A"],
        [
            "chr1\t100\t.\tA\tG\t50\tPASS\tDP=9\tGT:DP:AD:PL:PS"
            "\t0|1:70000:80000,5:0,70000,99:200000000"
        ],
    )
    _, rows = _rows(vcf)
    gt = _gt_dict(rows[0][1])["A"]
    assert (gt["dp"], gt["ad_ref"], gt["ad_alt"]) == (65535, 65535, 5)
    assert (gt["pl_ref_ref"], gt["pl_ref_alt"], gt["pl_alt_alt"]) == (0, 65535, 99)
    assert gt["ps"] == 200000000  # ps is UInt32: a phase set keeps its value
    write_parquet(rows[0][1], GENOTYPES_ARROW_SCHEMA, tmp_path / "g.parquet")


# ─── finding 40: samples but no GT ──────────────────────────────────────


def test_record_without_gt_builds_no_genotype_rows(tmp_path):
    """FORMAT=DP only: called counts are "not recorded" and there are no
    genotype rows, instead of a raw cyvcf2 exception."""
    vcf = _write_vcf(
        tmp_path / "nogt.vcf",
        ["A"],
        ["chr1\t200\t.\tA\tG\t50\tPASS\tDP=9\tDP\t10"],
    )
    _, rows = _rows(vcf)
    variant_row, genotype_rows, missing_rows = rows[0]
    variant = _variant_dict(variant_row)
    assert (variant["n_called"], variant["an_called"], variant["ac_called"]) == (
        None,
        None,
        None,
    )
    assert genotype_rows == []
    assert missing_rows == []


# ─── finding 7: partial calls under --keep-reference ────────────────────


def test_keep_reference_stores_no_row_for_a_partial_reference_call(tmp_path):
    """A call with a missing allele is never a proven 0/0: './0' and '0/.'
    get no row under --keep-reference. './1' stays a het (gt = 1), a full
    '0/0' and a haploid '0' are stored as gt = 0, './.' is never stored."""
    vcf = _write_vcf(
        tmp_path / "partial.vcf",
        ["A", "B", "C", "D", "E", "F"],
        ["chr1\t100\t.\tA\tG\t50\tPASS\tDP=9\tGT\t./0\t0/.\t./1\t0/0\t0\t./."],
    )
    _, rows = _rows(vcf, keep_reference=True)
    gts = _gt_dict(rows[0][1])
    assert {s: g["gt"] for s, g in gts.items()} == {"C": 1, "D": 0, "E": 0}
    variant = _variant_dict(rows[0][0])
    assert variant["n_called"] == 2  # D and E; A, B and C miss an allele
    assert [r[-1] for r in rows[0][2]] == ["F"]  # missing_genotypes: ./. only


class _NoAlleleArray:
    """A cyvcf2 record whose allele array cannot be read; the rest of the
    record (gt_types, FORMAT fields) is the real one."""

    def __init__(self, variant):
        self._variant = variant

    def __getattr__(self, name):
        if name == "genotype":
            raise Exception("couldn't get genotypes for variant")
        return getattr(self._variant, name)


def test_unreadable_allele_array_stores_no_gt0_row(tmp_path, caplog):
    """Without the allele array a './0' cannot be told from a '0/0', so
    under --keep-reference no gt=0 row is written for the record (and a
    warning names it), where a fallback to gt_types alone would store the
    half-call as a provable 0/0. Non-reference calls are kept."""
    import logging

    from ingest.routing import classify_header
    from ingest.vcf_rows import build_genotype_rows

    vcf = _write_vcf(
        tmp_path / "mask.vcf",
        ["A", "B", "C"],
        ["chr1\t100\t.\tA\tG\t50\tPASS\tDP=9\tGT\t./0\t0/0\t0/1"],
    )
    samples, (variant,) = _records(vcf)
    extra = classify_header(VCF(str(vcf)))["extra_format"]
    with caplog.at_level(logging.WARNING):
        rows = build_genotype_rows(
            _NoAlleleArray(variant), samples, extra, "i1", keep_reference=True
        )
    assert {s: g["gt"] for s, g in _gt_dict(rows).items()} == {"C": 1}
    assert "chr1:100" in caplog.text


@pytest.mark.parametrize("backend", ["duckdb", "chdb"])
def test_half_called_parent_is_no_denovo_evidence(vcfclick_home, tmp_path, backend):
    """Father './0' at a site where the child is 0/1: with --keep-reference
    the father has no stored row, so de novo does not report the site."""
    src = _write_vcf(
        tmp_path / "trio.vcf",
        ["kid", "dad", "mom"],
        [
            "chr1\t100\t.\tA\tG\t50\tPASS\tDP=9\tGT\t0/1\t./0\t0/0",
            "chr1\t200\t.\tC\tT\t50\tPASS\tDP=9\tGT\t0/1\t0/0\t0/0",
        ],
    )
    gz = bgzip_vcf(src, tmp_path / "trio.vcf.gz")
    run_cli(vcfclick_home, backend, "db", "create", "c2")
    run_cli(
        vcfclick_home,
        backend,
        "db",
        "ingest",
        "c2",
        str(gz),
        "--serial",
        "--keep-reference",
        "--ingest-id",
        "fam",
    )
    run_cli(
        vcfclick_home,
        backend,
        "db",
        "ped",
        "c2",
        "--proband",
        "kid",
        "--father",
        "dad",
        "--mother",
        "mom",
    )
    r = run_cli(
        vcfclick_home,
        backend,
        "db",
        "trio",
        "c2",
        "--proband",
        "kid",
        "--category",
        "denovo",
        "--format",
        "json",
    )
    doc = json.loads(r.stdout)
    positions = sorted(c["pos"] for c in doc["models"]["denovo"]["candidates"])
    assert positions == [200]


# ─── finding 39: a record with no ALT allele ────────────────────────────


def test_no_alt_record_names_the_remedy(tmp_path):
    """ALT '.' is not a multi-allelic site: the error says the record has
    no ALT allele and names the bcftools filter that drops such records.
    A real multi-allelic record keeps the norm remedy."""
    from ingest.vcf_rows import require_single_alt

    vcf = _write_vcf(
        tmp_path / "alt.vcf",
        ["A"],
        [
            "chr1\t100\t.\tA\t.\t50\tPASS\tDP=9\tGT\t0/0",
            "chr1\t200\t.\tA\tG,T\t50\tPASS\tDP=9\tGT\t1/2",
        ],
    )
    _, (ref_only, multi) = _records(vcf)
    with pytest.raises(ValueError) as e:
        require_single_alt(ref_only, "in.vcf.gz")
    msg = str(e.value)
    assert "chr1:100" in msg and "no ALT allele" in msg
    assert "bcftools view -e 'ALT=\".\"'" in msg
    assert "Multi-allelic" not in msg
    with pytest.raises(ValueError, match="Multi-allelic site at chr1:200"):
        require_single_alt(multi, "in.vcf.gz")
    assert "bcftools norm -m -" in str(
        pytest.raises(ValueError, require_single_alt, multi, "in.vcf.gz").value
    )


# ─── finding 38: the serial stage directory ─────────────────────────────


def test_list_dbs_skips_a_leftover_stage_dir(vcfclick_home):
    """A SIGKILLed serial ingest on an older version left
    <home>/dbs/vcfclick_ingest_xxxx behind; `db list` must not show it."""
    from storage import list_dbs

    dbs = vcfclick_home / "dbs"
    (dbs / "cohort").mkdir(parents=True)
    (dbs / "vcfclick_ingest_abc123").mkdir()
    assert list_dbs() == ["cohort"]


_SERIAL_SPY = """
import collections, json, sys
import ingest.vcf_load as vl
import storage

vl.BATCH_SIZE = 1
inserts = collections.Counter()
original = vl._commit_staged_ingest


def spy(staging_path, *args):
    print("STAGE", staging_path)
    print("DBS", json.dumps(storage.list_dbs()))
    sess = storage.get_session()
    query = sess.query

    def counted(sql, *a, **k):
        words = sql.split()
        if words[:2] == ["INSERT", "INTO"]:
            inserts[words[2]] += 1
        return query(sql, *a, **k)

    sess.query = counted
    try:
        return original(staging_path, *args)
    finally:
        sess.query = query


vl._commit_staged_ingest = spy
vl.ingest(sys.argv[1], "c", ingest_id="i1")
print("INSERTS", json.dumps(inserts))
"""


def _serial_spy(home, run_python, vcf) -> dict[str, str]:
    out = run_python(
        home,
        _SERIAL_SPY,
        str(vcf),
        VCFCLICK_BACKEND="duckdb",
        VCFCLICK_DB_NAME="b",
    )
    return dict(line.split(" ", 1) for line in out.splitlines() if " " in line)


def test_serial_stage_dir_lives_inside_the_db(vcfclick_home, run_python, tiny_vcf):
    """The serial stage directory is under the DB's own directory, so an
    in-flight or killed ingest never appears in `db list`."""
    lines = _serial_spy(vcfclick_home, run_python, tiny_vcf)
    assert Path(lines["STAGE"]).is_relative_to(vcfclick_home / "dbs" / "b")
    assert json.loads(lines["DBS"]) == ["b"]


def test_serial_ingest_wipes_a_stale_stage_dir(vcfclick_home, tiny_vcf):
    """A SIGKILLed serial ingest leaves <db>/staging/serial-<ingest_id>.
    The next ingest under that id removes it before staging, so its stale
    batches are neither imported nor left on disk."""
    run_cli(vcfclick_home, "duckdb", "db", "create", "b")
    stale = vcfclick_home / "dbs" / "b" / "staging" / "serial-i1"
    stale.mkdir(parents=True)
    (stale / "v_9999999999.parquet").write_bytes(b"not parquet")
    run_cli(
        vcfclick_home,
        "duckdb",
        "db",
        "ingest",
        "b",
        str(tiny_vcf),
        "--serial",
        "--ingest-id",
        "i1",
    )
    assert not stale.exists()
    assert _count_variants(vcfclick_home, "duckdb", "b", "i1") == 5


def _count_variants(home, backend, db, ingest_id) -> int:
    r = run_cli(
        home,
        backend,
        "db",
        "query",
        db,
        f"SELECT COUNT(*) FROM variants WHERE ingest_id = '{ingest_id}'",
        "--format",
        "JSONCompact",
    )
    return int(json.loads(r.stdout)["data"][0][0])


# ─── finding 43: one INSERT per staged table ────────────────────────────


def test_serial_commit_inserts_each_table_once(vcfclick_home, run_python, tiny_vcf):
    """Five one-site batches commit with one INSERT per staged table, as
    the parallel path does, not one INSERT per batch file."""
    inserts = json.loads(_serial_spy(vcfclick_home, run_python, tiny_vcf)["INSERTS"])
    assert inserts.get("variants") == 1, inserts
    assert inserts.get("genotypes") == 1, inserts
    assert inserts.get("missing_genotypes", 0) <= 1, inserts


# ─── findings 41/42: the per-batch row cap ──────────────────────────────


def test_batch_row_cap_bounds_worker_memory():
    """2,000,000 buffered rows measured 1.5 to 1.75 GB per process, and
    the 4 default workers hold a batch each. 250,000 rows keeps one
    process near 200 MB."""
    from ingest.vcf_rows import MAX_BATCH_ROWS

    assert MAX_BATCH_ROWS <= 250_000


# ─── finding 45: one rollback-on-failure helper ─────────────────────────


@pytest.mark.parametrize("tag", ["ingest", "parallel-ingest", "parquet-ingest"])
def test_ingest_failure_helper_rolls_back_only_after_commit(monkeypatch, caplog, tag):
    """The serial, parallel and Parquet loaders share one failure path:
    roll back only once engine writes began, log a rollback that fails,
    and say 'no engine writes' (not 'chDB') for a parse failure."""
    import logging

    import storage
    from ingest._rollback import handle_ingest_failure

    calls: list[str] = []
    monkeypatch.setattr(storage, "rollback_ingest", calls.append)
    with caplog.at_level(logging.WARNING):
        handle_ingest_failure(False, "i1", tag)
    assert calls == []
    assert f"[{tag}] failed during parse" in caplog.text
    assert "no engine writes" in caplog.text

    handle_ingest_failure(True, "i1", tag)
    assert calls == ["i1"]

    def broken(ingest_id):
        raise RuntimeError("disk full")

    monkeypatch.setattr(storage, "rollback_ingest", broken)
    caplog.clear()
    with caplog.at_level(logging.WARNING):
        handle_ingest_failure(True, "i1", tag)
    assert f"[{tag}] rollback FAILED: disk full" in caplog.text


def test_loaders_share_the_failure_helper():
    import inspect

    import ingest.parallel
    import ingest.parquet_load
    import ingest.vcf_load

    for module in (ingest.vcf_load, ingest.parallel, ingest.parquet_load):
        source = inspect.getsource(module)
        assert "handle_ingest_failure(" in source, module.__name__
        assert "rollback FAILED" not in source, module.__name__
