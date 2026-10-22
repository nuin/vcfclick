"""Byte-for-byte goldens of `db trio` output, on both backends.

The goldens under tests/goldens/trio/<backend>/<store>/ were captured from
the CLI at v0.14.3, BEFORE the trio-output refactor. They pin what the
desktop apps, scripts and the MCP briefing read:

  * text output (stdout and stderr) is compared verbatim;
  * JSON output is compared after removing the keys added since, additively
    (see ADDED_KEYS). Everything else must be byte-identical:
    same keys, same order, same values, same indent.

The JSON compatibility policy (docs/TRIO.md) is "additive only": a new
key may appear, an old one may not move, change type or disappear. When a
later step legitimately adds a key, it is added to ADDED_KEYS here in the
same commit; nothing else in this file may change for it.

Regenerate (only when an OLD key must change on purpose, and say why in
the commit):

    VCFCLICK_UPDATE_GOLDENS=1 .venv/bin/python -m pytest \\
        tests/test_trio_goldens.py -p no:cacheprovider
"""

from __future__ import annotations

import difflib
import gzip
import json
import os
import re
import shutil
from dataclasses import dataclass, field
from pathlib import Path

import pytest

from tests.conftest import FIXTURES, run_cli

GOLDENS = Path(__file__).parent / "goldens" / "trio"
UPDATE = os.environ.get("VCFCLICK_UPDATE_GOLDENS") == "1"
BACKENDS = ("chdb", "duckdb")

# Keys that may be added to a trio document, by location. Anything else
# that is not in the pre-refactor golden fails the comparison.
ADDED_KEYS = {
    "model": {"blocked_reason"},  # models.<name>.blocked_reason
    "gates": {"af_gate"},  # gates.af_gate
    "trio": {"family_id", "kind", "siblings"},  # trio.<key> (pedigree P1)
    "top": {"families"},  # the families of a no-`--proband` {"trios": [...]}
}

GIAB = FIXTURES / "giab"
CATEGORIES = ("denovo", "recessive", "dominant", "comphet")
# Goldens above this size are stored gzipped (mtime 0, so regenerating is
# byte-stable) to keep the public repo small; the test diffs decoded text.
GZIP_ABOVE = 32 * 1024


@dataclass(frozen=True)
class Store:
    """One fixture cohort and the `db trio` invocations pinned for it."""

    key: str
    vcf: Path
    ped_args: tuple[str, ...]
    runs: tuple[tuple[str, ...], ...]
    keep_reference: bool = True
    gene: tuple | None = None  # (symbol, chrom, start, end) seeded as an annotation
    gnomad: Path | None = None
    ingest_id: str = "g1"
    sample_note: str = field(default="", compare=False)


def _ped(path: Path) -> tuple[str, ...]:
    return (str(path),)


def _runs(
    proband: tuple[str, ...],
    *,
    extra: tuple[str, ...] = (),
    cats=CATEGORIES,
    explicit_all: bool = False,
) -> tuple[tuple[str, ...], ...]:
    """Every JSON and text run for one proband: the default (category all),
    optionally an explicit `--category all`, and each model on its own."""
    out = []
    for fmt in ("json", "text"):
        out.append((*proband, *extra, "--format", fmt))
        if explicit_all:
            out.append((*proband, *extra, "--category", "all", "--format", fmt))
        for cat in cats:
            out.append((*proband, *extra, "--category", cat, "--format", fmt))
    return tuple(out)


C1 = ("--proband", "C1")
STORES = (
    Store(
        "related_c1",
        FIXTURES / "related.vcf",
        ("--proband", "C1", "--father", "F", "--mother", "M"),
        _runs(C1, explicit_all=True),
    ),
    Store(
        # No --proband: every child in the pedigree, one report each.
        "related_siblings",
        FIXTURES / "related.vcf",
        ("--proband", "C1", "--proband", "C2", "--father", "F", "--mother", "M"),
        (
            ("--format", "json"),
            ("--format", "text"),
            ("--category", "denovo", "--format", "json"),
            ("--category", "denovo", "--format", "text"),
            ("--category", "comphet", "--format", "json"),
        ),
    ),
    Store(
        # A sparse ingest: de novo / dominant / comphet are blocked and the
        # text mode explains why on stderr.
        "related_sparse",
        FIXTURES / "related.vcf",
        ("--proband", "C1", "--father", "F", "--mother", "M"),
        (
            ("--proband", "C1", "--format", "json"),
            ("--proband", "C1", "--format", "text"),
            ("--proband", "C1", "--category", "denovo", "--format", "text"),
            ("--proband", "C1", "--category", "recessive", "--format", "json"),
        ),
        keep_reference=False,
    ),
    Store(
        "trio_default",
        FIXTURES / "trio.vcf.gz",
        _ped(FIXTURES / "trio.ped"),
        _runs(("--proband", "CHILD"), cats=("denovo",)),
        gene=("GENEX", "chr1", 1, 100000),
        ingest_id="fam1",
    ),
    Store(
        # Relaxed rarity so every model has candidates (AF 0.01 .. 0.2 pass).
        "trio_relaxed",
        FIXTURES / "trio.vcf.gz",
        _ped(FIXTURES / "trio.ped"),
        _runs(("--proband", "CHILD"), extra=("--max-af", "0.3")),
        gene=("GENEX", "chr1", 1, 100000),
        ingest_id="fam1",
    ),
    Store(
        "trio_sparse",
        FIXTURES / "trio.vcf.gz",
        _ped(FIXTURES / "trio.ped"),
        (
            ("--proband", "CHILD", "--format", "json"),
            ("--proband", "CHILD", "--format", "text"),
            ("--proband", "CHILD", "--max-af", "0.3", "--format", "json"),
        ),
        keep_reference=False,
        ingest_id="fam1",
    ),
    Store(
        # Real GIAB genotypes at CFTR, with a gene and gnomAD loaded so the
        # candidate enrichment keys (genes, gnomad_popmax) are populated.
        "giab_cftr",
        GIAB / "cftr_trio.vcf.gz",
        _ped(GIAB / "cftr_trio.ped"),
        _runs(("--proband", "HG002"))
        + (
            ("--proband", "HG002", "--gnomad-max-af", "0.01", "--format", "json"),
            ("--proband", "HG002", "--gnomad-max-af", "0.01", "--format", "text"),
        ),
        gene=("CFTR", "chr7", 117480025, 117668665),
        gnomad=FIXTURES / "gnomad_cftr.vcf.gz",
    ),
    Store(
        "giab_denovo",
        GIAB / "denovo_trio.vcf.gz",
        _ped(GIAB / "denovo_trio.ped"),
        _runs(("--proband", "HG002"), cats=("denovo",)),
    ),
)


def _golden_name(run: tuple[str, ...]) -> str:
    """A stable file stem for one invocation's arguments."""
    parts = [a.lstrip("-") for a in run if a not in ("--proband",)]
    return "__".join(p.replace("/", "_") for p in parts)


def _build(store: Store, home: Path, backend: str, monkeypatch, tmp_path) -> None:
    """Create the db, ingest, load the pedigree and the annotations."""
    monkeypatch.setenv("VCFCLICK_ANNOTATIONS_DB", str(tmp_path / "ann.duckdb"))
    if store.gene is not None:
        # In-process DuckDB (annotation store), never a chDB session.
        from annotations.db import get_connection

        conn = get_connection()
        sym, chrom, start, end = store.gene
        conn.execute(
            "INSERT INTO refseq_genes VALUES (?, ?, ?, ?, '+', 'id', 'protein_coding')",
            [sym, chrom, start, end],
        )
        conn.close()
    if store.gnomad is not None:
        run_cli(home, backend, "annotations", "load-gnomad", str(store.gnomad))
    run_cli(home, backend, "db", "create", "t")
    ingest = [
        "db", "ingest", "t", str(store.vcf),
        "--cohort", "fam", "--ingest-id", store.ingest_id, "--serial",
    ]  # fmt: skip
    if store.keep_reference:
        ingest.append("--keep-reference")
    run_cli(home, backend, *ingest)
    run_cli(home, backend, "db", "ped", "t", *store.ped_args)


def _strip_added(doc: dict) -> dict:
    """Remove the additive keys from a trio document (or a `trios` list)."""
    docs = doc["trios"] if "trios" in doc else [doc]
    if "trios" in doc:
        for key in ADDED_KEYS["top"]:
            doc.pop(key, None)
    for d in docs:
        for key in ADDED_KEYS["trio"]:
            d["trio"].pop(key, None)
        for key in ADDED_KEYS["gates"]:
            d["gates"].pop(key, None)
        for model in d["models"].values():
            for key in ADDED_KEYS["model"]:
                model.pop(key, None)
    return doc


# Log lines the embedded engine prints on its own, unrelated to vcfclick's
# output: the Linux chDB build (26.3.x) announces its allocator setting on
# stderr; the macOS build does not.
ENGINE_NOISE = re.compile(r"^Default dirty pages decay period: \d+ms\n", re.M)


def _engine_quiet(stderr: str) -> str:
    return ENGINE_NOISE.sub("", stderr)


def _normalise(run: tuple[str, ...], stdout: str) -> str:
    if "json" in run:
        return json.dumps(_strip_added(json.loads(stdout)), indent=2) + "\n"
    return stdout


def _read(path: Path) -> str | None:
    """The golden at `path`, or at `path`.gz when it was large; None if
    neither exists."""
    gz = path.with_name(path.name + ".gz")
    if path.exists():
        return path.read_text()
    if gz.exists():
        return gzip.decompress(gz.read_bytes()).decode()
    return None


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if len(text) <= GZIP_ABOVE:
        path.write_text(text)
        return
    gz = path.with_name(path.name + ".gz")
    with (
        open(gz, "wb") as raw,
        gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as f,
    ):
        f.write(text.encode())


def _check(path: Path, got: str | None, label: str, failures: list[str]) -> None:
    """Compare one captured stream against its golden file (None = no file,
    i.e. the stream must be empty)."""
    want = _read(path)
    if UPDATE:
        if got:
            _write(path, got)
        return
    if (want or "") == (got or ""):
        return
    diff = "".join(
        difflib.unified_diff(
            (want or "").splitlines(True),
            (got or "").splitlines(True),
            f"golden/{path.name}",
            f"actual/{label}",
            n=2,
        )
    )
    failures.append(f"{label}: differs from {path}\n{diff[:3000]}")


@pytest.mark.parametrize("backend", BACKENDS)
@pytest.mark.parametrize("store", STORES, ids=lambda s: s.key)
def test_trio_output_matches_goldens(
    store, backend, vcfclick_home, monkeypatch, tmp_path
):
    root = GOLDENS / backend / store.key
    if UPDATE and root.exists():
        shutil.rmtree(root)
    _build(store, vcfclick_home, backend, monkeypatch, tmp_path)

    failures: list[str] = []
    for run in store.runs:
        r = run_cli(vcfclick_home, backend, "db", "trio", "t", *run)
        stem = _golden_name(run)
        ext = "json" if "json" in run else "txt"
        _check(root / f"{stem}.{ext}", _normalise(run, r.stdout), stem, failures)
        _check(
            root / f"{stem}.err", _engine_quiet(r.stderr), f"{stem} (stderr)", failures
        )
    assert not failures, "\n\n".join(failures)


def test_goldens_have_no_orphans():
    """Every golden file belongs to a pinned run (a renamed run must not
    leave a stale golden that silently stops being checked)."""
    if UPDATE:
        pytest.skip("regenerating")
    expected: dict[str, set[str]] = {}
    for store in STORES:
        names = set()
        for run in store.runs:
            stem = _golden_name(run)
            for ext in ("json", "txt", "err"):
                names |= {f"{stem}.{ext}", f"{stem}.{ext}.gz"}
        expected[store.key] = names
    for backend in BACKENDS:
        for d in (GOLDENS / backend).iterdir():
            assert d.name in expected, f"golden dir for unknown store: {d}"
            stray = {p.name for p in d.iterdir()} - expected[d.name]
            assert not stray, f"orphan goldens in {d}: {sorted(stray)}"


def _round_af(text: str) -> str:
    """Round the AF values in a golden to 6 places. This is the ONE known
    backend divergence: info_AF is Float32, which chDB prints in its
    shortest form (0.01) and DuckDB expands as a double
    (0.009999999776482582). Old output stays byte-stable per backend, so the
    goldens differ there; everything else must match exactly."""
    return re.sub(
        r'("af": |AF=)(-?[0-9.]+(?:e-?[0-9]+)?)',
        lambda m: m.group(1) + repr(round(float(m.group(2)), 6)),
        text,
    )


def test_goldens_agree_across_backends():
    """The two backends print the same trio documents, up to the Float32
    AF representation (see _round_af). Any other difference is a real
    divergence and is never papered over by a per-backend golden."""
    if UPDATE:
        pytest.skip("regenerating")
    diffs = []
    for store in STORES:
        a, b = (GOLDENS / be / store.key for be in BACKENDS)
        names = {p.name for p in a.iterdir()} | {p.name for p in b.iterdir()}
        for name in sorted(names):
            ta, tb = (
                _read(a / name.removesuffix(".gz")),
                _read(b / name.removesuffix(".gz")),
            )
            if ta is None or tb is None or _round_af(ta) != _round_af(tb):
                diffs.append(f"{store.key}/{name}")
    assert not diffs, f"chDB and DuckDB goldens differ beyond AF rounding: {diffs}"


def test_no_golden_candidate_sits_on_an_af_gate_boundary():
    """The AF epsilon commit (1e-7 on `info_AF <= --max-af`) must not move
    a golden. That holds when no fixture site has AF strictly between a
    pinned --max-af and --max-af + 1e-7, or exactly 0.001 (the value whose
    Float32 round-up the epsilon exists for). AF == 0.01 is present in
    trio.vcf.gz; Float32(0.01) < 0.01, so it passes the gate both with and
    without the epsilon."""
    afs: set[float] = set()
    for store in STORES:
        if store.vcf.suffix == ".gz":
            import gzip

            text = gzip.open(store.vcf, "rt").read()
        else:
            text = store.vcf.read_text()
        for line in text.splitlines():
            if line.startswith("#"):
                continue
            m = re.search(r"(?:^|;)AF=([0-9.eE+-]+)", line.split("\t")[7])
            if m:
                afs.add(float(m.group(1)))
    assert afs, "fixtures carry INFO/AF"
    for max_af in (0.01, 0.001, 0.3):
        assert not [a for a in afs if max_af < a <= max_af + 1e-7], max_af
    assert 0.001 not in afs


def test_engine_noise_is_the_only_thing_stripped_from_stderr():
    note = "  note: this database has no stored hom-reference calls\n"
    noisy = "Default dirty pages decay period: 5000ms\n" + note
    assert _engine_quiet(noisy) == note
    assert _engine_quiet(note) == note
    # Only a whole line in that exact shape is engine noise.
    kept = "warning: Default dirty pages decay period: 5000ms\n"
    assert _engine_quiet(kept) == kept
