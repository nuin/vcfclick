"""The scope of a run: which genes and regions it looks at.

    vcfclick db trio fam1 --proband C1 --genes C9orf72,C19orf12 --region chr7:1-5000
    vcfclick db trio fam1 --proband C1 --gene-file panel.txt --exclude-region chr6:28510120-33480577

One module for every command that takes a gene list or a region (`db trio`
today; `db segregate` and `db case` take the same flags from the registry). A
`Scope` is built once per run from the resolved values (flag > preset >
default, cli/filter_registry.py) and is then pure data: intervals per contig,
in the annotation store's contig names, and the echo that says what was asked
and what it resolved to.

GENES. A gene list is resolved in ONE batched lookup on one connection of the
annotation store (`annotations.db.resolve_gene_tokens`): a current symbol in
any case (`c9orf72` is `C9orf72`), an Ensembl gene id (version ignored) and,
after `vcfclick annotations load-hgnc`, a previous symbol, an alias, an HGNC id
or an Entrez id. A token that matches nothing, or two genes at once, is
unresolved. A list of 20 tokens or fewer with any unresolved token, or a longer
one with more than 5% unresolved, is an ERROR (`--allow-unresolved-genes` runs
anyway): a panel that silently lost genes is the narrower search a clinician
did not ask for. A gene is its GENCODE span start..end (introns and UTRs, so a
splice-region variant is inside); promoters need `--gene-flank`. The store is
GRCh38 with one locus per symbol (a PAR gene is matched on chrX only); a
GRCh37 database cannot take a gene list.

REGIONS. `chrN:FROM-TO`, 1-based and inclusive, or a BED file (0-based
half-open, converted when read). Genes and regions together: a variant in
either is in scope. Exclusions (`--exclude-genes`, `--exclude-region`) are
subtracted.

APPLICATION, three layers, so no model can ignore a scope:

  (a) SQL pushdown. Every model builder takes `scope=None` and puts
      `scope_sql(alias, scope)` in its WHERE, contig spellings from
      storage.contigs and numbers through int(). Up to SCOPE_SQL_EXACT_MAX
      merged intervals it is exact: `(chrom IN (aliases) AND (pos BETWEEN a
      AND b OR ...)) OR ...`. Up to SCOPE_SQL_MAX_INTERVALS it is a cheap
      superset: the position's bucket (`intDiv(pos, B)`) in the set of
      buckets the intervals touch (an exclusion: the buckets it covers
      whole). Above that the SQL carries no scope. With no scope the SQL text
      is byte-identical to before. (Measured on chDB, 500 intervals as
      BETWEENs over a 200,000-site trio cost 0.5 s per statement, the bucket
      set 0.07 s against 0.1 s for no predicate: the engine evaluates every
      BETWEEN on every row of each joined side.)
  (b) A Python backstop: TrioAnalysis re-applies the scope to every row it
      returns (bisect over the merged intervals), so a builder that forgot to
      push it down is still right.
  (c) Counts come from the same rows (or, when the SQL is exact, the same
      predicate), so a count is never wider than its listing.

COMPOUND HET. Exclusions act BEFORE the gene grouping. Inclusions act AFTER it:
a gene is kept when it is a listed gene, or when at least one of its variants
lies in a listed region; each variant then carries `in_scope`, so a partner
outside a capture edge is flagged, not lost. A variant inside a listed and an
unlisted gene is reported under the listed one only.

`SCOPE_AWARE_BUILDERS` lists, per model, the SQL its builder emits under a
scope; trio models are added by `cli.trio_output.register_model`, and
`tests/test_scope_models_drift.py` fails for a model whose builder does not
carry the predicate or whose rows outside a planted scope come back.
"""

from __future__ import annotations

import hashlib
import json
from bisect import bisect_right
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

from cli import filter_registry as fr
from cli.warning_codes import warning
from storage.contigs import canonical_chrom, chrom_aliases

# Merged intervals up to which the scope is pushed into the SQL at all (as a
# bucket superset above SCOPE_SQL_EXACT_MAX); above it only the Python backstop
# applies it. Benchmarked in the timing job (tests/test_trio_timing.py). Both
# are read at call time, so a test can change them.
SCOPE_SQL_MAX_INTERVALS = 500
# Up to this many merged intervals the pushed predicate is exact BETWEENs (a
# count can be a SQL count); above it, a bucket set.
SCOPE_SQL_EXACT_MAX = 64
# The bucket set of a superset predicate: the smallest power-of-two bucket
# (from MIN_BUCKET_BP) that keeps it at most MAX_BUCKETS entries.
MIN_BUCKET_BP = 256
MAX_BUCKETS = 10_000
# The unresolved-genes rule: a short list may lose no gene, a long one at most
# this share (both unvalidated starting points, like the built-in presets).
UNRESOLVED_SMALL_LIST = 20
UNRESOLVED_MAX_FRACTION = 0.05
# A list longer than this is echoed as its count and hash.
ECHO_LIST_MAX = 500
GENE_BUILD = "GRCh38"  # the annotation store's assembly
SPAN_NOTE = (
    "gene spans only (start..end, introns and UTRs); use --gene-flank for promoters"
)
# PAR1 and PAR2 of GRCh38 chrX: a gene there has one locus in the store (chrX)
_PAR_GRCH38 = ((10_001, 2_781_479), (155_701_383, 156_030_895))

SCOPE_KEYS: tuple[str, ...] = (
    "genes",
    "exclude_genes",
    "gene_flank",
    "regions",
    "exclude_regions",
    "genes_label",
    "regions_build",
)

NEEDS_GENES = (
    "a gene list needs gene annotations: run `vcfclick annotations load` "
    "(GENCODE genes) first"
)


class ScopeError(Exception):
    """A scope the run must refuse (unresolved genes, no gene table, a build
    the gene coordinates do not fit). The message names the remedy."""


# ───────────────────────────── what was asked ─────────────────────────────


@dataclass(frozen=True)
class ScopeRequest:
    """The resolved scope values of one run (flag > preset > default)."""

    genes: tuple[str, ...] = ()
    exclude_genes: tuple[str, ...] = ()
    gene_flank: int = 0
    regions: tuple[str, ...] = ()
    exclude_regions: tuple[str, ...] = ()
    genes_label: str | None = None
    regions_build: str | None = None
    allow_unresolved: bool = False

    @property
    def any(self) -> bool:
        return bool(
            self.genes or self.exclude_genes or self.regions or self.exclude_regions
        )

    @classmethod
    def from_values(
        cls, values: Mapping[str, Any], *, allow_unresolved: bool = False
    ) -> ScopeRequest:
        """From {file key: value}. Raises ScopeError for a token that is not a
        gene or a region (a typed flag is checked by the registry's rules, the
        same as a preset value)."""
        out: dict[str, Any] = {}
        for key in SCOPE_KEYS:
            value = values.get(key)
            if value is None:
                continue
            spec = fr.REGISTRY[key]
            if spec.is_list and not value:
                continue
            checked, why = fr.check_value(spec, value)
            if why:
                raise ScopeError(f"{spec.flag or key}: {why}")
            out[key] = tuple(checked) if spec.is_list else checked
        return cls(**out, allow_unresolved=allow_unresolved)


# ───────────────────────────── intervals ─────────────────────────────


class IntervalSet:
    """Closed 1-based intervals, merged per contig under the annotation store's
    name (storage.contigs.canonical_chrom), so `7:100-200` and `chr7:150-300`
    are one interval and a row on `7` or `chr7` is looked up the same way."""

    def __init__(self, intervals: Iterable[tuple[str, int, int]] = ()) -> None:
        by: dict[str, list[tuple[int, int]]] = {}
        for chrom, lo, hi in intervals:
            by.setdefault(canonical_chrom(chrom), []).append((int(lo), int(hi)))
        self._by: dict[str, list[tuple[int, int]]] = {}
        for chrom in sorted(by):
            merged: list[tuple[int, int]] = []
            for lo, hi in sorted(by[chrom]):
                if merged and lo <= merged[-1][1] + 1:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
                else:
                    merged.append((lo, hi))
            self._by[chrom] = merged
        self._starts = {c: [lo for lo, _ in ivs] for c, ivs in self._by.items()}

    def __len__(self) -> int:
        return sum(len(v) for v in self._by.values())

    def __bool__(self) -> bool:
        return bool(self._by)

    def items(self) -> list[tuple[str, int, int]]:
        return [(c, lo, hi) for c, ivs in self._by.items() for lo, hi in ivs]

    def contains(self, chrom: str, pos: int) -> bool:
        c = canonical_chrom(chrom)
        ivs = self._by.get(c)
        if not ivs:
            return False
        i = bisect_right(self._starts[c], int(pos)) - 1
        return i >= 0 and ivs[i][1] >= int(pos)

    def sql(self, alias: str) -> str:
        """Exact: `((a.chrom IN (...) AND (a.pos BETWEEN lo AND hi OR ...)) OR
        ...)`."""
        from storage.sql import in_intervals_sql

        return self._per_contig(
            alias,
            {c: in_intervals_sql(f"{alias}.pos", ivs) for c, ivs in self._by.items()},
        )

    def bucket_size(self) -> int:
        """The smallest power of two from MIN_BUCKET_BP whose bucket set has at
        most MAX_BUCKETS entries."""
        size = MIN_BUCKET_BP
        while self._bucket_count(size) > MAX_BUCKETS:
            size *= 2
        return size

    def _bucket_count(self, size: int) -> int:
        return sum(hi // size - lo // size + 1 for _, lo, hi in self.items())

    def buckets_sql(self, alias: str, *, whole: bool = False) -> str:
        """The bucket form: `a.pos`'s bucket (`intDiv(pos, B)`) in the buckets
        the intervals touch, a superset of the exact predicate; `whole`: only
        the buckets they cover entirely, a subset (what an exclusion needs, so
        that NOT of it never removes a row outside the exclusion)."""
        from storage.sql import int_div_sql

        size = self.bucket_size()
        bucket = int_div_sql(f"{alias}.pos", size)
        by: dict[str, str] = {}
        for chrom, ivs in self._by.items():
            ids: set[int] = set()
            for lo, hi in ivs:
                first, last = lo // size, hi // size
                if whole:
                    first += lo % size != 0
                    last -= (hi + 1) % size != 0
                ids.update(range(first, last + 1))
            if ids:
                listed = ", ".join(str(i) for i in sorted(ids))
                by[chrom] = f"{bucket} IN ({listed})"
        return self._per_contig(alias, by)

    @staticmethod
    def _per_contig(alias: str, by: Mapping[str, str]) -> str:
        from storage.sql import chrom_in_sql

        parts = [
            f"({chrom_in_sql(f'{alias}.chrom', chrom_aliases(chrom))} AND {pred})"
            for chrom, pred in by.items()
        ]
        return "(" + " OR ".join(parts) + ")" if parts else "(1 = 0)"


# ───────────────────────────── the scope ─────────────────────────────


@dataclass(frozen=True)
class Scope:
    """Where a run looks. `include` None: everywhere (only exclusions)."""

    include: IntervalSet | None
    exclude: IntervalSet
    # the include intervals that came from regions, and the upper-cased listed
    # symbols: the compound-het rule keeps a gene by either
    regions_include: IntervalSet | None = None
    gene_symbols: frozenset[str] | None = None
    echo: dict = field(default_factory=dict, compare=False)
    warnings: tuple[dict, ...] = ()
    notes: tuple[str, ...] = ()

    def included(self, chrom: str, pos: int) -> bool:
        return self.include is None or self.include.contains(chrom, pos)

    def excluded(self, chrom: str, pos: int) -> bool:
        return bool(self.exclude) and self.exclude.contains(chrom, pos)

    def contains(self, chrom: str, pos: int) -> bool:
        return self.included(chrom, pos) and not self.excluded(chrom, pos)

    @property
    def interval_count(self) -> int:
        return len(self.include or ()) + len(self.exclude)

    @property
    def sql_mode(self) -> str:
        """How the scope reaches the SQL: `sql` (exact, up to
        SCOPE_SQL_EXACT_MAX intervals), `sql_buckets` (a superset, up to
        SCOPE_SQL_MAX_INTERVALS) or `python` (not at all). The Python
        backstop runs in every mode."""
        n = self.interval_count
        if n <= min(SCOPE_SQL_EXACT_MAX, SCOPE_SQL_MAX_INTERVALS):
            return "sql"
        if n <= SCOPE_SQL_MAX_INTERVALS:
            return "sql_buckets"
        return "python"

    @property
    def pushdown(self) -> bool:
        """Is the scope in the SQL at all?"""
        return self.sql_mode != "python"

    @property
    def exact_sql(self) -> bool:
        """Is the SQL predicate the scope itself (so a SQL count is right)?"""
        return self.sql_mode == "sql"

    @property
    def flags_variants(self) -> bool:
        """Do compound-het variants carry `in_scope`? Only when something is
        included: with exclusions alone every reported variant is in scope."""
        return self.include is not None

    def keep_gene(self, symbol: str, variants: Iterable[Sequence]) -> bool:
        """The compound-het rule, after grouping: a listed gene, or a gene with
        at least one variant (chrom, pos, ...) in a listed region."""
        if self.include is None:
            return True
        if self.gene_symbols is not None and symbol.upper() in self.gene_symbols:
            return True
        return self.regions_include is not None and any(
            self.regions_include.contains(v[0], int(v[1])) for v in variants
        )


def scope_sql(alias: str, scope: Scope | None, *, part: str = "all") -> str | None:
    """The scope as a SQL boolean on `alias.chrom` / `alias.pos`, or None (no
    scope, or more intervals than SCOPE_SQL_MAX_INTERVALS, or nothing in the
    part asked for). `part="exclude"`: the exclusions only, which is what the
    compound-het candidate query takes (inclusions act after the grouping).
    In `sql_buckets` mode it keeps a superset of the scope's rows, which the
    Python backstop then trims."""
    if scope is None or not scope.pushdown:
        return None
    exact = scope.exact_sql
    parts = []
    if part == "all" and scope.include is not None:
        inc = scope.include
        parts.append(inc.sql(alias) if exact else inc.buckets_sql(alias))
    if scope.exclude:
        exc = scope.exclude
        shown = exc.sql(alias) if exact else exc.buckets_sql(alias, whole=True)
        parts.append(f"NOT {shown}")
    return " AND ".join(parts) if parts else None


def filter_rows(rows: Sequence[Sequence], scope: Scope | None) -> list:
    """The rows (chrom, pos, ...) in scope: the Python backstop."""
    if scope is None:
        return list(rows)
    return [r for r in rows if scope.contains(r[0], int(r[1]))]


def drop_excluded(rows: Sequence[Sequence], scope: Scope | None) -> list:
    """The rows (chrom, pos, ...) not in an exclusion (before a gene grouping)."""
    if scope is None or not scope.exclude:
        return list(rows)
    return [r for r in rows if not scope.excluded(r[0], int(r[1]))]


def genes_in_scope(genes: dict, scope: Scope | None) -> tuple[dict, int]:
    """({gene: {paternal, maternal}} the scope keeps, how many it dropped)."""
    if scope is None or scope.include is None:
        return genes, 0
    kept = {
        sym: e
        for sym, e in genes.items()
        if scope.keep_gene(sym, [*e["paternal"], *e["maternal"]])
    }
    return kept, len(genes) - len(kept)


# ───────────────────────────── building one ─────────────────────────────


def _hash(tokens: Iterable[str], *, upper: bool) -> str:
    items = sorted({t.upper() if upper else t for t in tokens})
    return hashlib.sha256(json.dumps(items).encode()).hexdigest()[:12]


def compact_list(values: Sequence[str], *, upper: bool = False) -> Any:
    """A list as the echo shows it: itself, or {count, sha256_12} above
    ECHO_LIST_MAX entries."""
    if len(values) <= ECHO_LIST_MAX:
        return list(values)
    return {"count": len(values), "sha256_12": _hash(values, upper=upper)}


@dataclass
class _Resolution:
    genes: list = field(default_factory=list)  # GeneRange, one per symbol
    resolved_by: dict = field(default_factory=dict)
    unresolved: list = field(default_factory=list)
    ambiguous: list = field(default_factory=list)
    resolved: int = 0  # tokens resolved


def _resolve_list(tokens: Sequence[str], found: Mapping) -> _Resolution:
    from annotations.db import gene_token_key

    out = _Resolution()
    seen: set[str] = set()
    for token in tokens:
        matches = found.get(gene_token_key(token), [])
        if not matches:
            out.unresolved.append(token)
            continue
        best = matches[0][0]
        symbols = {g.gene_symbol for kind, g in matches if kind == best}
        if len(symbols) > 1:
            out.ambiguous.append(token)
            continue
        gene = next(g for kind, g in matches if kind == best)
        out.resolved += 1
        out.resolved_by[best] = out.resolved_by.get(best, 0) + 1
        if gene.gene_symbol not in seen:
            seen.add(gene.gene_symbol)
            out.genes.append(gene)
    return out


def unresolved_is_error(requested: int, missed: int) -> bool:
    """The rule: a list of UNRESOLVED_SMALL_LIST tokens or fewer may miss none;
    a longer one at most UNRESOLVED_MAX_FRACTION of them."""
    if not missed:
        return False
    if requested <= UNRESOLVED_SMALL_LIST:
        return True
    return missed / requested > UNRESOLVED_MAX_FRACTION


def _shown(tokens: Sequence[str], n: int = 10) -> str:
    head = ", ".join(t[:40] for t in tokens[:n])
    return head + (f" and {len(tokens) - n} more" if len(tokens) > n else "")


def _par_note(gene) -> str | None:
    if canonical_chrom(gene.chrom) not in ("chrX", "chrY"):
        return None
    if any(gene.start_pos <= hi and gene.end_pos >= lo for lo, hi in _PAR_GRCH38):
        return (
            f"{gene.gene_symbol} is a pseudo-autosomal gene; the store keeps one "
            f"locus per symbol, so it is matched on {gene.chrom} only"
        )
    return None


def _gene_echo(
    flag: str,
    tokens: Sequence[str],
    res: _Resolution,
    *,
    build_verified: bool,
    aliases_loaded: bool,
) -> dict:
    return {
        "requested": len(tokens),
        "resolved": res.resolved,
        "genes": len(res.genes),
        "resolved_by": dict(res.resolved_by),
        "unresolved": compact_list(res.unresolved),
        "ambiguous": compact_list(res.ambiguous),
        "sha256_12": _hash(tokens, upper=True),
        "build": GENE_BUILD,
        "build_verified": build_verified,
        "aliases_loaded": aliases_loaded,
    }


def _check_unresolved(flag: str, tokens, res: _Resolution, allow: bool) -> dict | None:
    missed = res.unresolved + res.ambiguous
    if not missed:
        return None
    what = f"{len(missed)} of {len(tokens)} {flag} tokens did not resolve: {_shown(missed)}"
    if res.ambiguous:
        what += f" ({len(res.ambiguous)} name more than one gene)"
    if unresolved_is_error(len(tokens), len(missed)) and not allow:
        raise ScopeError(
            f"{what}. Check the spelling; previous symbols, aliases, HGNC and "
            "Entrez ids need `vcfclick annotations load-hgnc`; "
            "--allow-unresolved-genes runs without them"
        )
    key = "genes" if flag == "--genes" else "exclude_genes"
    return warning("genes_unresolved", key, what + "; they are not in the scope")


def build_scope(request: ScopeRequest, *, build: str | None = None) -> Scope | None:
    """The Scope of a run, or None when nothing is asked. `build` is the
    database's assembly when known (GRCh37 / GRCh38), else None or
    `unknown`. Raises ScopeError."""
    if not request.any:
        return None
    warnings: list[dict] = []
    notes: list[str] = []
    include: list[tuple[str, int, int]] = []
    exclude: list[tuple[str, int, int]] = []
    region_ivs: list[tuple[str, int, int]] = []
    echo: dict[str, Any] = {}
    gene_symbols: frozenset[str] | None = None
    known_build = build if build in ("GRCh37", "GRCh38") else None

    if request.genes or request.exclude_genes:
        from annotations.db import (
            gene_count,
            get_connection,
            has_gene_aliases,
            resolve_gene_tokens,
        )

        if known_build == "GRCh37":
            raise ScopeError(
                "a gene list needs a GRCh38 database: the gene coordinates in the "
                "annotation store are GRCh38, and this database is GRCh37 "
                "(use --region with GRCh37 coordinates instead)"
            )
        conn = get_connection()
        try:
            if gene_count(conn) == 0:
                raise ScopeError(NEEDS_GENES)
            aliases_loaded = has_gene_aliases(conn)
            found = resolve_gene_tokens(
                [*request.genes, *request.exclude_genes], conn=conn
            )
        finally:
            conn.close()
        verified = known_build == GENE_BUILD
        if not verified:
            notes.append(
                "gene coordinates are GRCh38 (the annotation store); this "
                "database's build is not verified (pass --build GRCh38)"
            )
        flank = int(request.gene_flank)
        for flag, tokens, target in (
            ("--genes", request.genes, include),
            ("--exclude-genes", request.exclude_genes, exclude),
        ):
            if not tokens:
                continue
            res = _resolve_list(tokens, found)
            w = _check_unresolved(flag, tokens, res, request.allow_unresolved)
            if w:
                warnings.append(w)
            for gene in res.genes:
                target.append(
                    (
                        gene.chrom,
                        max(1, int(gene.start_pos) - flank),
                        int(gene.end_pos) + flank,
                    )
                )
                note = _par_note(gene)
                if note:
                    notes.append(note)
            block = _gene_echo(
                flag,
                tokens,
                res,
                build_verified=verified,
                aliases_loaded=aliases_loaded,
            )
            if flag == "--genes":
                block.update(flank=flank, label=request.genes_label)
                gene_symbols = frozenset(g.gene_symbol.upper() for g in res.genes)
                echo["genes"] = block
            else:
                echo["exclude_genes"] = block
        if "genes" in echo:
            echo["genes"]["excluded"] = echo.get("exclude_genes", {}).get("genes", 0)
            echo["note"] = SPAN_NOTE

    for text in request.regions:
        contig, lo, hi = fr.parse_region(text)
        region_ivs.append((contig, lo, hi))
    for text in request.exclude_regions:
        contig, lo, hi = fr.parse_region(text)
        exclude.append((contig, lo, hi))
    include += region_ivs
    if request.regions or request.exclude_regions:
        echo["regions"] = compact_list(request.regions)
        echo["exclude_regions"] = compact_list(request.exclude_regions)
        echo["regions_build"] = request.regions_build
        if (
            request.regions_build
            and known_build
            and request.regions_build != known_build
        ):
            warnings.append(
                warning(
                    "scope_build_mismatch",
                    "regions_build",
                    f"the regions are written in {request.regions_build}, and this "
                    f"database is {known_build}: they name other loci",
                )
            )

    has_include = bool(request.genes or request.regions)
    scope = Scope(
        include=IntervalSet(include) if has_include else None,
        exclude=IntervalSet(exclude),
        regions_include=IntervalSet(region_ivs) if request.regions else None,
        gene_symbols=gene_symbols if request.genes else None,
        echo=echo,
        warnings=tuple(warnings),
        notes=tuple(dict.fromkeys(notes)),
    )
    echo["intervals"] = scope.interval_count
    return scope


# ───────────────────────────── what a run says ─────────────────────────────


def scope_document(scope: Scope, stats: Mapping[str, Any]) -> dict:
    """The JSON `scope` object of a finished run: what was asked and resolved,
    how it was applied, and what the compound-het rule dropped and flagged."""
    doc = json.loads(json.dumps(scope.echo))
    doc["applied"] = scope.sql_mode
    doc["comphet_genes_dropped_out_of_scope"] = sum(
        stats.get("scope_genes_dropped", {}).values()
    )
    doc["n_outside"] = sum(stats.get("scope_outside", {}).values())
    doc["notes"] = list(scope.notes)
    return doc


def _plural(n: int, word: str) -> str:
    return f"{n} {word}" + ("" if n == 1 else "s")


def scope_line(scope: Scope) -> str:
    """`  scope: 3 genes (1 not found: FAM123X), 2 regions, 1 excluded region`."""
    echo = scope.echo
    bits = []
    genes = echo.get("genes")
    if genes:
        bit = _plural(genes["genes"], "gene")
        missed = genes["unresolved"] if isinstance(genes["unresolved"], list) else []
        missed = missed + (
            genes["ambiguous"] if isinstance(genes["ambiguous"], list) else []
        )
        if missed:
            bit += f" ({len(missed)} not found: {_shown(missed, 5)})"
        if genes.get("flank"):
            bit += f" +/- {genes['flank']} bp"
        if genes.get("label"):
            bit += f" [{genes['label']}]"
        bits.append(bit)
    excluded = echo.get("exclude_genes")
    if excluded:
        bits.append(_plural(excluded["genes"], "excluded gene"))
    regions = echo.get("regions")
    if regions:
        n = regions["count"] if isinstance(regions, dict) else len(regions)
        bits.append(_plural(n, "region"))
    ex_regions = echo.get("exclude_regions")
    if ex_regions:
        n = ex_regions["count"] if isinstance(ex_regions, dict) else len(ex_regions)
        bits.append(_plural(n, "excluded region"))
    line = "  scope: " + ", ".join(bits)
    if genes:
        line += f"; {SPAN_NOTE}"
    return line


def text_lines(scope: Scope) -> list[str]:
    """The lines a text run prints under its header: the scope line and one
    `WARNING:` line per scope warning (stdout)."""
    return [scope_line(scope)] + [
        f"  WARNING: {w['message']}" for w in scope.warnings if w["level"] == "warning"
    ]


# ───────────────────────────── builders ─────────────────────────────

# A builder returns the SQL statements a model runs under a scope, for a trio
# and its gates. The drift test calls each one with a planted scope and
# requires the scope predicate in every statement it returns.
ScopeBuilder = Callable[[Any, Any, Scope], list[str]]
SCOPE_AWARE_BUILDERS: dict[str, ScopeBuilder] = {}


def register_scope_builder(name: str, builder: ScopeBuilder) -> ScopeBuilder:
    """Declare that model (or phase) `name` takes a scope. Trio models are
    registered by `cli.trio_output.register_model`; another command registers
    its own builders (db segregate, the phases of db case)."""
    if name in SCOPE_AWARE_BUILDERS:
        raise ValueError(f"scope builder {name!r} is already registered")
    SCOPE_AWARE_BUILDERS[name] = builder
    return builder
