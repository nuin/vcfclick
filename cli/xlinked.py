"""The X-linked model of `db trio`: patterns, SQL, context and counters.

A haploid male call is stored like a diploid one (a lone `1` as gt 2, a
lone `0` as gt 0; schema/02_genotypes.sql), so hemizygosity is never read
from the data: it comes from the chromosome, the position (cli/xchrom.py)
and the sex (cli/sex.py). On sure non-PAR X (NPX):

  * a son's X is his mother's: his patterns read the proband and the
    mother only. The father is annotated (`father_flag`), never a filter.
  * a daughter has her father's X: the father must be hemizygous (stored
    0 or 2). A het father there is an artefact; the site is excluded and
    counted (`father_het_calls`).
  * a male het (gt 1, including a partial `./1`) is an artefact; excluded
    and counted (`male_het_calls`).

`PATTERNS` is the one table: the SQL is generated from it and `classify` is
its Python twin (the tests hold them equal on every candidate). Every
genotype combination of a resolved proband is either a pattern, an
artefact, `consistent_unpatterned` (Mendelian, but no model claims it, such
as a hom-alt mother) or `mendel_error_other`; the last two are counted.

Rarity on X uses min(--max-af, --max-af-x): a hemizygous male is affected
at the allele frequency q, not q squared.

No import of cli.main.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import NamedTuple

from cli.sex import SexEvidence, TrioSexes, chrx_counts, infer_sex, transmission
from cli.sex import resolve_trio_sexes as resolve_trio_sexes
from cli.warning_codes import warning
from cli.xchrom import (
    AMBIGUOUS_BAND,
    UNKNOWN,
    UNSUPPORTED_SEX_CONTIGS,
    X_NAMES,
    band_sql,
    npx_sql,
    resolve_build,
)


class _Any:
    def __repr__(self) -> str:
        return "ANY"


ANY = _Any()  # a role the pattern does not constrain

# The modes a candidate can carry, and those that apply per proband sex.
MODES: tuple[str, ...] = ("recessive", "dominant", "female_het", "denovo", "comphet")
APPLICABLE: dict[str, tuple[str, ...]] = {
    "male": ("recessive", "denovo"),
    "female": MODES,
}
# Modes that rest on a parent's stored 0/0 (a sparse ingest cannot prove it).
NEEDS_REFERENCE: frozenset[str] = frozenset(
    {"dominant", "female_het", "denovo", "comphet"}
)

# The counters of `context.excluded`, in order.
EXCLUDED_KEYS: tuple[str, ...] = (
    "male_het_calls",
    "father_het_calls",
    "ambiguous_band",
    "mendel_consistent_unpatterned",
    "mendel_error_other",
)


@dataclass(frozen=True)
class Pattern:
    """One genotype pattern: stored gt per role (a tuple of allowed values,
    None meaning no stored row, or ANY), the mode it reports, the parents
    whose 0/0 must carry at most --max-parent-alt-reads alt reads (`cap`) and
    the parents that transmit the allele (`transmits`)."""

    name: str
    sex: str
    mode: str
    p: tuple | _Any
    f: tuple | _Any
    m: tuple | _Any
    cap: tuple[str, ...] = ()
    transmits: tuple[str, ...] = ()


PATTERNS: tuple[Pattern, ...] = (
    Pattern("xlr_son", "male", "recessive", (2,), ANY, (1,), (), ("m",)),
    Pattern("x_denovo_son", "male", "denovo", (2,), ANY, (0,), ("m",)),
    Pattern(
        "xlr_hom_daughter", "female", "recessive", (2,), (2,), (1,), (), ("f", "m")
    ),
    Pattern("xhet_paternal", "female", "female_het", (1,), (2,), (0,), (), ("f",)),
    Pattern("xhet_maternal", "female", "female_het", (1,), (0,), (1, 2), (), ("m",)),
    Pattern("x_denovo_daughter", "female", "denovo", (1,), (0,), (0,), ("f", "m")),
    Pattern("x_denovo_daughter_mat", "female", "denovo", (2,), (2,), (0,), ("m",)),
    # the mirror: the maternal ALT inherited, the paternal one new
    Pattern("x_denovo_daughter_pat", "female", "denovo", (2,), (0,), (1, 2), ("f",)),
)
BY_NAME: dict[str, Pattern] = {p.name: p for p in PATTERNS}


def patterns_for(sex: str) -> tuple[Pattern, ...]:
    return tuple(p for p in PATTERNS if p.sex == sex)


def _matches(cell, value) -> bool:
    return cell is ANY or value in cell


def classify(sex: str | None, p, f, m) -> str | None:
    """The pattern name, `male_het`, `father_het`,
    `consistent_unpatterned`, `mendel_error_other`, or None (proband not
    carrying, sex unknown, or a parent with no stored call) for one NPX site.
    The Python twin of the SQL; `PATTERNS` drives both."""
    if sex not in APPLICABLE or p in (None, 0):
        return None
    if sex == "male":
        if p == 1:
            return "male_het"
    elif f == 1:
        return "father_het"
    for pat in patterns_for(sex):
        if _matches(pat.p, p) and _matches(pat.f, f) and _matches(pat.m, m):
            return pat.name
    if sex == "male":
        return "consistent_unpatterned" if m == 2 else None
    if f is None or m is None:
        return None
    if (p, f, m) in ((1, 2, 1), (2, 2, 2)):
        return "consistent_unpatterned"
    return "mendel_error_other"


_EXCLUDED_CLASS = {
    "male_het": "male_het_calls",
    "father_het": "father_het_calls",
    "consistent_unpatterned": "mendel_consistent_unpatterned",
    "mendel_error_other": "mendel_error_other",
}


# ─────────────────────────────── SQL ───────────────────────────────

_FIELDS = ("gt", "gq", "dp", "ad_ref", "ad_alt")


def _in(col: str, cell) -> str:
    """`col IN (...)`, with `OR col IS NULL` when None is allowed; `1 = 1`
    for ANY."""
    if cell is ANY:
        return "1 = 1"
    values = [str(int(v)) for v in cell if v is not None]
    parts = [f"{col} IN ({', '.join(values)})"] if values else []
    if None in cell:
        parts.append(f"{col} IS NULL")
    return "(" + " OR ".join(parts) + ")"


def pattern_sql(pat: Pattern, cap: int) -> str:
    """SQL boolean on the pivot `t`: the pattern's genotypes and, for each
    capped parent, at most `cap` alt reads (missing AD passes)."""
    parts = [_in("t.p_gt", pat.p), _in("t.f_gt", pat.f), _in("t.m_gt", pat.m)]
    parts += [f"(t.{r}_ad_alt IS NULL OR t.{r}_ad_alt <= {int(cap)})" for r in pat.cap]
    return "(" + " AND ".join(parts) + ")"


def _quality(role: str, min_gq: int, min_dp: int, gates) -> str:
    def one(column: str, minimum: int, strict: bool) -> str:
        c = f"t.{role}_{column}"
        if strict and minimum > 0:
            return f"({c} IS NOT NULL AND {c} >= {int(minimum)})"
        return f"({c} IS NULL OR {c} >= {int(minimum)})"

    return (
        f"{one('gq', min_gq, gates.strict_for('min_gq'))} AND "
        f"{one('dp', min_dp, gates.strict_for('min_dp'))}"
    )


def _fraction(role: str) -> tuple[str, str]:
    """(alt fraction, `depth is zero`) of a role's allele depths."""
    from storage.sql import allele_fraction_sql, to_int64

    alt, ref = f"t.{role}_ad_alt", f"t.{role}_ad_ref"
    zero = f"({to_int64(alt)} + {to_int64(ref)}) = 0"
    return allele_fraction_sql(alt, ref), zero


def _het_ab(role: str, lo: float, hi: float) -> str:
    frac, zero = _fraction(role)
    return (
        f"(t.{role}_gt IS NULL OR t.{role}_gt != 1 OR t.{role}_ad_ref IS NULL "
        f"OR t.{role}_ad_alt IS NULL OR {zero} OR {frac} BETWEEN {lo!r} AND {hi!r})"
    )


def _hemi_ab(role: str, minimum: float) -> str:
    frac, zero = _fraction(role)
    return (
        f"(t.{role}_gt IS NULL OR t.{role}_gt != 2 OR t.{role}_ad_ref IS NULL "
        f"OR t.{role}_ad_alt IS NULL OR {zero} OR {frac} >= {minimum!r})"
    )


def _ad_present(role: str) -> str:
    return (
        f"(t.{role}_gt IS NULL OR t.{role}_gt != 1 OR "
        f"(t.{role}_ad_ref IS NOT NULL AND t.{role}_ad_alt IS NOT NULL))"
    )


def male_x_dp(gates) -> int:
    """The DP floor of a hemizygous call: --min-male-x-dp, else half --min-dp."""
    if gates.min_male_x_dp is not None:
        return int(gates.min_male_x_dp)
    return int(gates.min_dp) // 2


def x_max_af(gates) -> float:
    """The X-model rarity cap: min(--max-af, --max-af-x)."""
    return min(float(gates.max_af), float(gates.max_af_x))


def gate_sql(sex: str, gates) -> list[str]:
    """The quality, allele-balance, rarity and FILTER gates of the X model for
    a proband of `sex`, applied only to the roles its patterns read: a son's
    proband (male DP floor) and mother; a daughter's proband, father (male DP
    floor) and mother."""
    from cli.trio_analysis import AF_GATE_SKIPPED, PASS_ONLY_SQL
    from storage.sql import af_gate_sql

    hemi_dp = male_x_dp(gates)
    if sex == "male":
        roles = {"p": hemi_dp, "m": gates.min_dp}
    else:
        roles = {"p": gates.min_dp, "f": hemi_dp, "m": gates.min_dp}
    parts = [_quality(r, gates.min_gq, dp, gates) for r, dp in roles.items()]
    parts.append(_het_ab("p", gates.min_ab, gates.max_ab))
    parts.append(_het_ab("m", gates.min_ab, gates.max_ab))
    if sex == "male":
        parts.append(_hemi_ab("p", gates.min_hemi_ab))
    if gates.ab_requested and gates.strict_for("min_ab"):
        parts += [_ad_present(r) for r in roles]
    if gates.af_gate != AF_GATE_SKIPPED:
        parts.append(af_gate_sql("v.info_AF", x_max_af(gates)))
    if gates.pass_only:
        parts.append(PASS_ONLY_SQL)
    return parts


def _pivot(trio, build: str) -> str:
    from storage.sql import genotype_pivot_sql

    return genotype_pivot_sql(
        trio.ingest_id,
        {"p": trio.proband, "f": trio.father, "m": trio.mother},
        _FIELDS,
        where=npx_sql(build),
    )


def xlinked_sql(trio, gates, sex: str, build: str) -> str:
    """The candidate rows of the X model for a proband of `sex`: one row per
    NPX site, (chrom, pos, ref, alt, p_gt, f_gt, m_gt, af, pattern).

    A pivot, not joins: on chDB a LEFT JOIN against the non-Nullable `gt`
    reads an absent parent as 0, a proven reference; `max(CASE ...)` is
    NULL on both backends (storage.sql.genotype_pivot_sql)."""
    cap = gates.max_parent_alt_reads
    pats = patterns_for(sex)
    case = " ".join(f"WHEN {pattern_sql(p, cap)} THEN '{p.name}'" for p in pats)
    where = " AND ".join(
        [
            "t.p_gt IS NOT NULL",
            "(" + " OR ".join(pattern_sql(p, cap) for p in pats) + ")",
            *gate_sql(sex, gates),
        ]
    )
    return (
        "SELECT t.chrom, t.pos, t.ref, t.alt, t.p_gt, t.f_gt, t.m_gt, "
        f"v.info_AF AS af, CASE {case} END AS pattern "
        f"FROM ({_pivot(trio, build)}) t "
        "LEFT JOIN variants v ON v.ingest_id = t.ingest_id AND v.chrom = t.chrom "
        "AND v.pos = t.pos AND v.ref = t.ref AND v.alt = t.alt "
        f"WHERE {where} ORDER BY t.chrom, t.pos, t.ref, t.alt"
    )


def genotype_classes_sql(trio, build: str) -> str:
    """Every NPX genotype combination where the proband carries, ungated:
    (p_gt, f_gt, m_gt, rows). `classify` turns it into the counters."""
    from storage.sql import count_expr

    return (
        f"SELECT t.p_gt, t.f_gt, t.m_gt, {count_expr()} "
        f"FROM ({_pivot(trio, build)}) t WHERE t.p_gt > 0 "
        "GROUP BY t.p_gt, t.f_gt, t.m_gt"
    )


def x_guard_sql(build: str, *, father: bool, proband_male: bool) -> str | None:
    """The X guard of the classic models (proband `g`, father `f`): on sure
    NPX drop a het father and, for a male proband, a het proband, both
    physically impossible genotypes. None when neither applies."""
    impossible = (["f.gt = 1"] if father else []) + (
        ["g.gt = 1"] if proband_male else []
    )
    if not impossible:
        return None
    where = npx_sql(build, chrom="g.chrom", pos="g.pos")
    return f"NOT ({where} AND ({' OR '.join(impossible)}))"


def comphet_x_sql(trio, gates, build: str) -> str:
    """The female compound-het candidates on NPX, (chrom, pos, ref, alt, af,
    origin): a proband het from a hemizygous-ALT father with a 0/0 mother, or
    from a carrier mother with a 0/0 father (a het father cannot be tagged
    paternal), under the X model's own gates (`gate_sql("female", ...)`): the
    hemizygous father at the male DP floor, the het-AB band on the proband and
    on a het mother, the X rarity cap. Grouped by gene in Python
    (cli/trio_analysis._comphet_genes)."""
    paternal = "t.f_gt = 2 AND t.m_gt = 0"
    maternal = "t.f_gt = 0 AND t.m_gt > 0"
    where = " AND ".join(
        [f"t.p_gt = 1 AND (({paternal}) OR ({maternal}))", *gate_sql("female", gates)]
    )
    return (
        "SELECT t.chrom, t.pos, t.ref, t.alt, v.info_AF AS af, "
        f"CASE WHEN {paternal} THEN 'paternal' WHEN {maternal} THEN 'maternal' END "
        f"AS origin FROM ({_pivot(trio, build)}) t "
        "LEFT JOIN variants v ON v.ingest_id = t.ingest_id AND v.chrom = t.chrom "
        "AND v.pos = t.pos AND v.ref = t.ref AND v.alt = t.alt "
        f"WHERE {where} ORDER BY t.chrom, t.pos, t.ref, t.alt"
    )


# ───────────────────────────── context ─────────────────────────────


@dataclass
class XContext:
    """What the X model knows about one trio before it looks for candidates."""

    has_x: bool
    build: str | None = None
    build_source: str | None = None
    sexes: TrioSexes = field(default_factory=TrioSexes)
    blocked_reason: str | None = None
    warnings: list[dict] = field(default_factory=list)
    conflict: str | None = None  # the text of a blocking sex conflict
    affected: dict[str, bool] = field(default_factory=dict)  # role -> affected
    # the NPX intervals the guard and the patterns use (unknown under a
    # build conflict)
    region_build: str = UNKNOWN

    @property
    def proband_sex(self) -> str | None:
        return self.sexes.proband_sex


def _query(session, sql: str) -> list:
    return json.loads(session.query(sql, "JSONCompact").bytes().decode())["data"]


def _has_x(session, trio) -> bool:
    from storage.sql import chrom_in_sql, sql_quote_str, str_list

    sql = (
        "SELECT 1 FROM genotypes "
        f"WHERE ingest_id = {sql_quote_str(trio.ingest_id)} "
        f"AND {chrom_in_sql('chrom', X_NAMES)} "
        f"AND sample_id IN ({str_list(trio.parents + (trio.proband,))}) LIMIT 1"
    )
    return bool(_query(session, sql))


def _unsupported_contigs(session, ingest_id: str) -> list[str]:
    from storage.sql import chrom_in_sql, sql_quote_str

    sql = (
        "SELECT DISTINCT chrom FROM variants "
        f"WHERE ingest_id = {sql_quote_str(ingest_id)} "
        f"AND {chrom_in_sql('chrom', UNSUPPORTED_SEX_CONTIGS)} ORDER BY chrom"
    )
    return [r[0] for r in _query(session, sql)]


def _max_x_pos(session, ingest_id: str) -> int | None:
    from storage.sql import chrom_in_sql, sql_quote_str

    sql = (
        "SELECT max(pos) FROM variants "
        f"WHERE ingest_id = {sql_quote_str(ingest_id)} "
        f"AND {chrom_in_sql('chrom', X_NAMES)}"
    )
    data = _query(session, sql)
    return None if not data or data[0][0] is None else int(data[0][0])


def _pedigree(session, trio) -> dict[str, tuple[str | None, str | None]]:
    """{sample: (sex, affected)} from the pedigree for the trio's members."""
    from storage import table_exists
    from storage.sql import sql_quote_str, str_list

    if not table_exists("pedigree"):
        return {}
    sql = (
        "SELECT sample_id, sex, affected FROM pedigree "
        f"WHERE ingest_id = {sql_quote_str(trio.ingest_id)} "
        f"AND sample_id IN ({str_list(trio.parents + (trio.proband,))})"
    )
    return {s: (sex, aff) for s, sex, aff in _query(session, sql)}


def build_context(session, trio, build_option: str = "auto") -> XContext:
    """Resolve the build and the sexes of one trio (both parents declared)."""
    if not _has_x(session, trio):
        ctx = XContext(has_x=False, blocked_reason="no_x_data")
        others = _unsupported_contigs(session, trio.ingest_id)
        if others:
            ctx.warnings.append(
                warning(
                    "unsupported_contig",
                    None,
                    f"{', '.join(others)}: sex contigs vcfclick does not "
                    "interpret; only chrX (chrX, X, chr23, 23) is analysed",
                )
            )
        return ctx
    build, source, conflict = resolve_build(
        build_option, _max_x_pos(session, trio.ingest_id)
    )
    ctx = XContext(has_x=True, build=build, build_source=source)
    ctx.region_build = UNKNOWN if conflict else build
    ped = _pedigree(session, trio)
    members = (trio.proband, trio.father, trio.mother)
    counts = chrx_counts(session, trio.ingest_id, members, ctx.region_build)
    evidence = {}
    for s in members:
        declared, _aff = ped.get(s, (None, None))
        sex, frac = infer_sex(*counts.get(s, (0, 0)))
        evidence[s] = SexEvidence(
            declared if declared in ("male", "female") else None, sex, frac
        )
    ev = evidence[trio.proband]
    if ev.declared is None and ev.inferred in ("male", "female"):
        from cli.sex import transmission_call

        call, tfrac = transmission_call(
            *transmission(
                session, trio.ingest_id, trio.proband, trio.father, trio.mother,
                ctx.region_build,
            )
        )  # fmt: skip
        evidence[trio.proband] = SexEvidence(
            None, ev.inferred, ev.het_frac, call, tfrac
        )
    ctx.sexes = resolve_trio_sexes(trio.proband, trio.father, trio.mother, evidence)
    ctx.warnings = list(ctx.sexes.warnings)
    ctx.conflict = ctx.sexes.conflict
    ctx.affected = {
        "proband": ped.get(trio.proband, (None, None))[1] == "affected",
        "f": ped.get(trio.father, (None, None))[1] == "affected",
        "m": ped.get(trio.mother, (None, None))[1] == "affected",
    }
    if conflict:
        ctx.blocked_reason = "build_conflict"
        ctx.conflict = (
            "--build GRCh37, but chrX has positions beyond the GRCh37 chrX "
            "length (155,270,560)"
        )
    else:
        ctx.blocked_reason = ctx.sexes.blocked_reason
    if source == "unknown":
        n = band_count(session, trio, build)
        bands = " / ".join(f"{lo:,}-{hi:,}" for lo, hi in AMBIGUOUS_BAND[UNKNOWN])
        ctx.warnings.append(
            warning(
                "build_unknown",
                None,
                f"{n} proband alt calls in {bands} left unclassified; pass --build",
            )
        )
    return ctx


def band_count(session, trio, build: str) -> int:
    """Proband non-reference calls in the ambiguous bands of `build`."""
    from storage.sql import count_expr, sql_quote_str

    if not AMBIGUOUS_BAND[build]:
        return 0
    sql = (
        f"SELECT {count_expr()} FROM genotypes "
        f"WHERE ingest_id = {sql_quote_str(trio.ingest_id)} "
        f"AND sample_id = {sql_quote_str(trio.proband)} AND gt > 0 "
        f"AND {band_sql(build)}"
    )
    return int(_query(session, sql)[0][0])


def excluded_counts(session, trio, ctx: XContext) -> dict[str, int]:
    """`context.excluded`: the NPX combinations no pattern reports, by class,
    and the proband's calls in the ambiguous bands (unknown build only)."""
    out = dict.fromkeys(EXCLUDED_KEYS, 0)
    if not ctx.has_x or ctx.proband_sex is None:
        return out
    for p, f, m, n in _query(session, genotype_classes_sql(trio, ctx.region_build)):
        key = _EXCLUDED_CLASS.get(
            classify(ctx.proband_sex, _int(p), _int(f), _int(m)) or ""
        )
        if key:
            out[key] += int(n)
    if ctx.build_source == "unknown":
        out["ambiguous_band"] = band_count(session, trio, ctx.region_build)
    return out


def _int(v):
    return None if v is None else int(v)


# ─────────────────────────────── rows ───────────────────────────────


class XRow(NamedTuple):
    """One X-linked candidate. The first four fields are the variant key, so
    the gnomAD filter and the annotation columns apply as to any row."""

    chrom: str
    pos: int
    ref: str
    alt: str
    proband_gt: int
    father_gt: int | None
    mother_gt: int | None
    af: float | None
    pattern: str
    mode: str


def mode_of(pattern: str, affected: dict[str, bool]) -> str:
    """The pattern's mode, `female_het` relabelled `dominant` only when the
    pedigree marks the proband AND the transmitting parent affected."""
    pat = BY_NAME[pattern]
    if pat.mode == "female_het" and affected.get("proband"):
        if all(affected.get(r) for r in pat.transmits):
            return "dominant"
    return pat.mode


def decode_rows(data: list, affected: dict[str, bool]) -> list[XRow]:
    return [
        XRow(
            c,
            int(p),
            r,
            a,
            int(pg),
            _int(fg),
            _int(mg),
            af,
            pat,
            mode_of(pat, affected),
        )
        for c, p, r, a, pg, fg, mg, af, pat in data
    ]


def father_flag(row: XRow, sex: str | None) -> str | None:
    """A son's father, who cannot pass him an X: hemizygous ALT or a het
    artefact, shown and never filtered on."""
    if sex != "male":
        return None
    return {2: "hemizygous_alt", 1: "anomalous_het"}.get(row.father_gt)


def hemizygous_roles(sex: str | None) -> list[str]:
    return ["proband", "father"] if sex == "male" else ["father"]


def parent_affected(row: XRow, affected: dict[str, bool]) -> bool:
    return any(affected.get(r) for r in BY_NAME[row.pattern].transmits)
