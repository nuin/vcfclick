"""Site filtering, polarisation and the `db popgen` computations.

`prepare()` turns the raw per-site counts of one ingestion into the set
of retained sites with per-group (n, k, het, called) arrays and a report
of what was dropped and why. `summary`, `sfs`, `fst` and `windows` work
only on that prepared data, so every subcommand sees the same sites.

Chromosomes are fetched and filtered one at a time; only retained sites
are kept, in compact integer types, so peak memory is one chromosome's
raw counts plus the retained arrays.

Choices (documented in docs/POPGEN.md):

  * one ingestion at a time (sample identity is (ingest_id, sample_id),
    and two ingestions do not share a site list);
  * autosomes only (sex chromosomes, MT and unplaced/decoy contigs are
    excluded; including sex chromosomes needs per-sample ploidy and is
    refused in this version);
  * the retained site set is shared by all groups: a site failing the
    call-rate threshold in any group is dropped for all, so per-group
    statistics and F_ST are computed over the same sites;
  * sites whose per-group called counts cannot be derived exactly
    (partial ./1, haploid or polyploid calls) are dropped and reported;
  * π uses each site's own sample size; θ_W sums 1/a_{n_i} over
    segregating sites; Tajima's D and Fay & Wu's H are computed on the
    SFS projected to a common sample size (they need a fixed n).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

import numpy as np

from popgen import estimators as est
from popgen.ancestral import (
    ANCESTRAL_IS_ALT,
    ANCESTRAL_IS_REF,
    aa_code,
    orientation_from_codes,
)
from popgen.counts import (
    ALL_GROUP,
    UNLABELLED,
    ChromCounts,
    Region,
    fetch_chrom,
    schema_features,
)

# Chromosomes excluded by "autosomes only", compared on the part of the
# name before the first "_" (so hg38 alt/random/decoy contigs such as
# chrX_KI270880v1_alt follow their chromosome) without a `chr` prefix,
# case-insensitively: sex chromosomes of XY and ZW systems and their
# pseudo-autosomal pseudo-contigs, mitochondria, and unplaced or decoy
# sequence (chrUn_*, EBV, hs37d5, RefSeq NT_/NW_ scaffolds, GRCh37
# GL/NC unplaced contigs).
NON_AUTOSOMAL = frozenset(
    {"X", "Y", "XY", "W", "Z", "M", "MT", "PAR1", "PAR2", "UN", "EBV", "NT", "NW"}
)
_NON_AUTOSOMAL_PREFIXES = ("HLA", "HS37D5", "HS38D1", "GL0", "KI2")
# Human RefSeq chromosome accessions: NC_000001..NC_000022 are autosomes,
# NC_000023/24 are X/Y and NC_012920 is the mitochondrion. Other naming
# schemes are the user's responsibility (restrict with --region).
_REFSEQ_HUMAN = re.compile(r"^NC_0000(\d\d)(\.\d+)?$", re.IGNORECASE)
_REFSEQ_NON_AUTOSOMAL = ("NC_012920", "NC_007605")  # MT, EBV

_BASES = frozenset("ACGT")
_SEQUENCE = frozenset("ACGTN")
_PASS = (None, "PASS", ".")

DROP_REASONS = (
    "non_autosomal",
    "not_biallelic",
    "variant_type",
    "filter",
    "inexact_group_counts",
    "call_rate",
    "maf",
)


class PopgenError(ValueError):
    """A request popgen cannot answer (bad scope, no data, ...)."""


def bare_chrom(chrom: str) -> str:
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


def is_autosome(chrom: str) -> bool:
    """True for an autosome name (human-style conventions; see
    NON_AUTOSOMAL). Anything unrecognised is treated as an autosome."""
    m = _REFSEQ_HUMAN.match(chrom)
    if m:
        return 1 <= int(m.group(1)) <= 22
    if chrom.upper().startswith(_REFSEQ_NON_AUTOSOMAL):
        return False
    base = bare_chrom(chrom.split("_", 1)[0]).upper()
    if base in NON_AUTOSOMAL:
        return False
    return not bare_chrom(chrom).upper().startswith(_NON_AUTOSOMAL_PREFIXES)


def chrom_order(chrom: str) -> tuple:
    bare = bare_chrom(chrom)
    if bare.isdigit():
        return (0, int(bare), "")
    return (1, 0, bare)


@dataclass(frozen=True)
class SiteFilters:
    include_indels: bool = False
    pass_only: bool = True
    min_call_rate: float = 0.9
    maf: float = 0.0


@dataclass
class GroupData:
    """Per-site arrays of one group over the retained sites, stored in a
    compact unsigned type (uint16 when every count fits)."""

    size: int
    called: np.ndarray  # called individuals
    k: np.ndarray  # ALT alleles
    het: np.ndarray  # heterozygous individuals

    @property
    def n(self) -> np.ndarray:
        """Called haplotypes (diploid): 2 × called individuals."""
        return 2 * self.called.astype(np.int64)

    def at(self, sel) -> tuple[np.ndarray, np.ndarray]:
        """(k, n) as int64 for the sites selected by `sel` (slice or mask)."""
        return self.k[sel].astype(np.int64), 2 * self.called[sel].astype(np.int64)


@dataclass
class Prepared:
    ingest_id: str
    grouping: str  # population | super_population | all
    ancestral: str  # aa | aa-high | ref | none
    missing_data_tracked: bool
    regions: list[Region]
    stored_chroms: list[str]  # every in-scope chromosome, as stored
    chrom_slices: list[tuple[str, int, int]]  # (stored name, start, stop)
    pos: np.ndarray
    orientation: np.ndarray  # -1 unpolarised, else ANCESTRAL_IS_REF/ALT
    groups: dict[str, GroupData]
    report: dict
    warnings: list[str] = field(default_factory=list)

    @property
    def n_sites(self) -> int:
        return len(self.pos)

    def derived(self, group: str, sel=slice(None)):
        """(polarised mask, derived count, n) over `sel` for `group`."""
        k, n = self.groups[group].at(sel)
        o = self.orientation[sel]
        return o >= 0, np.where(o == ANCESTRAL_IS_ALT, n - k, k), n


def _variant_type_ok(ref: str, alt: str, include_indels: bool) -> bool:
    ref, alt = ref.upper(), alt.upper()
    if len(ref) == 1 and len(alt) == 1:
        return ref in _BASES and alt in _BASES and ref != alt
    if not include_indels:
        return False
    return set(ref) <= _SEQUENCE and set(alt) <= _SEQUENCE and ref != alt


def _sequential_drop(keep: np.ndarray, test: np.ndarray, dropped: dict, why: str):
    dropped[why] += int((keep & ~test).sum())
    return keep & test


@dataclass(frozen=True)
class _Groups:
    """Group layout shared by every chromosome of one run."""

    names: list[str]  # sorted group labels (["all"] when ungrouped)
    sizes: np.ndarray  # samples per group
    n_samples: int
    grouped: bool  # False → the single group "all" = every sample


def _resolve_groups(sess, ingest_id: str, by: str, features, warnings) -> tuple:
    from popgen import counts as cq

    sample_ids = cq.samples(sess, ingest_id)
    if not sample_ids:
        raise PopgenError(f"ingestion {ingest_id!r} has no samples")
    group_labels: dict[str, str] = {}
    if by != "all" and features.populations_table:
        group_labels = cq.labels(sess, ingest_id, by)
    if by != "all" and not group_labels:
        if by == "super_population":
            raise PopgenError(
                f"no super_population labels for ingestion {ingest_id!r}; load "
                "a panel with a super-population column (`vcfclick db panel`)"
            )
        warnings.append(
            f"no population panel is loaded for ingestion {ingest_id!r}: "
            "computing for the whole cohort as one group 'all' (load one "
            "with `vcfclick db panel`, or pass --by all to silence this)"
        )
    stray = set(group_labels) - set(sample_ids)
    if stray:  # labels() joins to samples, so this is a programming error
        raise PopgenError(f"panel labels for samples not in {ingest_id!r}: {stray}")
    if not group_labels:
        return _Groups(
            [ALL_GROUP], np.array([len(sample_ids)]), len(sample_ids), False
        ), 0
    names = sorted(set(group_labels.values()))
    sizes = np.array([sum(1 for v in group_labels.values() if v == g) for g in names])
    unlabelled = len(sample_ids) - int(sizes.sum())
    if unlabelled < 0 or int(sizes.sum()) + unlabelled != len(sample_ids):
        raise PopgenError(
            f"inconsistent panel for {ingest_id!r}: {int(sizes.sum())} labelled "
            f"of {len(sample_ids)} samples"
        )
    return _Groups(names, sizes, len(sample_ids), True), unlabelled


def _filter_chunk(c: ChromCounts, groups: _Groups, filters: SiteFilters, features):
    """keep mask, per-reason drops and tracked mask for one chromosome."""
    n = len(c)
    dropped = dict.fromkeys(DROP_REASONS, 0)
    keep = np.ones(n, dtype=bool)

    # Sites are ordered by (pos, ref, alt): records sharing a position are
    # adjacent, so "more than one record here" is a neighbour comparison.
    same_prev = np.zeros(n, dtype=bool)
    same_prev[1:] = c.pos[1:] == c.pos[:-1]
    shared = same_prev.copy()
    shared[:-1] |= same_prev[1:]
    keep = _sequential_drop(keep, ~shared, dropped, "not_biallelic")

    vtype = np.array(
        [
            _variant_type_ok(r, a, filters.include_indels)
            for r, a in zip(c.ref, c.alt, strict=True)
        ],
        dtype=bool,
    )
    keep = _sequential_drop(keep, vtype, dropped, "variant_type")

    passing = np.array(
        [(not filters.pass_only) or f in _PASS for f in c.filter], dtype=bool
    )
    keep = _sequential_drop(keep, passing, dropped, "filter")

    missing_total = c.missing_total
    alt_total = c.alt_total
    tracked = (c.an_called >= 0) & features.called_columns
    derived_an = 2 * (groups.n_samples - missing_total)
    exact = ~tracked | ((derived_an == c.an_called) & (alt_total == c.ac_called))
    keep = _sequential_drop(keep, exact, dropped, "inexact_group_counts")

    call_ok = np.ones(n, dtype=bool)
    if groups.grouped:
        for gi, size in enumerate(groups.sizes):
            rate = (size - c.missing_by_bucket[:, gi]) / size
            call_ok &= rate >= filters.min_call_rate - 1e-12
    else:
        rate = (groups.n_samples - missing_total) / groups.n_samples
        call_ok &= rate >= filters.min_call_rate - 1e-12
    keep = _sequential_drop(keep, call_ok, dropped, "call_rate")

    if filters.maf > 0:
        p = np.divide(alt_total, derived_an, out=np.zeros(n), where=derived_an > 0)
        maf = np.minimum(p, 1 - p)
        keep = _sequential_drop(keep, maf >= filters.maf - 1e-12, dropped, "maf")

    evidence = (
        tracked
        & (c.n_called < groups.n_samples)
        & (c.an_called == 2 * c.n_called)
        & (missing_total == 0)
    )
    return keep, dropped, tracked, bool(evidence.any())


def _aa_codes(c: ChromCounts, idx: np.ndarray) -> np.ndarray:
    """`aa_code` of the sites `idx`, memoised on (REF, ALT, AA): the same
    few combinations repeat across a chromosome."""
    memo: dict[tuple, int] = {}
    out = np.empty(len(idx), dtype=np.int8)
    for j, i in enumerate(idx.tolist()):
        key = (c.ref[i], c.alt[i], c.aa[i])
        code = memo.get(key)
        if code is None:
            code = memo[key] = aa_code(*key)
        out[j] = code
    return out


def prepare(
    sess,
    ingest_id: str,
    regions: list[Region],
    by: str,
    filters: SiteFilters,
    ancestral: str | None,
    allow_untracked: bool = False,
) -> Prepared:
    """Fetch, filter and polarise the sites of one ingestion."""
    from popgen import counts as cq

    features = schema_features()
    warnings: list[str] = []
    groups, unlabelled = _resolve_groups(sess, ingest_id, by, features, warnings)
    grouping = by if groups.grouped else "all"
    buckets = (groups.names if groups.grouped else []) + [UNLABELLED]
    dtype = np.uint16 if 2 * int(groups.sizes.max()) < 2**16 else np.uint32

    dropped = dict.fromkeys(DROP_REASONS, 0)
    chrom_list = sorted(
        cq.chromosomes(sess, ingest_id, regions), key=lambda c: chrom_order(c[0])
    )
    in_scope = sum(n for _, n in chrom_list)
    pieces: dict[str, list] = {"pos": [], "aa": [], "called": [], "k": [], "het": []}
    slices: list[tuple[str, int, int]] = []
    tracked_all, missing_unrecorded, any_aa = True, False, False
    offset = 0
    for chrom, count in chrom_list:
        if not is_autosome(chrom):
            dropped["non_autosomal"] += count
            continue
        c = cq.fetch_chrom(
            sess,
            ingest_id,
            chrom,
            regions,
            grouping if groups.grouped else None,
            buckets,
            features,
        )
        keep, d, tracked, unrecorded = _filter_chunk(c, groups, filters, features)
        for reason, value in d.items():
            dropped[reason] += value
        tracked_all &= bool(tracked.all())
        missing_unrecorded |= unrecorded
        any_aa |= any(a is not None for a in c.aa)
        idx = np.flatnonzero(keep)
        if not len(idx):
            continue
        pieces["aa"].append(_aa_codes(c, idx))
        pieces["pos"].append(c.pos[idx])
        if groups.grouped:
            called = groups.sizes[None, :] - c.missing_by_bucket[idx, :-1]
            k, het = c.alt_by_bucket[idx, :-1], c.het_by_bucket[idx, :-1]
        else:
            called = (groups.n_samples - c.missing_total[idx])[:, None]
            k, het = c.alt_total[idx][:, None], c.het_total[idx][:, None]
        pieces["called"].append(called.astype(dtype))
        pieces["k"].append(k.astype(dtype))
        pieces["het"].append(het.astype(dtype))
        slices.append((chrom, offset, offset + len(idx)))
        offset += len(idx)

    missing_data_tracked = (
        features.called_columns
        and features.missing_table
        and tracked_all
        and not missing_unrecorded
    )
    if not missing_data_tracked and not allow_untracked:
        why = (
            "was loaded with --no-record-missing (per-sample missing calls "
            "unknown; sites with missing calls would be dropped)"
            if features.called_columns and features.missing_table and tracked_all
            else "predates called-genotype tracking (missing calls cannot be "
            "told from 0/0 and would be counted as 0/0)"
        )
        raise PopgenError(
            f"ingestion {ingest_id!r} {why}. Re-ingest the VCF for exact "
            "results, or pass --allow-untracked to compute anyway."
        )
    if not tracked_all:
        warnings.append(
            "this database (or ingestion) predates called-genotype tracking: "
            "missing calls cannot be told from 0/0 and are counted as 0/0 "
            "(missing_data_tracked: false). Re-ingest the VCF to fix."
        )
    elif missing_unrecorded:
        warnings.append(
            "this ingestion was loaded with --no-record-missing: per-group "
            "called counts are unknown, so sites with missing calls are "
            "dropped (missing_data_tracked: false)"
        )
    if dropped["inexact_group_counts"]:
        warnings.append(
            f"{dropped['inexact_group_counts']} site(s) dropped: per-group "
            "called counts cannot be derived exactly there (partially missing "
            "./1, haploid or polyploid calls, or missing calls not recorded)"
        )

    def cat(name: str, width: int, kind) -> np.ndarray:
        if pieces[name]:
            return np.concatenate(pieces[name])
        return np.zeros((0, width) if width else 0, dtype=kind)

    ncols = len(groups.names)
    columns: dict[str, list[np.ndarray]] = {}
    for x in ("called", "k", "het"):
        # One count type at a time, freeing its pieces first, so the
        # retained data is never held twice over.
        mat = cat(x, ncols, dtype)
        pieces[x].clear()
        columns[x] = [np.ascontiguousarray(mat[:, gi]) for gi in range(ncols)]
        del mat
    pos = cat("pos", 0, np.int64)
    if ancestral is None:
        ancestral = "aa" if any_aa else "none"
    orientation = orientation_from_codes(cat("aa", 0, np.int8), ancestral)

    group_data = {
        name: GroupData(
            size=int(groups.sizes[gi]),
            called=columns["called"][gi],
            k=columns["k"][gi],
            het=columns["het"][gi],
        )
        for gi, name in enumerate(groups.names)
    }

    n_retained = len(pos)
    n_polarised = int((orientation >= 0).sum())
    report = {
        "in_scope": in_scope,
        "retained": n_retained,
        "dropped": dropped,
        "polarised": n_polarised if ancestral != "none" else 0,
        "unpolarised": n_retained - n_polarised if ancestral != "none" else n_retained,
        "samples": groups.n_samples,
        "unlabelled_samples": unlabelled,
    }
    if ancestral in ("aa", "aa-high") and n_retained and n_polarised == 0:
        warnings.append(
            f"--ancestral {ancestral}: no retained site could be polarised "
            "(INFO/AA missing or matching neither REF nor ALT)"
        )
    return Prepared(
        ingest_id=ingest_id,
        grouping=grouping,
        ancestral=ancestral,
        missing_data_tracked=bool(missing_data_tracked),
        regions=regions,
        stored_chroms=[c for c, _ in chrom_list],
        chrom_slices=slices,
        pos=pos,
        orientation=orientation,
        groups=group_data,
        report=report,
        warnings=warnings,
    )


# ───────────────────────────── statistics ───────────────────────────────


def projection_size(prep: Prepared, group: str, project: int | None) -> int | None:
    """The SFS projection size: `project` if given, else the smallest
    number of called haplotypes at any retained site in the group."""
    if project is not None:
        return project
    called = prep.groups[group].called
    return 2 * int(called.min()) if len(called) else None


def projection_warnings(prep: Prepared, project: int | None) -> list[str]:
    """Explain a projection size too small for the statistics that use it
    (typically a site with almost no calls kept by --min-call-rate 0), or an
    explicit --project larger than some sites' called haplotypes, which
    leaves Tajima's D, Fay & Wu's H and the SFS on fewer sites than θ_W and
    π in the same row."""
    out = []
    for name in prep.groups:
        m = projection_size(prep, name, project)
        if project is not None and m is not None and prep.n_sites:
            _, n = prep.groups[name].at(slice(None))
            used = int(np.count_nonzero(n >= m))
            if used < len(n):
                most = int(n.max()) if len(n) else 0
                detail = (
                    f"no site has that many (at most {most}), so they are NA"
                    if used == 0
                    else f"they use only the {used} of {len(n)} sites with at least {m}"
                )
                out.append(
                    f"group {name}: --project {m} exceeds the called haplotypes "
                    f"at some sites; Tajima's D, Fay & Wu's H and the SFS need "
                    f"{m} per site, and {detail} (θ_W and π use every site)"
                )
        if m is None or m >= 4:
            continue
        what = "the SFS" if m < 2 else "Tajima's D and Fay & Wu's H"
        source = (
            "--project"
            if project is not None
            else "the smallest number of called haplotypes at a retained site"
        )
        out.append(
            f"group {name}: projection size is {m} ({source}), so {what} "
            f"cannot be computed (needs at least {2 if m < 2 else 4}); raise "
            "--min-call-rate or pass a larger --project"
        )
    return out


def _projected(k: np.ndarray, n: np.ndarray, m: int | None):
    """θ estimators of the projected SFS (closed form; see
    estimators.projected_thetas), or (None, 0) when m is unusable."""
    if m is None or m < 2:
        return None, 0
    return est.projected_thetas(k, n, m)


def _n_selected(prep: Prepared, sel) -> int:
    if isinstance(sel, slice):
        return len(range(*sel.indices(prep.n_sites)))
    return int(np.count_nonzero(sel))


def _diversity(prep: Prepared, group: str, sel, m: int | None) -> dict:
    """S, θ_W, π and Tajima's D for `group` over the sites in `sel`."""
    k, n = prep.groups[group].at(sel)
    theta_w = float(est.watterson_per_site(k, n).sum())
    pi = float(est.pi_per_site(k, n).sum())
    t, used = _projected(k, n, m)
    d = est.tajima_d(t.pi, t.s, m) if t is not None and used else None
    return {
        "sites": len(k),
        "projected_sites": int(used),
        "segregating_sites": int(est.segregating(k, n).sum()),
        "theta_w": theta_w,
        "pi": pi,
        "tajima_d": d,
    }


def _fay_wu(prep: Prepared, group: str, sel, m: int | None) -> dict:
    if prep.ancestral == "none":
        return {"fay_wu_h": None, "fay_wu_h_raw": None, "polarised_sites": 0}
    pol, d, n = prep.derived(group, sel)
    t, used = _projected(d[pol], n[pol], m)
    h = h_raw = None
    if t is not None and used:
        h = est.fay_wu_h_normalised(t.pi, t.theta_l, t.s, m)
        h_raw = t.pi - t.theta_h if t.s > 0 else None
    return {"fay_wu_h": h, "fay_wu_h_raw": h_raw, "polarised_sites": int(pol.sum())}


def summary(prep: Prepared, project: int | None = None) -> list[dict]:
    """Per-group diversity summary."""
    out = []
    every = slice(None)
    for name, g in prep.groups.items():
        m = projection_size(prep, name, project)
        div = _diversity(prep, name, every, m)
        k, n = g.at(every)
        het = est.heterozygosity(g.het, g.called, k, n)
        n_sites = div["sites"]
        out.append(
            {
                "group": name,
                "n_samples": g.size,
                "sites": n_sites,
                "segregating_sites": div["segregating_sites"],
                "theta_w": div["theta_w"],
                "theta_w_per_site": div["theta_w"] / n_sites if n_sites else None,
                "pi": div["pi"],
                "pi_per_site": div["pi"] / n_sites if n_sites else None,
                "tajima_d": div["tajima_d"],
                **_fay_wu(prep, name, every, m),
                "projection_n": m,
                "projected_sites": div["projected_sites"],
                "ho": het.ho,
                "he": het.he,
                "f": het.f,
            }
        )
    return out


def sfs(prep: Prepared, project: int | None = None) -> list[dict]:
    """Per-group projected SFS: folded always, unfolded when polarised."""
    out = []
    for name, g in prep.groups.items():
        m = projection_size(prep, name, project)
        entry: dict = {"group": name, "projection_n": m}
        if m is None or m < 2:
            entry.update(
                sites_used=0,
                sites_dropped=prep.n_sites,
                folded=None,
                unfolded=None,
                polarised_sites_used=0,
            )
            out.append(entry)
            continue
        k, n = g.at(slice(None))
        xi, used = est.project_sfs(k, n, m)
        entry["sites_used"] = used
        entry["sites_dropped"] = prep.n_sites - used
        entry["folded"] = est.fold(xi).tolist()
        if prep.ancestral == "none":
            entry["unfolded"] = None
            entry["polarised_sites_used"] = 0
        else:
            pol, d, n = prep.derived(name)
            uxi, uused = est.project_sfs(d[pol], n[pol], m)
            # No polarised site (every INFO/AA unknown): no spectrum, not zeros.
            entry["unfolded"] = uxi.tolist() if uused else None
            entry["polarised_sites_used"] = uused
        out.append(entry)
    return out


def group_pairs(prep: Prepared) -> list[tuple[str, str]]:
    names = list(prep.groups)
    return [(a, b) for i, a in enumerate(names) for b in names[i + 1 :]]


def _fst_pair(prep: Prepared, a: str, b: str, sel) -> dict:
    k1, n1 = prep.groups[a].at(sel)
    k2, n2 = prep.groups[b].at(sel)
    ok = (n1 >= 2) & (n2 >= 2)
    return {
        "fst": est.hudson_fst(k1[ok], n1[ok], k2[ok], n2[ok]),
        "sites": int(ok.sum()),
    }


def fst(prep: Prepared, window: int | None = None, step: int | None = None) -> dict:
    """Pairwise Hudson F_ST (ratio of averages), overall and per window."""
    pairs = group_pairs(prep)
    every = slice(None)
    result: dict = {
        "pairs": [
            {"group1": a, "group2": b, **_fst_pair(prep, a, b, every)} for a, b in pairs
        ]
    }
    if window:
        rows = []
        for chrom, start, end, sel in iter_windows(prep, window, step or window):
            row = {
                "chrom": chrom,
                "start": start,
                "end": end,
                "sites": _n_selected(prep, sel),
            }
            for a, b in pairs:
                row[f"fst_{a}_{b}"] = _fst_pair(prep, a, b, sel)["fst"]
            rows.append(row)
        result["windows"] = rows
    return result


def _stored_name(prep: Prepared, chrom: str) -> str:
    """The stored spelling of `chrom` (chr1 vs 1), when it is in scope."""
    if chrom in prep.stored_chroms:
        return chrom
    for stored in prep.stored_chroms:
        if bare_chrom(stored) == bare_chrom(chrom):
            return stored
    return chrom


def _window_spans(prep: Prepared) -> list[tuple[str, int, int, int, int]]:
    """(stored chrom, start bp, end bp, first site, stop site) spans to
    tile: the requested regions, or each chromosome with retained sites
    from 1 to its last retained site."""
    by_name = {c: (a, b) for c, a, b in prep.chrom_slices}
    spans = []
    if prep.regions:
        for r in prep.regions:
            name = _stored_name(prep, r.chrom)
            a, b = by_name.get(name, (0, 0))
            if r.start is None:
                if b > a:
                    spans.append((name, 1, int(prep.pos[b - 1]), a, b))
            else:
                spans.append((name, r.start, r.end if r.end else r.start, a, b))
        return spans
    for name, a, b in prep.chrom_slices:
        spans.append((name, 1, int(prep.pos[b - 1]), a, b))
    return spans


def iter_windows(prep: Prepared, window: int, step: int):
    """Yield (stored chrom, start, end, site slice) for sliding windows.

    Each chromosome's retained sites are a contiguous, position-sorted
    slice, so a window is two binary searches and a slice — no full-length
    masks. Tiling stops at the first window that reaches the span's end.
    """
    if window < 1 or step < 1:
        raise PopgenError("--window and --step must be positive")
    for chrom, lo, hi, a, b in _window_spans(prep):
        starts = np.arange(lo, hi + 1, step, dtype=np.int64)
        ends = np.minimum(starts + window - 1, hi)
        # Stop at the first window that reaches the span's end; if none does
        # (step > window leaves gaps), keep every window that starts in span.
        reach = ends >= hi
        last = int(np.argmax(reach)) if reach.any() else len(ends) - 1
        starts, ends = starts[: last + 1], ends[: last + 1]
        pos = prep.pos[a:b]
        firsts = a + np.searchsorted(pos, starts, side="left")
        stops = a + np.searchsorted(pos, ends, side="right")
        for s, e, f, t in zip(
            starts.tolist(), ends.tolist(), firsts.tolist(), stops.tolist(), strict=True
        ):
            yield chrom, s, e, slice(f, t)


def windows(
    prep: Prepared,
    window: int,
    step: int | None = None,
    with_fst: bool = False,
    project: int | None = None,
) -> list[dict]:
    """Sliding-window S, θ_W, π, Tajima's D per group (+ pairwise F_ST).

    θ_W and π are sums over the window's retained sites; the `_per_bp`
    columns divide by the window length, which assumes every position
    without a VCF record is callable and monomorphic (see POPGEN.md)."""
    sizes = {name: projection_size(prep, name, project) for name in prep.groups}
    pairs = group_pairs(prep) if with_fst else []
    rows = []
    for chrom, start, end, sel in iter_windows(prep, window, step or window):
        length = end - start + 1
        row: dict = {
            "chrom": chrom,
            "start": start,
            "end": end,
            "sites": _n_selected(prep, sel),
        }
        for name in prep.groups:
            div = _diversity(prep, name, sel, sizes[name])
            row[f"{name}_S"] = div["segregating_sites"]
            row[f"{name}_theta_w"] = div["theta_w"]
            row[f"{name}_pi"] = div["pi"]
            row[f"{name}_theta_w_per_bp"] = div["theta_w"] / length
            row[f"{name}_pi_per_bp"] = div["pi"] / length
            row[f"{name}_tajima_d"] = div["tajima_d"]
        for a, b in pairs:
            row[f"fst_{a}_{b}"] = _fst_pair(prep, a, b, sel)["fst"]
        rows.append(row)
    return rows


__all__ = [
    "ALL_GROUP",
    "ANCESTRAL_IS_REF",
    "PopgenError",
    "Prepared",
    "SiteFilters",
    "fetch_chrom",
    "fst",
    "prepare",
    "projection_warnings",
    "sfs",
    "summary",
    "windows",
]
