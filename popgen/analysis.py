"""Site filtering, polarisation and the `db popgen` computations.

`prepare()` turns the raw per-site counts of one ingestion into the set
of retained sites with per-group (n, k, het, called) arrays and a report
of what was dropped and why. `summary`, `sfs`, `fst` and `windows` work
only on that prepared data, so every subcommand sees the same sites.

Choices (documented in docs/POPGEN.md):

  * one ingestion at a time (sample identity is (ingest_id, sample_id),
    and two ingestions do not share a site list);
  * autosomes only (sex chromosomes and MT are excluded; including them
    needs per-sample ploidy and is refused in this version);
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

from dataclasses import dataclass, field

import numpy as np

from popgen import estimators as est
from popgen.ancestral import ANCESTRAL_IS_ALT, ANCESTRAL_IS_REF, polarise
from popgen.counts import (
    ALL_GROUP,
    Region,
    SiteCounts,
    fetch_counts,
    schema_features,
)

# Chromosomes excluded by "autosomes only": sex chromosomes of XY and ZW
# systems, the PAR pseudo-contig and mitochondria (with or without a
# `chr` prefix, any case).
NON_AUTOSOMAL = frozenset({"X", "Y", "XY", "W", "Z", "M", "MT"})

_BASES = frozenset("ACGT")
_SEQUENCE = frozenset("ACGTN")
_PASS = (None, "PASS", ".")


class PopgenError(ValueError):
    """A request popgen cannot answer (bad scope, no data, ...)."""


def bare_chrom(chrom: str) -> str:
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


def is_autosome(chrom: str) -> bool:
    return bare_chrom(chrom).upper() not in NON_AUTOSOMAL


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
    """Per-site arrays of one group over the retained sites."""

    size: int
    n: np.ndarray  # called haplotypes
    k: np.ndarray  # ALT alleles
    het: np.ndarray  # heterozygous individuals
    called: np.ndarray  # called individuals


@dataclass
class Prepared:
    ingest_id: str
    grouping: str  # population | super_population | all
    ancestral: str  # aa | aa-high | ref | none
    missing_data_tracked: bool
    regions: list[Region]
    chrom: list[str]
    pos: np.ndarray
    orientation: np.ndarray  # -1 unpolarised, else ANCESTRAL_IS_REF/ALT
    groups: dict[str, GroupData]
    report: dict
    warnings: list[str] = field(default_factory=list)

    @property
    def n_sites(self) -> int:
        return len(self.chrom)

    def derived(self, group: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """(mask of polarised sites, derived count, n) for `group`."""
        g = self.groups[group]
        mask = self.orientation >= 0
        d = np.where(self.orientation == ANCESTRAL_IS_ALT, g.n - g.k, g.k)
        return mask, d, g.n


def _variant_type_ok(ref: str, alt: str, include_indels: bool) -> bool:
    ref, alt = ref.upper(), alt.upper()
    if len(ref) == 1 and len(alt) == 1:
        return ref in _BASES and alt in _BASES and ref != alt
    if not include_indels:
        return False
    return set(ref) <= _SEQUENCE and set(alt) <= _SEQUENCE and ref != alt


def _sequential_drop(keep: np.ndarray, test: np.ndarray, dropped: dict, why: str):
    newly = keep & ~test
    dropped[why] = int(newly.sum())
    return keep & test


def prepare(
    sess,
    ingest_id: str,
    regions: list[Region],
    by: str,
    filters: SiteFilters,
    ancestral: str | None,
) -> Prepared:
    """Fetch, filter and polarise the sites of one ingestion."""
    from popgen import counts as cq

    features = schema_features()
    warnings: list[str] = []

    group_labels: dict[str, str] = {}
    grouping = "all"
    if by != "all" and features.populations_table:
        group_labels = cq.labels(sess, ingest_id, by)
        if group_labels:
            grouping = by
    if by == "super_population" and grouping == "all":
        raise PopgenError(
            f"no super_population labels for ingestion {ingest_id!r}; load a "
            "panel with a super-population column (`vcfclick db panel`)"
        )

    counts = fetch_counts(
        sess,
        ingest_id,
        regions,
        None if grouping == "all" else grouping,
        features,
        group_labels,
    )
    if counts.n_samples == 0:
        raise PopgenError(f"ingestion {ingest_id!r} has no samples")
    unlabelled = counts.n_samples - len(group_labels) if group_labels else 0

    keep, dropped, tracked = _filter_sites(counts, filters, features, warnings)
    tracked_all = bool(tracked.all()) if len(counts) else features.called_columns
    missing_unrecorded = _missing_unrecorded(counts, tracked)
    missing_data_tracked = (
        features.called_columns
        and features.missing_table
        and tracked_all
        and not missing_unrecorded
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

    if ancestral is None:
        ancestral = "aa" if any(a is not None for a in counts.aa) else "none"
    idx = np.flatnonzero(keep)
    orientation = np.array(
        [
            -1
            if (o := polarise(counts.ref[i], counts.alt[i], counts.aa[i], ancestral))
            is None
            else o
            for i in idx
        ],
        dtype=np.int64,
    )

    groups = {}
    for name in sorted(counts.group_size):
        size = counts.group_size[name]
        missing = counts.group_missing[name][idx]
        called = size - missing
        groups[name] = GroupData(
            size=size,
            n=2 * called,
            k=counts.group_alt[name][idx],
            het=counts.group_het[name][idx],
            called=called,
        )

    n_polarised = int((orientation >= 0).sum())
    report = {
        "in_scope": len(counts),
        "retained": int(len(idx)),
        "dropped": dropped,
        "polarised": n_polarised if ancestral != "none" else 0,
        "unpolarised": int(len(idx)) - n_polarised
        if ancestral != "none"
        else int(len(idx)),
        "samples": counts.n_samples,
        "unlabelled_samples": unlabelled,
    }
    if ancestral in ("aa", "aa-high") and len(idx) and n_polarised == 0:
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
        chrom=[counts.chrom[i] for i in idx],
        pos=counts.pos[idx],
        orientation=orientation,
        groups=groups,
        report=report,
        warnings=warnings,
    )


def _filter_sites(
    counts: SiteCounts, filters: SiteFilters, features, warnings: list[str]
):
    n = len(counts)
    keep = np.ones(n, dtype=bool)
    dropped: dict[str, int] = {}

    keep = _sequential_drop(
        keep,
        np.array([is_autosome(c) for c in counts.chrom], dtype=bool),
        dropped,
        "non_autosomal",
    )

    per_position: dict[tuple[str, int], int] = {}
    for c, p in zip(counts.chrom, counts.pos.tolist(), strict=True):
        per_position[(c, p)] = per_position.get((c, p), 0) + 1
    biallelic = np.array(
        [
            per_position[(c, p)] == 1
            for c, p in zip(counts.chrom, counts.pos.tolist(), strict=True)
        ],
        dtype=bool,
    )
    keep = _sequential_drop(keep, biallelic, dropped, "not_biallelic")

    vtype = np.array(
        [
            _variant_type_ok(r, a, filters.include_indels)
            for r, a in zip(counts.ref, counts.alt, strict=True)
        ],
        dtype=bool,
    )
    keep = _sequential_drop(keep, vtype, dropped, "variant_type")

    passing = np.array(
        [(not filters.pass_only) or f in _PASS for f in counts.filter], dtype=bool
    )
    keep = _sequential_drop(keep, passing, dropped, "filter")

    tracked = (counts.an_called >= 0) & features.called_columns
    derived_an = 2 * (counts.n_samples - counts.missing_total)
    exact = ~tracked | (
        (derived_an == counts.an_called) & (counts.alt_total == counts.ac_called)
    )
    keep = _sequential_drop(keep, exact, dropped, "inexact_group_counts")
    if dropped["inexact_group_counts"]:
        warnings.append(
            f"{dropped['inexact_group_counts']} site(s) dropped: per-group "
            "called counts cannot be derived exactly there (partially missing "
            "./1, haploid or polyploid calls, or missing calls not recorded)"
        )

    call_ok = np.ones(n, dtype=bool)
    for name, size in counts.group_size.items():
        rate = (size - counts.group_missing[name]) / size
        call_ok &= rate >= filters.min_call_rate - 1e-12
    keep = _sequential_drop(keep, call_ok, dropped, "call_rate")

    if filters.maf > 0:
        an = 2 * (counts.n_samples - counts.missing_total)
        p = np.divide(counts.alt_total, an, out=np.zeros(n), where=an > 0)
        maf = np.minimum(p, 1 - p)
        keep = _sequential_drop(keep, maf >= filters.maf - 1e-12, dropped, "maf")
    else:
        dropped["maf"] = 0
    return keep, dropped, tracked


def _missing_unrecorded(counts: SiteCounts, tracked: np.ndarray) -> bool:
    """True when the ingestion evidently skipped `missing_genotypes`: a
    site has fully missing diploid calls (n_called < N with exactly two
    called alleles per called sample) but no missing rows."""
    if not len(counts):
        return False
    evidence = (
        tracked
        & (counts.n_called < counts.n_samples)
        & (counts.an_called == 2 * counts.n_called)
        & (counts.missing_total == 0)
    )
    return bool(evidence.any())


# ───────────────────────────── statistics ───────────────────────────────


def projection_size(prep: Prepared, group: str, project: int | None) -> int | None:
    """The SFS projection size: `project` if given, else the smallest
    number of called haplotypes at any retained site in the group."""
    if project is not None:
        return project
    n = prep.groups[group].n
    return int(n.min()) if len(n) else None


def _projected(k: np.ndarray, n: np.ndarray, m: int | None):
    if m is None or m < 2:
        return None, 0
    return est.project_sfs(k, n, m)


def _diversity(prep: Prepared, group: str, mask: np.ndarray, m: int | None) -> dict:
    """S, θ_W, π and Tajima's D for `group` over the sites in `mask`."""
    g = prep.groups[group]
    k, n = g.k[mask], g.n[mask]
    n_sites = int(mask.sum())
    theta_w = float(est.watterson_per_site(k, n).sum())
    pi = float(est.pi_per_site(k, n).sum())
    xi, used = _projected(k, n, m)
    d = None
    if xi is not None and used:
        t = est.thetas_from_sfs(xi)
        d = est.tajima_d(t.pi, t.s, m)
    return {
        "sites": n_sites,
        "segregating_sites": int(est.segregating(k, n).sum()),
        "theta_w": theta_w,
        "pi": pi,
        "tajima_d": d,
    }


def _fay_wu(prep: Prepared, group: str, mask: np.ndarray, m: int | None) -> dict:
    if prep.ancestral == "none":
        return {"fay_wu_h": None, "fay_wu_h_raw": None, "polarised_sites": 0}
    pol, d, n = prep.derived(group)
    sel = pol & mask
    xi, used = _projected(d[sel], n[sel], m)
    h = h_raw = None
    if xi is not None and used:
        t = est.thetas_from_sfs(xi)
        h = est.fay_wu_h_normalised(t.pi, t.theta_l, t.s, m)
        h_raw = t.pi - t.theta_h if t.s > 0 else None
    return {"fay_wu_h": h, "fay_wu_h_raw": h_raw, "polarised_sites": int(sel.sum())}


def summary(prep: Prepared, project: int | None = None) -> list[dict]:
    """Per-group diversity summary."""
    out = []
    every = np.ones(prep.n_sites, dtype=bool)
    for name, g in prep.groups.items():
        m = projection_size(prep, name, project)
        div = _diversity(prep, name, every, m)
        het = est.heterozygosity(g.het, g.called, g.k, g.n)
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
        xi, used = est.project_sfs(g.k, g.n, m)
        entry["sites_used"] = used
        entry["sites_dropped"] = prep.n_sites - used
        entry["folded"] = est.fold(xi).tolist()
        if prep.ancestral == "none":
            entry["unfolded"] = None
            entry["polarised_sites_used"] = 0
        else:
            pol, d, n = prep.derived(name)
            uxi, uused = est.project_sfs(d[pol], n[pol], m)
            entry["unfolded"] = uxi.tolist()
            entry["polarised_sites_used"] = uused
        out.append(entry)
    return out


def group_pairs(prep: Prepared) -> list[tuple[str, str]]:
    names = list(prep.groups)
    return [(a, b) for i, a in enumerate(names) for b in names[i + 1 :]]


def _fst_pair(prep: Prepared, a: str, b: str, mask: np.ndarray) -> dict:
    ga, gb = prep.groups[a], prep.groups[b]
    sel = mask & (ga.n >= 2) & (gb.n >= 2)
    return {
        "fst": est.hudson_fst(ga.k[sel], ga.n[sel], gb.k[sel], gb.n[sel]),
        "sites": int(sel.sum()),
    }


def fst(prep: Prepared, window: int | None = None, step: int | None = None) -> dict:
    """Pairwise Hudson F_ST (ratio of averages), overall and per window."""
    pairs = group_pairs(prep)
    every = np.ones(prep.n_sites, dtype=bool)
    result: dict = {
        "pairs": [
            {"group1": a, "group2": b, **_fst_pair(prep, a, b, every)} for a, b in pairs
        ]
    }
    if window:
        rows = []
        for chrom, start, end, mask in iter_windows(prep, window, step or window):
            row = {"chrom": chrom, "start": start, "end": end, "sites": int(mask.sum())}
            for a, b in pairs:
                row[f"fst_{a}_{b}"] = _fst_pair(prep, a, b, mask)["fst"]
            rows.append(row)
        result["windows"] = rows
    return result


def _window_bounds(prep: Prepared) -> list[tuple[str, int, int]]:
    """(chrom, start, end) spans to tile: the requested regions, or each
    chromosome with retained sites from 1 to its last retained site."""
    spans = []
    if prep.regions:
        last: dict[str, int] = {}
        for c, p in zip(prep.chrom, prep.pos.tolist(), strict=True):
            last[bare_chrom(c)] = max(last.get(bare_chrom(c), 0), p)
        for r in prep.regions:
            if r.start is None:
                end = last.get(bare_chrom(r.chrom))
                if end:
                    spans.append((r.chrom, 1, end))
            else:
                spans.append((r.chrom, r.start, r.end if r.end else r.start))
        return spans
    by_chrom: dict[str, int] = {}
    for c, p in zip(prep.chrom, prep.pos.tolist(), strict=True):
        by_chrom[c] = max(by_chrom.get(c, 0), p)
    for c in sorted(by_chrom, key=chrom_order):
        spans.append((c, 1, by_chrom[c]))
    return spans


def iter_windows(prep: Prepared, window: int, step: int):
    """Yield (chrom, start, end, site mask) for sliding windows."""
    if window < 1 or step < 1:
        raise PopgenError("--window and --step must be positive")
    chroms = np.array([bare_chrom(c) for c in prep.chrom], dtype=object)
    for chrom, lo, hi in _window_bounds(prep):
        on_chrom = chroms == bare_chrom(chrom)
        start = lo
        while start <= hi:
            end = min(start + window - 1, hi)
            mask = on_chrom & (prep.pos >= start) & (prep.pos <= end)
            yield chrom, start, end, mask
            if end >= hi:
                break
            start += step


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
    for chrom, start, end, mask in iter_windows(prep, window, step or window):
        length = end - start + 1
        row: dict = {
            "chrom": chrom,
            "start": start,
            "end": end,
            "sites": int(mask.sum()),
        }
        for name in prep.groups:
            div = _diversity(prep, name, mask, sizes[name])
            row[f"{name}_S"] = div["segregating_sites"]
            row[f"{name}_theta_w"] = div["theta_w"]
            row[f"{name}_pi"] = div["pi"]
            row[f"{name}_theta_w_per_bp"] = div["theta_w"] / length
            row[f"{name}_pi_per_bp"] = div["pi"] / length
            row[f"{name}_tajima_d"] = div["tajima_d"]
        for a, b in pairs:
            row[f"fst_{a}_{b}"] = _fst_pair(prep, a, b, mask)["fst"]
        rows.append(row)
    return rows


__all__ = [
    "ALL_GROUP",
    "ANCESTRAL_IS_REF",
    "PopgenError",
    "Prepared",
    "SiteFilters",
    "fst",
    "prepare",
    "sfs",
    "summary",
    "windows",
]
