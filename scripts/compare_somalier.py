#!/usr/bin/env python3
"""Cross-check `vcfclick db relatedness` / `db qc` sex against somalier.

Optional validation helper (not run in CI; somalier ships a Linux binary).
See "Relatedness and sex cross-check with somalier" in docs/VALIDATION.md.

Inputs:
  --vcfclick-relatedness  JSON from `vcfclick db relatedness NAME --all --format json`
  --somalier-pairs        somalier.pairs.tsv from `somalier relate`
  --vcfclick-qc           JSON from `vcfclick db qc NAME --format json` (optional)
  --somalier-samples      somalier.samples.tsv from `somalier relate` (optional)

somalier reports *relatedness* (≈ 2 x kinship: 0.5 for first degree) and
IBS0 as a count over `n` sites; both are rescaled to vcfclick's KING
kinship and IBS0 fraction, then classified with the same thresholds
vcfclick uses (`cli.db_relatedness.classify`). Sex is inferred from
somalier's X_het / (X_het + X_hom_alt) with `cli.db_qc._infer_sex`, so
the two sides differ only in the markers, not the decision rules.

Exit status is 1 when any pair's relationship or any sample's sex
disagrees, or somalier finds relatives in a pair vcfclick did not report,
so the script can gate a validation run.

The somalier column names used here (`sample_a`, `sample_b`,
`relatedness`, `ibs0`, `n`; `sample_id`, `X_het`, `X_hom_alt`) and the
2 x kinship scaling were checked against somalier's relate.nim.

Usage:
    uv run python scripts/compare_somalier.py \\
        --vcfclick-relatedness rel.json --somalier-pairs somalier.pairs.tsv \\
        --vcfclick-qc qc.json --somalier-samples somalier.samples.tsv
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from cli.db_qc import _infer_sex  # noqa: E402
from cli.db_relatedness import classify  # noqa: E402


def _read_tsv(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        lines = list(fh)
    if lines and lines[0].startswith("#"):
        lines[0] = lines[0][1:]
    return list(csv.DictReader(lines, delimiter="\t"))


def _key(a: str, b: str) -> tuple[str, str]:
    return (a, b) if a <= b else (b, a)


def somalier_pairs(path: Path) -> dict[tuple[str, str], dict]:
    out = {}
    for r in _read_tsv(path):
        n = int(r["n"])
        rel = float(r["relatedness"])
        kinship = rel / 2
        ibs0 = int(r["ibs0"]) / n if n else float("nan")
        out[_key(r["sample_a"], r["sample_b"])] = {
            "kinship": round(kinship, 4),
            "relationship": classify(kinship if n else float("nan"), ibs0),
        }
    return out


def vcfclick_pairs(path: Path) -> dict[tuple[str, str], dict]:
    out = {}
    for res in json.loads(Path(path).read_text()):
        for p in res["pairs"]:
            out[_key(p["sample_a"], p["sample_b"])] = p
    return out


def somalier_sex(path: Path) -> dict[str, str]:
    out = {}
    for r in _read_tsv(path):
        het, hom = int(r["X_het"]), int(r["X_hom_alt"])
        out[r["sample_id"]] = _infer_sex(het, het + hom)[0]
    return out


def vcfclick_sex(path: Path) -> dict[str, str]:
    return {s["sample_id"]: s["inferred_sex"] for s in json.loads(path.read_text())}


def compare(args: argparse.Namespace) -> list[str]:
    """Print a comparison table; return a list of disagreement messages."""
    problems: list[str] = []
    vc = vcfclick_pairs(args.vcfclick_relatedness)
    so = somalier_pairs(args.somalier_pairs)
    shared = sorted(set(vc) & set(so))
    print(f"pairs: {len(shared)} shared, {len(vc)} vcfclick, {len(so)} somalier")
    print(
        f"{'sample_a':<16}{'sample_b':<16}{'vc_kin':>8}{'so_kin':>8}  vcfclick / somalier"
    )
    for k in shared:
        a, b = vc[k], so[k]
        vk = "n/a" if a["kinship"] is None else f"{a['kinship']:.3f}"
        flag = ""
        if "insufficient-data" in (a["relationship"], b["relationship"]):
            flag = "  (not compared)"
        elif a["relationship"] != b["relationship"]:
            flag = "  <-- differs"
            problems.append(
                f"{k[0]}/{k[1]}: vcfclick {a['relationship']}, somalier {b['relationship']}"
            )
        print(
            f"{k[0]:<16}{k[1]:<16}{vk:>8}{b['kinship']:>8.3f}  "
            f"{a['relationship']} / {b['relationship']}{flag}"
        )

    # `db relatedness` without --all lists related pairs only, so a pair it
    # omits is one it calls unrelated. somalier finding relatives there is a
    # disagreement, not a pair to skip.
    for k in sorted(set(so) - set(vc)):
        if so[k]["relationship"] not in ("unrelated", "insufficient-data"):
            problems.append(
                f"{k[0]}/{k[1]}: somalier {so[k]['relationship']}, not reported by vcfclick"
            )
            print(
                f"{k[0]:<16}{k[1]:<16}{'-':>8}{so[k]['kinship']:>8.3f}  "
                f"(none) / {so[k]['relationship']}  <-- differs"
            )

    if args.vcfclick_qc and args.somalier_samples:
        vs = vcfclick_sex(args.vcfclick_qc)
        ss = somalier_sex(args.somalier_samples)
        print(f"\nsex: {len(set(vs) & set(ss))} shared samples")
        for s in sorted(set(vs) & set(ss)):
            flag = ""
            if "unknown" in (vs[s], ss[s]):
                flag = "  (not compared)"
            elif vs[s] != ss[s]:
                flag = "  <-- differs"
                problems.append(f"{s}: sex vcfclick {vs[s]}, somalier {ss[s]}")
            print(f"{s:<16}{vs[s]:>10}{ss[s]:>10}{flag}")
    return problems


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--vcfclick-relatedness", type=Path, required=True)
    ap.add_argument("--somalier-pairs", type=Path, required=True)
    ap.add_argument("--vcfclick-qc", type=Path)
    ap.add_argument("--somalier-samples", type=Path)
    args = ap.parse_args(argv)
    problems = compare(args)
    if problems:
        print(f"\n{len(problems)} disagreement(s):", file=sys.stderr)
        for p in problems:
            print(f"  {p}", file=sys.stderr)
        return 1
    print("\nno disagreements among the compared pairs and samples")
    return 0


if __name__ == "__main__":
    sys.exit(main())
