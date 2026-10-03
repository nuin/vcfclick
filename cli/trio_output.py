"""Text and JSON presentation of trio inheritance candidates."""

from __future__ import annotations

import click

from cli.trio_analysis import TrioAnalysis

_NEEDS_REF = {"denovo", "dominant", "comphet"}


def _annotation_status() -> dict:
    """Which annotation sources hold data (each lookup is skipped if not)."""
    try:
        from annotations.db import get_connection

        conn = get_connection()
        have = {}
        for key, table in (
            ("genes", "refseq_genes"),
            ("clinvar", "clinvar_variants"),
            ("gnomad", "gnomad_af"),
        ):
            try:
                have[key] = (
                    conn.execute(f"SELECT count(*) FROM {table}").fetchone()[0] > 0
                )
            except Exception:  # noqa: BLE001 - table absent in older stores
                have[key] = False
        return have
    except Exception:  # noqa: BLE001 - no annotation store at all
        return {"genes": False, "clinvar": False, "gnomad": False}


def _enrich(chrom: str, pos, ref: str, alt: str, have: dict) -> dict:
    from annotations import clinvar_lookup, gene_at, gnomad_af

    out = {"genes": [], "gnomad_popmax": None, "clinvar": None}
    if have["genes"]:
        out["genes"] = [g.gene_symbol for g in gene_at(chrom, int(pos))]
    if have["gnomad"]:
        g = gnomad_af(chrom, int(pos), ref, alt)
        out["gnomad_popmax"] = None if g is None else g.popmax
    if have["clinvar"]:
        c = clinvar_lookup(chrom, int(pos), ref, alt)
        out["clinvar"] = None if c is None else c.clin_sig
    return out


def _num(v):
    if v is None:
        return None
    f = float(v)
    return int(f) if f.is_integer() else f


def _trio_json(analysis: TrioAnalysis, category: str, limit: int) -> dict:
    trio, gates = analysis.trio, analysis.gates
    has_ref = analysis.has_reference
    have = _annotation_status()

    def candidate(row: list) -> dict:
        chrom, pos, ref, alt, pgt, fgt, mgt, af = row
        return {
            "chrom": chrom,
            "pos": int(pos),
            "ref": ref,
            "alt": alt,
            "proband_gt": int(pgt),
            "father_gt": int(fgt),
            "mother_gt": int(mgt),
            "af": _num(af),
            **_enrich(chrom, pos, ref, alt, have),
        }

    models = {}
    for cat in ("denovo", "recessive", "dominant"):
        if category not in ("all", cat):
            continue
        rows = analysis.detail_rows(cat)
        models[cat] = {
            "count": len(rows),
            "blocked": cat in _NEEDS_REF and not has_ref,
            "truncated": len(rows) > limit,
            "candidates": [candidate(r) for r in rows[:limit]],
        }
    if category in ("all", "comphet"):
        genes = analysis.comphet_genes()
        entries = []
        for sym in sorted(genes)[:limit]:
            entry = {"gene": sym}
            for origin in ("paternal", "maternal"):
                entry[origin] = [
                    {"chrom": c, "pos": int(p), "ref": r, "alt": a, "af": _num(af)}
                    for c, p, r, a, af in genes[sym][origin]
                ]
            entries.append(entry)
        models["comphet"] = {
            "count": len(genes),
            "blocked": not has_ref,
            "needs_genes": not have["genes"],
            "truncated": len(genes) > limit,
            "genes": entries,
        }
    return {
        "trio": {
            "ingest_id": trio.ingest_id,
            "proband": trio.proband,
            "father": trio.father,
            "mother": trio.mother,
        },
        "keep_reference": has_ref,
        "gates": {
            "min_gq": gates.min_gq,
            "min_dp": gates.min_dp,
            "max_af": gates.max_af,
            "min_ab": gates.min_ab,
            "max_ab": gates.max_ab,
            "gnomad_max_af": analysis.gnomad_max_af,
        },
        "annotations": have,
        "models": models,
    }


def trio_text(analysis: TrioAnalysis, category: str) -> None:
    trio = analysis.trio
    proband, father, mother = trio.proband, trio.father, trio.mother
    has_ref = analysis.has_reference
    click.echo(f"trio: proband={proband} father={father} mother={mother}")
    if not has_ref and (category in _NEEDS_REF or category == "all"):
        click.echo(
            "  note: this database has no stored hom-reference calls, so "
            "de-novo/dominant/comphet cannot prove a parent is 0/0 (vs "
            "no-call). Re-ingest with `--keep-reference` for those models.",
            err=True,
        )

    if category == "all":
        for cat in ["denovo", "recessive", "dominant"]:
            n = analysis.count(cat)
            blocked = (
                "" if has_ref or cat not in _NEEDS_REF else "  (needs --keep-reference)"
            )
            click.echo(f"  {cat:10s} {n:>6}{blocked}")
        genes = analysis.comphet_genes()
        blocked = "" if has_ref else "  (needs --keep-reference)"
        click.echo(f"  {'comphet':10s} {len(genes):>6}{blocked}  genes")
        return

    if category == "comphet":
        genes = analysis.comphet_genes()
        click.echo(f"\ncomphet candidate genes: {len(genes)}")
        for sym in sorted(genes):
            entry = genes[sym]
            click.echo(f"  {sym}")
            for origin in ("paternal", "maternal"):
                for chrom, pos, ref, alt, af in entry[origin]:
                    af_s = "NA" if af is None else f"{af}"
                    click.echo(f"    {origin:8s} {chrom}:{pos} {ref}>{alt}  AF={af_s}")
        return

    data = analysis.detail_rows(category)
    click.echo(f"\n{category} candidates: {len(data)}")
    for row in data:
        chrom, pos, ref, alt, pgt, fgt, mgt, af = row
        af_s = "NA" if af is None else f"{af}"
        click.echo(
            f"  {chrom}:{pos} {ref}>{alt}  proband_gt={pgt} "
            f"father_gt={fgt} mother_gt={mgt}  AF={af_s}"
        )
