"""`vcfclick discover` — preview which VCF fields will land where.

For any new VCF (different lab, different pipeline, different annotator),
running

    vcfclick discover input.vcf.gz

reports — without ingesting anything — which INFO/FORMAT fields will
end up in typed columns and which will land in the info_extra /
format_extra overflow Maps, plus suggested ClickHouse types for any
overflow field that could be promoted.

This makes the "adapts to any VCF" claim concrete: you can see, in
advance, the schema-routing decisions the ingester will make.
"""

from __future__ import annotations

import click
from cyvcf2 import VCF

from cli.main import cli
from ingest.vcf_load import (
    FORMAT_PAIR,
    FORMAT_SCALAR,
    FORMAT_TRIPLE,
    INFO_FLAG,
    INFO_PAIR,
    INFO_SCALAR,
)


# Map of VCF Type → suggested ClickHouse base type. The Integer mapping
# uses Int32 (signed) because VCF spec allows negatives (e.g., some
# imputation scores). Float → Float32 is wide enough for QC stats.
_BASE_TYPE = {
    "Integer": "Int32",
    "Float": "Float32",
    "String": "String",
    "Character": "String",
}


def _promotion_hint(vcf_type: str, number: str, prefix: str, field: str) -> str | None:
    """Return a suggested DDL fragment for a promotable overflow field,
    or None if the field can't be cleanly promoted (variable-length lists,
    unknown types). `prefix` is "info_" for variants, "" for genotypes.

    Column-naming convention follows the existing schema: INFO-derived
    columns keep VCF-style casing (info_AC, info_MQ) while FORMAT-derived
    columns are lowercased (gq, dp, ad_ref) — match that here.
    """
    name = field if prefix else field.lower()

    if vcf_type == "Flag":
        return f"{prefix}{name} UInt8 DEFAULT 0"

    base = _BASE_TYPE.get(vcf_type)
    if base is None:
        return None

    # Number=. → variable-length list. Can't typed-column it cleanly;
    # leave in the overflow Map as a stringified value.
    if number == ".":
        return None

    nullable = f"Nullable({base})"
    if number in ("1", "0", "A"):
        # A scalar per row (we decompose multi-allelic, so Number=A is a scalar).
        return f"{prefix}{name} {nullable}"
    if number == "R":
        # ref + alt pair.
        return f"{prefix}{name}_ref {nullable}, {prefix}{name}_alt {nullable}"
    if number == "G":
        # ref/ref, ref/alt, alt/alt triple (biallelic).
        return (
            f"{prefix}{name}_ref_ref {nullable}, "
            f"{prefix}{name}_ref_alt {nullable}, "
            f"{prefix}{name}_alt_alt {nullable}"
        )
    # Fixed-size lists with k != 1,2,3 — uncommon, skip.
    return None


def _iter_header(vcf: VCF, kind: str):
    """Yield (id, type, number, description) for INFO or FORMAT entries."""
    for h in vcf.header_iter():
        try:
            d = h.info(extra=True)
        except Exception:
            continue
        if str(d.get("HeaderType", "")).upper() != kind:
            continue
        fid = d.get("ID")
        if not fid:
            continue
        yield (
            fid,
            d.get("Type", "String"),
            str(d.get("Number", ".")),
            d.get("Description", "").strip('"').strip(),
        )


def _print_section(
    section_name: str,
    container: str,
    typed_ids: set[str],
    fields: list[tuple[str, str, str, str]],
    prefix: str,
) -> int:
    """Print one INFO or FORMAT section. Returns # overflow fields found."""
    typed = sorted({f[0] for f in fields if f[0] in typed_ids})
    overflow = sorted([f for f in fields if f[0] not in typed_ids])

    click.echo(f"{section_name} fields:")
    click.echo(f"  typed   ({len(typed):>3}) → {container}.<typed columns>")
    if typed:
        click.echo(f"    {', '.join(typed)}")
    click.echo(f"  overflow ({len(overflow):>3}) → {container}.{section_name.lower()}_extra Map")
    if not overflow:
        click.echo()
        return 0

    promotable = 0
    for fid, vtype, number, desc in overflow:
        hint = _promotion_hint(vtype, number, prefix, fid)
        if hint:
            promotable += 1
            click.echo(f"    {fid:25s} {vtype:8s} Number={number:<3s} → promote: {hint}")
        else:
            click.echo(f"    {fid:25s} {vtype:8s} Number={number:<3s}   (variable-length / unknown — keep in Map)")
        if desc:
            short = desc if len(desc) <= 80 else desc[:77] + "..."
            click.echo(f"      └ {short}")
    click.echo()
    return promotable


@cli.command()
@click.argument("vcf_path", type=click.Path(exists=True, dir_okay=False))
def discover(vcf_path: str) -> None:
    """Report which fields of VCF_PATH will land in typed columns vs the
    info_extra / format_extra overflow Maps. Suggests promotion DDL for
    overflow fields that have a fixed shape."""
    vcf = VCF(vcf_path)

    click.echo(f"VCF:       {vcf_path}")
    click.echo(f"Samples:   {len(vcf.samples):,}")
    click.echo()

    info_typed_ids = set(INFO_SCALAR) | set(INFO_PAIR) | set(INFO_FLAG)
    info_fields = list(_iter_header(vcf, "INFO"))
    n_info_prom = _print_section(
        "INFO", "variants", info_typed_ids, info_fields, prefix="info_"
    )

    format_typed_ids = {"GT"} | set(FORMAT_SCALAR) | set(FORMAT_PAIR) | set(FORMAT_TRIPLE)
    format_fields = list(_iter_header(vcf, "FORMAT"))
    n_fmt_prom = _print_section(
        "FORMAT", "genotypes", format_typed_ids, format_fields, prefix=""
    )

    if n_info_prom + n_fmt_prom > 0:
        click.echo("To promote an overflow field to a typed column:")
        click.echo("  1. Add the column(s) to schema/01_variants.sql or 02_genotypes.sql.")
        click.echo("  2. Add a routing entry to ingest/vcf_load.py:")
        click.echo("       INFO_SCALAR['<FIELD>'] = 'info_<field>'")
        click.echo("       (or INFO_PAIR / INFO_FLAG / FORMAT_SCALAR / FORMAT_PAIR / FORMAT_TRIPLE)")
        click.echo("  3. Re-create the DB and re-ingest. (In-place ALTER works but back-fills NULL.)")
