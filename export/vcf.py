"""Export a database, or a slice of it, back to VCF (`vcfclick db export`).

One ingestion per file: different ingestions carry different samples and
site lists, so merging them would invent genotypes for samples that were
never genotyped at a site.

Genotypes live in a sparse table, so reconstructing a VCF needs one rule
for samples with no row at a site:
- ingestion made with --keep-reference: 0/0 calls are stored, so absent
  means a no-call and is written ./.;
- otherwise only non-reference calls were stored and a hom-ref call cannot
  be told apart from a no-call. Absent is written 0/0 by default
  (--absent-as nocall writes ./. instead), and the header says so.

Exported per genotype: GT, GQ, DP, AD, FT. Phase is not preserved: the
store keeps a phased flag but not which haplotype carries the ALT, so all
calls are written unphased. PL/GL are not exported. INFO is rebuilt from
the typed info_* columns plus the info_extra map. FILTER PASS and unset (.)
are both stored as NULL at ingest, so both come back as ".".

.vcf.gz output is BGZF (tabix-indexable); .vcf is plain text; "-" is stdout.
"""

from __future__ import annotations

import tempfile

from .vcf_context import FORMAT_FIELDS, FORMAT_PLACEHOLDER
from .vcf_context import ExportError as ExportError
from .vcf_context import ExportRequest
from .vcf_format import header, record
from .vcf_query import iter_records
from .vcf_query import parse_region as parse_region
from .vcf_query import prepare_export
from .vcf_query import read_bed as read_bed
from .vcf_writer import BgzfWriter as BgzfWriter
from .vcf_writer import write_vcf


def export_vcf(
    sess,
    db_name: str,
    out: str,
    *,
    ingest_id: str | None = None,
    samples: list[str] | None = None,
    regions: list[tuple[str, int | None, int | None]] | None = None,
    where: str | None = None,
    pass_only: bool = False,
    sites_only: bool = False,
    absent_as: str = "ref",
    version: str = "",
) -> dict:
    """Write the VCF. Returns a summary dict (variants, samples, path)."""
    request = ExportRequest(
        db_name=db_name,
        ingest_id=ingest_id,
        samples=samples,
        regions=regions,
        where=where,
        pass_only=pass_only,
        sites_only=sites_only,
        absent_as=absent_as,
        version=version,
    )
    context = prepare_export(sess, request)
    n_variants = 0
    with tempfile.TemporaryFile("w+", encoding="utf-8") as body:
        for rec, calls in iter_records(sess, context):
            body.write(record(rec, calls, context))
            n_variants += 1
        used_fmt = [
            field
            for field in FORMAT_FIELDS
            if field == "GT" or field in context.fmt_seen
        ]
        text_header = header(context, used_fmt if context.samples else [])
        write_vcf(out, text_header, body, FORMAT_PLACEHOLDER, used_fmt)
    return {
        "path": out,
        "ingest_id": context.ingest_id,
        "variants": n_variants,
        "samples": len(context.samples),
        "absent_written_as": context.absent_gt,
    }
