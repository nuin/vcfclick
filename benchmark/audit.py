"""Auditability layer — surface per-error (FN/FP) variants with their annotation
context so a disagreement can be understood, not just counted. hap.py gives an
opaque intermediate VCF; this joins each error to gene / ClinVar / gnomAD.

Query layer (DuckDB); the benchmark core stays SQL-free.
"""

from __future__ import annotations

import duckdb

# An FN is a truth-side miss; an FP is a query-side spurious call.
_SIDE = {"FN": "truth", "FP": "query"}


def annotated_errors(concordance, ann_path: str, kind: str = "FN") -> list[dict]:
    """Return the `kind` (FN or FP) errors from the concordance frame, each joined
    to its gene, ClinVar significance, and gnomAD AF (NULL where unannotated)."""
    side = _SIDE.get(kind)
    if side is None:
        raise ValueError(f"kind must be 'FN' or 'FP', got {kind!r}")

    # `side`/`kind` are validated against _SIDE above (not injectable). Open the
    # annotation store directly read-only rather than ATTACH an interpolated path.
    sql = f"""
        SELECT c.chrom, c.pos, c.ref, c.alt, c.vtype,
               gn.gene_symbol AS gene,
               v.clin_sig     AS clin_sig,
               g.af           AS af
        FROM conc c
        LEFT JOIN refseq_genes gn
            ON c.chrom=gn.chrom AND c.pos BETWEEN gn.start_pos AND gn.end_pos
        LEFT JOIN clinvar_variants v
            ON c.chrom=v.chrom AND c.pos=v.pos AND c.ref=v.ref AND c.alt=v.alt
        LEFT JOIN gnomad_af g
            ON c.chrom=g.chrom AND c.pos=g.pos AND c.ref=g.ref AND c.alt=g.alt
        WHERE c.side='{side}' AND c.bd='{kind}' AND c.filter_view='ALL'
        ORDER BY c.chrom, c.pos
    """
    con = duckdb.connect(ann_path, read_only=True)
    try:
        con.register("conc", concordance)
        cur = con.execute(sql)
        cols = [d[0] for d in cur.description]
        rows = cur.fetchall()
    finally:
        con.close()
    return [dict(zip(cols, r)) for r in rows]
