"""Generate the scope fixtures (deterministic, hand-planted, no randomness).

    python tests/fixtures/make_scope_vcf.py

Writes, next to this file:

  scope.vcf.gz (+ .tbi)       one trio (CHILD male, FATHER, MOTHER) on chr7 and
                              chrX, with candidates of every model inside and
                              outside a small gene panel, and region edges
  scope_bare.vcf.gz (+ .tbi)  the same sites on contigs `7` and `X`
  scope.ped                   the trio
  scope_genes.gff3            a tiny GENCODE-style GFF3: the panel genes (with
                              the lower-case letters HGNC gives them), a gene
                              that overlaps one of them, genes outside, a PAR
                              gene on chrX and chrY, an X gene
  hgnc_aliases.tsv            a small HGNC complete-set extract for the panel
  scope_truth.json            what each site is, written by hand below

The sites pass every default gate (GQ 99, DP 40, allele balance 0.5, INFO/AF
0.0005), so the only thing that removes one is the scope. Ingest with
--keep-reference (de novo, dominant and compound het need the parents' 0/0).
The gene names and ids are those of real genes; the coordinates are not.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent.parent))

from tests.fixtures._vcfgen import Site, write_vcf  # noqa: E402

SAMPLES = ["CHILD", "FATHER", "MOTHER"]
PED = (
    "FAM1\tCHILD\tFATHER\tMOTHER\t1\t2\n"
    "FAM1\tFATHER\t0\t0\t1\t1\n"
    "FAM1\tMOTHER\t0\t0\t2\t1\n"
)

HET = "0/1:20,20:99:40"
REF = "0/0:40,0:99:40"
HOM = "1/1:0,40:99:40"
HEMI = "1/1:0,30:99:30"  # a son's hemizygous ALT on non-PAR X

# (name, chrom, pos, child, father, mother, model)
SITES = [
    ("dn_c9", "chr7", 15_000, HET, REF, REF, "denovo"),
    ("dn_before_edge", "chr7", 29_999, HET, REF, REF, "denovo"),
    ("dn_at_start", "chr7", 30_000, HET, REF, REF, "denovo"),
    ("dn_at_end", "chr7", 40_000, HET, REF, REF, "denovo"),
    ("dn_after_edge", "chr7", 40_001, HET, REF, REF, "denovo"),
    ("dn_other", "chr7", 75_000, HET, REF, REF, "denovo"),
    ("dom_c9", "chr7", 18_000, HET, HET, REF, "dominant"),
    ("ch_c19_pat", "chr7", 32_000, HET, HET, REF, "comphet"),
    ("ch_ovl_pat", "chr7", 35_500, HET, HET, REF, "comphet"),
    ("ch_ovl_mat", "chr7", 35_800, HET, REF, HET, "comphet"),
    ("ch_c19_mat", "chr7", 38_000, HET, REF, HET, "comphet"),
    ("rec_c12", "chr7", 52_000, HOM, HET, HET, "recessive"),
    ("ch_c12_pat", "chr7", 55_000, HET, HET, REF, "comphet"),
    ("ch_c12_mat", "chr7", 59_000, HET, REF, HET, "comphet"),
    ("ch_other_pat", "chr7", 72_000, HET, HET, REF, "comphet"),
    ("ch_other_mat", "chr7", 78_000, HET, REF, HET, "comphet"),
    ("dom_out", "chr7", 85_000, HET, HET, REF, "dominant"),
    ("rec_out", "chr7", 90_000, HOM, HET, HET, "recessive"),
    ("x_in", "chrX", 50_050_000, HEMI, REF, HET, "xlinked"),
    ("x_out", "chrX", 60_000_000, HEMI, REF, HET, "xlinked"),
]

# (symbol, chrom, start, end, strand, Ensembl gene id)
GENES = [
    ("C9orf72", "chr7", 10_000, 20_000, "-", "ENSG00000147894.16"),
    ("C19orf12", "chr7", 30_000, 40_000, "-", "ENSG00000131943.19"),
    ("OVL1", "chr7", 35_000, 36_000, "+", "ENSG00000000101.1"),
    ("C12orf57", "chr7", 50_000, 60_000, "+", "ENSG00000111678.11"),
    ("OTHER1", "chr7", 70_000, 80_000, "+", "ENSG00000000102.1"),
    ("XGENE1", "chrX", 50_000_000, 50_100_000, "+", "ENSG00000000103.1"),
    ("XOTHER1", "chrX", 59_950_000, 60_050_000, "+", "ENSG00000000104.1"),
    # a PAR1 gene, on chrX first and chrY second: the loader keeps the first
    ("XG", "chrX", 2_700_000, 2_760_000, "+", "ENSG00000124343.14"),
    ("XG", "chrY", 2_700_000, 2_760_000, "+", "ENSG00000124343.14_PAR_Y"),
]

HGNC_COLUMNS = (
    "hgnc_id", "symbol", "name", "status", "prev_symbol", "alias_symbol",
    "entrez_id", "ensembl_gene_id",
)  # fmt: skip
HGNC_ROWS = [
    ("HGNC:28337", "C9orf72", "C9orf72-SMCR8 complex subunit", "Approved",
     "FTDALS1", "DENNL72|ALSFTD", "203228", "ENSG00000147894"),
    ("HGNC:25443", "C19orf12", "chromosome 19 open reading frame 12", "Approved",
     '"NBIA4|SPG43"', "SHARED1", "83636", "ENSG00000131943"),
    ("HGNC:29150", "C12orf57", "chromosome 12 open reading frame 57", "Approved",
     "", "C10|SHARED1", "113246", "ENSG00000111678"),
    ("HGNC:99999", "WITHDRAWN1", "a withdrawn entry", "Entry Withdrawn",
     "C9orf72", "", "", ""),
    # approved in HGNC, absent from the gene store, and also listed as a
    # previous symbol of C19orf12 above: it must not resolve to C19orf12
    ("HGNC:90001", "SPG43", "an approved gene the store lacks", "Approved",
     "", "", "", ""),
]  # fmt: skip

PANEL = ["C9orf72", "C19orf12", "C12orf57"]
# the region the edge tests use (inclusive, 1-based)
EDGE_REGION = "chr7:30000-40000"
# what the panel keeps, by hand (the genes span their whole start..end)
IN_PANEL = {
    "denovo": ["dn_c9", "dn_at_start", "dn_at_end"],
    "recessive": ["rec_c12"],
    "dominant": ["dom_c9", "ch_c19_pat", "ch_ovl_pat", "ch_ovl_mat", "ch_c19_mat",
                 "ch_c12_pat", "ch_c12_mat"],
    "comphet": ["C12orf57", "C19orf12"],  # sorted, as the report lists genes
    "xlinked": [],
}  # fmt: skip


def _sites(bare: bool) -> list[Site]:
    out = []
    for _, chrom, pos, child, father, mother, _ in sorted(
        SITES, key=lambda s: (s[1], s[2])
    ):
        name = chrom[3:] if bare else chrom
        out.append(
            Site(name, pos, "A", "G", (child, father, mother), info="AF=0.0005",
                 fmt="GT:AD:GQ:DP")
        )  # fmt: skip
    return out


def _gff() -> str:
    lines = ["##gff-version 3"]
    for symbol, chrom, start, end, strand, gene_id in GENES:
        attrs = f"ID={gene_id};gene_id={gene_id};gene_type=protein_coding;gene_name={symbol}"
        lines.append(
            "\t".join(
                [chrom, "HAVANA", "gene", str(start), str(end), ".", strand, ".", attrs]
            )
        )
    return "\n".join(lines) + "\n"


def _bgzip(path: Path) -> None:
    gz = path.with_name(path.name + ".gz")
    with open(gz, "wb") as fh:
        subprocess.run(["bgzip", "-c", str(path)], stdout=fh, check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", str(gz)], check=True)
    path.unlink()


def main() -> None:
    for name, bare in (("scope", False), ("scope_bare", True)):
        contigs = ("7", "X") if bare else ("chr7", "chrX")
        _bgzip(write_vcf(HERE / f"{name}.vcf", SAMPLES, _sites(bare), contigs=contigs))
    (HERE / "scope.ped").write_text(PED)
    (HERE / "scope_genes.gff3").write_text(_gff())
    tsv = ["\t".join(HGNC_COLUMNS)] + ["\t".join(r) for r in HGNC_ROWS]
    (HERE / "hgnc_aliases.tsv").write_text("\n".join(tsv) + "\n")
    truth = {
        "sites": {n: {"chrom": c, "pos": p, "model": m} for n, c, p, *_, m in SITES},
        "panel": PANEL,
        "edge_region": EDGE_REGION,
        "in_panel": IN_PANEL,
        "note": "hand-written; see make_scope_vcf.py for the reason of each site",
    }
    (HERE / "scope_truth.json").write_text(json.dumps(truth, indent=2) + "\n")
    print(f"wrote {len(SITES)} sites x {len(SAMPLES)} samples, {len(GENES)} genes")


if __name__ == "__main__":
    main()
