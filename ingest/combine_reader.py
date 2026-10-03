"""Read and normalize prioritized call sets into a shared allele union."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from ingest.combine_records import (
    _PASSTHROUGH,
    CombineError,
    _format_arr,
    _gt_str,
    _hdr_fields,
    _hdr_lines_by_id,
    _is_pass,
    _project_info,
    _remap_gt,
    _sample_cell,
    _Union,
)


class AlleleNormalizer:
    """Optional reference normalization, with contig order from the FASTA index."""

    def __init__(self, reference, atomize, contigs):
        self.reference = reference
        self.atomize = atomize
        if reference is None:
            return
        if not Path(reference).exists():
            raise CombineError(f"reference not found: {reference}")
        try:
            from benchmark.normalize import atomize as atomize_fn
            from benchmark.normalize import left_align, trim
            from benchmark.reference import Reference

            self.ref = Reference(reference)
        except ImportError as e:
            raise CombineError(
                "--reference requires pyfaidx; install 'vcfclick[benchmark]'."
            ) from e
        self.left_align = left_align
        self.trim = trim
        self.atomize_fn = atomize_fn
        fai = Path(f"{reference}.fai")
        if fai.exists():
            for line in fai.read_text().splitlines():
                if line.strip():
                    contigs.setdefault(line.split("\t")[0], len(contigs))

    def rows(self, variant, path):
        alts = list(variant.ALT)
        if self.reference is None:
            if len(alts) > 1:
                raise CombineError(
                    f"{path} has a multi-allelic site at "
                    f"{variant.CHROM}:{variant.POS}. Decompose first "
                    f"(bcftools norm -m -), or pass --reference to split "
                    f"internally."
                )
            alt = alts[0] if alts else "."
            return [((variant.CHROM, variant.POS, variant.REF, alt), None)]
        if not alts:
            # Monomorphic records still count at their position in site mode.
            return [((variant.CHROM, variant.POS, variant.REF, "."), None)]
        rows = []
        for j, alt in enumerate(alts, start=1):
            npos, nref, nalt = self.normalize(variant, alt)
            alleles = (
                self.atomize_fn(npos, nref, nalt)
                if self.atomize
                else [(npos, nref, nalt)]
            )
            rows.extend(
                ((variant.CHROM, pos, ref, alt), j) for pos, ref, alt in alleles
            )
        return rows

    def normalize(self, variant, alt):
        try:
            return self.left_align(
                self.ref.fetch, variant.CHROM, variant.POS, variant.REF, alt
            )
        except ValueError:
            # Unshiftable alleles retain the minimal-trim representation.
            return self.trim(variant.POS, variant.REF, alt)


@dataclass
class CombinedInputs:
    """Accumulated records, sample order, and carried header definitions."""

    carry_info: bool
    union: _Union = field(
        default_factory=lambda: _Union({}, {}, {}, {}, {}, {}, {}, set())
    )
    samples: list[str] = field(default_factory=list)
    seen_samples: set[str] = field(default_factory=set)
    info_headers: dict[str, str] = field(default_factory=dict)
    filter_headers: dict[str, str] = field(default_factory=dict)
    info_numbers: dict[str, str] = field(default_factory=dict)

    def read(self, paths, normalizer):
        from cyvcf2 import VCF

        for idx, path in enumerate(paths):
            vcf = VCF(str(path))
            try:
                samples = list(vcf.samples)
                self.add_header(vcf, samples)
                for variant in vcf:
                    rows = normalizer.rows(variant, path)
                    self.add_record(variant, rows, samples, idx)
            finally:
                vcf.close()

    def add_header(self, vcf, samples):
        if self.carry_info:
            self.add_info_headers(vcf.raw_header)
            for name, line in _hdr_lines_by_id(vcf.raw_header, "FILTER").items():
                self.filter_headers.setdefault(name, line)
        for sample in samples:
            if sample not in self.seen_samples:
                self.seen_samples.add(sample)
                self.samples.append(sample)
        # Header sequence order precedes the first appearance of variants.
        for contig in vcf.seqnames:
            self.union.contigs.setdefault(contig, len(self.union.contigs))

    def add_info_headers(self, raw_header):
        for name, line in _hdr_lines_by_id(raw_header, "INFO").items():
            if name not in self.info_headers:
                self.info_headers[name] = line
                self.info_numbers[name] = _hdr_fields(line)[0]
            elif _hdr_fields(self.info_headers[name]) != _hdr_fields(line):
                raise CombineError(
                    f"incompatible INFO header for ID={name!r} across "
                    f"inputs: '{self.info_headers[name]}' vs '{line}'"
                )

    def add_record(self, variant, rows, samples, idx):
        genotypes = variant.genotypes
        fmt_arrs = {f: _format_arr(variant, f) for f in _PASSTHROUGH}
        self.union.fields.update(f for f, arr in fmt_arrs.items() if arr is not None)
        record_pass = _is_pass(variant)
        for key, alt_index in rows:
            self.add_presence(key, idx, record_pass)
            if self.carry_info and key not in self.union.meta:
                cols = str(variant).rstrip("\n").split("\t")
                self.union.meta[key] = {
                    "qual": cols[5],
                    "filter": cols[6],
                    "info": _project_info(cols[7], alt_index, self.info_numbers),
                }
            gts = self.union.gts.setdefault(key, {})
            for sample_idx, sample in enumerate(samples):
                if sample in gts or sample_idx >= len(genotypes):
                    continue
                gt = (
                    _gt_str(genotypes[sample_idx])
                    if alt_index is None
                    else _remap_gt(genotypes[sample_idx], alt_index)
                )
                if gt is not None:
                    gts[sample] = _sample_cell(
                        gt, variant, fmt_arrs, sample_idx, alt_index
                    )

    def add_presence(self, key, idx, record_pass):
        union = self.union
        union.contigs.setdefault(key[0], len(union.contigs))
        spos = (key[0], key[1])
        union.inputs.setdefault(key, set()).add(idx)
        union.site_inputs.setdefault(spos, set()).add(idx)
        if record_pass:
            union.passes.setdefault(key, set()).add(idx)
            union.site_passes.setdefault(spos, set()).add(idx)
