"""Record-for-record emulation of GATK 3.8-1 CombineVariants.

`vcfclick combine --gatk3` runs this engine instead of the normalising one.
It follows GATKVariantContextUtils.simpleMerge and CombineVariants.map from
GATK 3.8-1 (and htsjdk 2.12.0, the VCF library inside that jar) step by
step, so raw caller output merges the way the Java tool merges it:

* No normalisation. Records merge only when they start at the same
  position. Multi-allelic records stay multi-allelic. Differently spelled
  indels stay separate records.
* Records at one position are split by variant type (SNP, MNP, INDEL,
  SYMBOLIC, MIXED), with GATK's allele-subset rule for records of different
  types. Each type group becomes one output record.
* `--min-callsets N` is GATK's `--minimumN`: records at the position minus
  records with a FILTER other than PASS or `.` must reach N. It counts
  records, not inputs, so one input can satisfy it with a SNP and an indel.
* The genotype of a sample is the one in the first record (input order, or
  `--priority` order) that has the sample. A `./.` there wins over a later
  called genotype. UNSORTED and PRIORITIZE in input order are the same.
* QUAL is the first non-missing QUAL. FILTER is PASS when any record is not
  filtered, else the sorted union. ID is the comma-joined union.
* With `--merge-info-with-max-ac` INFO is that of the record with the
  highest INFO AC (first record wins a tie), then AC/AF/AN are recomputed
  from the merged genotypes. When no record has an AC, INFO holds only the
  recomputed AC/AF/AN and the `set` key is absent.
* Header lines are merged in Java HashMap order of the input names, as
  htsjdk's smartMergeHeaders does, and sorted like a TreeSet of lines.

Differences from the Java tool that are not bugs: a malformed numeric field
raises CombineError where Java raises NumberFormatException, and the
`##GATKCommandLine` and `##reference` lines are not written.
"""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass, field
from decimal import ROUND_HALF_UP, Decimal
from math import floor
from pathlib import Path

from ingest.combine_records import CombineError

# VariantContext.Type enum order (htsjdk 2.12.0).
TYPES = ("NO_VARIATION", "SNP", "MNP", "INDEL", "SYMBOLIC", "MIXED")

SET_KEY = "set"
MAX_GENOTYPE_QUAL = 99

# htsjdk VCFStandardHeaderLines (2.12.0): id -> (count, type, description).
# count is an int, "A", "R", "G" or ".".
STANDARD_FORMAT = {
    "GT": ("1", "String", "Genotype"),
    "GQ": ("1", "Integer", "Genotype Quality"),
    "DP": (
        "1",
        "Integer",
        "Approximate read depth (reads with MQ=255 or with bad mates are filtered)",
    ),
    "PL": (
        "G",
        "Integer",
        "Normalized, Phred-scaled likelihoods for genotypes as defined in the "
        "VCF specification",
    ),
    "AD": (
        "R",
        "Integer",
        "Allelic depths for the ref and alt alleles in the order listed",
    ),
    "FT": (".", "String", "Genotype-level filter"),
    "PQ": ("1", "Float", "Read-backed phasing quality"),
}
STANDARD_INFO = {
    "END": ("1", "Integer", "Stop position of the interval"),
    "DB": ("0", "Flag", "dbSNP Membership"),
    "DP": ("1", "Integer", "Approximate read depth; some reads may have been filtered"),
    "SB": ("1", "Float", "Strand Bias"),
    "AF": (
        "A",
        "Float",
        "Allele Frequency, for each ALT allele, in the same order as listed",
    ),
    "AC": (
        "A",
        "Integer",
        "Allele count in genotypes, for each ALT allele, in the same order as listed",
    ),
    "AN": ("1", "Integer", "Total number of alleles in called genotypes"),
    "MQ0": ("1", "Integer", "Total Mapping Quality Zero Reads"),
    "MQ": ("1", "Float", "RMS Mapping Quality"),
    "SOMATIC": ("0", "Flag", "Somatic event"),
}
_STANDARD = {"FORMAT": STANDARD_FORMAT, "INFO": STANDARD_INFO}


# --- Java behaviours -------------------------------------------------------


def java_string_hash(s: str) -> int:
    h = 0
    for ch in s:
        h = (h * 31 + ord(ch)) & 0xFFFFFFFF
    return h


def java_hashmap_order(hashes: list[int]) -> list[int]:
    """Indices of keys (given in insertion order, as unsigned 32-bit hash
    codes) in the iteration order of a java.util.HashMap holding them."""
    cap = 16
    while len(hashes) > cap * 0.75:
        cap *= 2

    def bucket(h: int) -> int:
        return (h ^ (h >> 16)) & (cap - 1)

    return sorted(range(len(hashes)), key=lambda i: bucket(hashes[i]))


def java_round(x: float) -> int:
    return floor(x + 0.5)


def java_format(x: float, places: int) -> str:
    """String.format("%.<places>f", x): HALF_UP on the shortest repr digits."""
    q = Decimal(1).scaleb(-places)
    return str(Decimal(repr(float(x))).quantize(q, rounding=ROUND_HALF_UP))


def format_vcf_double(d: float) -> str:
    """htsjdk VCFEncoder.formatVCFDouble."""
    if d < 1:
        if d < 0.01:
            if abs(d) >= 1e-20:
                return "%.3e" % d
            return "0.00"
        return java_format(d, 3)
    return java_format(d, 2)


def format_qual(log10_p_error: float) -> str:
    s = java_format(-10.0 * log10_p_error, 2)
    return s[:-3] if s.endswith(".00") else s


def allele_hash(bases: str) -> int:
    h = 1
    for i, ch in enumerate(bases):
        h += (i + 1) * ord(ch)
    return h & 0xFFFFFFFF


def _is_symbolic(allele: str) -> bool:
    """htsjdk Allele.wouldBeSymbolicAllele: `*` (length 1) is not symbolic, so
    it is typed by length like any other allele."""
    return len(allele) > 1 and (
        allele.startswith("<") or allele.endswith(">") or "[" in allele or "]" in allele
    )


# --- records ---------------------------------------------------------------

Allele = tuple  # (bases, is_ref)


@dataclass
class Genotype:
    alleles: list  # list of bases str or None (no-call)
    phased: bool
    gq: int | None = None
    dp: str | None = None
    ad: str | None = None
    pl: str | None = None
    attrs: dict = field(default_factory=dict)


@dataclass
class Vc:
    source: str
    src_idx: int
    chrom: str
    pos: int
    id: str
    ref: str
    alts: list
    qual: float | None  # as log10 P(error), None when missing
    filters: list | None  # None: not applied (FILTER '.'); [] PASS
    info: dict
    genotypes: dict  # sample -> Genotype

    @property
    def alleles(self) -> list:
        return [(self.ref, True)] + [(a, False) for a in self.alts]

    @property
    def filtered(self) -> bool:
        return bool(self.filters)

    @property
    def type(self) -> str:
        if not self.alts:
            return "NO_VARIATION"
        t = None
        for a in self.alts:
            if _is_symbolic(a):
                bt = "SYMBOLIC"
            elif len(a) == len(self.ref):
                bt = "SNP" if len(a) == 1 else "MNP"
            else:
                bt = "INDEL"
            if t is None:
                t = bt
            elif bt != t:
                return "MIXED"
        return t


# --- header ----------------------------------------------------------------

_KV = re.compile(r'(\w+)=("(?:[^"\\]|\\.)*"|[^,>]*)')


def _unquote(value: str) -> str:
    return value[1:-1] if len(value) >= 2 and value[0] == value[-1] == '"' else value


def _parse_structured(line: str) -> dict:
    inner = line[line.index("<") + 1 : line.rindex(">")]
    return {k: v for k, v in _KV.findall(inner)}


@dataclass
class HeaderLine:
    kind: str  # INFO, FORMAT, FILTER, contig, other
    key: str  # smartMerge map key
    text: str
    id: str = ""
    count: str = ""
    type: str = ""
    desc: str = ""

    def compound_text(self) -> str:
        return (
            f"##{self.kind}=<ID={self.id},Number={self.count},Type={self.type},"
            f'Description="{self.desc}">'
        )

    def same_excluding_description(self, o: HeaderLine) -> bool:
        return (self.count, self.type, self.kind) == (o.count, o.type, o.kind)


def _repair_standard(line: HeaderLine) -> HeaderLine:
    std = _STANDARD.get(line.kind, {}).get(line.id)
    if std is None:
        return line
    count, typ, desc = std
    fixed = count.isdigit()
    bad_count_type = (line.count.isdigit() != fixed) or (
        not fixed and line.count != count
    )
    bad_count = line.count.isdigit() and not bad_count_type and line.count != count
    if bad_count_type or bad_count or line.type != typ:
        out = HeaderLine(line.kind, line.key, "", line.id, count, typ, desc)
        out.text = out.compound_text()
        return out
    return line


def _parse_header_line(raw: str) -> HeaderLine | None:
    if not raw.startswith("##") or raw.startswith("##fileformat"):
        return None
    body = raw[2:]
    name, _, value = body.partition("=")
    if name in ("INFO", "FORMAT") and value.startswith("<"):
        f = _parse_structured(raw)
        if "ID" not in f:
            return None
        line = HeaderLine(
            name,
            f"{name}-{f['ID']}",
            "",
            f["ID"],
            f.get("Number", "."),
            f.get("Type", "String"),
            _unquote(f.get("Description", '""')),
        )
        line.text = line.compound_text()
        return _repair_standard(line)
    if name == "FILTER" and value.startswith("<"):
        f = _parse_structured(raw)
        desc = _unquote(f.get("Description", '""'))
        line = HeaderLine("FILTER", f"FILTER-{f['ID']}", "", f["ID"], "", "", desc)
        line.text = f'##FILTER=<ID={f["ID"]},Description="{desc}">'
        return line
    if name == "contig" and value.startswith("<"):
        f = _parse_structured(raw)
        return HeaderLine("contig", f"contig-{f['ID']}", raw, f["ID"])
    if value.startswith("<") and "ID=" in value:
        f = _parse_structured(raw)
        return HeaderLine("other", f"{name}-{f.get('ID', '')}", raw)
    return HeaderLine("other", name, raw)


def _merge_headers(per_input: list[list[HeaderLine]]) -> dict[str, HeaderLine]:
    """htsjdk VCFUtils.smartMergeHeaders over headers in the given order."""
    merged: dict[str, HeaderLine] = {}
    for lines in per_input:
        for line in sorted(lines, key=lambda ln: ln.text):
            other = merged.get(line.key)
            if other is None:
                merged[line.key] = line
                continue
            if line.text == other.text:
                continue
            if line.kind != other.kind:
                raise CombineError(
                    f"Incompatible header types: {line.text} {other.text}"
                )
            if line.kind in ("INFO", "FORMAT") and not line.same_excluding_description(
                other
            ):
                if line.type == other.type:
                    other.count = "."
                    other.text = other.compound_text()
                elif {line.type, other.type} <= {"Integer", "Float"}:
                    pass  # htsjdk keeps the first line, whichever type it has
                else:
                    raise CombineError(
                        "Incompatible header types, collision between these two "
                        f"types: {line.text} {other.text}"
                    )
    return merged


# --- reading ---------------------------------------------------------------


def _open_text(path: Path):
    if str(path).endswith((".gz", ".bgz")):
        return gzip.open(path, "rt")
    return open(path)


def _parse_qual(q: str) -> float | None:
    if q == ".":
        return None
    return float(q) / -10.0


def _parse_genotype(
    keys: list[str], cell: str, where: str, bases: list[str]
) -> Genotype:
    vals = cell.split(":")
    gt = Genotype([None, None], False)
    for key, val in zip(keys, vals):
        if key == "GT":
            gt.phased = "|" in val
            try:
                gt.alleles = [
                    None if a == "." else bases[int(a)] for a in re.split(r"[/|]", val)
                ]
            except (ValueError, IndexError) as e:
                raise CombineError(f"{where}: bad genotype {val!r}") from e
            continue
        if val == ".":
            continue
        try:
            if key == "GQ":
                gt.gq = java_round(float(val))
            elif key == "DP":
                gt.dp = str(int(val))
            elif key == "AD":
                gt.ad = ",".join(str(int(x)) for x in val.split(","))
            elif key == "PL":
                gt.pl = ",".join(str(int(x)) for x in val.split(","))
            elif key == "GL":
                gls = [float(x) for x in val.split(",")]
                top = max(gls)
                gt.pl = ",".join(str(java_round(-10.0 * (g - top))) for g in gls)
            else:
                gt.attrs[key] = val
        except ValueError as e:
            raise CombineError(f"{where}: cannot parse FORMAT {key}={val!r}") from e
    return gt


def _read_vcf(path: Path, src_idx: int, name: str):
    """Return (header lines, samples, records, contig order)."""
    header: list[HeaderLine] = []
    samples: list[str] = []
    records: list[Vc] = []
    contigs: list[str] = []
    with _open_text(path) as fh:
        for raw in fh:
            raw = raw.rstrip("\n")
            if raw.startswith("##"):
                hl = _parse_header_line(raw)
                if hl is not None:
                    header.append(hl)
                    if hl.kind == "contig" and hl.id not in contigs:
                        contigs.append(hl.id)
                continue
            if raw.startswith("#"):
                samples = raw.split("\t")[9:]
                continue
            if not raw:
                continue
            c = raw.split("\t")
            where = f"{path}:{c[0]}:{c[1]}"
            if len(c) < 8:
                raise CombineError(f"{where}: fewer than 8 columns")
            alts = [] if c[4] in (".", "") else [a.upper() for a in c[4].split(",")]
            flt = None if c[6] == "." else ([] if c[6] == "PASS" else c[6].split(";"))
            info: dict = {}
            if c[7] not in (".", ""):
                for item in c[7].split(";"):
                    k, eq, v = item.partition("=")
                    info[k] = v if eq else True
            gts: dict = {}
            if len(c) > 9:
                keys = c[8].split(":")
                for s, cell in zip(samples, c[9:]):
                    gts[s] = _parse_genotype(keys, cell, where, [c[3].upper(), *alts])
            records.append(
                Vc(
                    name,
                    src_idx,
                    c[0],
                    int(c[1]),
                    c[2],
                    c[3].upper(),
                    alts,
                    _parse_qual(c[5]),
                    flt,
                    info,
                    gts,
                )
            )
    return header, samples, records, contigs


# --- merging ---------------------------------------------------------------


def _subset(vc1: Vc, vc2: Vc) -> bool:
    """GATKVariantContextUtils.allelesAreSubset."""
    if vc1.ref != vc2.ref:
        return False
    return all(a in vc2.alts for a in vc1.alts)


def _separate_by_type(vcs: list[Vc]) -> dict[str, list[Vc]]:
    """GATKVariantContextUtils.separateVariantContextsByType."""
    mapped: dict[str, list[Vc]] = {}
    for vc in vcs:
        vc_type = vc.type
        add_to_own = True
        for t in TYPES:
            if t == vc_type or t not in mapped:
                continue
            lst = mapped[t]
            for k, other in enumerate(lst):
                if _subset(other, vc):
                    lst.pop(k)
                    if not lst:
                        del mapped[t]
                    mapped.setdefault(vc_type, []).append(other)
                    break
                if _subset(vc, other):
                    mapped[t].append(vc)
                    add_to_own = False
                    break
        if add_to_own:
            mapped.setdefault(vc_type, []).append(vc)
    return mapped


def _is_usable_alt(a: str) -> bool:
    return not (a == "*" or _is_symbolic(a))


def _resolve_alleles(ref_allele: str, vc: Vc, current: list[Allele]):
    """resolveIncompatibleAlleles: (mapping or None, values in Java order)."""
    if vc.ref == ref_allele:
        return None, vc.alleles
    extra = ref_allele[len(vc.ref) :]
    entries: list[tuple[Allele, Allele]] = []
    for a in vc.alts:
        if _is_usable_alt(a):
            ext = (a + extra, False)
            for b in current:
                if ext == b:
                    ext = b
            entries.append(((a, False), ext))
        else:
            entries.append(((a, False), (a, False)))
    entries.append(((vc.ref, True), (ref_allele, True)))
    order = java_hashmap_order([allele_hash(k[0]) for k, _ in entries])
    mapping = dict(entries)
    return mapping, [entries[i][1] for i in order]


def _incompatible(alleles: list[Allele], vc_alleles: list[Allele]) -> bool:
    for a, b in zip(alleles, vc_alleles):
        if a != b:
            return True
    return len(alleles) != len(vc_alleles)


@dataclass
class Merged:
    chrom: str
    pos: int
    id: str
    alleles: list  # ref first
    qual: float | None
    filters: list | None
    info: dict
    genotypes: dict


def _max_ac(vc: Vc) -> list[int] | None:
    raw = vc.info.get("AC")
    if raw is None or raw is True:
        return None
    try:
        return [int(x.strip()) for x in raw.split(",")]
    except ValueError as e:
        raise CombineError(
            f"{vc.chrom}:{vc.pos} in {vc.source}: INFO AC={raw!r} is not an integer list"
        ) from e


def _simple_merge(
    vcs: list[Vc],
    n_inputs: int,
    merge_info_with_max_ac: bool,
    filtered_are_uncalled: bool,
) -> Merged | None:
    vcs = [v for v in vcs if not filtered_are_uncalled or not v.filtered]
    if not vcs:
        return None
    first = vcs[0]
    ref_allele = first.ref
    for v in vcs:
        if len(v.ref) > len(ref_allele):
            ref_allele = v.ref
        elif len(v.ref) == len(ref_allele) and v.ref != ref_allele:
            raise CombineError(
                f"inconsistent references for the same position {v.chrom}:{v.pos}, "
                f"{ref_allele} vs. {v.ref}"
            )

    alleles: list[Allele] = []  # LinkedHashSet semantics
    filters: set[str] = set()
    attributes: dict = {}
    inconsistent: set[str] = set()
    variant_sources: set[str] = set()
    rs_ids: list[str] = []
    non_bool_counts: dict[str, int] = {}
    depth = 0
    max_ac = -1
    vc_with_max_ac: Vc | None = None
    log10_p_error: float | None = None
    any_filters_applied = False
    n_filtered = 0
    genotypes: dict[str, Genotype] = {}

    for vc in vcs:
        n_filtered += 1 if vc.filtered else 0
        if vc.alts:
            variant_sources.add(vc.source)
        mapping, values = _resolve_alleles(ref_allele, vc, alleles)
        for a in values:
            if a not in alleles:
                alleles.append(a)
        for sample, g in vc.genotypes.items():
            if sample in genotypes:
                continue
            if mapping is not None:
                g = _remap_genotype(g, vc, mapping, ref_allele)
            genotypes[sample] = g
        if log10_p_error is None and vc.qual is not None:
            log10_p_error = vc.qual
        filters.update(vc.filters or [])
        any_filters_applied |= vc.filters is not None

        if "DP" in vc.info:
            try:
                depth += int(vc.info["DP"])
            except (TypeError, ValueError) as e:
                raise CombineError(
                    f"{vc.chrom}:{vc.pos} in {vc.source}: INFO DP={vc.info['DP']!r} is not an integer"
                ) from e
        if vc.id not in (".", "") and vc.id not in rs_ids:
            rs_ids.append(vc.id)
        if merge_info_with_max_ac:
            counts = _max_ac(vc)
            for ac in counts or []:
                if ac > max_ac:
                    max_ac = ac
                    vc_with_max_ac = vc
        for key, value in vc.info.items():
            if value is not True:
                non_bool_counts[key] = non_bool_counts.get(key, 0) + 1
            if key not in inconsistent:
                found = key in attributes
                bound = attributes.get(key)
                bound_missing = found and bound == "."
                if found and bound != value and not bound_missing:
                    inconsistent.add(key)
                    del attributes[key]
                elif not found or bound_missing:
                    attributes[key] = value

    for key, n in non_bool_counts.items():
        if n < len(vcs):
            attributes.pop(key, None)

    # Allele-dependent genotype fields are stale when any record's allele
    # list differs from the merged one.
    for vc in vcs:
        if len(vc.alleles) == 1:
            continue
        if _incompatible(alleles, vc.alleles):
            for g in genotypes.values():
                g.pl = None
                g.ad = None
            break

    if merge_info_with_max_ac:
        out_attrs = dict(vc_with_max_ac.info) if vc_with_max_ac is not None else {}
    else:
        out_attrs = attributes

    if _clear_filters(n_filtered, len(vcs)):
        filters.clear()

    if n_filtered == 0 and len(variant_sources) == n_inputs:
        set_value = "Intersection"
    elif n_filtered == len(vcs):
        set_value = "FilteredInAll"
    elif not variant_sources:
        set_value = "ReferenceInAll"
    else:
        names: list[str] = []
        for vc in vcs:
            if vc.alts:
                n = ("filterIn" + vc.source) if vc.filtered else vc.source
                if n not in names:
                    names.append(n)
        set_value = "-".join(names)
    if not merge_info_with_max_ac or vc_with_max_ac is not None:
        out_attrs[SET_KEY] = set_value
    if depth > 0 and not merge_info_with_max_ac:
        out_attrs["DP"] = str(depth)

    longest = first
    # ref allele first, the other alleles in set order
    ordered = [(ref_allele, True)] + [a for a in alleles if not a[1]]
    return Merged(
        first.chrom,
        longest.pos,
        ",".join(rs_ids) if rs_ids else ".",
        ordered,
        log10_p_error,
        sorted(filters) if any_filters_applied else None,
        out_attrs,
        genotypes,
    )


def _clear_filters(n_filtered: int, n: int) -> bool:
    """KEEP_IF_ANY_UNFILTERED: drop the filters when any record is unfiltered."""
    return n_filtered != n


def _remap_genotype(g: Genotype, vc: Vc, mapping: dict, ref_allele: str) -> Genotype:
    def remap(base):
        if base is None:
            return None
        if base == vc.ref:
            return ref_allele
        m = mapping.get((base, False))
        return m[0] if m else base

    return Genotype(
        [remap(a) for a in g.alleles], g.phased, g.gq, g.dp, g.ad, g.pl, dict(g.attrs)
    )


def _genotype_index_alleles(g: Genotype, bases_to_idx: dict, where: str) -> str:
    out = []
    for a in g.alleles:
        if a is None:
            out.append(".")
        elif a in bases_to_idx:
            out.append(str(bases_to_idx[a]))
        else:
            raise CombineError(
                f"{where}: genotype allele {a} is not an allele of the merged record"
            )
    return ("|" if g.phased else "/").join(out)


def _chromosome_counts(merged: Merged, samples: list[str]) -> None:
    """VariantContextUtils.calculateChromosomeCounts(builder, removeStale=false)."""
    gts = [merged.genotypes[s] for s in samples if s in merged.genotypes]
    if not gts:
        return
    alts = [a[0] for a in merged.alleles[1:]]
    called = [a for g in gts for a in g.alleles if a is not None]
    an = len(called)
    merged.info["AN"] = str(an)
    if not alts:
        merged.info.pop("AC", None)
        merged.info.pop("AF", None)
        return
    counts = [called.count(a) for a in alts]
    freqs = [0.0 if an == 0 else c / an for c in counts]
    merged.info["AC"] = ",".join(str(c) for c in counts)
    merged.info["AF"] = ",".join(format_vcf_double(f) for f in freqs)


def _format_cell(
    g: Genotype | None, keys: list[str], bases_to_idx: dict, where: str
) -> str:
    if g is None:
        return "./."
    vals = {"GT": _genotype_index_alleles(g, bases_to_idx, where)}
    if g.gq is not None:
        vals["GQ"] = str(min(g.gq, MAX_GENOTYPE_QUAL))
    if g.dp is not None:
        vals["DP"] = g.dp
    if g.ad is not None:
        vals["AD"] = g.ad
    if g.pl is not None:
        vals["PL"] = g.pl
    vals.update(g.attrs)
    out = [vals.get(k, ".") for k in keys]
    while len(out) > 1 and out[-1] == ".":
        out.pop()
    return ":".join(out)


def _format_keys(genotypes: list[Genotype]) -> list[str]:
    keys: set[str] = set()
    for g in genotypes:
        keys.update(g.attrs)
        if g.gq is not None:
            keys.add("GQ")
        if g.dp is not None:
            keys.add("DP")
        if g.ad is not None:
            keys.add("AD")
        if g.pl is not None:
            keys.add("PL")
    return ["GT"] + sorted(keys)


# --- driver ----------------------------------------------------------------


def _read_fai(reference: Path) -> list[tuple[str, int]]:
    fai = Path(f"{reference}.fai")
    if not fai.exists():
        raise CombineError(f"reference index not found: {fai} (run samtools faidx)")
    out = []
    for line in fai.read_text().splitlines():
        if line.strip():
            f = line.split("\t")
            out.append((f[0], int(f[1])))
    return out


def combine_gatk3(
    inputs: list[Path],
    output: Path,
    names: list[str],
    *,
    min_callsets: int = 1,
    merge_info_with_max_ac: bool = False,
    filtered_are_uncalled: bool = False,
    priority: list[str] | None = None,
    reference: str | Path | None = None,
) -> int:
    """Write the GATK3 CombineVariants merge of `inputs` to plain-text `output`.

    Returns the number of records written."""
    n_inputs = len(inputs)
    if priority is not None:
        if sorted(priority) != sorted(names):
            raise CombineError("--priority must list every input name exactly once.")
    per_input: list[list[HeaderLine]] = []
    all_samples: set[str] = set()
    by_pos: dict[tuple[str, int], list[Vc]] = {}
    header_contigs: list[str] = []
    for idx, (path, name) in enumerate(zip(inputs, names)):
        header, samples, records, contigs = _read_vcf(path, idx, name)
        per_input.append(header)
        all_samples.update(samples)
        for c in contigs:
            if c not in header_contigs:
                header_contigs.append(c)
        for vc in records:
            by_pos.setdefault((vc.chrom, vc.pos), []).append(vc)

    # Headers merge in the iteration order of a HashMap keyed by input name.
    order = java_hashmap_order([java_string_hash(n) for n in names])
    merged_hdr = _merge_headers([per_input[i] for i in order])
    samples_out = sorted(all_samples)

    contig_lines: list[str] = []
    if reference is not None:
        fai = _read_fai(Path(reference))
        contig_order = [c for c, _ in fai]
        contig_lines = [f"##contig=<ID={c},length={n}>" for c, n in fai]
    else:
        contig_order = list(header_contigs)
        contig_lines = [hl.text for hl in merged_hdr.values() if hl.kind == "contig"]
        contig_lines.sort(
            key=lambda t: contig_order.index(re.search(r"ID=([^,>]+)", t).group(1))
        )
    for c, _p in by_pos:
        if c not in contig_order:
            if reference is not None:
                raise CombineError(
                    f"contig {c} is not in the reference index {reference}.fai."
                )
            contig_order.append(c)
    cidx = {c: i for i, c in enumerate(contig_order)}

    # AC/AF/AN and set lines added by CombineVariants.initialize.
    extra = [
        HeaderLine("INFO", "", "", i, *STANDARD_INFO[i][:2], STANDARD_INFO[i][2])
        for i in ("AC", "AF", "AN")
    ]
    for hl in extra:
        hl.text = hl.compound_text()
    set_line = HeaderLine(
        "INFO",
        "",
        "##INFO=<ID=set,Number=1,Type=String,Description="
        '"Source VCF for the merged record in CombineVariants">',
    )
    info_ids = {hl.id for hl in merged_hdr.values() if hl.kind == "INFO"} | {
        "AC",
        "AF",
        "AN",
        SET_KEY,
    }
    format_ids = {hl.id for hl in merged_hdr.values() if hl.kind == "FORMAT"}
    filter_ids = {hl.id for hl in merged_hdr.values() if hl.kind == "FILTER"}

    prio_idx = {n: i for i, n in enumerate(priority)} if priority else None
    out_lines: list[str] = []
    used_pl = False
    n_written = 0
    for key in sorted(by_pos, key=lambda k: (cidx[k[0]], k[1])):
        vcs = by_pos[key]
        n_filtered = sum(1 for v in vcs if v.filtered)
        if min_callsets > 1 and len(vcs) - n_filtered < min_callsets:
            continue
        groups = _separate_by_type(vcs)
        if "NO_VARIATION" in groups and len(groups) > 1:
            refs = groups.pop("NO_VARIATION")
            for t in TYPES:
                if t in groups:
                    groups[t].extend(refs)
                    break
        for t in TYPES:
            if t not in groups:
                continue
            members = groups[t]
            if prio_idx is not None:  # simpleMerge sorts after the type split
                members = sorted(members, key=lambda v: prio_idx[v.source])
            m = _simple_merge(
                members, n_inputs, merge_info_with_max_ac, filtered_are_uncalled
            )
            if m is None:
                continue
            _chromosome_counts(m, samples_out)
            line, pl = _render(m, samples_out, info_ids, format_ids, filter_ids)
            used_pl |= pl
            out_lines.append(line)
            n_written += 1

    header = _output_header(
        merged_hdr, extra, set_line, contig_lines, used_pl, format_ids
    )
    with open(output, "w") as fh:
        fh.write("\n".join(header) + "\n")
        fh.write(
            "\t".join(
                [
                    "#CHROM",
                    "POS",
                    "ID",
                    "REF",
                    "ALT",
                    "QUAL",
                    "FILTER",
                    "INFO",
                    "FORMAT",
                ]
                + samples_out
            )
            + "\n"
        )
        for line in out_lines:
            fh.write(line)
    return n_written


def _render(m: Merged, samples: list[str], info_ids, format_ids, filter_ids):
    where = f"{m.chrom}:{m.pos}"
    for k in m.info:
        if k not in info_ids:
            raise CombineError(
                f"Key {k} found in VariantContext field INFO at {where} but this key "
                "isn't defined in the VCF header of any input."
            )
    for f in m.filters or []:
        if f not in filter_ids:
            raise CombineError(
                f"Key {f} found in VariantContext field FILTER at {where} but this key "
                "isn't defined in the VCF header of any input."
            )
    gts = [m.genotypes[s] for s in samples if s in m.genotypes]
    keys = _format_keys(gts)
    used_pl = any(g.pl is not None for g in gts)
    for k in keys:
        if k != "GT" and k not in format_ids and not (k == "PL" and used_pl):
            raise CombineError(
                f"Key {k} found in VariantContext field FORMAT at {where} but this key "
                "isn't defined in the VCF header of any input."
            )
    bases_to_idx = {a[0]: i for i, a in enumerate(m.alleles)}
    cells = [
        _format_cell(m.genotypes.get(s), keys, bases_to_idx, where) for s in samples
    ]
    if m.filters is None:
        flt = "."
    elif not m.filters:
        flt = "PASS"
    else:
        flt = ";".join(m.filters)
    info_items = []
    for k in sorted(m.info):
        v = m.info[k]
        info_items.append(k if v is True or v == "" else f"{k}={v}")
    cols = [
        m.chrom,
        str(m.pos),
        m.id,
        m.alleles[0][0],
        ",".join(a[0] for a in m.alleles[1:]) or ".",
        "." if m.qual is None else format_qual(m.qual),
        flt,
        ";".join(info_items) or ".",
    ]
    if samples:
        cols += [":".join(keys), *cells]
    return "\t".join(cols) + "\n", used_pl


def _output_header(merged_hdr, extra, set_line, contig_lines, used_pl, format_ids):
    texts: set[str] = {"##fileformat=VCFv4.2"}
    for hl in merged_hdr.values():
        if hl.kind != "contig":
            texts.add(hl.text)
    for hl in extra:
        texts.add(hl.text)
    texts.add(set_line.text)
    if used_pl and "PL" not in format_ids:
        pl = STANDARD_FORMAT["PL"]
        texts.add(f'##FORMAT=<ID=PL,Number={pl[0]},Type={pl[1]},Description="{pl[2]}">')
    texts.add("##vcfclickCombine=gatk3-compatible")
    texts.discard("##fileformat=VCFv4.2")
    body = sorted(texts)
    # contig lines sort as a group, in reference order, at the text position
    # of "##contig".
    marker = "##contig="
    pos = next((i for i, t in enumerate(body) if t > marker), len(body))
    return ["##fileformat=VCFv4.2", *body[:pos], *contig_lines, *body[pos:]]
