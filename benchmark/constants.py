"""GA4GH-style annotation vocabularies for `vcfclick benchmark`.

Casing follows hap.py convention. These strings are a P1 approximation and must
be pinned against a real hap.py intermediate-VCF fixture before any hap.py-parity
claim (see design doc, open questions).
"""

from __future__ import annotations

# BD — benchmarking decision.
BD_TP = "TP"
BD_FP = "FP"
BD_FN = "FN"
BD_N = "N"  # outside confident region / unknown

# BK — match kind.
BK_GM = "gm"  # genotype match
BK_AM = "am"  # allele match, genotype mismatch
BK_LM = "lm"  # local haplotype match (P2 only)
BK_NONE = "."

# BVT — variant type.
VT_SNP = "SNP"
VT_INDEL = "INDEL"
VT_COMPLEX = "COMPLEX"
VT_NOCALL = "NOCALL"

# BLT — biallelic locus (zygosity) type.
BLT_HET = "het"
BLT_HOMALT = "homalt"
BLT_HETALT = "hetalt"
BLT_HOMREF = "homref"
BLT_HAPLOID = "haploid"
BLT_NOCALL = "nocall"

# Filter views. PASS/ALL are classified independently (a query call failing FILTER
# is *absent* in the PASS view, so its matched truth call becomes FN there).
FILTER_ALL = "ALL"
FILTER_PASS = "PASS"
