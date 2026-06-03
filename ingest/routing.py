"""VCF field → table column routing.

Single source of truth for which INFO/FORMAT fields land in typed
columns versus the info_extra / format_extra Map overflow. Extend
the routing tables below when promoting a field from overflow to typed.

Pairs with ingest.vcf_load (which consumes these tables in the
row builders) and cli.discover (which uses them to suggest new
typed columns from VCF headers it has seen in the wild).
"""

from __future__ import annotations

from cyvcf2 import VCF

INFO_SCALAR = {
    "AC": "info_AC",
    "AF": "info_AF",
    "AN": "info_AN",
    "DP": "info_DP",
    "MQ": "info_MQ",
    "MQ0": "info_MQ0",
    "NS": "info_NS",
    "BQ": "info_BQ",
    "SB": "info_SB",
    "END": "info_END",
    "CIGAR": "info_CIGAR",
    "AA": "info_AA",
    "QD": "info_QD",
    "FS": "info_FS",
    "SOR": "info_SOR",
    "MQRankSum": "info_MQRankSum",
    "ReadPosRankSum": "info_ReadPosRankSum",
    "ExcessHet": "info_ExcessHet",
    "InbreedingCoeff": "info_InbreedingCoeff",
    "MLEAC": "info_MLEAC",
    "MLEAF": "info_MLEAF",
    "BaseQRankSum": "info_BaseQRankSum",
    "ClippingRankSum": "info_ClippingRankSum",
}

INFO_PAIR = {
    "AD": ("info_AD_ref", "info_AD_alt"),
}

INFO_FLAG = {
    "SOMATIC": "info_SOMATIC",
    "VALIDATED": "info_VALIDATED",
    "DB": "info_DB",
    "H2": "info_H2",
    "H3": "info_H3",
    "1000G": "info_1000G",
}

FORMAT_SCALAR = {
    "GQ": "gq",
    "DP": "dp",
    "MQ": "mq",
    "FT": "ft",
    "PS": "ps",
    "PQ": "pq",
}

FORMAT_PAIR = {
    "AD": ("ad_ref", "ad_alt"),
    "ADF": ("adf_ref", "adf_alt"),
    "ADR": ("adr_ref", "adr_alt"),
}

FORMAT_TRIPLE = {
    "PL": ("pl_ref_ref", "pl_ref_alt", "pl_alt_alt"),
    "GL": ("gl_ref_ref", "gl_ref_alt", "gl_alt_alt"),
}


def _typed_format_fields() -> set[str]:
    return {"GT"} | FORMAT_SCALAR.keys() | FORMAT_PAIR.keys() | FORMAT_TRIPLE.keys()


def classify_header(vcf: VCF) -> dict[str, list[str]]:
    typed_info, extra_info = [], []
    typed_format, extra_format = [], []
    typed_format_ids = _typed_format_fields()

    for h in vcf.header_iter():
        try:
            d = h.info(extra=True)
        except Exception:
            continue
        kind = str(d.get("HeaderType", "")).upper()
        field = d.get("ID")
        if not field:
            continue
        if kind == "INFO":
            if field in INFO_SCALAR or field in INFO_PAIR or field in INFO_FLAG:
                typed_info.append(field)
            else:
                extra_info.append(field)
        elif kind == "FORMAT":
            if field in typed_format_ids:
                typed_format.append(field)
            else:
                extra_format.append(field)

    return {
        "typed_info": sorted(typed_info),
        "extra_info": sorted(extra_info),
        "typed_format": sorted(typed_format),
        "extra_format": sorted(extra_format),
    }
