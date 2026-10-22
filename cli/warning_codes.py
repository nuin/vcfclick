"""The closed list of warnings a filtered run can carry.

One shape everywhere (`{"code", "level", "key", "message"}`), so the text
header, the JSON `effective.warnings` and a desktop Checks row read the same
thing. A code is registered here, with its level, by the change that first
emits it; `warning()` refuses a code that is not, so a typo cannot invent a
vocabulary the apps do not know. A consumer must still treat an unknown code
as a plain string (a newer vcfclick may emit more).

`level` is `warning` (printed as a `WARNING:` line on stdout, because a run
that is quieter than the clinician believes is the failure these exist for)
or `info` (a note: the run did what was asked, and here is what that meant).
"""

from __future__ import annotations

WARNING = "warning"
INFO = "info"

CODES: dict[str, str] = {
    # a requested gate cannot be enforced: its field is missing on too many rows
    "gate_unbound": WARNING,
    # a gate is enforced, but not on every row
    "gate_partial": INFO,
    # the proband's median depth is close to min_dp: the preset may be too
    # strict for this assay
    "thin_coverage": WARNING,
    # neither an INFO/AF gate nor gnomAD narrows the candidates by frequency
    "no_population_frequency_filter": WARNING,
    # --gnomad-max-af was asked for and no gnomAD data is loaded
    "gnomad_not_loaded": WARNING,
    # gnomAD is loaded but has records on chromosomes that hold too little of
    # the cohort's variants (another naming, a slice elsewhere): it judges nothing
    "gnomad_no_coverage": WARNING,
    # after the run: the gnomAD filter found a record for none of the candidates
    "gnomad_no_record": WARNING,
    # a preset key this command does not take was skipped, and skipping it
    # makes the search wider
    "skipped_widening": WARNING,
    # the preset's max_af was not applied because INFO/AF is the cohort's own
    "preset_max_af_not_applied": INFO,
}


def warning(code: str, key: str | None, message: str) -> dict:
    """One warning object. `key` is the filter key it is about, or None."""
    if code not in CODES:
        raise ValueError(f"unregistered warning code {code!r}")
    return {"code": code, "level": CODES[code], "key": key, "message": message}
