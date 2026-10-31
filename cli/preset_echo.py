"""What a preset run says about itself: the `preset` object, the text header.

Pure functions, so the same echo can be built for `db trio` today and for the
other commands that take a preset later, and so the web and MCP resolvers can
produce the identical `preset` object from the same inputs.

The JSON policy (docs/TRIO.md) is additive only: `preset` and `effective` are
top-level keys present only when a preset (or a strictness flag) is in use, so
a plain run's document is byte-identical to what it was. The text is the same:
a `preset:` and a `gates:` line and any `WARNING:` lines follow the `trio:`
line, and only for a run with a preset.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from cli import filter_registry as fr
from cli.preset_schema import EXTENSION_PREFIX
from cli.presets import ApplyPlan, Located

# gates always listed in the text `gates:` line, whatever their source
_ALWAYS = ("min_gq", "min_dp", "min_ab", "max_ab", "max_af")


def compact(key: str, value: Any) -> Any:
    """A value as the echo shows it: a gene or region list above
    cli.scope.ECHO_LIST_MAX entries as {count, sha256_12}, so a 20,000-gene
    panel is not repeated in every result document."""
    spec = fr.REGISTRY.get(key)
    if spec is None or not spec.is_list or not isinstance(value, (list, tuple)):
        return value
    from cli.scope import compact_list

    return compact_list(list(value), upper=spec.type == "genes")


def _compact_overridden(overridden: Sequence[dict]) -> list[dict]:
    return [
        {
            **o,
            "preset": compact(o["key"], o["preset"]),
            "flag": compact(o["key"], o["flag"]),
        }
        for o in overridden
    ]


def skipped_items(
    plan: ApplyPlan, command: str, *, skip_unavailable: bool
) -> list[dict]:
    """The preset keys that were not applied: a known key the command does not
    take (`not_applicable`) and, when `--preset-skip-unavailable` was given, a
    key this vcfclick does not know (`unavailable`)."""
    out = []
    for key in plan.not_applicable:
        out.append(
            {
                "key": key,
                "reason": "not_applicable",
                "widening": fr.REGISTRY[key].widening,
                "message": f"{key} is not a filter of `{command}`; the preset's "
                "value was not applied",
            }
        )
    if skip_unavailable:
        for key in plan.unavailable:
            out.append(
                {
                    "key": key,
                    "reason": "unavailable",
                    "widening": True,
                    "message": f"{key[:40]} is unknown to this vcfclick; the "
                    "preset's value was not applied",
                }
            )
    return out


def resolved_block(
    resolved: Mapping[str, fr.Resolved], bound: Mapping[str, bool] | None
) -> dict:
    """{key: {"value", "from", "bound"?}} in registry order. `bound` appears
    only for a gate the effective-gate report judged."""
    out = {}
    for key, r in resolved.items():
        entry: dict[str, Any] = {"value": compact(key, r.value), "from": r.source}
        if bound is not None and key in bound:
            entry["bound"] = bound[key]
        out[key] = entry
    return out


def preset_block(
    located: Located,
    plan: ApplyPlan,
    resolved: Mapping[str, fr.Resolved],
    overridden: Sequence[dict],
    *,
    command: str,
    skip_unavailable: bool = False,
    bound: Mapping[str, bool] | None = None,
    notes: Sequence[str] = (),
    models_not_run: Sequence[str] = (),
) -> dict:
    """The JSON `preset` object. `requested` is the file's filters,
    `resolved` the value each filter took and where it came from, `overridden`
    the preset keys a flag or --preset-ignore changed, `skipped` the keys not
    applied, `models_not_run` the models the preset's category leaves out
    (their counts need --funnel)."""
    p = located.preset
    meta = p.meta
    return {
        "name": p.name,
        "source": p.source,
        "revision": p.revision,
        "status": p.status,
        "schema_version": p.doc.get("schema_version"),
        "fingerprint": p.fingerprint,
        "pinned": located.pinned,
        "meta": {
            "approved_by": meta.get("approved_by"),
            "approved_on": meta.get("approved_on"),
            "locked": bool(meta.get("locked", False)),
        },
        "path": p.path.name if p.path is not None and p.source != "builtin" else None,
        "requested": {
            k: compact(k, v)
            for k, v in p.filters.items()
            if not k.startswith(EXTENSION_PREFIX)
        },
        "resolved": resolved_block(resolved, bound),
        "overridden": _compact_overridden(overridden),
        "skipped": skipped_items(plan, command, skip_unavailable=skip_unavailable),
        "models_not_run": [{"model": m, "count": None} for m in models_not_run],
        "covers": p.builtin.covers if p.builtin else None,
        "does_not_cover": p.builtin.does_not_cover if p.builtin else None,
        "notes": list(notes),
    }


# ───────────────────────────────── text ─────────────────────────────────


def _show(value: Any) -> str:
    if value is None:
        return "none"
    if isinstance(value, bool):
        return "on" if value else "off"
    if isinstance(value, (list, tuple)):
        return ",".join(str(v) for v in value)
    if isinstance(value, float):
        return f"{value:g}"
    return str(value)


def preset_line(block: Mapping[str, Any]) -> str:
    """`preset: rare-recessive (builtin, rev 1, a1b2c3d4e5f6, starting point)`."""
    bits = [block["source"]]
    if block["revision"] is not None:
        bits.append(f"rev {block['revision']}")
    bits.append(str(block["fingerprint"]))
    if block["status"] == "starting_point":
        bits.append("starting point")
    if block["pinned"]:
        bits.append("pinned")
    return f"preset: {block['name']} ({', '.join(bits)})"


def gates_line(
    resolved: Mapping[str, fr.Resolved], overridden: Sequence[dict], af_gate: str
) -> str:
    """`  gates: min_gq=20 [preset]  max_af=0.005 [flag; preset 0.01] ...`: every
    value that matters and where it came from."""
    over = {o["key"]: o for o in overridden}
    screen = resolved.get("parent_alt_screen")  # trio only
    parts = []
    for key, r in resolved.items():
        if fr.REGISTRY[key].is_list:
            continue  # genes and regions have their own `scope:` line
        # a derived default (min_male_x_dp: half of min_dp) is a default too
        default = r.source in (fr.FROM_DEFAULT, fr.FROM_DERIVED)
        listed = key in _ALWAYS or not default or key in over
        if key == "max_parent_alt_frac" and screen is not None:
            listed = bool(screen.value)  # only meaningful with the screen
        if not listed:
            continue
        source = r.source
        if key in over:
            was = _show(over[key]["preset"])
            source = f"{'off' if over[key]['flag'] == fr.OFF else 'flag'}; preset {was}"
        if key == "max_af" and af_gate == "skipped_cohort_derived":
            source += "; skipped: INFO/AF is cohort-derived"
        parts.append(f"{key}={_show(r.value)} [{source}]")
    return "  gates: " + "  ".join(parts)


def header_lines(
    block: Mapping[str, Any],
    resolved: Mapping[str, fr.Resolved],
    af_gate: str,
    warnings: Sequence[Mapping[str, Any]],
) -> list[str]:
    """The lines that follow `trio:` for a run with a preset: the preset, its
    gates with their sources, and one WARNING line per warning (stdout, so a
    pipe that drops stderr still shows them)."""
    lines = [
        preset_line(block),
        gates_line(resolved, block["overridden"], af_gate),
    ]
    lines += [f"  WARNING: {w['message']}" for w in warnings if w["level"] == "warning"]
    return lines
