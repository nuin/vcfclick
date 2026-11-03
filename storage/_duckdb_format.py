"""Render DuckDB result rows in chDB's output formats, byte for byte.

chDB's renderer is the reference: `db query --format JSON` is read by the
desktop apps, and every caller that parses JSONCompact, TSV or CSV must
see the same text on both backends. Each value is rendered from its
DuckDB result type, so a REAL prints as the shortest Float32 text and a
DOUBLE as the shortest Float64 text, both in ClickHouse's notation
(fixed point for exponents -6..20, otherwise `1e21` / `1.5e-8`).

What cannot match: type names carry no `Nullable(...)` or
`LowCardinality(...)` (DuckDB's result type has neither), and the JSON
`statistics` block holds `elapsed` only (rows and bytes read are chDB
counters), with no `rows_before_limit_at_least`.
"""

from __future__ import annotations

import datetime as _dt
import decimal
import math
import struct
import uuid
from json.encoder import encode_basestring as _c_json_string

# Format name (lower case, as chDB accepts any case) -> canonical name.
FORMAT_ALIASES = {
    "json": "JSON",
    "jsoncompact": "JSONCompact",
    "jsoneachrow": "JSONEachRow",
    "tsv": "TSV",
    "tabseparated": "TSV",
    "tsvwithnames": "TSVWithNames",
    "tabseparatedwithnames": "TSVWithNames",
    "csv": "CSV",
    "csvwithnames": "CSVWithNames",
}


# ---------------------------------------------------------------- types


class _Type:
    """A result column type: DuckDB id plus children for nested types."""

    __slots__ = ("id", "children", "name")

    def __init__(self, id_: str, children: list, name: str) -> None:
        self.id = id_
        self.children = children
        self.name = name


def column_type(t) -> _Type:
    """Wrap a DuckDBPyType from `cursor.description`."""
    tid = t.id
    if tid in ("list", "array"):
        return _Type("list", [column_type(t.children[0][1])], str(t))
    if tid == "map":
        return _Type(
            "map",
            [column_type(t.children[0][1]), column_type(t.children[1][1])],
            str(t),
        )
    if tid == "struct":
        return _Type("struct", [(n, column_type(c)) for n, c in t.children], str(t))
    if tid == "decimal":
        return _Type("decimal", [t.children[0][1], t.children[1][1]], str(t))
    return _Type(tid, [], str(t))


_CH_NAMES = {
    "utinyint": "UInt8",
    "usmallint": "UInt16",
    "uinteger": "UInt32",
    "ubigint": "UInt64",
    "uhugeint": "UInt128",
    "tinyint": "Int8",
    "smallint": "Int16",
    "integer": "Int32",
    "bigint": "Int64",
    "hugeint": "Int128",
    "float": "Float32",
    "double": "Float64",
    "varchar": "String",
    "blob": "String",
    "boolean": "Bool",
    "date": "Date",
    "timestamp": "DateTime",
    "uuid": "UUID",
}


def clickhouse_type_name(t: _Type) -> str:
    """The ClickHouse name of a DuckDB result type, without nullability."""
    if t.id == "list":
        return f"Array({clickhouse_type_name(t.children[0])})"
    if t.id == "map":
        k, v = t.children
        return f"Map({clickhouse_type_name(k)}, {clickhouse_type_name(v)})"
    if t.id == "struct":
        inner = ", ".join(f"{n} {clickhouse_type_name(c)}" for n, c in t.children)
        return f"Tuple({inner})"
    if t.id == "decimal":
        return f"Decimal({t.children[0]}, {t.children[1]})"
    return _CH_NAMES.get(t.id, t.name)


# -------------------------------------------------------------- scalars


def _shortest_float32(v: float) -> str:
    """Fewest significant digits that read back as the same Float32."""
    target = struct.unpack("f", struct.pack("f", v))[0]
    for p in range(1, 10):
        s = f"{target:.{p - 1}e}"
        if struct.unpack("f", struct.pack("f", float(s)))[0] == target:
            return s
    return repr(target)


def float_text(v: float, single: bool = False) -> str:
    """ClickHouse's text for a float: shortest round-trip digits, fixed
    point when the decimal exponent is in -6..20, else `d.ddde<exp>`."""
    if math.isnan(v):
        return "nan"
    if math.isinf(v):
        return "inf" if v > 0 else "-inf"
    sign = "-" if math.copysign(1.0, v) < 0 else ""
    if v == 0:
        return sign + "0"
    shortest = _shortest_float32(v) if single else repr(abs(v))
    tup = decimal.Decimal(shortest).copy_abs().normalize().as_tuple()
    digits = "".join(map(str, tup.digits))
    point = len(digits) + tup.exponent  # value = 0.digits * 10**point
    exponent = point - 1
    if -6 <= exponent < 21:
        if point <= 0:
            body = "0." + "0" * (-point) + digits
        elif point >= len(digits):
            body = digits + "0" * (point - len(digits))
        else:
            body = digits[:point] + "." + digits[point:]
    else:
        body = digits[0] + ("." + digits[1:] if len(digits) > 1 else "")
        body += f"e{exponent}"
    return sign + body


def _double_text(v: float) -> str:
    """`float_text` of a Float64. Between 1e-4 and 1e16 Python's repr is
    already fixed point with the shortest digits, so only a trailing `.0`
    goes; elsewhere the slow path decides the notation."""
    if 1e-4 <= abs(v) < 1e16:
        r = repr(v)
        return r[:-2] if r.endswith(".0") else r
    return float_text(v)


def _float32_text(v: float) -> str:
    """`float_text(v, single=True)`. Between 1e-4 and 1e16 numpy prints a
    float32 fixed point with its shortest round-trip digits."""
    if 1e-4 <= abs(v) < 1e16:
        import numpy as np

        return np.format_float_positional(np.float32(v), unique=True, trim="-")
    return float_text(v, single=True)


def _bool_text(v: bool) -> str:
    return "true" if v else "false"


_INT_IDS = frozenset(
    (
        "utinyint",
        "usmallint",
        "uinteger",
        "ubigint",
        "uhugeint",
        "tinyint",
        "smallint",
        "integer",
        "bigint",
        "hugeint",
    )
)


def _scalar_number_fn(t: _Type):
    """Text of a non-null value of a number or bool column, chosen once per
    column; None for any other column type."""
    if t.id in _INT_IDS:
        return str
    if t.id == "boolean":
        return _bool_text
    if t.id == "double":
        return _double_text
    if t.id == "float":
        return _float32_text
    return None


def _nullable(fn, null: str):
    def text(v):
        return null if v is None else fn(v)

    return text


def _decimal_text(v: decimal.Decimal) -> str:
    """Decimal without trailing zeros, as ClickHouse prints it."""
    s = format(v, "f")
    if "." in s:
        s = s.rstrip("0").rstrip(".")
    return "0" if s in ("-0", "") else s


def _datetime_text(v) -> str:
    if isinstance(v, _dt.datetime):
        return v.strftime("%Y-%m-%d %H:%M:%S")
    if isinstance(v, _dt.date):
        return v.isoformat()
    return str(v)


def _as_str(v) -> str:
    if isinstance(v, bytes | bytearray | memoryview):
        return bytes(v).decode("utf-8", errors="replace")
    if isinstance(v, uuid.UUID):
        return str(v)
    return str(v)


def _number_text(v, t: _Type) -> str | None:
    """Text of a numeric/bool value, or None when `v` is not a number."""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        return float_text(v, single=t.id == "float")
    if isinstance(v, decimal.Decimal):
        return _decimal_text(v)
    return None


def _is_quoted_scalar(v) -> bool:
    return not isinstance(v, bool | int | float | decimal.Decimal)


def _map_items(v) -> list:
    if isinstance(v, dict):
        if set(v) == {"key", "value"} and isinstance(v["key"], list):
            return list(zip(v["key"], v["value"], strict=False))
        return list(v.items())
    return list(v)


# ----------------------------------------------------------------- JSON


def json_string(s: str) -> str:
    """A JSON string as chDB writes it: the C encoder's escapes (control
    characters as \\b \\f \\n \\r \\t or \\u00xx, quote, backslash), plus
    `/`, U+2028 and U+2029, which chDB also escapes."""
    out = _c_json_string(s)
    if "/" in out:
        out = out.replace("/", "\\/")
    if "\u2028" in out:
        out = out.replace("\u2028", "\\u2028")
    if "\u2029" in out:
        out = out.replace("\u2029", "\\u2029")
    return out


def json_value(v, t: _Type, indent: int | None = None) -> str:
    """One JSON value. `indent` (a tab depth) lays maps and structs out one
    entry per line, as chDB's JSON format does; None is compact."""
    if v is None:
        return "null"
    if t.id == "list":
        return "[" + ",".join(json_value(x, t.children[0]) for x in v) + "]"
    if t.id in ("map", "struct"):
        if t.id == "map":
            kt, vt = t.children
            items = [(_key_text(k, kt), x, vt) for k, x in _map_items(v)]
        else:
            items = [(n, v.get(n), ct) for n, ct in t.children]
        if not items:
            return "{}"
        if indent is None:
            return (
                "{"
                + ",".join(
                    f"{json_string(k)}:{json_value(x, ct)}" for k, x, ct in items
                )
                + "}"
            )
        pad = "\t" * (indent + 1)
        body = ",\n".join(
            f"{pad}{json_string(k)}: {json_value(x, ct, indent + 1)}"
            for k, x, ct in items
        )
        return "{\n" + body + "\n" + "\t" * indent + "}"
    if isinstance(v, float) and not math.isfinite(v):
        return "null"
    num = _number_text(v, t)
    if num is not None:
        return num
    return json_string(_datetime_text(v) if _is_temporal(v) else _as_str(v))


def _is_temporal(v) -> bool:
    return isinstance(v, _dt.date | _dt.datetime | _dt.time)


def _key_text(k, t: _Type) -> str:
    """A map key as the string a JSON object key holds."""
    num = _number_text(k, t)
    if num is not None:
        return num
    return _datetime_text(k) if _is_temporal(k) else _as_str(k)


def _json_meta(columns: list[str], types: list[_Type]) -> str:
    entries = [
        "\t\t{\n"
        f'\t\t\t"name": {json_string(c)},\n'
        f'\t\t\t"type": {json_string(clickhouse_type_name(t))}\n'
        "\t\t}"
        for c, t in zip(columns, types, strict=True)
    ]
    return '\t"meta":\n\t[\n' + ",\n".join(entries) + "\n\t],\n\n"


def _json_tail(n_rows: int, elapsed: float) -> str:
    return (
        f'\t"rows": {n_rows},\n\n'
        '\t"statistics":\n\t{\n'
        f'\t\t"elapsed": {float_text(elapsed)}\n'
        "\t}\n}\n"
    )


# ------------------------------------------------------------ TSV / CSV

_TSV_ESCAPES = {
    "\\": "\\\\",
    "'": "\\'",
    "\b": "\\b",
    "\f": "\\f",
    "\n": "\\n",
    "\r": "\\r",
    "\t": "\\t",
    "\0": "\\0",
}


_TSV_TABLE = str.maketrans(_TSV_ESCAPES)


def _escape(s: str) -> str:
    return s.translate(_TSV_TABLE)


def _quoted_text(v, t: _Type) -> str:
    """A value inside an Array / Map / Tuple text form: strings in single
    quotes with backslash escapes, NULL as NULL."""
    if v is None:
        return "NULL"
    if t.id == "list":
        return "[" + ",".join(_quoted_text(x, t.children[0]) for x in v) + "]"
    if t.id == "map":
        kt, vt = t.children
        return (
            "{"
            + ",".join(
                f"{_quoted_text(k, kt)}:{_quoted_text(x, vt)}" for k, x in _map_items(v)
            )
            + "}"
        )
    if t.id == "struct":
        return "(" + ",".join(_quoted_text(v.get(n), ct) for n, ct in t.children) + ")"
    num = _number_text(v, t)
    if num is not None:
        return num
    return "'" + _escape(_datetime_text(v) if _is_temporal(v) else _as_str(v)) + "'"


def _tsv_field(v, t: _Type) -> str:
    if v is None:
        return "\\N"
    if t.id in ("list", "map", "struct"):
        return _quoted_text(v, t)
    num = _number_text(v, t)
    if num is not None:
        return num
    return _escape(_datetime_text(v) if _is_temporal(v) else _as_str(v))


def _csv_quote(s: str) -> str:
    return '"' + s.replace('"', '""') + '"'


def _csv_fields(v, t: _Type) -> list[str]:
    """One value as CSV fields; a Tuple spreads over one field per element."""
    if v is None:
        return ["\\N"]
    if t.id == "struct":
        out = []
        for n, ct in t.children:
            out.extend(_csv_fields(v.get(n), ct))
        return out
    if t.id in ("list", "map"):
        return [_csv_quote(_quoted_text(v, t))]
    num = _number_text(v, t)
    if num is not None:
        return [num]
    return [_csv_quote(_datetime_text(v) if _is_temporal(v) else _as_str(v))]


# --------------------------------------------------- per-column formatters
#
# Each formatter is chosen once per column from its result type. Integer,
# bool, float and string columns get a direct function; every other type
# goes through the general per-value path above.


def _json_fn(t: _Type, indent: int | None = None):
    num = _scalar_number_fn(t)
    if num is not None:
        if t.id in ("double", "float"):

            def finite(v):
                return num(v) if math.isfinite(v) else "null"

            return _nullable(finite, "null")
        return _nullable(num, "null")
    if t.id == "varchar":
        return _nullable(json_string, "null")
    return lambda v: json_value(v, t, indent)


def _tsv_fn(t: _Type):
    num = _scalar_number_fn(t)
    if num is not None:
        return _nullable(num, "\\N")
    if t.id == "varchar":
        return _nullable(_escape, "\\N")
    return lambda v: _tsv_field(v, t)


def _csv_fn(t: _Type):
    num = _scalar_number_fn(t)
    if num is not None:
        return _nullable(num, "\\N")
    if t.id == "varchar":
        return _nullable(_csv_quote, "\\N")
    return lambda v: ",".join(_csv_fields(v, t))


def render_json(columns, types, rows, elapsed: float) -> str:
    """chDB `JSON`: meta with types, data as one object per row."""
    keys = [f"\t\t\t{json_string(c)}: " for c in columns]
    fns = [_json_fn(t, 3) for t in types]
    data = [
        "\t\t{\n"
        + ",\n".join([k + f(v) for k, f, v in zip(keys, fns, r, strict=True)])
        + "\n\t\t}"
        for r in rows
    ]
    return (
        "{\n"
        + _json_meta(columns, types)
        + '\t"data":\n\t[\n'
        + ",\n".join(data)
        + "\n\t],\n\n"
        + _json_tail(len(rows), elapsed)
    )


def render_jsoncompact(columns, types, rows, elapsed: float) -> str:
    """chDB `JSONCompact`: meta with types, data as one array per row."""
    fns = [_json_fn(t) for t in types]
    data = [
        "\t\t[" + ", ".join([f(v) for f, v in zip(fns, r, strict=True)]) + "]"
        for r in rows
    ]
    return (
        "{\n"
        + _json_meta(columns, types)
        + '\t"data":\n\t[\n'
        + ",\n".join(data)
        + "\n\t],\n\n"
        + _json_tail(len(rows), elapsed)
    )


def render_jsoneachrow(columns, types, rows) -> str:
    keys = [f"{json_string(c)}:" for c in columns]
    fns = [_json_fn(t) for t in types]
    return "".join(
        [
            "{"
            + ",".join([k + f(v) for k, f, v in zip(keys, fns, r, strict=True)])
            + "}\n"
            for r in rows
        ]
    )


def render_tsv(columns, types, rows, names: bool = False) -> str:
    fns = [_tsv_fn(t) for t in types]
    lines = ["\t".join(_escape(c) for c in columns)] if names else []
    lines += ["\t".join([f(v) for f, v in zip(fns, r, strict=True)]) for r in rows]
    return "".join(line + "\n" for line in lines)


def render_csv(columns, types, rows, names: bool = False) -> str:
    fns = [_csv_fn(t) for t in types]
    lines = [",".join(_csv_quote(c) for c in columns)] if names else []
    lines += [",".join([f(v) for f, v in zip(fns, r, strict=True)]) for r in rows]
    return "".join(line + "\n" for line in lines)


def render(fmt: str, columns, duck_types, rows, elapsed: float) -> bytes:
    """Render rows in canonical format `fmt` (a FORMAT_ALIASES value)."""
    types = [column_type(t) for t in duck_types]
    if fmt == "JSON":
        text = render_json(columns, types, rows, elapsed)
    elif fmt == "JSONCompact":
        text = render_jsoncompact(columns, types, rows, elapsed)
    elif fmt == "JSONEachRow":
        text = render_jsoneachrow(columns, types, rows)
    elif fmt in ("TSV", "TSVWithNames"):
        text = render_tsv(columns, types, rows, names=fmt == "TSVWithNames")
    elif fmt in ("CSV", "CSVWithNames"):
        text = render_csv(columns, types, rows, names=fmt == "CSVWithNames")
    else:  # pragma: no cover - callers pass FORMAT_ALIASES values only
        raise ValueError(f"Unknown format {fmt}")
    return text.encode()
