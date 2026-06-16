"""Server-side natural-language → SQL for the web UI.

The browser sends the question plus a bring-your-own API key (never
stored server-side); this turns it into a single read-only SQL query
using the same SCHEMA_DESCRIPTION briefing the MCP server hands an LLM.
Uses stdlib urllib so the `[web]` extra needs no HTTP client dependency.
"""

from __future__ import annotations

import json
import urllib.error
import urllib.request

DEFAULT_MODEL = {
    "gemini": "gemini-flash-latest",
    "anthropic": "claude-haiku-4-5-20251001",
}

_INSTRUCTION = (
    "You translate the user's question into ONE read-only SQL query for the "
    "schema described above, to run on an embedded ClickHouse (chDB) / DuckDB "
    "backend. Output ONLY the SQL — no prose, no explanation, no markdown "
    "fences. Prefer a LIMIT on row-returning queries. Never write INSERT, "
    "UPDATE, DELETE, DROP, ALTER, or CREATE."
)


class LLMError(RuntimeError):
    """Raised when the provider call fails or returns no usable SQL."""


def _post(url: str, payload: dict, headers: dict, timeout: float = 45.0) -> dict:
    body = json.dumps(payload).encode()
    req = urllib.request.Request(url, data=body, headers=headers, method="POST")
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return json.loads(resp.read().decode())
    except urllib.error.HTTPError as e:
        detail = e.read().decode(errors="replace")[:400]
        raise LLMError(f"{e.code} {e.reason}: {detail}") from e
    except urllib.error.URLError as e:
        raise LLMError(f"could not reach the LLM provider: {e.reason}") from e


def _strip_sql(text: str) -> str:
    """Pull the bare SQL out of a model reply that may wrap it in fences."""
    t = text.strip()
    if t.startswith("```"):
        t = t.split("\n", 1)[1] if "\n" in t else t[3:]
        if t.endswith("```"):
            t = t[:-3]
        if t.lstrip().lower().startswith("sql"):
            t = t.lstrip()[3:]
    return t.strip().rstrip(";").strip()


def _gemini(key: str, model: str, system: str, question: str) -> str:
    url = (
        f"https://generativelanguage.googleapis.com/v1beta/models/"
        f"{model}:generateContent?key={key}"
    )
    payload = {
        "systemInstruction": {"parts": [{"text": system}]},
        "contents": [{"role": "user", "parts": [{"text": question}]}],
        "generationConfig": {
            "temperature": 0,
            "maxOutputTokens": 2048,
            "thinkingConfig": {"thinkingBudget": 0},
        },
    }
    data = _post(url, payload, {"content-type": "application/json"})
    candidates = data.get("candidates") or []
    if not candidates:
        raise LLMError("Gemini returned no candidates (check the API key / model).")
    content = candidates[0].get("content") or {}
    parts = content.get("parts") or []
    text = "".join(p.get("text", "") for p in parts)
    if not text.strip():
        raise LLMError("Gemini returned an empty response.")
    return _strip_sql(text)


def _anthropic(key: str, model: str, system: str, question: str) -> str:
    payload = {
        "model": model,
        "max_tokens": 1024,
        "temperature": 0,
        "system": system,
        "messages": [{"role": "user", "content": question}],
    }
    headers = {
        "content-type": "application/json",
        "x-api-key": key,
        "anthropic-version": "2023-06-01",
    }
    data = _post("https://api.anthropic.com/v1/messages", payload, headers)
    blocks = data.get("content") or []
    text = "".join(b.get("text", "") for b in blocks if b.get("type") == "text")
    if not text.strip():
        raise LLMError("Anthropic returned an empty response.")
    return _strip_sql(text)


def generate_sql(
    provider: str, key: str, model: str, question: str, schema: str
) -> str:
    """Turn `question` into a SQL string via the chosen provider.

    `schema` is the SCHEMA_DESCRIPTION briefing. Raises LLMError on any
    provider/auth/parse failure so the route can surface a clean message.
    """
    if not key:
        raise LLMError("No API key provided.")
    provider = provider.lower()
    model = model or DEFAULT_MODEL.get(provider, "")
    system = f"{schema}\n\n{_INSTRUCTION}"
    if provider == "gemini":
        return _gemini(key, model, system, question)
    if provider == "anthropic":
        return _anthropic(key, model, system, question)
    raise LLMError(f"unknown provider: {provider!r} (expected 'gemini' or 'anthropic')")
