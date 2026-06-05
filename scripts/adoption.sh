#!/usr/bin/env bash
# Adoption snapshot for vcfclick.
#
# Pulls current numbers from PyPI (pypistats.org), GitHub (gh api), and
# the GitHub release-asset endpoint. No auth needed for pypistats; `gh`
# reuses your existing `gh auth` login for the GitHub side.
#
# Run as needed:
#     scripts/adoption.sh
#
# Or piped to a dated log:
#     mkdir -p logs && scripts/adoption.sh > "logs/adoption-$(date +%Y-%m-%d).txt"

set -euo pipefail

OWNER=nuin
REPO=vcfclick
PYPI_PKG=vcfclick

for tool in curl python3 gh; do
    command -v "$tool" >/dev/null 2>&1 || { echo "missing: $tool" >&2; exit 1; }
done

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

echo "=== vcfclick adoption snapshot ==="
echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo

# ─── PyPI ────────────────────────────────────────────────────────────
#
# Approach: fetch each API response into a temp file, then heredoc-feed
# the Python script and read the file by name. Lets the Python script
# use f-strings + double-quoted dict keys without any bash-quote escape
# nightmares.

curl -s "https://pypistats.org/api/packages/${PYPI_PKG}/recent" > "$WORK/pypi_recent.json"
curl -s "https://api.pepy.tech/api/v2/projects/${PYPI_PKG}"     > "$WORK/pepy.json"

echo "PyPI ($PYPI_PKG):"
python3 - "$WORK" <<'PY'
import json, sys, pathlib
work = pathlib.Path(sys.argv[1])

try:
    d = json.loads((work / "pypi_recent.json").read_text())["data"]
    print(f"  last 24h:   {d.get('last_day', 0):>6,}")
    print(f"  last week:  {d.get('last_week', 0):>6,}")
    print(f"  last month: {d.get('last_month', 0):>6,}")
except Exception as e:
    print(f"  pypistats.org error: {e}")

try:
    d = json.loads((work / "pepy.json").read_text())
    total = d.get("total_downloads")
    if total:
        print(f"  total (all-time): {total:,}")
except Exception:
    pass  # pepy may not have indexed this package yet
PY
echo

# ─── GitHub repo ────────────────────────────────────────────────────

gh api "repos/${OWNER}/${REPO}" > "$WORK/repo.json"
gh api "repos/${OWNER}/${REPO}/releases" > "$WORK/releases.json"

echo "GitHub ($OWNER/$REPO):"
python3 - "$WORK" <<'PY'
import json, sys, pathlib
work = pathlib.Path(sys.argv[1])
d = json.loads((work / "repo.json").read_text())
print(f"  stars:        {d['stargazers_count']:>6,}")
print(f"  forks:        {d['forks_count']:>6,}")
print(f"  watchers:     {d['subscribers_count']:>6,}")
print(f"  open issues:  {d['open_issues_count']:>6,}")
PY
echo

# ─── Release assets ─────────────────────────────────────────────────

echo "Release assets:"
python3 - "$WORK" <<'PY'
import json, sys, pathlib
work = pathlib.Path(sys.argv[1])
releases = json.loads((work / "releases.json").read_text())
for rel in releases:
    tag = rel.get("tag_name", "?")
    assets = rel.get("assets", [])
    if not assets:
        continue
    print(f"  {tag}:")
    for a in assets:
        print(f"    {a['name']:<40s} {a['download_count']:>6,} downloads")
PY
echo

# ─── Traffic (admin-only) ───────────────────────────────────────────

echo "Traffic (last 14 days; visible only to repo admins):"
if gh api "repos/${OWNER}/${REPO}/traffic/views"  > "$WORK/views.json"  2>/dev/null &&
   gh api "repos/${OWNER}/${REPO}/traffic/clones" > "$WORK/clones.json" 2>/dev/null; then
    python3 - "$WORK" <<'PY'
import json, sys, pathlib
work = pathlib.Path(sys.argv[1])
v = json.loads((work / "views.json").read_text())
c = json.loads((work / "clones.json").read_text())
print(f"  views:         {v.get('count', 0):>6,}")
print(f"  unique views:  {v.get('uniques', 0):>6,}")
print(f"  clones:        {c.get('count', 0):>6,}")
print(f"  unique clones: {c.get('uniques', 0):>6,}")
PY
else
    echo "  (skipped — no push access to ${OWNER}/${REPO})"
fi
