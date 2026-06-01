# Licensing

vcfclick ships under a **single open-source license** (AGPL-3.0 or BSL —
to be selected before first public release). All code, all features, no
tiers. Academic, commercial, and clinical use are all permitted by the
license itself.

This is a deliberate choice. The project is built as a portfolio /
distribution piece for [Bioinformat](#) — its job is
to be technically credible and widely used, not to maximise per-seat
revenue. Revenue, when it comes, comes from the operational side
described below — not from gating features.

## How Bioinformat makes money around vcfclick

The OSS itself stays free. Three revenue surfaces sit *next to* it:

1. **Hosted SaaS** — a cheap managed instance (Plausible-style pricing).
   Covers infrastructure plus a modest margin. The point is convenience
   for users who don't want to operate ClickHouse + DuckDB + monthly
   ClinVar refreshes themselves, not feature parity gating.

2. **Self-hosted support contracts** — install, customise, and support
   contracts for research labs and bioinformatics companies that want
   to run vcfclick on their own infrastructure. Annual flat fee.

3. **Consulting engagements** — Bioinformat's existing line of work.
   vcfclick as a portfolio piece draws inbound conversations; some of
   those become engagements that have nothing to do with vcfclick.

None of these require the OSS to be feature-limited. The hosted SaaS
sells operational ease; the support contracts sell time and expertise;
consulting sells everything else.

## Status

The open-source license itself (AGPL-3.0 vs BSL) has not been chosen.
This is the next strategic decision before the first public commit.

- **AGPL-3.0** — strong copyleft; network-deployment changes must be
  shared back. Protects against a competitor wrapping vcfclick as a
  closed SaaS. Some commercial users avoid AGPL on principle.
- **BSL** — time-delayed conversion to a permissive license (e.g.,
  Apache 2 after 3 years). Lets Bioinformat restrict commercial
  competition during the active window without scaring research users
  for normal use. Pioneered by MariaDB, used by Sentry, CockroachDB.

For a portfolio-first project with a Bioinformat-operated hosted SaaS as
the eventual revenue lane, BSL is the more aligned choice. Decide
before first public commit.
