# License

vcfclick is licensed under the **Apache License, Version 2.0**. The full
licence text is in [`LICENSE`](LICENSE) at the repository root.

## Why Apache 2.0

Three specific reasons, in order of weight for this project:

1. **It matches the bioinformatics convention.** GATK, Hail, most
   Broad / Sanger tools, and Polars are Apache 2.0 or MIT. Researchers
   see "Apache 2.0" and have no institutional friction; some
   commercial labs have policies that refuse AGPL outright. For a
   research-tooling project, adoption friction is the bigger cost
   than the protection a restrictive licence would buy.

2. **The patent grant matters here.** Apache 2.0 includes an explicit
   patent licence that MIT does not. For a tool that may eventually
   sit alongside closed-source variant callers (DRAGEN and similar)
   and that downstream users may want to deploy in clinical-adjacent
   contexts, the explicit patent grant is what legal reviewers look
   for.

3. **It does not preclude future relicensing.** If vcfclick ever
   becomes attractive enough to be worth cloning as a competing
   managed service, future versions can be relicensed to BSL or
   dual-licensed with a commercial offering. Several established
   projects have done exactly this (Sentry, Elastic, Grafana,
   CockroachDB). Going *from* permissive *to* restrictive is socially
   contentious but legally clean; the reverse is much harder.

## What this means for users

You can do essentially anything with vcfclick:

- use it in private or commercial projects, modified or unmodified,
- redistribute it,
- include it in larger works,
- build commercial products on top of it.

The conditions are minimal: retain the copyright notice, state any
modifications, and include the licence with redistributions. The full
list is in `LICENSE`, section 4. There is no obligation to share back
your modifications.

## What this means for contributors

Contributions are accepted under the same Apache 2.0 terms (see
`LICENSE`, section 5). You retain the copyright on your contributions;
you grant the project an Apache-licensed copy. No CLA is required at
this time.
