# Citing vcfclick

vcfclick is archived on [Zenodo](https://zenodo.org/), which mints a
DOI for every tagged GitHub release plus a stable "concept DOI" that
always resolves to the latest release. The concept DOI is the one
you cite from a paper; version-specific DOIs are for reproducibility
when you need to pin to an exact release.

## Quick form

If you use vcfclick in published research, please cite the concept
DOI. Once Zenodo mints the first DOI for this project, replace the
placeholder below with the real one.

```bibtex
@software{nuin_vcfclick,
  author       = {Nuin, Paulo A. F.},
  title        = {vcfclick: embedded VCF databases with
                  auditable natural-language queries},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.PLACEHOLDER},
  url          = {https://github.com/nuin/vcfclick}
}
```

## Setup notes

The integration is one-time and works as follows:

1. **Authorize Zenodo on GitHub** — at
   <https://zenodo.org/account/settings/github/>, log in with the
   GitHub account that owns the repo and grant access.

2. **Toggle the `nuin/vcfclick` repo "on"** in that same page's
   repository list. The toggle adds a webhook to the repo that fires
   on each future GitHub release.

3. **Cut a new GitHub release** AFTER toggling. Zenodo does NOT
   retroactively DOI releases that existed before integration was
   enabled — see [Zenodo's docs][zenodo-gh] for the current rule.
   The smallest acceptable trigger is a no-code-change patch bump
   (e.g. v0.3.2) created via `gh release create`.

4. **Replace the ORCID placeholder.** [`.zenodo.json`](../.zenodo.json)
   currently lists the author ORCID as `0000-0000-0000-0000`. Replace
   it with a real ORCID before the first release-with-Zenodo so the
   DOI metadata ships clean.

5. **First DOI mints in ~5 minutes** of the GitHub release event.
   Update this file's BibTeX block with the real DOI once it appears
   at <https://zenodo.org/account/settings/github/repository/nuin/vcfclick>.

[zenodo-gh]: https://docs.github.com/en/repositories/archiving-a-github-repository/referencing-and-citing-content
