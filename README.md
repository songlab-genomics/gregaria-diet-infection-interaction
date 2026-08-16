# Effect of Diet on Metarhizium Infection in S. Gregaria

A [workflowr][] project.

[workflowr]: https://github.com/workflowr/workflowr

## Data and full results

GitHub contains the analysis code, sample metadata, rRNA exclusion list, and
rendered website. Large analysis inputs and full-resolution result trees are
stored separately in the Dryad deposit working directory.

On the project computer, the current data package is located at:

```text
/Users/maevatecher/Dropbox/3. Research Papers/Locusts/Ongoing/gregaria-diet-infection-interaction-dryad
```

After cloning the GitHub repository or downloading the Dryad package, set the
data location and create the project links once:

```bash
export GREGARIA_DIET_DRYAD_DIR=/path/to/gregaria-diet-infection-interaction-dryad
bash scripts/setup_dryad_links.sh
```

The R Markdown pages keep project-relative paths such as
`data/reference/...` and `output/rmd_runs/...`. The setup script links those
paths to the external data package, so the analysis code remains portable and
does not depend on files committed to GitHub.

## Zenodo code archive

Create the Zenodo code-and-website zip from a clean commit rather than zipping
the working directory directly:

```bash
bash scripts/build_zenodo_release.sh v1.0.0
```

The script uses `git archive`, validates the resulting zip, and reports its
SHA-256 checksum. It includes only committed code, metadata, and rendered site
files. Local R history files, `.git`, and the absolute links to the separately
deposited Dryad data are intentionally excluded.

## GitHub Pages deployment

The rendered workflowR website is committed under `docs/`. In the repository
Pages settings, use **Deploy from a branch**, with branch `main` and folder
`/docs`. The included Pages workflow is retained as a manual fallback only; use
it only after intentionally changing the repository Pages source to
**GitHub Actions**. Keeping both deployment modes active on every push can
create competing Pages deployments.
