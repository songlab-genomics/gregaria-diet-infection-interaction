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
