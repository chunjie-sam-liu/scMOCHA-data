#!/usr/bin/env bash
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-09-17
# @DESCRIPTION: Run the whole stage: extract the archives, run steps 01 to 05
#               for every sample, then the cross-sample step and the workbook.
#               Usage: bash newplots/compare-with-mgatk/run-all.sh [sample_id ...]
# @VERSION: v0.1.0

set -euo pipefail

stagedir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repodir="$(cd "${stagedir}/../.." && pwd)"

# Read from the SAMPLES registry in config.R, never a copy kept here. A second
# list silently stops running whatever was added to the registry but not to the
# copy; this file held five samples while the registry held ten.
mapfile -t samples < <(
  pixi run --manifest-path "${repodir}/pixi.toml" Rscript -e '
    source(file.path(
      path.expand(Sys.getenv("REPODIR")),
      "newplots", "compare-with-mgatk", "config.R"
    ))
    cat(SAMPLE_IDS, sep = "\n")
  ' 2> /dev/null
)

if [[ "${#samples[@]}" -eq 0 ]]; then
  echo "could not read SAMPLE_IDS from config.R" >&2
  exit 1
fi

if [[ $# -gt 0 ]]; then
  samples=("$@")
fi

echo "running ${#samples[@]} samples"

steps=(
  01-load-harmonize
  02-cell-inclusion
  03-variant-funnel-overlap
  04-heteroplasmy-spectrum
  05-vmr-strand
  08-call-af-depth
)

bash "${stagedir}/00-extract-archives.sh"

cd "${repodir}"

for s in "${samples[@]}"; do
  for step in "${steps[@]}"; do
    echo "=== ${s} :: ${step} ==="
    pixi run Rscript "newplots/compare-with-mgatk/${step}.R" --sample="${s}"
  done
done

# 07 reads every sample's step 03 cache and the step 02 and 08 tables, and
# writes the cross-sample tables that 06 binds into the workbook, so 07 runs
# after all per-sample steps and 06 runs last.
echo "=== cross-sample ==="
pixi run Rscript newplots/compare-with-mgatk/07-cross-sample.R
pixi run Rscript newplots/compare-with-mgatk/06-summary-workbook.R
