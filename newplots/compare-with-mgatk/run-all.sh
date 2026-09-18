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

samples=(
  GSE149689_GSM4509019_3PV3
  GSE163314_GSM4976997_3PV2
  GSE163668_GSM4995445_5PR2
  GSE181279_GSM5494116_5PPE
  GSE271107_GSM8369876_3PV3
)

if [[ $# -gt 0 ]]; then
  samples=("$@")
fi

steps=(
  01-load-harmonize
  02-cell-inclusion
  03-variant-funnel-overlap
  04-heteroplasmy-spectrum
  05-vmr-strand
)

bash "${stagedir}/00-extract-archives.sh"

cd "${repodir}"

for s in "${samples[@]}"; do
  for step in "${steps[@]}"; do
    echo "=== ${s} :: ${step} ==="
    pixi run Rscript "newplots/compare-with-mgatk/${step}.R" --sample="${s}"
  done
done

# 07 writes the cross-sample tables that 06 binds into the workbook, so it runs
# first and 06 runs last.
echo "=== cross-sample ==="
pixi run Rscript newplots/compare-with-mgatk/07-cross-sample.R
pixi run Rscript newplots/compare-with-mgatk/06-summary-workbook.R
