#!/usr/bin/env bash
# Metainfo ----------------------------------------------------------------
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2026-09-17
# @DESCRIPTION: Unpack each sample archive into its own directory under the
#               stage input root, extracting only the seven files this stage
#               reads. Idempotent: a sample whose files are already present is
#               skipped.
#               Usage: bash newplots/compare-with-mgatk/00-extract-archives.sh
# @VERSION: v0.1.0

set -euo pipefail

stagedir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repodir="$(cd "${stagedir}/../.." && pwd)"

set -o allexport
# shellcheck source=/dev/null
source "${repodir}/.env"
set +o allexport

root="${ISILON_BASE}/compare-with-mgatk"
outroot="${root}/samples"

if ! command -v unzip > /dev/null 2>&1; then
  echo "unzip not found on PATH" >&2
  exit 1
fi

# Keep in step with ARCHIVE_MEMBERS in config.R. The allele-count matrices,
# cell.rds, cell.signac.rds and the PNGs are never read by this stage, so
# skipping them extracts about 1.3 GB instead of 3.6 GB.
members=(
  cell.variant_stats.tsv.gz
  cell.variant_stats_mgatk_original.tsv.gz
  cell.cell_heteroplasmic_df.tsv.gz
  cell.cell_heteroplasmic_df_mgatk_original.tsv.gz
  cell.cell_heteroplasmic_df_raw.tsv.gz
  cell.depthTable.txt
  cell.coverage.txt.gz
)

# sample_id:archive. Two archives separate GSE from GSM with "-"; sample_id
# normalises that so one token is safe as a path and as a factor level.
archives=(
  "GSE149689_GSM4509019_3PV3:GSE149689_GSM4509019_3PV3.zip"
  "GSE163314_GSM4976997_3PV2:GSE163314_GSM4976997_3PV2.zip"
  "GSE163668_GSM4995445_5PR2:GSE163668-GSM4995445_5PR2.zip"
  "GSE181279_GSM5494116_5PPE:GSE181279-GSM5494116_5PPE.zip"
  "GSE271107_GSM8369876_3PV3:GSE271107_GSM8369876_3PV3.zip"
)

mkdir -p "${outroot}"

for entry in "${archives[@]}"; do
  sample_id="${entry%%:*}"
  archive="${entry#*:}"
  src="${root}/${archive}"
  dest="${outroot}/${sample_id}"

  if [[ ! -f "${src}" ]]; then
    echo "missing archive: ${src}" >&2
    exit 1
  fi

  mkdir -p "${dest}"

  missing=0
  for m in "${members[@]}"; do
    if [[ ! -s "${dest}/${m}" ]]; then
      missing=1
    fi
  done

  if [[ "${missing}" -eq 0 ]]; then
    echo "[skip]    ${sample_id}: ${#members[@]} files already present"
    continue
  fi

  echo "[extract] ${sample_id} from ${archive}"
  unzip -o -q -j "${src}" "${members[@]}" -d "${dest}"

  for m in "${members[@]}"; do
    if [[ ! -s "${dest}/${m}" ]]; then
      echo "extraction produced no ${m} for ${sample_id}" >&2
      exit 1
    fi
  done

  echo "[done]    ${sample_id}: $(du -sh "${dest}" | cut -f1)"
done

echo "input root: ${outroot}"
