#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-01-02 15:14:02
# @DESCRIPTION:
# @VERSION: v0.0.1

gses=(
  "GSE140881"
  "GSE142595"
  "GSE149313"
  "GSE154386"
  "GSE159117"
  "GSE162117"
  "GSE167825"
  "GSE179566"
  "GSE188632"
  "GSE192391"
)

done=(
  "GSE140881"
  "GSE142595"
  "GSE149313"
  "GSE154386"
  "GSE159117"
  "GSE162117"
  "GSE167825"
  "GSE179566"
  "GSE188632"
  "GSE192391"
)

basedir="/home/liuc9/github/scMOCHA-data/data/scfoundation"

# /home/liuc9/github/scMOCHA-data/src/01-sra-metadata.R
sra_metadata() {
  for gse in "${gses[@]}"; do
    echo "Rscript /home/liuc9/github/scMOCHA-data/src/01-sra-metadata.R -g ${gse} -b ${basedir}"
    Rscript /home/liuc9/github/scMOCHA-data/src/01-sra-metadata.R -g ${gse} -b ${basedir} &
  done
}
# sra_metadata

# /home/liuc9/github/scMOCHA-data/src/02-sra-download-dump.R
sra_download_dump() {
  for gse in "${gses[@]}"; do
    echo "Rscript /home/liuc9/github/scMOCHA-data/src/02-sra-download-dump.R -g ${gse} -b ${basedir}"
    Rscript /home/liuc9/github/scMOCHA-data/src/02-sra-download-dump.R -g ${gse} -b ${basedir} &
  done
}
# sra_download_dump

# /home/liuc9/github/scMOCHA-data/data/scfoundation/GSE140881/00.${gseid}.prefetch.sh
prefetch() {
  for gse in "${gses[@]}"; do
    echo "bash /home/liuc9/github/scMOCHA-data/data/scfoundation/${gse}/00.${gse}.prefetch.sh"
    cd /home/liuc9/github/scMOCHA-data/data/scfoundation/${gse}
    bash 00.${gse}.prefetch.sh
  done
}
# prefetch

# /home/liuc9/github/scMOCHA-data/data/scfoundation/GSE140881/01.${gseid}.prefetch.check.sh
prefetch_check() {
  for gse in "${gses[@]}"; do
    echo "bash /home/liuc9/github/scMOCHA-data/data/scfoundation/${gse}/01.${gse}.prefetch.check.sh"
    cd /home/liuc9/github/scMOCHA-data/data/scfoundation/${gse}
    bash 01.${gse}.prefetch.check.sh
  done
}
# prefetch_check

# rm all GL KI and NC_ files
rm_gl_ki_nc() {
  for gse in "${gses[@]}"; do
    # gse="GSE179566"
    echo "rm /home/liuc9/github/scMOCHA-data/data/scfoundation/${gse}/sra/*/GL*"
    rm /home/liuc9/github/scMOCHA-data/data/scfoundation/${gse}/sra/*/GL*
    rm /home/liuc9/github/scMOCHA-data/data/scfoundation/${gse}/sra/*/KI*
    rm /home/liuc9/github/scMOCHA-data/data/scfoundation/${gse}/sra/*/NC_*
  done
}

# rm_gl_ki_nc

# /home/liuc9/github/scMOCHA-data/data/scfoundation/GSE140881/02.${gseid}.dump.slrm
dump_slrm() {
  for gse in "${gses[@]}"; do
    echo "bash /home/liuc9/github/scMOCHA-data/data/scfoundation/${gse}/02.${gse}.dump.slrm"
    cd /home/liuc9/github/scMOCHA-data/data/scfoundation/${gse}
    sbatch 02.${gse}.dump.slrm
  done
}

# dump_slrm

# /home/liuc9/github/scMOCHA-data/src/03-sra-rename-gsm-merge.R
sra_rename_gsm_merge() {
  for gse in "${gses[@]}"; do
    echo "Rscript /home/liuc9/github/scMOCHA-data/src/03-sra-rename-gsm-merge.R -g ${gse} -b ${basedir}"
    Rscript /home/liuc9/github/scMOCHA-data/src/03-sra-rename-gsm-merge.R -g ${gse} -b ${basedir} &
  done
}
sra_rename_gsm_merge

# /home/liuc9/github/scMOCHA-data/src/04-scmocha-conf.R
scmocha_conf() {
  for gse in "${gses[@]}"; do
    echo "Rscript /home/liuc9/github/scMOCHA-data/src/04-scmocha-conf.R -g ${gse} -b ${basedir}"
    Rscript /home/liuc9/github/scMOCHA-data/src/04-scmocha-conf.R -g ${gse} -b ${basedir} &
  done
}
# scmocha_conf

# run 04.{gse}.batch.sh
scmocha_batch_run() {
  for gse in "${gses[@]}"; do
    echo "bash /home/liuc9/github/scMOCHA-data/data/04.${gse}.batch.sh"
    cd /home/liuc9/github/scMOCHA-data/data/${gse}
    bash 04.${gse}.batch.sh
  done
}
# scmocha_batch_run

# /home/liuc9/github/scMOCHA-data/src/05-parse-log.R
parse_log() {
  for gse in "${gses[@]}"; do
    echo "Rscript /home/liuc9/github/scMOCHA-data/src/05-parse-log.R -g ${gse} -b ${basedir}"
    Rscript /home/liuc9/github/scMOCHA-data/src/05-parse-log.R -g ${gse} -b ${basedir} &
  done
}
# parse_log

# /home/liuc9/github/scMOCHA-data/data/GSE279945/05.{gseid}.scmocha.cptargz.sh

cptargz() {
  for gse in "${gses[@]}"; do
    echo "bash /home/liuc9/github/scMOCHA-data/data/${gse}/05.${gse}.scmocha.cptargz.sh"
    cd /home/liuc9/github/scMOCHA-data/data/${gse}
    bash 05.${gse}.scmocha.cptargz.sh
  done
}
# cptargz

# /home/liuc9/github/scMOCHA-data/data/GSE279945/07.GSE279945.scmocha.untargz.sh
untargz() {
  for gse in "${gses[@]}"; do
    echo "bash /home/liuc9/github/scMOCHA-data/data/${gse}/07.${gse}.scmocha.untargz.sh"
    cd /home/liuc9/github/scMOCHA-data/data/${gse}
    bash 07.${gse}.scmocha.untargz.sh
  done
}
# untargz

# wait
# /home/liuc9/github/scMOCHA-data/src/06-collect-variants.R
collect_variants() {
  for gse in "${gses[@]}"; do
    echo "Rscript /home/liuc9/github/scMOCHA-data/src/06-collect-variants.R -g ${gse} -b ${basedir}"
    Rscript /home/liuc9/github/scMOCHA-data/src/06-collect-variants.R -g ${gse} -b ${basedir} &
  done
}
# collect_variants
