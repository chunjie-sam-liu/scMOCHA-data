#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2024-12-19 15:21:45
# @DESCRIPTION:
# @VERSION: v0.0.1

gses=(
  GSE161354
  GSE175524
  GSE206283
  GSE226598
  GSE235050
  GSE261140
  GSE279945
)

done=(
  GSE161354 # done
  GSE175524 # done
  GSE206283 # done
  GSE226598 # done
  GSE235050 # done
  GSE261140 # some errors
  GSE279945 # some errors
)

# /home/liuc9/github/scMOCHA-data/src/01-sra-metadata.R
sra_metadata() {
  for gse in "${gses[@]}"; do
    echo "Rscript /home/liuc9/github/scMOCHA-data/src/01-sra-metadata.R -g ${gse}"
    Rscript /home/liuc9/github/scMOCHA-data/src/01-sra-metadata.R -g ${gse} &
  done
}

# /home/liuc9/github/scMOCHA-data/src/02-sra-download-dump.R
sra_download_dump() {
  for gse in "${gses[@]}"; do
    echo "Rscript /home/liuc9/github/scMOCHA-data/src/02-sra-download-dump.R -g ${gse}"
    Rscript /home/liuc9/github/scMOCHA-data/src/02-sra-download-dump.R -g ${gse} &
  done
}
# sra_download_dump

# /home/liuc9/github/scMOCHA-data/src/03-sra-rename-gsm-merge.R
sra_rename_gsm_merge() {
  for gse in "${gses[@]}"; do
    echo "Rscript /home/liuc9/github/scMOCHA-data/src/03-sra-rename-gsm-merge.R -g ${gse}"
    Rscript /home/liuc9/github/scMOCHA-data/src/03-sra-rename-gsm-merge.R -g ${gse} &
  done
}
# sra_rename_gsm_merge

# /home/liuc9/github/scMOCHA-data/src/04-scmocha-conf.R
scmocha_conf() {
  for gse in "${gses[@]}"; do
    echo "Rscript /home/liuc9/github/scMOCHA-data/src/04-scmocha-conf.R -g ${gse}"
    Rscript /home/liuc9/github/scMOCHA-data/src/04-scmocha-conf.R -g ${gse} &
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
    echo "Rscript /home/liuc9/github/scMOCHA-data/src/05-parse-log.R -g ${gse}"
    Rscript /home/liuc9/github/scMOCHA-data/src/05-parse-log.R -g ${gse} &
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
    echo "Rscript /home/liuc9/github/scMOCHA-data/src/06-collect-variants.R -g ${gse}"
    Rscript /home/liuc9/github/scMOCHA-data/src/06-collect-variants.R -g ${gse} &
  done
}
collect_variants
