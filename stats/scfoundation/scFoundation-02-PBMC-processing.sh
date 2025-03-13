#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-02-28 10:18:16
# @DESCRIPTION:
# @VERSION: v0.0.1

gses=(
  "GSE140881" # only one read
  "GSE142595" # under run, Seq-Well method
  "GSE143353" # done
  "GSE147794" # finished, some are not success, need rerun
  "GSE148215" # under run, done
  "GSE152981" # CITE-seq
  "GSE153098" # not run, CD8
  "GSE153421" # under run, done
  "GSE162528" # not run
  "GSE162708" # only on sample,
  "GSE163314" # under run, done
  "GSE163633" # under run, some are not success, need rerun, under run
  "GSE164690" # under run, some are not success, need rerun, done
  "GSE165087" # not run
  "GSE165496" # under run, not run.
  "GSE165822" # under run, not run.
  "GSE166638" # Seq-Well
  "GSE167825" # under run, some are not success, need rerun, under run
  "GSE168453" # under run, some are not success, need rerun, under run
  "GSE174125" # under run, done
  "GSE178756" # only one sample
  "GSE179566" # only one sample
  "GSE184703" # under run, some are not success, need rerun, under run
  "GSE190839" # not run
  "GSE192391" # not run
)

gses_bone_marrow=(
  "GSE183267"
  "GSE153056"
  "GSE163278"
  "GSE165087"
  "GSE173320"
  "GSE154109"
  "GSE182020"
  "GSE149136"
)

basedir="/home/liuc9/github/scMOCHA-data/data/scfoundation2/PBMC"
basedir=$(realpath "$basedir")

# /home/liuc9/github/scMOCHA-data/src/01-sra-metadata.R
sra_metadata() {
  for gse in "${gses[@]}"; do
    echo "Rscript /home/liuc9/github/scMOCHA-data/src/01-sra-metadata.R -g ${gse} -b ${basedir}"
    Rscript /home/liuc9/github/scMOCHA-data/src/01-sra-metadata.R -g ${gse} -b ${basedir} &
  done
}
# sra_metadata

# SraRunTable

sra_run_table_gseid() {
  local gse=$1
  # create 00.edirect.gds.${gse}.sh edirect gds download xml
  echo "esearch -db gds -query '${gse}[Accession]' | efetch -format docsum > ${basedir}/${gse}/${gse}.edirect.gds.xml" >${basedir}/${gse}/00.edirect.gds.${gse}.sh
  bash ${basedir}/${gse}/00.edirect.gds.${gse}.sh

  # xml2json
  # python /home/liuc9/github/scMOCHA-data/src/gds_xml2json.py -i ${basedir}/${gse}/${gse}.edirect.gds.xml -o ${basedir}/${gse}/${gse}.edirect.gds.json -p
  #
  if [[ ! -f "${basedir}/${gse}/${gse}.edirect.gds.json" ]] || [[ ! -s "${basedir}/${gse}/${gse}.edirect.gds.json" ]]; then
    python /home/liuc9/github/scMOCHA-data/src/gds_xml2json.py -i ${basedir}/${gse}/${gse}.edirect.gds.xml -o ${basedir}/${gse}/${gse}.edirect.gds.json -p
  else
    echo "File ${basedir}/${gse}/${gse}.edirect.gds.json already exists and is not empty, skipping xml2json conversion"
  fi
  # json2sraruntable
  # python /home/liuc9/github/scMOCHA-data/src/json2sraruntable.py -r ${basedir}/${gse}/${gse}.edirect.gds.json
  if [[ ! -f "${basedir}/${gse}/${gse}.SraRunTable" ]] || [[ ! -s "${basedir}/${gse}/${gse}.SraRunTable" ]]; then
    python /home/liuc9/github/scMOCHA-data/src/json2sraruntable.py -r ${basedir}/${gse}/${gse}.edirect.gds.json
  else
    echo "File ${basedir}/${gse}/${gse}.SraRunTable already exists and is not empty, skipping json2sraruntable conversion"
  fi
  # biosample_runinfo2csv.py
  # python /home/liuc9/github/scMOCHA-data/src/biosample_runinfo2csv.py -i ${basedir}/${gse}/${gse}.edirect.biosample.runinfo -o ${basedir}/${gse}/${gse}.edirect.biosample.csv
  if [[ ! -f "${basedir}/${gse}/${gse}.edirect.biosample.runinfo" ]] || [[ ! -s "${basedir}/${gse}/${gse}.edirect.biosample.runinfo" ]]; then
    python /home/liuc9/github/scMOCHA-data/src/biosample_runinfo2csv.py -i ${basedir}/${gse}/${gse}.edirect.biosample.runinfo -o ${basedir}/${gse}/${gse}.edirect.biosample.csv
  else
    echo "File ${basedir}/${gse}/${gse}.SraRunTable.csv already exists and is not empty, skipping biosample_runinfo2csv conversion"
  fi
}

sra_run_table_gseid_force() {
  local gse=$1
  # create 00.edirect.gds.${gse}.sh edirect gds download xml
  echo "esearch -db gds -query '${gse}[Accession]' | efetch -format docsum > ${basedir}/${gse}/${gse}.edirect.gds.xml" >${basedir}/${gse}/00.edirect.gds.${gse}.sh
  bash ${basedir}/${gse}/00.edirect.gds.${gse}.sh

  # xml2json
  python /home/liuc9/github/scMOCHA-data/src/gds_xml2json.py -i ${basedir}/${gse}/${gse}.edirect.gds.xml -o ${basedir}/${gse}/${gse}.edirect.gds.json -p

  # json2sraruntable
  python /home/liuc9/github/scMOCHA-data/src/json2sraruntable.py -r ${basedir}/${gse}/${gse}.edirect.gds.json

  # biosample_runinfo2csv
  python /home/liuc9/github/scMOCHA-data/src/biosample_runinfo2csv.py -i ${basedir}/${gse}/${gse}.edirect.biosample.runinfo -o ${basedir}/${gse}/${gse}.edirect.biosample.csv
}

sra_run_table() {
  for gse in "${gses[@]}"; do
    # sra_run_table_gseid "${gse}" #&
    sra_run_table_gseid_force "${gse}" #&
  done
}
# sra_run_table

# /home/liuc9/github/scMOCHA-data/src/02-sra-download-dump.R
sra_download_dump() {
  for gse in "${gses[@]}"; do
    echo "Rscript /home/liuc9/github/scMOCHA-data/src/02-sra-download-dump.R -g ${gse} -b ${basedir}"
    Rscript /home/liuc9/github/scMOCHA-data/src/02-sra-download-dump.R -g ${gse} -b ${basedir} &
  done
}
# sra_download_dump

# /home/liuc9/github/scMOCHA-data/data/scfoundation/GSE140881/00.${gseid}.prefetch.slrm
prefetch() {
  for gse in "${gses[@]}"; do
    echo "bash ${basedir}/${gse}/00.${gse}.prefetch.sh"
    cd ${basedir}/${gse}
    # bash 00.${gse}.prefetch.sh
    sbatch 00.${gse}.prefetch.slrm
  done
}
# prefetch

# /home/liuc9/github/scMOCHA-data/data/scfoundation/GSE140881/01.${gseid}.prefetch.check.sh
prefetch_check() {
  for gse in "${gses[@]}"; do
    # echo "bash ${basedir}/${gse}/01.${gse}.prefetch.check.sh"
    cd ${basedir}/${gse}
    # bash 01.${gse}.prefetch.check.sh
    sbatch 01.${gse}.prefetch.check.slrm
  done
}
# prefetch_check

# /home/liuc9/github/scMOCHA-data/data/scfoundation2/PBMC/GSE140881/02.GSE140881.dump.slrm
dump_slrm() {
  for gse in "${gses[@]}"; do
    echo "bash ${basedir}/${gse}/02.${gse}.dump.slrm"
    cd ${basedir}/${gse}
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
# sra_rename_gsm_merge

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
    echo "bash ${basedir}/${gse}/04.${gse}.batch.sh"
    cd ${basedir}/${gse}
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
    echo "bash ${basedir}/${gse}/05.${gse}.scmocha.cptargz.sh"
    cd ${basedir}/${gse}
    bash 05.${gse}.scmocha.cptargz.sh
  done
}
# cptargz

# /home/liuc9/github/scMOCHA-data/data/GSE279945/07.GSE279945.scmocha.untargz.sh
untargz() {
  for gse in "${gses[@]}"; do
    echo "bash ${basedir}/${gse}/07.${gse}.scmocha.untargz.sh"
    cd ${basedir}/${gse}
    bash 07.${gse}.scmocha.untargz.sh
  done
}
# untargz

# /home/liuc9/github/scMOCHA-data/src/06-collect-variants.R
collect_variants() {
  for gse in "${gses[@]}"; do
    echo "Rscript /home/liuc9/github/scMOCHA-data/src/06-collect-variants.R -g ${gse} -b ${basedir}"
    Rscript /home/liuc9/github/scMOCHA-data/src/06-collect-variants.R -g ${gse} -b ${basedir} &
  done
}
collect_variants
