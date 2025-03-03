#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-02-28 10:18:16
# @DESCRIPTION:
# @VERSION: v0.0.1

gses=(
  "GSE165496"
  "GSE165822"
  "GSE140881"
  "GSE162708"
  "GSE163314"
  "GSE162528"
  "GSE165087"
  "GSE178756"
  "GSE192391"
  "GSE168453"
  "GSE164690"
  "GSE167825"
  "GSE163633"
  "GSE166638"
  "GSE190839"
  "GSE142595"
  "GSE147794"
  "GSE153098"
  "GSE179566"
  "GSE184703"
  "GSE143353"
  "GSE148215"
  "GSE152981"
  "GSE153421"
  "GSE174125"
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

# /home/liuc9/github/scMOCHA-data/src/01-sra-metadata.R
sra_metadata() {
  for gse in "${gses[@]}"; do
    echo "Rscript /home/liuc9/github/scMOCHA-data/src/01-sra-metadata.R -g ${gse} -b ${basedir}"
    Rscript /home/liuc9/github/scMOCHA-data/src/01-sra-metadata.R -g ${gse} -b ${basedir}
  done
}
# sra_metadata

# SraRunTable
sra_run_table() {
  for gse in "${gses[@]}"; do
    # edirect gds download xml
    if [[ ! -f "${basedir}/${gse}/${gse}.edirect.gds.xml" ]]; then
      bash ${basedir}/${gse}/00.edirect.gds.${gse}.sh
    else
      echo "File ${basedir}/${gse}/${gse}.edirect.gds.xml already exists, skipping edirect gds download"
    fi
    # xml2json
    echo "python /home/liuc9/github/scMOCHA-data/src/xml2json.py -i ${basedir}/${gse}/${gse}.edirect.gds.xml -o ${basedir}/${gse}/${gse}.edirect.gds.json"
    if [[ ! -f "${basedir}/${gse}/${gse}.edirect.gds.json" ]]; then
      python /home/liuc9/github/scMOCHA-data/src/xml2json.py -i ${basedir}/${gse}/${gse}.edirect.gds.xml -o ${basedir}/${gse}/${gse}.edirect.gds.json -p
    else
      echo "File ${basedir}/${gse}/${gse}.edirect.gds.json already exists, skipping xml2json conversion"
    fi
    # json2sraruntable
    if [[ ! -f "${basedir}/${gse}/${gse}.SraRunTable" ]]; then
      python /home/liuc9/github/scMOCHA-data/src/json2sraruntable.py -r ${basedir}/${gse}/${gse}.edirect.gds.json
    else
      echo "File ${basedir}/${gse}/${gse}.SraRunTable already exists, skipping json2sraruntable conversion"
    fi
  done
}
# sra_run_table

# /home/liuc9/github/scMOCHA-data/src/02-sra-download-dump.R
sra_download_dump() {
  for gse in "${gses[@]}"; do
    echo "Rscript /home/liuc9/github/scMOCHA-data/src/02-sra-download-dump.R -g ${gse} -b ${basedir}"
    Rscript /home/liuc9/github/scMOCHA-data/src/02-sra-download-dump.R -g ${gse} -b ${basedir}
  done
}
sra_download_dump

# /home/liuc9/github/scMOCHA-data/data/scfoundation/GSE140881/00.${gseid}.prefetch.sh
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
    echo "bash ${basedir}/${gse}/01.${gse}.prefetch.check.sh"
    cd ${basedir}/${gse}
    bash 01.${gse}.prefetch.check.sh
  done
}
# prefetch_check
