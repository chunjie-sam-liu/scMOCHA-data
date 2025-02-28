#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-02-28 10:18:16
# @DESCRIPTION:
# @VERSION: v0.0.1

gses=(
  "GSE183267"
  "GSE153056"
  "GSE163278"
  "GSE165087"
  "GSE173320"
  "GSE154109"
  "GSE182020"
  "GSE149136"
)

basedir="/home/liuc9/github/scMOCHA-data/data/scfoundation2/Bone-Marrow"

# /home/liuc9/github/scMOCHA-data/src/01-sra-metadata.R
sra_metadata() {
  for gse in "${gses[@]}"; do
    echo "Rscript /home/liuc9/github/scMOCHA-data/src/01-sra-metadata.R -g ${gse} -b ${basedir}"
    Rscript /home/liuc9/github/scMOCHA-data/src/01-sra-metadata.R -g ${gse} -b ${basedir}
    bash ${basedir}/${gse}/00.edirect.gds.${gse}.sh
  done
}
sra_metadata

# /home/liuc9/github/scMOCHA-data/src/02-sra-download-dump.R
sra_download_dump() {
  for gse in "${gses[@]}"; do
    echo "python /home/liuc9/github/scMOCHA-data/src/xml2json.py -i ${basedir}/${gse}/${gse}.edirect.gds.xml -o ${basedir}/${gse}/${gse}.edirect.gds.json"
    python /home/liuc9/github/scMOCHA-data/src/xml2json.py -i ${basedir}/${gse}/${gse}.edirect.gds.xml -o ${basedir}/${gse}/${gse}.edirect.gds.json -p
    python /home/liuc9/github/scMOCHA-data/src/json2sraruntable.py -r ${basedir}/${gse}/${gse}.edirect.gds.json
    echo "Rscript /home/liuc9/github/scMOCHA-data/src/02-sra-download-dump.R -g ${gse} -b ${basedir}"
    Rscript /home/liuc9/github/scMOCHA-data/src/02-sra-download-dump.R -g ${gse} -b ${basedir}
  done
}
# sra_download_dump
