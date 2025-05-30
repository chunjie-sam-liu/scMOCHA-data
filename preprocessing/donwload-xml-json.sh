#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-03-03 14:57:47
# @DESCRIPTION:
# @VERSION: v0.0.1

gses=(
  "GSE149689"
  "GSE155223"
  "GSE159117"
  "GSE163668"
  "GSE171555"
  "GSE181279"
  "GSE214865"
  "GSE226602"
  "GSE261140"
  "GSE154386"
  "GSE155673"
  "GSE161354"
  "GSE166992"
  "GSE175499"
  "GSE188632"
  "GSE220189"
  "GSE233844"
  "GSE279945"
  "GSE149313"
  "GSE154567"
  "GSE157344"
  "GSE162117"
  "GSE167825"
  "GSE175524"
  "GSE206283"
  "GSE226598"
  "GSE235050"
)
basedir="/mnt/isilon/u01_project/large-scale/liuc9/raw"
basedir=$(realpath "$basedir")

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
    sra_run_table_gseid_force "${gse}" #&
  done
}
sra_run_table
# sra_run_table_gseid_force GSE233844
