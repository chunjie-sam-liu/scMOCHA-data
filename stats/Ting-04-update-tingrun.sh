#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2025-01-20 14:55:00
# @DESCRIPTION:
# @VERSION: v0.0.1

gseids=(
  GSE214865
  GSE220189
  GSE233844
  GSE175499
)

cj_dir=/home/liuc9/github/scMOCHA-data/data
# make ${gseid}.srrid.list

make_srrid_list() {
  for gse in "${gseids[@]}"; do
    echo "Making ${gse}.srrid.list"
    if [ -f ${cj_dir}/${gse}/${gse}.srrid.list ]; then
      rm ${cj_dir}/${gse}/${gse}.srrid.list
    fi
    srrids=$(find ${cj_dir}/${gse}/final/ -maxdepth 1 -mindepth 1 -type d -exec basename {} \;)
    # echo ${srrids[@]} >${cj_dir}/${gse}/${gse}.srrid.list
    for srrid in ${srrids[@]}; do
      echo ${srrid} >>${cj_dir}/${gse}/${gse}.srrid.list
    done
  done
}
# make_srrid_list

parse_variants_gseid() {
  local gseid=$1
  echo "Parsing variants for $gseid"
  # Update the target directory
  data_dir="$basedir/$gseid"
  final_dir="$data_dir/final"

  Rscript /home/liuc9/github/scMOCHA-data/src/06-collect-variants.R -g ${gseid}

}

parse_variants_gseids() {
  for gseid in ${gseids[@]}; do
    parse_variants_gseid ${gseid} &
  done
}

parse_variants_gseids
