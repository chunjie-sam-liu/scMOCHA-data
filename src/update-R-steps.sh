#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2024-11-26 15:53:38
# @DESCRIPTION:

# Number of input parameters
param=$#
basedir="/mnt/isilon/u01_project/large-scale/liuc9/raw"
# gseids=(GSE149689 GSE154567 GSE155223 GSE155673 GSE157344 GSE163668 GSE166992 GSE171555 GSE181279 GSE226602)
# gseids=(GSE149689 GSE155223 GSE155673 GSE157344 GSE163668 GSE166992 GSE171555 GSE181279 GSE226602)
gseids=(GSE149689 GSE155223 GSE157344 GSE163668 GSE166992 GSE171555 GSE181279 GSE226602)

cp_final_dir() {
  local gseid=$1
  echo "Updating $gseid"
  # Update the target directory
  data_dir="$basedir/$gseid"
  final_dir="$data_dir/final"

  # check if $final_dir has any tar.gz file, if not run ${data_dir}/05.${gseid}.scmocha.cpfinal.sh
  if [ -z "$(ls -A $final_dir)" ]; then
    echo "No tar.gz files in $final_dir"
    cmd="bash ${data_dir}/05.${gseid}.scmocha.cpfinal.sh"
    echo $cmd
    eval $cmd
  else
    echo "There are tar.gz files in $final_dir"
  fi
}

unzip_final_dir() {
  local gseid=$1
  echo "Updating $gseid"
  # Update the target directory
  data_dir="$basedir/$gseid"
  final_dir="$data_dir/final"
  # unzip_dir="$data_dir/unzip"

  # check if $unzip_dir has any tar.gz file, if not run ${data_dir}/06.${gseid}.scmocha.unzip.sh
  if [ -z "$(ls -A $unzip_dir)" ]; then
    # echo "No tar.gz files in $unzip_dir"
    bash ${data_dir}/06.${gseid}.scmocha.clear.sh
    bash ${data_dir}/07.${gseid}.scmocha.unfinal.sh

  else
    echo "There are tar.gz files in $unzip_dir"
  fi
}

# update_final_dir "GSE171555"
cp_final_dirs() {
  gseids_for_final=(GSE157344 GSE171555)
  for gseid in ${gseids_for_final[@]}; do
    update_final_dir $gseid
  done
}

# update_final_dirs

unzip_final_dirs() {
  gseids_for_unzip=(GSE157344 GSE171555)
  for gseid in ${gseids_for_unzip[@]}; do
    unzip_final_dir $gseid
  done
}

# unzip_final_dirs

parse_variants_gseid() {
  local gseid=$1
  echo "Parsing variants for $gseid"
  # Update the target directory
  data_dir="$basedir/$gseid"
  final_dir="$data_dir/final"

  Rscript /home/liuc9/github/scMOCHA-data/src/06-parse-variants.R -g ${gseid}

}

parse_variants_gseids() {
  for gseid in ${gseids[@]}; do
    parse_variants_gseid ${gseid} &
  done
}

parse_variants_gseids
