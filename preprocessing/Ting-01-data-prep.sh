#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2024-12-19 15:03:11
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

ting_dir=/home/liuc9/large-scale/ting/raw
cj_dir=/home/liuc9/github/scMOCHA-data/data

# ! soft link --------------------------------------------------------------------

make_soft_link() {
  for gse in "${gses[@]}"; do
    if [ -e "${cj_dir}/${gse}" ]; then
      echo "${cj_dir}/${gse} exists"
      rm "${cj_dir}/${gse}"
    fi
    ln -s "${ting_dir}/${gse}" "${cj_dir}/${gse}"
  done
}

# make_soft_link

# ! remove Running results in the gse folder --------------------------------------------------------------------
rmmove_running_result() {
  for gse in "${gses[@]}"; do
    if [ -e "${cj_dir}/${gse}/Running_results" ]; then
      echo "rm -rf ${cj_dir}/${gse}/Running_results"
      rm -rf "${cj_dir}/${gse}/Running_results" &
    fi
  done
}

# rmmove_running_result
gses=(
  GSE214865
  GSE220189
  GSE233844
  GSE175499
)
# ! make sra dir and put the SRR* dir into sra dir --------------------------------------------------------------------
make_sra_dir() {
  for gse in "${gses[@]}"; do
    if [ -e "${cj_dir}/${gse}/sra" ]; then
      echo "${cj_dir}/${gse}/sra exists"
      # rm -rf "${cj_dir}/${gse}/sra"
    fi
    mkdir "${cj_dir}/${gse}/sra"
    mv ${cj_dir}/${gse}/SRR* "${cj_dir}/${gse}/sra"
  done
}
make_sra_dir

# ! move fastq file to gsm --------------------------------------------------------------------

move_fastq_to_gsm() {
  for gse in "${gses[@]}"; do
    if [ -e "${cj_dir}/${gse}/fastq" ]; then
      echo "${cj_dir}/${gse}/fastq exists"
      # rm -rf "${cj_dir}/${gse}/gsm"
      mv ${cj_dir}/${gse}/fastq ${cj_dir}/${gse}/gsm
    fi
  done
}
# move_fastq_to_gsm
