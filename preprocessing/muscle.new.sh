#!/usr/bin/env bash
# Author: Chunjie Liu
# Contact: chunjie.sam.liu.at.gmail.com
# Date: 2025-01-20
# Description: Optimized parallel scMOCHA data processing pipeline with proper synchronization and resource management
# Version: 0.2

# Configuration Variables
MAX_CONCURRENT_PREFETCH=2
MAX_CONCURRENT_DUMPS=3
MAX_CONCURRENT_SCMOCHA=2
SLURM_CHECK_INTERVAL=60
FILE_CHECK_INTERVAL=30
MAX_RETRIES=3

gses=(
  GSE130646
  GSE143704
)

cj_dir=/mnt/isilon/u01_project/large-scale/liuc9/raw/Muscle
basedir=/mnt/isilon/u01_project/large-scale/liuc9/raw/Muscle
basedir=$(realpath "$basedir")

# Utility Functions
log_message() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

log_error() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $1" >&2
}

# Function to wait for SLURM jobs to complete
wait_for_slurm_jobs() {
  local job_pattern="$1"
  local timeout_minutes="${2:-1440}"  # Default 24 hours
  local start_time=$(date +%s)

  log_message "Waiting for SLURM jobs matching pattern: ${job_pattern}"

  while true; do
    local job_count=$(squeue -u $(whoami) --name="${job_pattern}" --noheader 2>/dev/null | wc -l)
    if [ "$job_count" -eq 0 ]; then
      log_message "All ${job_pattern} jobs completed"
      break
    fi

    # Check timeout
    local current_time=$(date +%s)
    local elapsed_minutes=$(( (current_time - start_time) / 60 ))
    if [ "$elapsed_minutes" -ge "$timeout_minutes" ]; then
      log_error "Timeout waiting for ${job_pattern} jobs after ${timeout_minutes} minutes"
      return 1
    fi

    log_message "Still waiting for ${job_count} ${job_pattern} jobs... (${elapsed_minutes}/${timeout_minutes} min)"
    sleep $SLURM_CHECK_INTERVAL
  done
  return 0
}

# Function to check if files exist and are not empty
check_files_exist() {
  local files=("$@")
  for file in "${files[@]}"; do
    local retries=0
    while [ ! -f "$file" ] || [ ! -s "$file" ]; do
      if [ $retries -ge $MAX_RETRIES ]; then
        log_error "File $file not found or empty after $MAX_RETRIES attempts"
        return 1
      fi
      log_message "Waiting for $file... (attempt $((retries + 1))/$MAX_RETRIES)"
      sleep $FILE_CHECK_INTERVAL
      ((retries++))
    done
    log_message "File verified: $file"
  done
  return 0
}

# Function to wait for background processes
wait_for_background_jobs() {
  local function_name="$1"
  log_message "Waiting for all ${function_name} background processes to complete..."
  wait
  local exit_code=$?
  if [ $exit_code -eq 0 ]; then
    log_message "${function_name} completed successfully"
  else
    log_error "${function_name} completed with errors (exit code: $exit_code)"
    return $exit_code
  fi
  return 0
}

# Function to run with retry logic
run_with_retry() {
  local command="$1"
  local max_attempts="${2:-$MAX_RETRIES}"
  local attempt=1

  while [ $attempt -le $max_attempts ]; do
    log_message "Executing: $command (attempt $attempt/$max_attempts)"
    if eval "$command"; then
      log_message "Command succeeded: $command"
      return 0
    else
      log_error "Command failed: $command (attempt $attempt/$max_attempts)"
      if [ $attempt -lt $max_attempts ]; then
        local wait_time=$((attempt * 30))  # Exponential backoff
        log_message "Retrying in ${wait_time} seconds..."
        sleep $wait_time
      fi
      ((attempt++))
    fi
  done

  log_error "Command failed after $max_attempts attempts: $command"
  return 1
}

# Stage 1: SRA Metadata Collection (Already parallelized)
sra_metadata() {
  log_message "Starting SRA metadata collection for ${#gses[@]} GSEs"
  for gse in "${gses[@]}"; do
    log_message "Starting metadata collection for $gse"
    Rscript /home/liuc9/github/scMOCHA-data/src/01-sra-metadata.R -g ${gse} -b ${basedir} &
  done
  wait_for_background_jobs "sra_metadata"
}

# Stage 2: SRA Run Table Generation (Now parallelized)
sra_run_table_gseid() {
  local gse=$1
  log_message "Processing SRA run table for $gse"

  # Create edirect script
  echo "esearch -db gds -query '${gse}[Accession]' | efetch -format docsum > ${basedir}/${gse}/${gse}.edirect.gds.xml" > ${basedir}/${gse}/00.edirect.gds.${gse}.sh

  # Execute edirect
  if ! run_with_retry "bash ${basedir}/${gse}/00.edirect.gds.${gse}.sh"; then
    log_error "Failed to download GDS XML for $gse"
    return 1
  fi

  # XML to JSON conversion
  if [[ ! -f "${basedir}/${gse}/${gse}.edirect.gds.json" ]] || [[ ! -s "${basedir}/${gse}/${gse}.edirect.gds.json" ]]; then
    if ! run_with_retry "python /home/liuc9/github/scMOCHA-data/src/gds_xml2json.py -i ${basedir}/${gse}/${gse}.edirect.gds.xml -o ${basedir}/${gse}/${gse}.edirect.gds.json -p"; then
      log_error "Failed to convert XML to JSON for $gse"
      return 1
    fi
  else
    log_message "JSON file already exists for $gse, skipping conversion"
  fi

  # JSON to SRA run table
  if [[ ! -f "${basedir}/${gse}/${gse}.SraRunTable" ]] || [[ ! -s "${basedir}/${gse}/${gse}.SraRunTable" ]]; then
    if ! run_with_retry "python /home/liuc9/github/scMOCHA-data/src/json2sraruntable.py -r ${basedir}/${gse}/${gse}.edirect.gds.json"; then
      log_error "Failed to generate SRA run table for $gse"
      return 1
    fi
  else
    log_message "SRA run table already exists for $gse, skipping generation"
  fi

  # Biosample runinfo to CSV
  if [[ -f "${basedir}/${gse}/${gse}.edirect.biosample.runinfo" ]] && [[ -s "${basedir}/${gse}/${gse}.edirect.biosample.runinfo" ]]; then
    if [[ ! -f "${basedir}/${gse}/${gse}.edirect.biosample.csv" ]] || [[ ! -s "${basedir}/${gse}/${gse}.edirect.biosample.csv" ]]; then
      if ! run_with_retry "python /home/liuc9/github/scMOCHA-data/src/biosample_runinfo2csv.py -i ${basedir}/${gse}/${gse}.edirect.biosample.runinfo -o ${basedir}/${gse}/${gse}.edirect.biosample.csv"; then
        log_error "Failed to convert biosample runinfo to CSV for $gse"
        return 1
      fi
    else
      log_message "Biosample CSV already exists for $gse, skipping conversion"
    fi
  else
    log_message "No biosample runinfo file found for $gse, skipping CSV conversion"
  fi

  log_message "Completed SRA run table processing for $gse"
  return 0
}

sra_run_table() {
  log_message "Starting SRA run table generation for ${#gses[@]} GSEs"
  for gse in "${gses[@]}"; do
    sra_run_table_gseid "${gse}" &
  done
  wait_for_background_jobs "sra_run_table"

  # Verify required files exist
  log_message "Verifying SRA run table files..."
  for gse in "${gses[@]}"; do
    local required_files=(
      "${basedir}/${gse}/${gse}.edirect.gds.json"
      "${basedir}/${gse}/${gse}.SraRunTable"
    )
    if ! check_files_exist "${required_files[@]}"; then
      log_error "Required files missing for $gse"
      return 1
    fi
  done
}

# Stage 3: SRA Download Setup (Already parallelized)
sra_download_dump() {
  log_message "Starting SRA download dump setup for ${#gses[@]} GSEs"
  for gse in "${gses[@]}"; do
    log_message "Setting up download scripts for $gse"
    Rscript /home/liuc9/github/scMOCHA-data/src/02-sra-download-dump.R -g ${gse} -b ${basedir} &
  done
  wait_for_background_jobs "sra_download_dump"

  # Verify setup files exist
  log_message "Verifying download setup files..."
  for gse in "${gses[@]}"; do
    local required_files=(
      "${basedir}/${gse}/00.${gse}.prefetch.sh"
      "${basedir}/${gse}/01.${gse}.prefetch.check.sh"
      "${basedir}/${gse}/02.${gse}.dump.slrm"
    )
    if ! check_files_exist "${required_files[@]}"; then
      log_error "Required setup files missing for $gse"
      return 1
    fi
  done
}

# Stage 4: Data Prefetching (Now resource-aware parallelized)
prefetch_sh() {
  log_message "Starting prefetch for ${#gses[@]} GSEs with max ${MAX_CONCURRENT_PREFETCH} concurrent"
  local active_jobs=0
  local pids=()

  for gse in "${gses[@]}"; do
    # Wait if we've reached the concurrent limit
    if [ $active_jobs -ge $MAX_CONCURRENT_PREFETCH ]; then
      log_message "Waiting for prefetch slot to become available..."
      wait ${pids[0]}  # Wait for the first job to complete
      pids=("${pids[@]:1}")  # Remove first PID from array
      ((active_jobs--))
    fi

    log_message "Starting prefetch for $gse"
    (
      cd ${basedir}/${gse}
      if ! run_with_retry "bash 00.${gse}.prefetch.sh"; then
        log_error "Prefetch failed for $gse"
        exit 1
      fi
      log_message "Prefetch completed for $gse"
    ) &

    local pid=$!
    pids+=($pid)
    ((active_jobs++))
  done

  # Wait for all remaining prefetch jobs
  log_message "Waiting for all prefetch jobs to complete..."
  for pid in "${pids[@]}"; do
    wait $pid
    if [ $? -ne 0 ]; then
      log_error "A prefetch job failed"
      return 1
    fi
  done
  log_message "All prefetch jobs completed successfully"
}

# Stage 5: Prefetch Validation (Now parallelized)
prefetch_check() {
  log_message "Starting prefetch validation for ${#gses[@]} GSEs"
  for gse in "${gses[@]}"; do
    (
      cd ${basedir}/${gse}
      log_message "Validating prefetch for $gse"
      if ! run_with_retry "bash 01.${gse}.prefetch.check.sh"; then
        log_error "Prefetch validation failed for $gse"
        exit 1
      fi
      log_message "Prefetch validation completed for $gse"
    ) &
  done
  wait_for_background_jobs "prefetch_check"
}

# Stage 6: Data Dumping (Resource-aware with SLURM monitoring)
dump_slrm() {
  log_message "Starting dump jobs for ${#gses[@]} GSEs with max ${MAX_CONCURRENT_DUMPS} concurrent"
  local active_jobs=0
  local submitted_jobs=()

  for gse in "${gses[@]}"; do
    # Wait if we've reached the concurrent limit
    if [ $active_jobs -ge $MAX_CONCURRENT_DUMPS ]; then
      log_message "Reached concurrent dump limit, waiting for jobs to complete..."
      wait_for_slurm_jobs "dump" 720  # 12 hour timeout for dump jobs
      active_jobs=0
      submitted_jobs=()
    fi

    log_message "Submitting dump job for $gse"
    cd ${basedir}/${gse}
    if sbatch 02.${gse}.dump.slrm; then
      submitted_jobs+=($gse)
      ((active_jobs++))
      log_message "Dump job submitted for $gse"
    else
      log_error "Failed to submit dump job for $gse"
      return 1
    fi
  done

  # Wait for all dump jobs to complete
  log_message "Waiting for all dump jobs to complete..."
  wait_for_slurm_jobs "dump" 1440  # 24 hour timeout
}

# Stage 7: GSM Renaming and Merging (Already parallelized, add dependency check)
sra_rename_gsm_merge() {
  log_message "Starting GSM renaming and merging for ${#gses[@]} GSEs"

  # First verify that dump stage completed successfully
  for gse in "${gses[@]}"; do
    log_message "Checking for FASTQ files from dump stage for $gse"
    # Check for existence of expected FASTQ files (pattern may vary)
    local fastq_dir="${basedir}/${gse}/fastq"
    if [ ! -d "$fastq_dir" ] || [ -z "$(find "$fastq_dir" -name "*.fastq.gz" -o -name "*.fastq" 2>/dev/null)" ]; then
      log_error "No FASTQ files found for $gse in $fastq_dir"
      return 1
    fi
  done

  for gse in "${gses[@]}"; do
    log_message "Starting GSM renaming for $gse"
    Rscript /home/liuc9/github/scMOCHA-data/src/03-sra-rename-gsm-merge.R -g ${gse} -b ${basedir} &
  done
  wait_for_background_jobs "sra_rename_gsm_merge"
}

# Stage 8: scMOCHA Configuration (Already parallelized, add dependency check)
scmocha_conf() {
  log_message "Starting scMOCHA configuration for ${#gses[@]} GSEs"

  # Verify that GSM renaming completed successfully
  for gse in "${gses[@]}"; do
    local expected_config="${basedir}/${gse}/04.${gse}.batch.sh"
    if [ -f "$expected_config" ]; then
      log_message "Configuration script already exists for $gse"
    fi
  done

  for gse in "${gses[@]}"; do
    log_message "Generating scMOCHA configuration for $gse"
    Rscript /home/liuc9/github/scMOCHA-data/src/04-scmocha-conf.R -g ${gse} -b ${basedir} &
  done
  wait_for_background_jobs "scmocha_conf"

  # Verify configuration files were created
  for gse in "${gses[@]}"; do
    local required_files=(
      "${basedir}/${gse}/04.${gse}.batch.sh"
    )
    if ! check_files_exist "${required_files[@]}"; then
      log_error "Required configuration files missing for $gse"
      return 1
    fi
  done
}

# Stage 9: scMOCHA Batch Execution (Resource-aware parallelization)
scmocha_batch_run() {
  log_message "Starting scMOCHA batch execution for ${#gses[@]} GSEs with max ${MAX_CONCURRENT_SCMOCHA} concurrent"
  local active_jobs=0
  local pids=()

  for gse in "${gses[@]}"; do
    # Wait if we've reached the concurrent limit
    if [ $active_jobs -ge $MAX_CONCURRENT_SCMOCHA ]; then
      log_message "Waiting for scMOCHA slot to become available..."
      wait ${pids[0]}  # Wait for the first job to complete
      local exit_code=$?
      if [ $exit_code -ne 0 ]; then
        log_error "A scMOCHA job failed"
        return 1
      fi
      pids=("${pids[@]:1}")  # Remove first PID from array
      ((active_jobs--))
    fi

    log_message "Starting scMOCHA batch execution for $gse"
    (
      cd ${basedir}/${gse}
      if ! run_with_retry "bash 04.${gse}.batch.sh"; then
        log_error "scMOCHA batch execution failed for $gse"
        exit 1
      fi
      log_message "scMOCHA batch execution completed for $gse"
    ) &

    local pid=$!
    pids+=($pid)
    ((active_jobs++))
  done

  # Wait for all remaining scMOCHA jobs
  log_message "Waiting for all scMOCHA batch jobs to complete..."
  for pid in "${pids[@]}"; do
    wait $pid
    if [ $? -ne 0 ]; then
      log_error "A scMOCHA batch job failed"
      return 1
    fi
  done
  log_message "All scMOCHA batch jobs completed successfully"
}

# Stage 10: Log Parsing (Already parallelized)
parse_log() {
  log_message "Starting log parsing for ${#gses[@]} GSEs"
  for gse in "${gses[@]}"; do
    log_message "Parsing logs for $gse"
    Rscript /home/liuc9/github/scMOCHA-data/src/05-parse-log.R -g ${gse} -b ${basedir} &
  done
  wait_for_background_jobs "parse_log"
}

# Stage 11: Archive Compression (Now parallelized)
cptargz() {
  log_message "Starting archive compression for ${#gses[@]} GSEs"
  for gse in "${gses[@]}"; do
    (
      cd ${basedir}/${gse}
      log_message "Compressing archives for $gse"
      if ! run_with_retry "bash 05.${gse}.scmocha.cptargz.sh"; then
        log_error "Archive compression failed for $gse"
        exit 1
      fi
      log_message "Archive compression completed for $gse"
    ) &
  done
  wait_for_background_jobs "cptargz"
}

# Stage 12: Archive Decompression (Now parallelized)
untargz() {
  log_message "Starting archive decompression for ${#gses[@]} GSEs"
  for gse in "${gses[@]}"; do
    (
      cd ${basedir}/${gse}
      log_message "Decompressing archives for $gse"
      if ! run_with_retry "bash 07.${gse}.scmocha.untargz.sh"; then
        log_error "Archive decompression failed for $gse"
        exit 1
      fi
      log_message "Archive decompression completed for $gse"
    ) &
  done
  wait_for_background_jobs "untargz"
}

# Stage 13: Variant Collection (Already parallelized)
collect_variants() {
  log_message "Starting variant collection for ${#gses[@]} GSEs"
  for gse in "${gses[@]}"; do
    log_message "Collecting variants for $gse"
    Rscript /home/liuc9/github/scMOCHA-data/src/06-collect-variants.R -g ${gse} -b ${basedir} &
  done
  wait_for_background_jobs "collect_variants"
}

# Main execution pipeline with proper error handling
main() {
  log_message "Starting scMOCHA data processing pipeline"
  log_message "Processing GSEs: ${gses[*]}"
  log_message "Base directory: $basedir"
  log_message "Configuration: MAX_CONCURRENT_PREFETCH=$MAX_CONCURRENT_PREFETCH, MAX_CONCURRENT_DUMPS=$MAX_CONCURRENT_DUMPS, MAX_CONCURRENT_SCMOCHA=$MAX_CONCURRENT_SCMOCHA"

  # Execute pipeline stages with proper dependencies
  local stages=(
    "sra_metadata"
    "sra_run_table"
    "sra_download_dump"
    "prefetch_sh"
    "prefetch_check"
    "dump_slrm"
    "sra_rename_gsm_merge"
    "scmocha_conf"
    "scmocha_batch_run"
    "parse_log"
    "cptargz"
    "untargz"
    "collect_variants"
  )

  for stage in "${stages[@]}"; do
    log_message "========================================"
    log_message "Starting stage: $stage"
    log_message "========================================"

    if ! $stage; then
      log_error "Pipeline failed at stage: $stage"
      exit 1
    fi

    log_message "Stage completed successfully: $stage"
  done

  log_message "========================================"
  log_message "scMOCHA data processing pipeline completed successfully!"
  log_message "========================================"
}

# Execute main pipeline if script is run directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  main "$@"
fi