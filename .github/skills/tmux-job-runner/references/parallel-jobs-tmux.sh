#!/usr/bin/env bash
# Author: Chunjie Liu
# Contact: chunjie.sam.liu.at.gmail.com
# Date: {TODAY}
# Description: Launch parallel tmux sessions for {SCRIPT_DESC} (one per {DIMENSION})
# Version: 0.1

set -euo pipefail

# --- Configurable variables (override via environment) ---
REPO_DIR="${REPO_DIR:-/home/liuc9/github/scMOCHA-data}"
SCRIPT_PATH="${SCRIPT_PATH:-${REPO_DIR}/high-res/{NN-script-name.R}}"
RSCRIPT="${RSCRIPT:-/scr1/users/liuc9/tools/miniforge3/envs/renv/bin/Rscript}"
CONDA_SH="${CONDA_SH:-/scr1/users/liuc9/tools/miniforge3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-renv}"
SESSION_PREFIX="${SESSION_PREFIX:-{PREFIX}}"
NTHREAD="${NTHREAD:-1}"
# {ADD SCRIPT-SPECIFIC ARGS HERE}
CELLTYPES="${CELLTYPES:-}"
VERBOSE="${VERBOSE:-true}"

# --- Derived paths ---
OUTPUT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TIMESTAMP="$(date '+%Y%m%d_%H%M%S')"

# --- Parallel dimension: one tmux session per item ---
# {CUSTOMIZE: change array name and values to match your parallel dimension}
LEVELS_LIST=("L1" "L2")

for LEVEL in "${LEVELS_LIST[@]}"; do
  SESSION_NAME="${SESSION_PREFIX}_${LEVEL}_${TIMESTAMP}"
  LOG_DIR="${OUTPUT_DIR}/logs-launch_${TIMESTAMP}"
  LOG_FILE="${LOG_DIR}/${SESSION_NAME}.log"
  STATUS_FILE="${LOG_DIR}/${SESSION_NAME}.status"
  META_FILE="${LOG_DIR}/${SESSION_NAME}.meta"
  JOB_SCRIPT="${LOG_DIR}/${SESSION_NAME}.sh"

  mkdir -p "${LOG_DIR}"

  # Skip if session already exists (parallel-safe)
  if tmux has-session -t "${SESSION_NAME}" 2>/dev/null; then
    echo "tmux session already exists: ${SESSION_NAME} - skipping" >&2
    continue
  fi

  # --- Build optional CLI args ---
  verbose_arg=""
  if [[ "${VERBOSE}" == "true" ]]; then
    verbose_arg="--verbose"
  fi

  celltypes_arg=""
  if [[ -n "${CELLTYPES}" ]]; then
    printf -v celltypes_arg ' --celltypes %q' "${CELLTYPES}"
  fi

  # --- Quote helper ---
  escape_sq() {
    printf "%s" "$1" | sed "s/'/'\\\\''/g"
  }

  repo_dir_sq="$(escape_sq "${REPO_DIR}")"
  script_path_sq="$(escape_sq "${SCRIPT_PATH}")"
  rscript_sq="$(escape_sq "${RSCRIPT}")"
  conda_sh_sq="$(escape_sq "${CONDA_SH}")"
  conda_env_sq="$(escape_sq "${CONDA_ENV}")"
  level_sq="$(escape_sq "${LEVEL}")"
  log_file_sq="$(escape_sq "${LOG_FILE}")"
  status_file_sq="$(escape_sq "${STATUS_FILE}")"
  output_dir_sq="$(escape_sq "${OUTPUT_DIR}")"
  session_name_sq="$(escape_sq "${SESSION_NAME}")"
  verbose_arg_sq="$(escape_sq "${verbose_arg}")"

  # --- Write metadata ---
  cat > "${META_FILE}" <<EOF
session_name=${SESSION_NAME}
log_file=${LOG_FILE}
status_file=${STATUS_FILE}
job_script=${JOB_SCRIPT}
repo_dir=${REPO_DIR}
script_path=${SCRIPT_PATH}
level=${LEVEL}
celltypes=${CELLTYPES}
created_at=$(date '+%F %T')
EOF

  # --- Generate job script ---
  cat > "${JOB_SCRIPT}" <<EOF
#!/usr/bin/env bash
set -euo pipefail

on_exit() {
  local status=\$?
  echo "[\$(date '+%F %T')] session='${session_name_sq}' finish status=\${status}"
  printf "%s\n" "\${status}" > '${status_file_sq}'
  exit \${status}
}

trap on_exit EXIT
exec > >(tee -a '${log_file_sq}') 2>&1

echo "[\$(date '+%F %T')] session='${session_name_sq}' start"
echo "[\$(date '+%F %T')] output_dir='${output_dir_sq}'"
echo "[\$(date '+%F %T')] level='${level_sq}'"

set +u
source '${conda_sh_sq}'
conda activate '${conda_env_sq}'
set -u

cd '${repo_dir_sq}'

'${rscript_sq}' '${script_path_sq}' \\
  --levels '${level_sq}' \\
  --nthread '${NTHREAD}' \\
  ${verbose_arg_sq}${celltypes_arg}
EOF

  chmod +x "${JOB_SCRIPT}"
  tmux new-session -d -s "${SESSION_NAME}" "${JOB_SCRIPT}"

  echo "Started tmux session: ${SESSION_NAME}"
  echo "  Log file: ${LOG_FILE}"
  echo "  Status file: ${STATUS_FILE}"
  echo ""
done

echo "All parallel sessions launched. Check with: tmux ls | grep ${SESSION_PREFIX}"
