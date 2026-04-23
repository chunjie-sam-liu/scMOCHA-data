#!/usr/bin/env bash
# Author: Chunjie Liu
# Contact: chunjie.sam.liu.at.gmail.com
# Date: {TODAY}
# Description: Launch {SCRIPT_DESC} job in a detached tmux session with logging
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
# {ADD SCRIPT-SPECIFIC ARGS HERE, e.g.:}
# LEVELS="${LEVELS:-L1,L2}"
# MIN_CELLS="${MIN_CELLS:-15}"
CELLTYPES="${CELLTYPES:-}"
VERBOSE="${VERBOSE:-true}"

# --- Derived paths (do not change) ---
OUTPUT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TIMESTAMP="$(date '+%Y%m%d_%H%M%S')"
SESSION_NAME="${SESSION_PREFIX}_${TIMESTAMP}"
LOG_DIR="${LOG_DIR:-${OUTPUT_DIR}/logs-launch_${TIMESTAMP}}"
LOG_FILE="${LOG_FILE:-${LOG_DIR}/${SESSION_NAME}.log}"
STATUS_FILE="${STATUS_FILE:-${LOG_DIR}/${SESSION_NAME}.status}"
META_FILE="${META_FILE:-${LOG_DIR}/${SESSION_NAME}.meta}"
JOB_SCRIPT="${JOB_SCRIPT:-${LOG_DIR}/${SESSION_NAME}.sh}"

mkdir -p "${LOG_DIR}"

# --- Preflight checks ---
if [[ ! -f "${SCRIPT_PATH}" ]]; then
  echo "Missing R script: ${SCRIPT_PATH}" >&2
  exit 1
fi

if [[ ! -f "${CONDA_SH}" ]]; then
  echo "Missing conda setup script: ${CONDA_SH}" >&2
  exit 1
fi

if ! command -v tmux >/dev/null 2>&1; then
  echo "tmux not found in PATH" >&2
  exit 1
fi

if tmux has-session -t "${SESSION_NAME}" 2>/dev/null; then
  echo "tmux session already exists: ${SESSION_NAME}" >&2
  exit 1
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

# --- Quote helper for heredoc embedding ---
escape_sq() {
  printf "%s" "$1" | sed "s/'/'\\\\''/g"
}

repo_dir_sq="$(escape_sq "${REPO_DIR}")"
script_path_sq="$(escape_sq "${SCRIPT_PATH}")"
rscript_sq="$(escape_sq "${RSCRIPT}")"
conda_sh_sq="$(escape_sq "${CONDA_SH}")"
conda_env_sq="$(escape_sq "${CONDA_ENV}")"
log_file_sq="$(escape_sq "${LOG_FILE}")"
status_file_sq="$(escape_sq "${STATUS_FILE}")"
output_dir_sq="$(escape_sq "${OUTPUT_DIR}")"
session_name_sq="$(escape_sq "${SESSION_NAME}")"
verbose_arg_sq="$(escape_sq "${verbose_arg}")"
# {ADD escape_sq CALLS FOR SCRIPT-SPECIFIC ARGS}

# --- Write metadata ---
cat > "${META_FILE}" <<EOF
session_name=${SESSION_NAME}
log_file=${LOG_FILE}
status_file=${STATUS_FILE}
job_script=${JOB_SCRIPT}
repo_dir=${REPO_DIR}
script_path=${SCRIPT_PATH}
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

set +u
source '${conda_sh_sq}'
conda activate '${conda_env_sq}'
set -u

cd '${repo_dir_sq}'

'${rscript_sq}' '${script_path_sq}' \\
  --nthread '${NTHREAD}' \\
  ${verbose_arg_sq}${celltypes_arg}
EOF

chmod +x "${JOB_SCRIPT}"

# --- Launch ---
tmux new-session -d -s "${SESSION_NAME}" "${JOB_SCRIPT}"

echo "Started tmux session: ${SESSION_NAME}"
echo "Log file: ${LOG_FILE}"
echo "Status file: ${STATUS_FILE}"
echo "Metadata file: ${META_FILE}"
