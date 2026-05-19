#!/usr/bin/env bash
set -euo pipefail

PROJECT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_NAME="descartes-r45-dada2-1.38-cutadapt-5.2"
LOG_DIR="${PROJECT}/logs/versions"
mkdir -p "${LOG_DIR}" "${PROJECT}/logs/stdout_stderr"

if conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
  echo "Environment already exists: ${ENV_NAME}"
else
  if command -v mamba >/dev/null 2>&1; then
    mamba env create -f "${PROJECT}/config/environment.yml"
  else
    conda env create -f "${PROJECT}/config/environment.yml"
  fi
fi

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${ENV_NAME}"
conda list > "${LOG_DIR}/conda_list_${ENV_NAME}.txt"
R --version > "${LOG_DIR}/R_version_${ENV_NAME}.txt"
cutadapt --version > "${LOG_DIR}/cutadapt_version_${ENV_NAME}.txt"
Rscript -e 'cat(as.character(packageVersion("dada2")), "\n")' > "${LOG_DIR}/dada2_version_${ENV_NAME}.txt"
Rscript -e 'sessioninfo::session_info()' > "${LOG_DIR}/session_info_${ENV_NAME}.txt"
echo "Environment ready: ${ENV_NAME}"
