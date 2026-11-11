#!/usr/bin/env bash
# Build the `sempgs` environment (Python + R) for this repository.
#
#   ./setup_env.sh            create (or update) the env, then build rpf and NPSOL-enabled OpenMx in it
#   ./setup_env.sh --check    only run check_env.py against an existing env
#
# Needs micromamba (brew install micromamba) or conda/mamba on PATH. Takes ~10 min, mostly compiling OpenMx.
# Afterwards:  micromamba activate sempgs   (or: micromamba run -n sempgs <cmd>)
set -euo pipefail
cd "$(dirname "$0")"
ENV=sempgs

if command -v micromamba >/dev/null; then MM=micromamba; export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/mamba}"
elif command -v mamba >/dev/null;    then MM=mamba
elif command -v conda >/dev/null;    then MM=conda
else echo "Install micromamba first: brew install micromamba" >&2; exit 1; fi
run() { "$MM" run -n "$ENV" "$@"; }

if [[ "${1:-}" != "--check" ]]; then
  if "$MM" env list | awk '{print $1}' | grep -qx "$ENV"; then
    echo "== updating env $ENV"; "$MM" env update -y -n "$ENV" -f environment.yml 2>/dev/null || "$MM" update -y -n "$ENV" -f environment.yml
  else
    echo "== creating env $ENV"; "$MM" create -y -f environment.yml
  fi

  # A previously interrupted R install leaves a lock dir that makes the next install fail.
  rm -rf "$(run Rscript -e 'cat(.libPaths()[1])')"/00LOCK-*

  echo "== building rpf (CRAN source) and NPSOL-enabled OpenMx (getOpenMx.R from the OpenMx team)"
  run Rscript -e '
    options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = parallel::detectCores())
    Sys.setenv(MAKEFLAGS = paste0("-j", parallel::detectCores()))
    if (!requireNamespace("rpf", quietly = TRUE)) install.packages("rpf", type = "source")
    if (!requireNamespace("OpenMx", quietly = TRUE) || !OpenMx::imxHasNPSOL())
      source("https://vipbg.vcu.edu/vipbg/OpenMx2/software/getOpenMx.R")
  '
fi

echo "== checking env"
run python check_env.py
