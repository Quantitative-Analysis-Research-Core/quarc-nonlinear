#!/bin/zsh
# Regenerate every committed characterization report, then the figures and
# the lab report. Restartable: a finished stage leaves a .done marker and is
# skipped on the next run; the long stages also resume inside themselves
# (see tests/matlab/rerun_reports.m).
#
#   nohup caffeinate -i tests/rerun_reports.sh > tests/artifacts/rerun/run.log 2>&1 &
#
# Stops at the first failing stage, so a fault is never papered over.

set -u
repo=${0:A:h:h}
cd "$repo"
matlab=/Applications/MATLAB_R2025b.app/bin/matlab
python=${QUARC_PYTHON:-$HOME/Projects/quarc-nonlinear/.venv/bin/python}
quarto=$HOME/.local/quarto/bin/quarto
marks=tests/artifacts/rerun
mkdir -p "$marks"

stamp() { date '+%Y-%m-%d %H:%M:%S'; }

for stage in fast curves evolve full dysts noise noiserqa; do
  if [[ -f $marks/$stage.done ]]; then
    echo "$(stamp) skip $stage (done)"
    continue
  fi
  echo "$(stamp) start $stage"
  if ! $matlab -batch "addpath('tests/matlab'); rerun_reports('$stage')" \
        > "$marks/$stage.log" 2>&1; then
    echo "$(stamp) FAILED $stage -- see $marks/$stage.log"
    exit 1
  fi
  touch "$marks/$stage.done"
  echo "$(stamp) done $stage"
done

echo "$(stamp) start post-processing"
python3 python/merge_noise_sweep_chunks.py || { echo "$(stamp) FAILED merge"; exit 1; }
$python python/rescore_ratios.py || { echo "$(stamp) FAILED rescore"; exit 1; }
Rscript R/make_figures.R || { echo "$(stamp) FAILED figures"; exit 1; }
$quarto render docs/lab-report/lyapunov_estimators_lab_report.qmd || { echo "$(stamp) FAILED render"; exit 1; }
echo "$(stamp) ALL DONE"
