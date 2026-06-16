#!/usr/bin/env bash
# Canonical figure-regeneration manifest for the v44 model.
# Run after ANY change to src/build_model_v44_heldt.py or the calibration so every
# figure reflects the current model. Writes .png (preview) + .pdf (paper) into simulations/.
#
#   ./venv/bin/... not needed -- this script calls the venv python directly.
#   bash simulations/regenerate_figures.sh
#
# Two tiers:
#   FAST  -- every simulations/fig_v44_*.py (per-condition; auto-globbed, so new
#            fig_v44_* scripts are picked up with no edit here).
#   SLOW  -- ensemble figures that don't match the fig_v44_* glob and take minutes
#            (N-cell heterogeneity draws). Listed explicitly below.
set -u
cd "$(dirname "$0")/.." || exit 1
PY=./venv/bin/python

# --- SLOW ensemble figures (explicit; not caught by the fig_v44_* glob) ---
SLOW=(
  simulations/sim_g0_bifurcation.py        # fig_v44_g0_bifurcation -- cyclin D1 / birth-p27 bifurcation (immediate vs transient-G0); ~10 min, N=140
)

ok=0; fail=0
run() {
  if $PY "$1" > "/tmp/regenfig_$(basename "$1").log" 2>&1; then
    echo "OK   $(basename "$1")"; ok=$((ok+1))
  else
    echo "FAIL $(basename "$1")  (see /tmp/regenfig_$(basename "$1").log)"; fail=$((fail+1))
  fi
}

echo "== FAST: per-condition fig_v44_*.py =="
for f in simulations/fig_v44_*.py; do run "$f"; done

echo "== SLOW: ensemble figures =="
for f in "${SLOW[@]}"; do run "$f"; done

echo "------------------------------------------------------------"
echo "regenerated: $ok ok, $fail failed"
[ "$fail" -eq 0 ]
