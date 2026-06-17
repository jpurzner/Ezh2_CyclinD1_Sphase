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

# --- SLOW ensemble / analysis figures (explicit; not caught by the fig_v44_* glob) ---
# Each entry is "script [args]"; cached scripts get --fresh so they actually recompute
# against the changed model (without it they would silently reuse a stale .npz cache).
SLOW=(
  "simulations/sim_g0_bifurcation.py"            # fig_v44_g0_bifurcation -- cyclin D1 / birth-p27 bifurcation (immediate vs transient-G0); ~10 min, N=140
  "simulations/sim_ezh2_phaseplane.py --fresh"   # fig_v44_ezh2_phaseplane -- EZH2-CyclinD1 nullcline portrait + bifurcation diagrams; ~35 min (cached)
  "simulations/sim_mitogen_withdrawal.py --fresh" # fig_v44_mitogen_withdrawal -- sudden/gradual mitogen withdrawal x depth x +/-feedback, post-withdrawal divisions; ~2.5 h (cached)
  "simulations/sim_vismo_withdrawal.py --fresh"   # fig_v44_vismo_withdrawal -- MB+vismodegib (partial-withdrawal analog) x dose x +/-feedback; EZH2i rescue dose/heterogeneity-resolved; ~2.5 h (cached)
)

ok=0; fail=0
run() {
  local script="$1"
  if $PY "$@" > "/tmp/regenfig_$(basename "$script").log" 2>&1; then
    echo "OK   $(basename "$script")"; ok=$((ok+1))
  else
    echo "FAIL $(basename "$script")  (see /tmp/regenfig_$(basename "$script").log)"; fail=$((fail+1))
  fi
}

echo "== FAST: per-condition fig_v44_*.py =="
for f in simulations/fig_v44_*.py; do run "$f"; done

echo "== SLOW: ensemble / analysis figures =="
for entry in "${SLOW[@]}"; do run $entry; done   # unquoted -> split script + args

echo "------------------------------------------------------------"
echo "regenerated: $ok ok, $fail failed"
[ "$fail" -eq 0 ]
