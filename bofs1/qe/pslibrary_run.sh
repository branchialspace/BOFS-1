#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BOFS1_DIR="$SCRIPT_DIR/../.."
PSLIBRARY_DIR="$BOFS1_DIR/pslibrary"
FUNCTIONAL="rel-pbesol"

# ============================================================================
# GLOBAL LD1.X PARAMETERS - Edit these to inject into all pseudopotentials
# ============================================================================

INPUTP_PARAMS=(
    "lsave_wfc=.true."
    "lgipaw_reconstruction=.true."
    "use_paw_as_gipaw=.true."
)

INPUT_PARAMS=(

)

# ============================================================================

inject_params() {
    local file="$1"
    [ ! -f "$file" ] && return
    
    awk -v np="${#INPUTP_PARAMS[@]}" -v ni="${#INPUT_PARAMS[@]}" '
    BEGIN {in_inputp=0; in_input=0; done_inputp=0; done_input=0}
    /^[[:space:]]*&inputp/ {in_inputp=1; print; next}
    /^[[:space:]]*&input/ {in_input=1; print; next}
    in_inputp && /^[[:space:]]*\// {
        if (!done_inputp && np>0) {
            for(i=0; i<np; i++) print "   " ARGV[ARGC-np+i] ","
            done_inputp=1
        }
        in_inputp=0; print; next
    }
    in_input && /^[[:space:]]*\// {
        if (!done_input && ni>0) {
            for(i=0; i<ni; i++) print "   " ARGV[ARGC-np-ni+i] ","
            done_input=1
        }
        in_input=0; print; next
    }
    {print}
    ' "$file" "${INPUTP_PARAMS[@]}" "${INPUT_PARAMS[@]}" > "$file.tmp" && mv "$file.tmp" "$file"
}

# Inject parameters
cd "$PSLIBRARY_DIR/$FUNCTIONAL"
for f in paw_ps_*.job us_ps_*.job; do
    inject_params "$f"
done

# Run PSLibrary
source "$BOFS1_DIR/miniforge3/etc/profile.d/conda.sh"
conda activate "$BOFS1_DIR/bofs1_env"
. ../make_ps
