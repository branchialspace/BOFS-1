#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BOFS1_DIR="$SCRIPT_DIR/../.."
PSLIBRARY_DIR="$BOFS1_DIR/pslibrary"
FUNCTIONAL="rel-pbe"

# ============================================================================
# GLOBAL LD1.X PARAMETERS - Edit these to inject into all pseudopotentials
# ============================================================================

INPUTP_PARAMS=(
    "lsave_wfc=.true."
    "lgipaw_reconstruction=.true."
    "use_paw_as_gipaw=.true."
)

INPUT_PARAMS=(
    # Add &input parameters here if needed
)

# ============================================================================

inject_params() {
    local file="$1"
    [ ! -f "$file" ] && return

    perl -i -pe '
        BEGIN {
            $in_inputp = 0;
            $in_input = 0;
            @inputp_params = @ARGV[1..$#ARGV];
            $#ARGV = 0;

            # Extract parameter names (everything before =)
            @param_names = ();
            foreach my $param (@inputp_params) {
                if ($param =~ /^([^=]+)=/) {
                    push @param_names, $1;
                }
            }
        }

        if (/^\s*&inputp/) {
            $in_inputp = 1;
            $in_input = 0;
        }
        elsif (/^\s*&input/ && !$in_inputp) {
            $in_input = 1;
        }
        elsif ($in_inputp && /^\s*\/\s*$/) {
            # Inject all parameters before closing /
            foreach my $param (@inputp_params) {
                print "   $param,\n";
            }
            $in_inputp = 0;
        }
        elsif ($in_input && /^\s*\/\s*$/) {
            $in_input = 0;
        }
        elsif ($in_inputp) {
            # Check if this line matches any of our parameter names
            my $skip = 0;
            foreach my $name (@param_names) {
                if (/^\s*$name\s*=/) {
                    $skip = 1;
                    last;
                }
            }
            # Skip this line if it matches one of our parameters
            $_ = "" if $skip;
        }
    ' "$file" "${INPUTP_PARAMS[@]}"
}

# Inject parameters into job files
cd "$PSLIBRARY_DIR"

for f in paw_ps_*.job us_ps_*.job; do
    [ -f "$f" ] && inject_params "$f"
done

# Run PSLibrary
cd "$FUNCTIONAL"
source "$BOFS1_DIR/miniforge3/etc/profile.d/conda.sh"
conda activate "$BOFS1_DIR/bofs1_env"
. ../make_ps
