#!/usr/bin/env bash
set -euo pipefail
cd "$HOME/sites/FEMBasicsLectureNotes"

export JULIA_NUM_THREADS=auto

exec $HOME/.juliaup/bin/julia --project="pluto-slider-server-environment" -e '
import Pkg; Pkg.instantiate(); 
import PlutoSliderServer; 
PlutoSliderServer.run_git_directory(".";)
'