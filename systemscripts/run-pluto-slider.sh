#!/usr/bin/env bash
set -euo pipefail
cd "$HOME/sites/FEMBasicsLectureNotes"

# Choose a data dir for exports & logs
export JULIA_NUM_THREADS=auto

# Note: tweak arguments to match your repo & needs
exec julia--project="pluto-slider-server-environment" -e '
import Pkg; Pkg.instantiate(); 
import PlutoSliderServer; 
PlutoSliderServer.run_git_directory(\".\";
)'