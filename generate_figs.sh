#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# echo "Processing trajectories..."
# bash "$SCRIPT_DIR/post_processing/process_trajectory.sh"

# echo "Running analysis..."
# bash "$SCRIPT_DIR/analysis/run_analysis.sh"

echo "Generating figures..."
bash "$SCRIPT_DIR/plotting/run_plotting.sh"

echo "Done."