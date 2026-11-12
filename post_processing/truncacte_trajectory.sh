#!/bin/bash
set -e

# This command truncates the trajectory of each BDNF prodomain variant to 2 us for consistent comparison. The trajectories are  slightly over 2 us.
# This command will not be required to be run by users. The truncated trajectory will be provided in the Zenodo
# Run: bash truncate_trajectory.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

ROOT_DIR="$(dirname "$SCRIPT_DIR")"

REPLICA="9"
TRAJ_DIR="$ROOT_DIR/trajectories/replica_$REPLICA/raw"

for seq in F M V L A Y I 
do
	echo "Protein" | gmx trjconv \
	-s ${TRAJ_DIR}/${seq}66/${seq}66.tpr \
	-f ${TRAJ_DIR}/${seq}66/${seq}66_full.xtc \
	-o ${TRAJ_DIR}/${seq}66/${seq}66.xtc \
	-e 2000000
done
