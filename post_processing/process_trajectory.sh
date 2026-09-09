#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

ROOT_DIR="$(dirname "$SCRIPT_DIR")"

TRAJ_DIR="$ROOT_DIR/trajectories/raw"
PROCESSED_TRAJ_DIR="$ROOT_DIR/trajectories/post_processed"
DATA_DIR="$ROOT_DIR/data"

SEQUENCES="${SEQUENCES:-F66 M66 V66 L66 A66 Y66 I66}"
EQUILIBRATION_TIME=8000

mkdir -p "$PROCESSED_TRAJ_DIR"
mkdir -p "$DATA_DIR"

# Make protein whole, wrap, and center
python3 manipulate_pbc.py \
  --raw_traj_path "$TRAJ_DIR" \
  --processed_traj_path "$PROCESSED_TRAJ_DIR" \
  --sequences "$SEQUENCES"

# Calculate the minimum distance between the protein and its periodic image. Write frame numbers to a list where distance is < cutoff
python3 calc_mindist.py \
  --raw_traj_path "$TRAJ_DIR" \
  --input_traj_path "$PROCESSED_TRAJ_DIR" \
  --output_data_path "$DATA_DIR" \
  --sequences "$SEQUENCES" 

# Remove periodic image frames from trajectory
for seq in $SEQUENCES; do
  vmd -dispdev none \
      -e run_PIF_script.tcl \
      -eofexit \
      -args \
      "$TRAJ_DIR/${seq}/${seq}.gro" \
      "$PROCESSED_TRAJ_DIR/${seq}_wrapped_centered.xtc" \
      "$DATA_DIR" \
      "$PROCESSED_TRAJ_DIR" \
      &> "${PROCESSED_TRAJ_DIR}/${seq}_wrapped_centered_PIF_del.log"
    gmx trjconv -f ${PROCESSED_TRAJ_DIR}/${seq}_wrapped_centered_PIF_del.gro -o ${PROCESSED_TRAJ_DIR}/${seq}_wrapped_centered_PIF_del.xtc
done

# Remove equilibration time 
for seq in $SEQUENCES; do
  vmd -dispdev none \
      -e run_equil_script.tcl \
      -eofexit \
      -args \
      "$EQUILIBRATION_TIME" \
      "$TRAJ_DIR/${seq}/${seq}.gro" \
      "$PROCESSED_TRAJ_DIR/${seq}_wrapped_centered_PIF_del.xtc" \
      "$PROCESSED_TRAJ_DIR" \
      &> "${PROCESSED_TRAJ_DIR}/${seq}_wrapped_centered_PIF_equil_del.log" 
  gmx trjconv \
  -f ${PROCESSED_TRAJ_DIR}/${seq}_wrapped_centered_PIF_equil_del.gro \
  -o ${PROCESSED_TRAJ_DIR}/${seq}_wrapped_centered_PIF_equil_del.xtc
done

