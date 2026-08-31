#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

ROOT_DIR="$(dirname "$SCRIPT_DIR")"

TRAJ_DIR="$ROOT_DIR/trajectories/raw"
PROCESSED_TRAJ_DIR="$ROOT_DIR/trajectories/post_processed5"
DATA_DIR="$ROOT_DIR/data5"

mkdir -p "$PROCESSED_TRAJ_DIR"
mkdir -p "$DATA_DIR"

# python3 calc_mindist.py --input_traj_path "$TRAJ_DIR" --output_data_path "$DATA_DIR"

for seq in F66 M66 V66 L66 A66 Y66 I66; do
  vmd -dispdev none \
      -e run_PIF_script.tcl \
      -eofexit \
      -args \
      "$TRAJ_DIR/${seq}/${seq}_resid_23-113-capped.gro" \
      "$TRAJ_DIR/${seq}/${seq}_trajout_cut.xtc" \
      "$DATA_DIR" \
      "$PROCESSED_TRAJ_DIR" \
      &> "${PROCESSED_TRAJ_DIR}/${seq}_PIF_del_full_traj.log"
    gmx trjconv -f ${PROCESSED_TRAJ_DIR}/${seq}_PIF_del_full_traj.gro -o ${PROCESSED_TRAJ_DIR}/${seq}_PIF_del_full_traj.xtc
done

python3 manipulate_pbc.py --raw_traj_path "$TRAJ_DIR" --processed_traj_path "$PROCESSED_TRAJ_DIR"

for seq in F66 M66 V66 L66 A66 Y66 I66; do
  vmd -dispdev none \
      -e run_equil_script.tcl \
      -eofexit \
      -args \
      "$TRAJ_DIR/${seq}/${seq}_resid_23-113-capped.gro" \
      "$PROCESSED_TRAJ_DIR/${seq}_PIF_del_full_traj_wrapped_centered.xtc" \
      "$PROCESSED_TRAJ_DIR" \
      &> "${PROCESSED_TRAJ_DIR}/${seq}_wrapped_PIF_equil_del.log" 
  gmx trjconv \
  -f ${PROCESSED_TRAJ_DIR}/${seq}_wrapped_PIF_equil_del.gro \
  -o ${PROCESSED_TRAJ_DIR}/${seq}_wrapped_PIF_equil_del.xtc
done

