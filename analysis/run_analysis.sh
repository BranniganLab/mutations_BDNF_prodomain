#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

ROOT_DIR="$(dirname "$SCRIPT_DIR")"

TRAJ_DIR="$ROOT_DIR/trajectories"
PROCESSED_TRAJ_DIR="$ROOT_DIR/trajectories/post_processed6"
DATA_DIR="$ROOT_DIR/data6"

for seq in F66 M66 V66 L66 A66 Y66 I66; do
  vmd -dispdev none \
      -e run_analysis.tcl \
      -eofexit \
      -args \
      "$TRAJ_DIR/raw/${seq}/${seq}_resid_23-113-capped.gro" \
      "$PROCESSED_TRAJ_DIR/${seq}_PIF_del_full_traj_wrapped_centered.xtc" \
      "$PROCESSED_TRAJ_DIR/${seq}_wrapped_PIF_equil_del.xtc" \
      "$DATA_DIR"
done