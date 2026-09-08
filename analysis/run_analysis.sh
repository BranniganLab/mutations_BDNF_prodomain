#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

ROOT_DIR="$(dirname "$SCRIPT_DIR")"

TRAJ_DIR="$ROOT_DIR/trajectories"
PROCESSED_TRAJ_DIR="$ROOT_DIR/trajectories/post_processed11"
DATA_DIR="$ROOT_DIR/data11"
BLOBULATOR_TRAJ=""$ROOT_DIR/blobulator/VMD_scripts"

SEQUENCES="${SEQUENCES:-F66 M66 V66 L66 A66 Y66 I66}"
L_MIN=4
H_STAR=0.37
HYDROPHOBICITY_SCALE="Kyte-Doolittle"
RESID_TO_REASSIGN=65
RESID_NEW_USER1=1
RESID_NEW_USER2=8.0
N_TERMINAL_BLOB_INDEX=2
C_TERMINAL_BLOB_INDEX=16
BLOB_CONTACT_CUTOFF=5.5
VARIANT_BLOB_INDEX=8
MEDIATOR_BLOB_INDEX=12
STATE1="VM_MN"
STATE2="MN_NC"
RESIDUE_CONTACT_CUTOFF=6

for seq in $SEQUENCES; do
  vmd -dispdev none \
      -e run_analysis.tcl \
      -eofexit \
      -args \
      "$TRAJ_DIR/raw/${seq}/${seq}.gro" \
      "$PROCESSED_TRAJ_DIR/${seq}_wrapped_centered_PIF_del.xtc" \
      "$PROCESSED_TRAJ_DIR/${seq}_wrapped_centered_PIF_equil_del.xtc" \
      "$DATA_DIR" \
      "$L_MIN" \
      "$H_STAR" \
      "$HYDROPHOBICITY_SCALE" \
      "$RESID_TO_REASSIGN" \
      "$RESID_NEW_USER1" \
      "$RESID_NEW_USER2" \
      "$N_TERMINAL_BLOB_INDEX" \
      "$C_TERMINAL_BLOB_INDEX" \
      "$BLOB_CONTACT_CUTOFF" \
      "$VARIANT_BLOB_INDEX" \
      "$MEDIATOR_BLOB_INDEX" \
      "$STATE1" \
      "$STATE2" \
      "$RESIDUE_CONTACT_CUTOFF" 
done