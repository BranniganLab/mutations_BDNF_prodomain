#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

ROOT_DIR="$(dirname "$SCRIPT_DIR")"

REPLICA="0" #1, 9
TRAJ_DIR="$ROOT_DIR/trajectories/replica_$REPLICA/raw"
PROCESSED_TRAJ_DIR="$ROOT_DIR/trajectories/replica_$REPLICA/post_processed"
DATA_DIR="$ROOT_DIR/data/replica_$REPLICA"
BLOBULATOR_DIR="${BLOBULATOR_DIR:-$ROOT_DIR/../blobulator}"
BLOB_GUI_DIR="$BLOBULATOR_DIR/VMD_scripts"

SEQUENCES="${SEQUENCES:-V66 F66 M66 L66 A66 Y66 I66}"
L_MIN=5 #4 5
H_STAR=0.33 #0.33, 0.40, 0.37
HYDROPHOBICITY_SCALE="Kyte-Doolittle"
N_TERMINAL_BLOB_INDEX=2
C_TERMINAL_BLOB_INDEX=12 #12, 15, 16
BLOB_CONTACT_CUTOFF=5.5 #5, 6 - in Angstroms
VARIANT_BLOB_INDEX=6 #6, 8
MEDIATOR_BLOB_INDEX=8 #8, 10 12
STATE1="VM_MN"
STATE2="MN_NC"
STATE3="VM_NC"
RESIDUE_CONTACT_CUTOFF=6 #in Angstroms
REF_SEQ="V66"

for SEQ in $SEQUENCES; do
  vmd -dispdev none \
      -e run_analysis.tcl \
      -eofexit \
      -args \
      "$TRAJ_DIR/${SEQ}/${SEQ}.gro" \
      "$PROCESSED_TRAJ_DIR/${SEQ}_wrapped_centered_PIF_del.xtc" \
      "$PROCESSED_TRAJ_DIR/${SEQ}_wrapped_centered_PIF_equil_del.xtc" \
      "$DATA_DIR" \
      "$BLOB_GUI_DIR" \
      "$L_MIN" \
      "$H_STAR" \
      "$HYDROPHOBICITY_SCALE" \
      "$N_TERMINAL_BLOB_INDEX" \
      "$C_TERMINAL_BLOB_INDEX" \
      "$BLOB_CONTACT_CUTOFF" \
      "$VARIANT_BLOB_INDEX" \
      "$MEDIATOR_BLOB_INDEX" \
      "$STATE1" \
      "$STATE2" \
      "$STATE3" \
      "$RESIDUE_CONTACT_CUTOFF" \
      "$SEQ" \
      "$REF_SEQ"
done