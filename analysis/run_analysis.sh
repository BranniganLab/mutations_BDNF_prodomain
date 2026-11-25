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
REF_SEQ="V66"

#Blobulation parameters
L_MIN=4 #3 4 5
H_STAR=0.37 #0.33, 0.40, 0.37
HYDROPHOBICITY_SCALE="Kyte-Doolittle"

MUTATION_RESIDUE=66
MEDIATOR_RESIDUE=95

STATE1="VM_MN"
STATE2="MN_NC"
STATE3="VM_NC"

BLOB_CONTACT_CUTOFF=5.5 #5, 6 in Angstroms
RESIDUE_CONTACT_CUTOFF=6 #in Angstroms

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
      "$MUTATION_RESIDUE" \
      "$MEDIATOR_RESIDUE" \
      "$BLOB_CONTACT_CUTOFF" \
      "$STATE1" \
      "$STATE2" \
      "$STATE3" \
      "$RESIDUE_CONTACT_CUTOFF" \
      "$SEQ" \
      "$REF_SEQ"
done
