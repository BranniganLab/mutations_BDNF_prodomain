# Mutations of the BDNF prodomain paper

This repository contains scripts for trajectory post-processing, analysis, and plotting scripts for the figures in the BDNF prodomain paper. Raw trajectory files can be found in [this Zenodo](link) repository.

## Requirements

- VMD
- GROMACS
- Python 3.11
- Blobulator
- Conda or Miniconda

The Python dependencies are listed in [`environment.yml`](environment.yml).

## Installation

Clone this repository:

```bash
git clone git@github.com:BranniganLab/mutations_BDNF_prodomain.git
cd mutations_BDNF_prodomain
```

Create and activate the Conda environment:

```bash
conda env create -f environment.yml
conda activate mutations-bdnf-prodomain
```

VMD and GROMACS must be installed separately and should be available from the command line.

## Blobulator dependency

The Blobulator is maintained in a separate repository and is required for the analysis workflow.

Clone the Blobulator repository separately:

```bash
cd ..
git clone git@github.com:BranniganLab/blobulator.git
```

The recommended directory structure is:

```text
GitHub/
├── mutations_BDNF_prodomain/
└── blobulator/
```

The analysis scripts expect the Blobulator VMD scripts to be located at:

```text
blobulator/VMD_scripts/
```

If the Blobulator repository is located elsewhere, set its path before running the analysis:

```bash
export BLOBULATOR_DIR="/path/to/blobulator"
```

## Repository organization

```text
analysis/
    VMD and Tcl analysis scripts

plotting/
    Python scripts for generating figures

post_processing/
    Trajectory processing and minimum-distance calculations

environment.yml
    Conda environment specification

generate_figs.sh
    Main script that runs the trajectory processing, analysis, and plotting
```

## Workflow

### 1. Download trajectories

Download the trajectory files from [this Zenodo](link) and place the folder into the mutations_BDNF_prodomain folder.

### 2. Generate figures

From the mutations_BDNF_prodomain directory, run:

```bash
./generate_figs.sh
```

First process_trajectory.sh is called, which performs trajectory preprocessing, periodic-boundary handling, minimum-distance analysis, PIF filtering, and equilibration-frame removal.

Then, run_analysis.sh is called, which performs the VMD/Tcl calculations, including:
- Radius of gyration
- Radius of gyration of Blobulator blobs
- Blob midpoint calculations
- Blob contact calculations
- Residue side-chain contacts
- State classification
- Odds-ratio calculations

Lastly, run_plotting.sh is called, which reads the generated data files and save figures in the configured output directory.

## Configuration

Important analysis values are defined near the top of the shell scripts, including:

- Sequence names
- Equilibration time
- Blob indices
- Contact cutoffs
- Residue ranges
- State names
- Input and output directories

Review these values before running the workflow.

The sequence list can be overridden without editing the script:

```bash
SEQUENCES="F66 M66 V66 L66 A66 Y66 I66" bash process_trajectory.sh
```

Similarly, the Blobulator location can be overridden with:

```bash
BLOBULATOR_DIR="/path/to/blobulator" bash run_analysis.sh
```

## Note
The trajectory files and external Blobulator repository are not included in this repository. Anyone reproducing the analysis must obtain those files separately.