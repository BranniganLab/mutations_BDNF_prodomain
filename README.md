# Mutations of the BDNF Prodomain

This repository contains scripts for trajectory post-processing, Blobulator analysis, contact analysis, odds-ratio calculations, and figure generation for the BDNF prodomain mutation study.

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
git clone <REPOSITORY_URL>
cd mutations_BDNF_prodomain
```

Create and activate the Conda environment:

```bash
conda env create -f environment.yml
conda activate mutations-bdnf-prodomain
```

VMD and GROMACS must be installed separately and available from the command line.

## Blobulator dependency

Blobulator is maintained in a separate repository and is required for the analysis workflow.

Clone Blobulator separately:

```bash
cd ..
git clone <BLOBULATOR_REPOSITORY_URL> blobulator
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

If Blobulator is located elsewhere, set its path before running the analysis:

```bash
export BLOBULATOR_DIR="/path/to/blobulator"
```

Record the Blobulator commit or version used for the analysis so that the results can be reproduced:

```bash
cd /path/to/blobulator
git rev-parse HEAD
```

## Repository organization

```text
analysis/
    VMD and Tcl analysis scripts

post_processing/
    Trajectory processing and minimum-distance calculations

plotting/
    Python scripts for generating figures

data/
    Input and output analysis data

trajectories/
    Input and processed trajectories

environment.yml
    Conda environment specification
```

## Workflow

### 1. Prepare trajectories

Place the required trajectory and topology files in the expected directories under `trajectories/`.

Run:

```bash
cd post_processing
bash process_trajectory.sh
```

This step performs trajectory preprocessing, periodic-boundary handling, minimum-distance analysis, PIF filtering, and equilibration-frame removal.

### 2. Run trajectory analysis

From the analysis directory, run:

```bash
cd ../analysis
bash run_analysis.sh
```

The analysis script runs the VMD/Tcl calculations, including:

- Radius of gyration
- Radius of gyration of Blobulator blobs
- Blob midpoint calculations
- Blob contact calculations
- Residue side-chain contacts
- State classification
- Odds-ratio calculations

### 3. Generate figures

Run:

```bash
cd ../plotting
bash run_plotting.sh
```

The plotting scripts read the generated data files and save figures in the configured output directory.

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

## Reproducibility

For each analysis, record:

- Git commit of this repository
- Git commit or version of Blobulator
- GROMACS version
- VMD version
- Python version
- Operating system
- Values of analysis parameters and cutoffs
- Input trajectory and topology files

The trajectory files and external Blobulator repository are not necessarily included in this repository. Anyone reproducing the analysis must obtain those files separately.