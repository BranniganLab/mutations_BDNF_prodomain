# Mutations of the BDNF prodomain

This repository contains scripts for trajectory post-processing, analysis, and plotting scripts for the figures in the BDNF prodomain paper. Raw trajectory files can be found in [this Zenodo](https://zenodo.org/records/22710682).

## Requirements

- VMD
- GROMACS
- Python 3.11
- Blobulator v1.1.0
- Conda or Miniconda

Note: VMD and GROMACS must be installed separately and should be available from the command line.

The Python dependencies are listed in [`environment.yml`](environment.yml).

## Workflow

### 1. Clone this repository

```bash
git clone git@github.com:BranniganLab/mutations_BDNF_prodomain.git
cd mutations_BDNF_prodomain
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

### 2. Create and activate the Conda environment

```bash
conda env create -f environment.yml
conda activate mutations-bdnf-prodomain
```

### 3. Clone the Blobulator repository

The Blobulator is maintained in a separate repository and is required for the analysis workflow. This workflow uses Blobulator v1.1.0. Clone the [Blobulator](https://github.com/BranniganLab/blobulator) repository separately:

```bash
cd ..
git clone --branch v1.1.0 --single-branch git@github.com:BranniganLab/blobulator.git blobulator-v1.1.0
```

The recommended directory structure is:

```text
GitHub/
├── mutations_BDNF_prodomain/
└── blobulator-v1.1.0/
```

The analysis scripts expect the blobulation script for VMD to be located at:

```text
blobulator-v1.1.0/VMD_scripts/
```

### 4. Generate figures

Generate the figures for the paper:

```bash
cd mutations_BDNF_prodomain/
./generate_figs.sh
```

Note: This will directly download the raw trajectories from Zenodo.

Once the trajectories are downloaded and unzipped, process_trajectory.sh is called, which performs trajectory preprocessing, periodic-boundary handling, minimum-distance analysis, periodic image frame filtering, and equilibration-frame removal. The processed trajectories are placed in a folder called "post_processed" in the trajectories directory.

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

- Input and output directories
- Sequence names
- Equilibration time
- Blobulation arguments
- Blob indices
- Contact cutoffs
- Residue ranges
- State names
