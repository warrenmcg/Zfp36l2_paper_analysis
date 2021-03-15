# Zfp36l2 Paper Analysis

This repository contains all of the necessary scripts and metadata to reproduce
the RNA-Seq analysis in the Zfp36l2 paper.

# Steps to reproduce the analysis

0. You must be using a Linux or MacOS X system.

1. Please make sure bioconda installed on your system. If you do not, go to
[this website](https://bioconda.github.io/user/install.html#install-conda) and
follow at least step 1.

2. Make sure you have Snakemake installed on your system. Please follow
the instructions at [this website](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

3. Clone this repository onto your computer.

4. Make sure the raw data is contained in the `data` directory.

5. Run the `run_analysis.sh` bash script. This will initialize the conda environment,
make sure the appropriate R packages are installed, initialize the annotations & metadata,
and run the analysis.

# Any problems? Questions?

If there are any problems running the pipeline or any questions about the analysis, please submit an issue.

# Remaining to-dos

- [ ] Clean up conda environment yaml file with only essential packages
