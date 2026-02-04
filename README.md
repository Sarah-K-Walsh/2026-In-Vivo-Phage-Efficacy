# Phylogenetic Patterns of Bacterial Virulence and Phage Efficacy

This repository contains data and analysis scripts from **Walsh et al., (2026) Investigating phylogenetic patterns of bacterial virulence and phage efficacy in vivo**.

## Using this repository

### Platform Support
![Windows](https://img.shields.io/badge/Windows-blue?logo=microsoftwindows)
![macOS](https://img.shields.io/badge/macOS-black?logo=apple)
![Linux](https://img.shields.io/badge/Linux-grey?logo=linux)

### Dependencies
![R Version](https://img.shields.io/badge/R-4.5.2-2980b9)
![tidyverse](https://img.shields.io/badge/tidyverse-2.0.0-bbd4f1)
![here](https://img.shields.io/badge/here-1.0.2-bbd4f1)
![patchwork](https://img.shields.io/badge/patchwork-1.3.0-bbd4f1)
![viridisLite](https://img.shields.io/badge/viridisLite-0.4.2-bbd4f1)
![scales](https://img.shields.io/badge/scales-1.4.0-bbd4f1)


![MCMCglmm](https://img.shields.io/badge/MCMCglmm-2.36-1abc9c)
![minpack.lm](https://img.shields.io/badge/minpack.lm-1.2.4-1abc9c)
![ggtree](https://img.shields.io/badge/ggtree-3.16.0-1abc9c)
![ape](https://img.shields.io/badge/ape-5.8.1-1abc9c)

### Running Scripts

The included script `scripts/00_setup.R` can be used to install all package dependencies at once.

Scripts in this repository use the `here` library to dynamically set paths. For this to work correctly, Rstudio must be opened by double-clicking on one of the files in `scripts/`. Path errors will appear if Rstudio was first opened using a shortcut or a script from a different location.


## Contents
| Item                        | Description                                                                              |
|-----------------------------|------------------------------------------------------------------------------------------|
| `data/`                     | Contains all data files used in this study                                               |
| └─ `alphafold/`             | Contains HA alphafold model CIF file and SVG render                                      |
| └─ `neut75.csv`             | Neutralisation assay 1:75 dilution raw data                                              |
| └─ `sample_metadata.csv`    | Serum sample metadata file                                                               |
| └─ `gisaid_260126pn.pdf`    | GISAID record and DOI link                                                               |
| `scripts/`                  | Contains all analysis scripts used in this study                                         |
| └─ `00_setup.R`             | Convenience R script to install all dependencies used in this repository                 |
| └─ `01_analysis.R`          | R script for forward model building and statistical analysis of neutralisation data      |
