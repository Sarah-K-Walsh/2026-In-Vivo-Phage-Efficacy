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

The included script `Scripts/00_setup.R` can be used to install all package dependencies at once.

Scripts in this repository use the `here` library to dynamically set paths. For this to work correctly, Rstudio must be opened by double-clicking on one of the files in `Scripts/`. Path errors will appear if Rstudio was first opened using a shortcut or a script from a different location.


## Contents
| Item                         | Description                                                                              |
|------------------------------|------------------------------------------------------------------------------------------|
| `Data/`                      | Contains all data files used in this study                                               |
| └─ `Individuals/`            | Example images of individual segmented larva                                             |
| └─ `Raw/`                    | Example unprocessed raw images                                                           |
| └─ `Validation/`             | Example validation images of well & larva segmentation                                   |
| └─ `Bacteria_order.csv/`     | Ladderised phylogenetic order of bacterial strains for figures                           |
| └─ `Bacteria_phylogeny.nwk/` | Maximum clade credibility Staphylococcaceae phylogeny                                    |
| └─ `Galleria_main.csv/`      | Mortality and melanisation dataset                                                       |
| └─ `Galleria_weights.csv/`   | Weight and pixel area pilot dataset                                                      |
| └─ `In_vitro.csv/`           | In vitro phage efficacy dataset                                                          |
| `Models/`                    | Contains all Rdata MCMCglmm model files used in analysis                                 |
| `Scripts/`                   | Contains all analysis scripts used in this study                                         |
| └─ `00_setup.R`              | Convenience R script to install all dependencies used in this repository                 |
| └─ `01_models.R`             | R script for data wrangling and model fitting                                            |
| `.here`                      | Empty text file used for dynamic pathing with the here() library                         |
