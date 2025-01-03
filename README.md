# Kernel-based entropy coefficient estimation

This Matlab repository contains the code and data used as presented in the paper W.D. Widanage et al. "A system identification approach to estimate lithium-ion battery entropy coefficients".


## Requirements

- MATLAB (version R2021 and above)
- System Identification Toolbox
- Signal Processing Toolbox

## Usage

Download or clone the Github repo and add it to your Matlab path.

### Generating reference temperautre profiles and obtaining the entropy coefficient 
1. To use the "EntropyCoeffEstimator.m" class see the notebooks in the Example_Notebooks folder

### Paper figure creation
To create the figures as presented in the paper, run the following Matlab script.
1. Run `Paper_Figures.m` found in the Create_Paper_Figures folder

## Data

### measurements_Aug2023 (as used in the paper)
- Contains the kernel-based and potentiometric based data from 100% to 0% performed every 5% SoC

### measurements_Jul2022
- Contains kernel-based and potentiometric based data from 100% to 0% performed every 10% SoC

## License
BSD 3-Clause License

## Citation
TBC

