<img src="logo.png" width="300px">  

# About

This repository contains the scripts and model outputs used to produce the results in the manuscript Ireson et al (2026, submitted to VZJ).

The simulations use [soilice v1.0](https://doi.org/10.5281/zenodo.19832264).

`soilice` v1.0 is a coupled mass and heat transport model for frozen soils. The code is written in Python and is designed to be concise, readable, and easy to customize, allowing alternative constitutive relationships and process representations to be readily tested. It is platform independent and uses just-in-time compilation and an ODE solver for computational efficiency while maintaining excellent mass and energy conservation.

## Reproducing the model runs

A recommended way to install the required packages and reproduce the model runs is:

```bash
# Clone this repository and navigate to the root folder
git clone XXX
cd XXX

# Create and activate a virtual environment (optional)
python -m venv .venv
source .venv/bin/activate

# Install the required packages
python -m pip install --upgrade pip
python -m pip install -r requirements.txt

# Run the model scripts (for example, the infiltration simulations)
cd Section4_Infiltration
make clean
make
```

Each folder contains a Makefile that runs the required model and plotting scripts and produces the corresponding figures and model output files.

The unfrozen flow benchmarks use MATLAB to evaluate the Mathias and Sander analytical solution. These calculations were performed using MATLAB R2025a Update 1 (25.1.0.2973910), 64-bit (maca64). The Makefile for these benchmarks runs the required MATLAB script in addition to the Python scripts.
