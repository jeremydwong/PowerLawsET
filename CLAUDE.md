# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

### Project Setup
```bash
# Activate the Julia package environment
julia --project=.

# Or from within Julia REPL
using Pkg; Pkg.activate(".")
Pkg.instantiate()
```

### Running Notebooks
```bash
# Start Pluto notebook server
julia -e "using Pluto; Pluto.run(host=\"0.0.0.0\", port=1234, launch_browser=false)"

# Or use the alias in devcontainer
pluto

# Start Jupyter notebook (if using IJulia)
julia -e "using IJulia; notebook()"
```

### GitHub Codespaces Setup
This repository includes devcontainer configuration for GitHub Codespaces:

1. Click "Code" → "Codespaces" → "Create codespace on main"
2. The environment automatically installs Julia 1.10.0 and all dependencies
3. Pluto server runs on port 1234 (forwarded automatically)
4. Use `pluto` command to start the notebook server
5. Open `src/notebooks/reach_notebook.jl` in Pluto for interactive analysis

### Running Core Functions
```julia
# Basic usage pattern
using PowerLawsEt
f, mainseq, fit = replicate_pub_whk_sim("reachgrid_ctmech.mat")

# Plot surfaces
f = plot_vt_surfaces(mainseq)

# Fit power laws to data
fit = fit_powerlaws_to_oc(mainseq)
```

## Architecture

### Package Structure
- **PowerLawsEt.jl**: Main module that exports core functions
- **reach_et.jl**: Contains the reach_et submodule with all implementation details

### Core Data Types
- **main_sequence**: Struct holding distance, duration, peak_speed, and time_valuation data
- **fit_**: Mutable struct containing fitted power law functions and parameters

### Key Functions
- `replicate_pub_whk_sim()`: Loads default reaching data and generates surface plots with power law fits
- `plot_vt_surfaces()`: Creates 3D surface plots of duration and peak speed vs time valuation and distance
- `fit_powerlaws_to_oc()`: Fits power law functions to velocity and duration data
- `simple_grid_fillmissing()`: Interpolates missing values in 2D grid data

### Data Organization
- **data/reaching/**: Contains .mat files with reaching movement simulation data
- **data/reaching/subjects/**: Individual subject data files (ms_[condition]_subject_[N].mat)
- Conditions include: "eday", "fast", "slow", "infp"

### Power Law Models
The code implements power law relationships:
- Velocity: `V = k_v * c_t^(1/4) * (1/M)^(1/4) * L^(3/4)`
- Duration: `T = k_t * (1/c_t)^(1/4) * M^(1/4) * L^(1/4)`

Where:
- `c_t`: time valuation cost
- `L`: distance
- `M`: mass parameter
- `k_v`, `k_t`: fitted scaling constants

### Visualization
- Uses GLMakie for 3D surface plotting with interactive camera controls
- Surfaces show duration and peak speed as functions of time valuation and distance
- Optional interactive sliders for azimuth and elevation control