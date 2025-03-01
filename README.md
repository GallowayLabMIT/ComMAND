# ComMAND

This repo contains all modeling, analysis, and code to create the figures in https://www.biorxiv.org/content/10.1101/2024.06.25.600629


## Modeling details

### Directory layout

- `modeling/mathematica` contains a Mathematica notebook with algebraic calculations for hte steady-state solution.
- `modeling/reference` contains the reference code for the [Equalizer](https://dx.doi.org/10.1038/s41467-021-23889-0) paper.
- `modeling/src` contains Julia code that simulates both the steady state and stochastic code. Generally, these simulations output simulation results as `.h5` or `.csv` files that are processed into figures via other Jupyter notebooks.
- `modeling/Project.toml` contains the Julia dependencies needed to run the simulations. This file is automatically used by Julia.
- `modeling/Manifest.toml` contains the specific "pinned" version of the dependencies.

### Julia setup

1. Open a terminal in the modeling figure folder (e.g. the folder with the `Manifest.toml` and `Project.toml` files)
2. Launch Julia 1.10 (other versions may work, but this is the version which we used) within this project folder (e.g. `julia --project=.`)
3. Using the Pkg interface, instantiate the pinned environment (`] instantiate`).
4. Run the desired Julia file directly using `julia --project=. src/steady_state.jl` or launch one of the Jupyter notebooks within the instantiated environment.

### Steady-state modeling

The file `modeling/src/miR_iFFL.jl` contains the library code implementing the steady state solution. The key files `steady_state.jl`, `steady_state_all_param_sweep.jl`, and `steady_state_distributions.jl` generate dataframes used for plotted figures.

### Stochastic modeling

The Jupyter notebooks within `modeling/src` both generate auxiliary plots and generate `.csv` files containing the results of stochastic simulations. Launch these notebooks using a Julia kernel within the instantiated environment.