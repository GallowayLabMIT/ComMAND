# ComMAND

This repo contains all modeling, data analysis, and code to create the figures in [Love et al. "Model-guided design of microRNA-based gene circuits supports precise dosage of transgenic cargoes into diverse primary cells"](https://www.biorxiv.org/content/10.1101/2024.06.25.600629).


## Modeling

### Directory layout

- `modeling/mathematica` contains a Mathematica notebook with algebraic calculations for the steady-state analytical solution.
- `modeling/reference` contains the reference code for the [Equalizer](https://dx.doi.org/10.1038/s41467-021-23889-0) paper, work that inspired this model.
- `modeling/src` contains Julia code that computes both the steady-state analytical solutions and the stochastic simulations. Generally, this modeling outputs results as `.gzip` files that are processed into figures via other Jupyter notebooks.
- `modeling/Project.toml` contains the Julia dependencies needed to perform the modeling. This file is automatically used by Julia.
- `modeling/Manifest.toml` contains the specific "pinned" versions of the dependencies.

### Julia setup

1. Open a terminal in the modeling figure folder (i.e., the folder with the `Manifest.toml` and `Project.toml` files)
2. Launch Julia 1.10 (other versions may work, but this is the version which we used) within this project folder (e.g. `julia --project=.`)
3. Using the Pkg interface, instantiate the pinned environment (`] instantiate`).
4. Run the desired Julia file directly using `julia --project=. src/steady_state.jl` or launch one of the Jupyter notebooks within the instantiated environment.

### Steady-state solutions

The file `modeling/src/miR_iFFL.jl` contains the library code implementing the steady-state solution. The key files `steady_state.jl`, `steady_state_all_param_sweep.jl`, and `steady_state_distributions.jl` generate dataframes (`.gzip` files) used for plotted figures.

### Stochastic simulations

The Jupyter notebooks within `modeling/src` both generate auxiliary plots and generate `.gzip` files containing the results of stochastic simulations. Launch these notebooks using a Julia kernel within the instantiated environment.

### Modeling results

The steady-state solution and stochastic simulation results used to generate the figures in the paper have also been deposited at Zenodo (DOI: [10.5281/zenodo.14946133](https://dx.doi.org/10.5281/zenodo.14946133)). Download `output.gzip`, unzip the directory, and move the directory to the root directory of this repository. Then, the figures with modeling results can be generated without running the Julia code.


## Figure generation

The directory `figures` contains all the code necessary to generate the plots in main and supplementary figures for the paper. This requires the raw and analyzed data deposited at Zenodo (DOI: [10.5281/zenodo.14946133](https://dx.doi.org/10.5281/zenodo.14946133)).
Note that this code generates data and modeling plots, but does not add the cartoons/other graphics in the paper figures (these were created in Adobe Illustrator).

### Directory contents

- Jupyter notebooks beginning `fig#` each generate a single output `.pdf` file containing all the plots for that figure.
- `base.py` contains functions that load subsets of the data and perform analysis (gating cells, computing summary statistics, etc.). It also contains helper functions for data processing and global figure style information. 
- `for-review.ipynb` generates plots presented to reviewers but not included in the paper. These plots can be found in the Transparent Peer Review file included in the Supplemental Information associated with the paper.

### Python setup

1. Create a virtual environment in the repository directory using `python -m venv env`
2. Activate the virtual environment using `source env/bin/activate` (MacOS, Linux) or `.\env\Scripts\activate` (Windows)
3. Install the current package versions for this project using `pip install -r requirements.txt`
4. Download the raw and analyzed data from Zenodo (DOI: [10.5281/zenodo.14946133](https://dx.doi.org/10.5281/zenodo.14946133)), specifically the `data.zip` file. Unzip this file.
5. Create a file in the root directory of the repo called `datadir.txt` that contains the absolute path to the `data` directory you just downloaded. This should be a single line.
6. Run the code to generate the modeling results (see "Modeling" section above). Alternatively, download the modeling results directly (see "Modeling results" above).


## Data analysis

This repo contains preliminary analysis in Python of the data used to create the figures in the paper. These Jupyter notebooks verify gating strategies, explore trends in the data, etc. but are not necessary to generate the paper figures.
Running this code requires the raw and analyzed data deposited at Zenodo (DOI: [10.5281/zenodo.14946133](https://dx.doi.org/10.5281/zenodo.14946133)), and some require additional data not included in the paper.
Note: Some notebooks may produce errors or not run correctly.

### Directory layout

- `flow` contains preliminary analysis of flow cytometry data
- `qPCR` contains preliminary analysis of RT-qPCR data

