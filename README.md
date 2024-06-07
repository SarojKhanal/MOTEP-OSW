# Multi-Objective Transmission Expansion: An Offshore Wind Power Integration Case Study (MOTEP-OSW)

This repository contains the source code used to implement the model described in our paper, "Multi-Objective Transmission Expansion: An Offshore Wind Power Integration Case Study." The paper has been accepted for publication in a future issue of the IEEE Transactions on Energy Markets, Policy and Regulation, and is available as an early access article on [IEEE Xplore](https://ieeexplore.ieee.org/abstract/document/10504955). A preprint version of the paper is also available on [arXiv](https://arxiv.org/abs/2311.09563).

## Citing
If you benefit from the theoretical frameworks, methodologies, and findings presented in this repository or the paper, please use the following citation in your work:

```bibtex
@article{khanal2024multi,
  title={Multi-Objective Transmission Expansion: An Offshore Wind Power Integration Case Study},
  author={Khanal, Saroj and Graf, Christoph and Liang, Zhirui and Dvorkin, Yury and {\"U}nel, Bur{\c{c}}in},
  journal={IEEE Transactions on Energy Markets, Policy and Regulation},
  year={2024},
  publisher={IEEE}
}
```

## Test Systems
### 8-Zone Test System Based on ISO-NE Data ("Nodal", DC-OPF Model)
  - Network topology based on [Krishnamurthy et. al, 2015](https://ieeexplore.ieee.org/abstract/document/7039273)
  - Can be utilized for any so-called nodal model
  - Example configurations:
    - `ISONE8busTS_N5X5/`: Indicates 5 normal and 5 extreme days.
    - `ISONE8busTS_osw1_17_N5X5/`: Fixed POI (Point of Interconnection) specification with one degree of freedom (potential landfall nodes/zones of offshore power) and 17 candidate lines.
    - `ISONE8busTS_osw6_51_N5X5/`: Optimized POI specification with six degrees of freedom and 51 candidate lines.

The new operational scenarios (corresponding specifications denoted by `_NS`) better capture the effects of variable renewables and energy storage than those based on net load. The scenarios were generated using clustering methods from [Scott et al. (2019)](https://www.sciencedirect.com/science/article/pii/S0306261919312772), which used individual time series of wind, solar, and load, resulting in more temporally coupled samples.

### 9-Zone Test System Based on PJM Data ("Zonal", Pipes-and-Bubbles or Transportation Model):
  - Network topology based on [the U.S. EPA's Integrated Planning Model](https://www.epa.gov/power-sector-modeling/integrated-planning-model-ipm-run-files-supporting-scenarios)
  - Can be utilized for any so-called zonal model
  - Configuration adaptations include:
    - `EPA_PJM9Zones_N5X5R/`: With retirement data.
    - `EPA_PJM9Zones_osw2_17_N5X5R/` and `EPA_PJM9Zones_osw3_52_N5X5R/`: Different POI specifications with specified numbers of candidate lines.
  - `_NS` suffix denotes specifications incorporating the new scenarios.

## Running the Simulations

   - Begin by configuring the `cases.csv`, which sets the parameters and paths for input files. Note that the specification names used for the output directories are derived (starting) from the 5th column header in `cases.csv`. For serial execution of multiple model specifications, start them from the 5th column. For example, two headers, `SO` and `MO` in the 5th and 6th columns, can be used to run the two specifications serially in a single run. For each columns starting from 5th, update only parameter values that differs from defaults as defaults will apply from the `Default Value` column. See `cases_isone.csv` and `cases_pjm.csv` for references.
   - Run `TEP.jl` for simulations.
   - For PJM Case Study, update the `model.jl` and `output.jl` to `model_pjm.jl` and `output_pjm.jl` respectively in `TEP.jl` file. E.g., `include(joinpath(@__DIR__, "src", "model.jl"))` should look like `include(joinpath(@__DIR__, "src", "model_pjm.jl"))`.

## Inputs

Several input files from the `inputs/` directory are required for the model to run successfully. Each file is described in detail below, focusing on the ISO-NE Test System. The PJM Test System data is also collected and curated using a methodology similar to that used for ISO-NE, as described in the paper. [`raw data/`](inputs/raw%20data/) and [`raw data epa/`](inputs/raw%20data%20epa/) contains the raw, intermediate and final files for ISO-NE and PJM Test Systems respectively.

### Network Topology: Nodes/Zones, Lines, and Generators

- **`mappings_assumptions.csv`**
  - Details on per unit (p.u.) conversion bases and node/zone numbering.
  - Important for ensuring that nodal/zonal data are properly ordered as per the specifications in individual files.

- **`lines.csv`**
  - Provides transmission line data from [Krishnamurthy et. al, 2015](https://ieeexplore.ieee.org/abstract/document/7039273), complemented by line limits from https://github.com/jipkim/Trilevel

- **`generators_loc_emission_damages.csv`**
  - Enhanced `generators.csv`, which is the base generator data, with average marginal emissions values computed using [InMAP](https://inmap.run/).

### Operations: Load and Renewable Generation Time Series
- **`weights_of_scenarios.csv`**
  - Weights of representative days (operational scenarios), used in conjunction with time series data.
  - Derived from clustering algorithms to prioritize certain days based on similarity and significance. See `inputs/raw data*/Wind Scenarios*/` for details.

- **`load.csv`**
  - Hourly load data sourced from the ISO-NE website.
  - **Source**: [ISO-NE Energy, Load, and Demand Reports](https://www.iso-ne.com/isoexpress/web/reports/load-and-demand/-/tree/zone-info).

- **`solar.csv` and `wind.csv`**
  - Actual hourly generation data for solar and wind respectively, for representative days (computed from clustering algorithms, see [Wind Scenaarios](inputs/raw%20data/Wind%20Scenarios) for ISO-NE, as an example) of the year 2022.
  - **Columns**: Nodes/zones in order, consistent with `nodes.csv`.
  - **Rows**: Stacked by representative days, 24 hours at a time.
  - **Source**: [Daily Generation by Fuel Type, ISO-NE](https://www.iso-ne.com/isoexpress/web/reports/operations/-/tree/daily-gen-fuel-type).

- **`solar_normalized.csv` and `wind_normalized.csv`**
  - Timeseries data normalized by installed capacity for solar and wind respectively, to standardize the scale of generation across different zones and conditions.

**_Notes:_** To run the model with data sets other than those found in `inputs/`, ensure that all files are stored in the designated directory and are properly structured to avoid errors. To run the model with the pre-existing dataset, simply run `TEP.jl` and edit `cases.csv` as desired.

## Outputs
The results of the simulation run will appear in the `outputs/` directory, which is organized in line with the specification names listed in `cases.csv`. Each `outputs/` subdirectory contains a number of files that are arranged to facilitate the extraction and analysis of the data. Below is a list of the primary output files along with their structures.

### Configuration
- **`case*.csv`**
  - Contains the run configuration for each specification.

### Investment Decisions
- **`new_build_gen_cap_*_[power/energy].csv`**
  - Details new buildout capacities in terms of power or energy (for batteries).
  - **Rows**: Nodes/zones in order.
  - **Columns**: Epochs in order, e.g., 1 through 4.

- **`new_line_decision.csv`**
  - New line investment decisions coded as binary (0 for no investment, 1 for investment).
  - **Rows**: Lines in order, consistent with `lines.csv`.
  - **Columns**: Categorized by cable types for each epoch in order.

- **`line_upgrade_decision.csv`**
  - Decisions on existing line upgrades.
  - **Rows**: Lines in order, matching `lines.csv`.
  - **Columns**: Epochs in order.

- **`line_upgrade_decision_size.csv`**
  - Size of line upgrade decisions (used only in pipes-and-bubbles model, in this study the PJM test system).

### Investment Costs
- **`annual_cost_*_investment.csv`**
  - Discounted annual costs related to various categories such as line investments and externalities.
  - **Rows**: Sequential epochs, e.g., 1 through 4.

### Operations
- **`*load*.csv`** and related files (`*line_flow`, `theta.csv`, `wind_spillage.csv`, `*_out`, `[ch/dis]_node.csv`)
  - Timeseries data for loads, line flows, bus voltage angles, renewable curtailment, power output, and charging/discharging at nodes.
  - **Rows**: Stacked by nodes/zones, then by representative days.
  - **Columns**: Stacked by 24-hour periods for each epoch in sequence.

- **`Emissions*.csv`**
  - Total annual emissions for various pollutants like CO2 and NOx.
  - **Rows**: Epochs in order, from 1 through 4.

## Post-Processing After Run
The following post-processing scripts help to parse output files and prepare them for visualization and analysis of simulation results. Please see the complete instructions below for particular post-processing steps:

### Processing Costs and Capacities
- **Script**: `src/utils/process_costs_caps.csv`
- This script processes output files to extract capacity additions and cost outcomes.
- **Setup**:
  - For the ISO-NE Test System, set `n_lines_existing` to `12`, and for the PJM, set it to `13`.
  - Set `n_buses_existing` to `8` for the ISO-NE and `9` for the PJM.
  - Specify the result directory path: `result_dir = "/Users/skhanal/Globus_local/MOTEP-OSW/outputs/paper_draft/"`
  - Use a base data directory for comparison across different runs: `data_dir = "/Users/skhanal/Globus_local/MOTEP-OSW/inputs/ISONE8busTS_osw1_17_N5X5/"`
  - Specify the base case (specification) directory name to compare different specifications, e.g., `base_case_dirname = "SO"`

### Processing Topologies
- **Script**: `process_topologies.py`
- This script is used for processing and visualizing transmission buildouts.
- **Setup**:
  - Ensure you have the necessary Python packages installed. Refer to the `pyproject.toml` and `poetry.lock` files for setup guidance.
  - Update the `result_dir` and `data_dir` as required.
  - For PJM specific processing, rename scripts to use `_pjm.py` suffix (e.g., replace `process_transmissions.py` with `process_transmissions_pjm.py` and `plots.py` with `plots_pjm.py`).  

  **_Notes_:** The IPM v6 Regions shapefiles can be downloaded from the [U.S. EPA website](https://www.epa.gov/sites/default/files/2019-08/ipm_v6_regions.zip) and should be kept in the `inputs/` directory. Despite being necessary for executing `plots_pjm.py`, this file is excluded because of file size restrictions.

### Processing Generation Dispatch
- **Script**: `process_gens.py`
- This script processes generation dispatch.
- **Setup**:
  - Update the `result_dir`.
  - For PJM runs, set `GSw_PJM = True` in the `process_generations.py` script. You can use the same script the compute air quality damage details by simply changing `csv_path = os.path.join(subdir, "Pg_reduced.csv")` to `csv_path = os.path.join(subdir, "Pg_aq.csv")`, and the `generations.csv` then will have combined `Pg_aq.csv` from all specifications. You can also rename `generations.csv` to something like `Pg_aq_all.csv` avoid confusions in `final_df.to_csv(os.path.join(result_dir, "generations.csv"), index=False)`

### Computing Differences
- **Script**: `process_difference.py`
- This script is used to compute the differences between identical files from two specifications.
- **Setup**:
  - Set the base directory: `dir_base = "/Users/skhanal/Globus_local/MOTEP-OSW/outputs/ISONE/MO"`
  - Set the comparison directory: `dir_compare = "/Users/skhanal/Globus_local/MOTEP-OSW/outputs/ISONE/MO OPOI"`
  - Results will be stored in a directory named `MO OPOI-MO` as CSV files with the differences, which can be used for further analysis such as plotting air quality impacts using the `plot_airquality_damage_costs.ipynb` Jupyter notebook.

## Requirements
See `Project.toml` for the required packages.
