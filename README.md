This repository contains the code and input data supporting the manuscript  
"Resource-Environmental Efficiency Trap in Swine Production: A Multi-Frontier SBM-DEA Analysis."

## Repository Structure

- `Codes/`: Core analysis scripts (`.R` and `.py`).
- `Data/`: Input data required to run the analyses (access to figshare).
- `LICENSE`: Apache-2.0 license.
- `README.md`: Unified documentation, including dependencies and usage.

## Code Modules

### Main Pipeline

- `Codes/run_pipeline.R`: End-to-end launcher for pre-flight checks, main model estimation, and post-analysis modules.
- `Codes/multi_frontier_sbm_dea.R`: Multi-frontier SBM-DEA estimation and core output generation.

### Frontier and Diagnostic Modules

- `Codes/meta_frontier_module.R`: Meta-frontier and technology gap ratio construction.
- `Codes/meta_frontier_diagnostics.R`: Meta-frontier diagnostic checks and validation plots.
- `Codes/slack_analysis_module.R`: Slack decomposition across frontiers and scale groups.
- `Codes/frontier_rank_comparison.R`: Cross-frontier rank consistency analysis.

### Scale, Mechanism, and Causal Modules

- `Codes/scale_efficiency_mechanism_analysis.R`: Scale-efficiency relationship and mechanism analysis.
- `Codes/scale_expansion_did_analysis.R`: R-based difference-in-differences analysis for expansion dynamics.
- `Codes/causal_inference_did.py`: Python-based DID and treatment-control matching workflow.
- `Codes/preexperiment_weight_ratio.R`: Pre-experiment weight-ratio calibration.

### Robustness Modules

- `Codes/temporal_effect_robustness.R`: Temporal control and period-level robustness analysis.
- `Codes/monte_carlo_sensitivity.R`: Monte Carlo sensitivity analysis.
- `Codes/underreporting_sensitivity.R`: Underreporting scenario sensitivity analysis.

## Dependencies

### R Environment

- Reference runtime: R 4.2 or higher.
- Core packages used across scripts:
  `readxl`, `writexl`, `dplyr`, `tidyr`, `stringr`, `purrr`, `tibble`, `data.table`,
  `ggplot2`, `ggrepel`, `scales`, `patchwork`, `cowplot`, `mgcv`, `AER`, `censReg`,
  `sandwich`, `lmtest`, `broom`, `MASS`, `Metrics`, and `reshape2`.
- Package installation uses standard CRAN provisioning.

### Python Environment

- Reference runtime: Python 3.10 or higher.
- Packages:
  `pandas`, `numpy`, `statsmodels`, `scipy`, `matplotlib`, `seaborn`, and `openpyxl`.
- Package installation uses standard PyPI provisioning.

## How to Run

1. Set the working directory to the repository root.
2. Provide the raw input dataset locally as `./swine_farm_data.xlsx` before running the full pipeline.
3. Run the main pipeline:
   - In R: `source("Codes/run_pipeline.R")`
4. (Optional) Run Python causal inference:
   - `python Codes/causal_inference_did.py`

## Reproducibility Scope

The package supports reproduction of the manuscript's primary frontier estimates, scale-mechanism analysis, causal inference outputs, and robustness checks, conditional on input data availability.
