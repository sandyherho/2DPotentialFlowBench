# Supporting Material for "Benchmarking Modern Scientific Computing Platforms for 2D Potential Flow Solver"

[![DOI](https://zenodo.org/badge/909530309.svg)](https://doi.org/10.5281/zenodo.14572734)
[![License: WTFPL](https://img.shields.io/badge/License-WTFPL-brightgreen.svg)](http://www.wtfpl.net/about/)
[![No Maintenance Intended](http://unmaintained.tech/badge.svg)](http://unmaintained.tech/)

This repository contains the code, data, and visualization outputs for benchmarking Python, Julia, and R implementations of a 2D potential flow solver around dual square obstacles.

## Repository Structure

```plaintext
.
├── README.md                   # This file
├── LICENSE.txt                 # License information
├── code/
│   ├── sim.py                 # Python implementation
│   ├── sim.jl                 # Julia implementation
│   ├── sim.R                  # R implementation
│   ├── matrix_plot.py         # Statistical significance matrix visualization
│   ├── violin_plot.py         # Performance violin plots
│   ├── performance.py         # Performance metrics computation
│   ├── stat_perf.py          # Statistical analysis of performance data
│   └── viz.py                # Flow visualization
├── outputs/
│   ├── data/                 # Performance measurements and analysis
│   │   ├── memory_usage_performance_stats.csv
│   │   ├── execution_time_performance_stats.csv
│   │   ├── memory_usage_comparison_stats.csv
│   │   └── execution_time_comparison_stats.csv
│   └── figs/                 # Generated figures
│       ├── memory_significance.png
│       ├── memory_violin.png
│       ├── runtime_significance.png
│       ├── runtime_violin.png
│       ├── convergence.png
│       ├── potential_surface.png
│       └── streamlines_pressure.png
```

## Citation

If you use this repository in your research, please cite:

```bibtex
@article{herhoKabanAnwar2024,
  author = {Herho, S. H. S. and Kaban, S. N. and Anwar, I. P.},
  title = {{Benchmarking Modern Scientific Computing Platforms for 2D Potential Flow Solver}},
  journal = {TBD},
  year = {2024},
  volume = {},
  number = {},
  pages = {},
  doi = {}
}
```

## Usage

### Prerequisites
Ensure you have the following dependencies installed:

- **Python**: `numpy`, `scipy`, `pandas`, `matplotlib`
- **Julia**: `LinearAlgebra`, `DataFrames`, `CSV`
- **R**: `tidyverse`

### Running the Simulations

1. **Python Implementation**
   ```bash
   python ./code/sim.py
   ```

2. **Julia Implementation**
   ```bash
   julia ./code/sim.jl
   ```

3. **R Implementation**
   ```bash
   Rscript ./code/sim.R
   ```

### Performance Analysis

Run the performance analysis scripts:
```bash
python ./code/performance.py
python ./code/stat_perf.py
```

### Visualization

Generate all visualization plots:
```bash
python ./code/matrix_plot.py
python ./code/violin_plot.py
python ./code/viz.py
```

### Outputs
- Performance measurement data is saved in `outputs/data/`
- Visualization figures are saved in `outputs/figs/`

### Available Figures

- `memory_significance.png`: Statistical significance matrix for memory usage
- `memory_violin.png`: Memory usage violin plots
- `runtime_significance.png`: Statistical significance matrix for runtime
- `runtime_violin.png`: Runtime violin plots
- `convergence.png`: L2 norm convergence plot
- `potential_surface.png`: Velocity potential surface plot
- `streamlines_pressure.png`: Streamlines and pressure coefficient plot

## License
This repository is licensed under the WTFPL License - see [WTFPL](http://www.wtfpl.net/about/) for details.
