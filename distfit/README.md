# Distribution Fitter (`distfit`)

`distfit` is a high-performance C++ tool to fit discrete probability distributions to empirical data. It supports multiple distributions, computes maximum likelihood estimates, and evaluates model fit using negative log-likelihood (NLL) and information-theoretic criteria such as AIC and Akaike weights. The project can handle single files or entire directories of datasets and generate plots comparing empirical and theoretical cumulative distribution functions (CDFs).

---

## Features

- Fits the following distributions:
  - Geometric
  - Poisson
  - Altmann–Zeta (optional)
  - Truncated Zeta (optional)
- Supports multiple optimizers:
  - RandomOptimizer
  - RandomStallOptimizer
- Computes:
  - Negative log-likelihood (NLL)
  - Akaike Information Criterion (AIC)
  - Akaike weights (model probabilities)
- Writes fitted CDFs to text files for each dataset
- Generates plots using an R script (`plot_cdf.R`)
- Parallelized computation via OpenMP and asynchronous tasks
- Supports single dataset files or directories

---

## Requirements

- **C++ compiler** with C++17 support (e.g., `g++`)
- **OpenMP** support for parallel processing
- **R** installed for plotting (`plot_cdf.R` must be present)
- Standard C++ libraries (`<filesystem>`, `<future>`, etc.)

---

## Building

1. Clone or download the repository.
2. Navigate to the project root.
3. Run:

```bash
make
```

This will compile the source files in `src/` and produce the executable `distfit`.

### Clean build

To remove compiled objects and the executable:

```bash
make clean
```

---

## Usage

```bash
./distfit [options] <file_or_directory>
```

### Options

- `-d <directory>`: Treat the argument as a directory containing dataset files.
- `-fast`: Use a faster optimization strategy (RandomStallOptimizer only).
- `-compare`: Compare results between RandomOptimizer and RandomStallOptimizer.
- `-alternative`: Enable additional distributions (Altmann–Zeta and Truncated Zeta).

### Examples

- Fit a single file:

```bash
./distfit mydata.txt
```

- Fit all datasets in a directory:

```bash
./distfit -d datasets/
```

- Run fast optimization and enable alternative distributions:

```bash
./distfit -fast -alternative -d datasets/
```

---

## Output

- **Console output**:
  - Negative log-likelihoods for each distribution and dataset
  - Optimized parameters
  - AIC table
  - Akaike weights
  - Best model per dataset

- **Fitted CDFs**: Written to `fitted_models/` as text files:
  - Columns: `K EmpiricalCDF TheoreticalCDF`

- **Plots**: Generated via `plot_cdf.R`, saved as `.png` files alongside the text output.

---

## Project Structure

```
.
├── Makefile             # Build configuration
├── README.md
├── src/                 # C++ source files
├── obj/                 # Compiled object files (generated)
├── fitted_models/       # Output fitted CDFs and plots
├── plot_cdf.R           # R script for plotting CDFs
```

---

## Notes

- The project uses **parallelization** (OpenMP + async) for faster computation.
- The log-likelihood provides a **global measure of fit**, but local deviations may still exist. Examine CDF plots to assess local fitting quality.
- Poisson fits are often poor for heavy-tailed distributions; only a few datasets may resemble a Geometric distribution.
- Altmann–Zeta and Truncated Zeta fits are usually very similar, except in the distribution tails.

