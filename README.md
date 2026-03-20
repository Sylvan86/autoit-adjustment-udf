# Adjustment UDF — Least Squares Adjustment for AutoIt

**You have measurements. They don't perfectly agree. You need the best answer.**

Least squares adjustment finds the optimal values from redundant, contradictory observations — and tells you exactly how reliable the results are. It's the mathematical foundation behind positioning, navigation, calibration, curve fitting, regression, signal processing, and virtually any domain where measurements meet reality.

This UDF brings that capability to AutoIt with a radically simple interface: describe your model as string formulas, add measurements with their uncertainties, and call `_adj_solve()`. No matrix algebra, no numerical recipes, no external tools — just your problem and your data.

Underneath, it's a complete adjustment engine powered by OpenBLAS, with nonlinear iterative solvers, variance component estimation, robust estimators, and full statistical diagnostics — the kind of tooling normally reserved for MATLAB or specialized scientific software.

## Quick Start

```autoit
#include "Adjustment.au3"

; 5 measurements of a distance — find the best estimate
Local $mSystem = _adj_createSystem()
_adj_addObsFunction($mSystem, "M1", "X", 10.02, 1.0)
_adj_addObsFunction($mSystem, "M2", "X",  9.98, 1.0)
_adj_addObsFunction($mSystem, "M3", "X", 10.01, 1.0)
_adj_addObsFunction($mSystem, "M4", "X", 10.03, 1.0)
_adj_addObsFunction($mSystem, "M5", "X",  9.99, 1.0)
_adj_setInitialValue($mSystem, "X", 10.0)

_adj_solve($mSystem)
ConsoleWrite(_adj_displayResults($mSystem))

; Result: X = 10.006 ± 0.009
```

This is just the simplest case — a mean value. The real strength lies in **simultaneously adjusting different measurements with different formulas**. For example, determining a point's position from distance *and* direction measurements, or fitting a curve to data while enforcing constraints.

→ **[Tutorial: Learn step by step](docs/tutorial/01_getting_started.md)**

## Feature Overview

For experts — the full scope of the UDF:

### Model Types

| Base Type | Weighted | Generalized (full covariance matrix) | Description |
|-----------|----------|--------------------------------------|-------------|
| OLS | WLS | GLS | Overdetermined system |
| LSE | WLSE | GLSE | With parameter constraints |
| CLS | WCLS | GCLS | Condition adjustment |
| GLM | WGLM | GGLM | Gauss-Helmert model (most general form) |

Model type detection is automatic based on the input structure.

### Solvers and Iteration Methods

| Component | Options |
|-----------|---------|
| Iteration method | Gauss-Newton, Levenberg-Marquardt (Nielsen damping) |
| Linear solver | QR decomposition (DGELSY), Singular Value Decomposition (DGELSD) |
| Jacobian matrices | Numerical (Central, Forward, Backward, Ridder, Higham) or analytical |
| Scaling | Jacobi equilibration (automatic column scaling) |

### Statistics and Diagnostics

| Category | Measures |
|----------|----------|
| Basic statistics | A posteriori variance factor s₀², degrees of freedom, vᵀPv |
| Accuracy | Cofactor matrix Qxx, standard deviations (parameters, observations) |
| Controllability | Redundancy numbers r_i, cofactor matrices Qvv and Qŷ |
| Model validation | Global test (χ²) |
| Outlier diagnostics | Baarda test (w-statistic), Pope test (τ-statistic), p-values, MDB |

### Variance Component Estimation (VCE)

Helmert method: Separate variance factors for different observation groups (e.g. distances vs. angles). Iterative until convergence.

### Robust Estimation (IRLS)

| Estimator | Tuning constant | Breakdown point |
|-----------|----------------|-----------------|
| L1 (Median) | — | 50% |
| Huber | c = 1.345 | ~5% |
| Hampel | a = 1.7, b = 3.4, c = 8.5 | ~25% |
| Biweight (Tukey) | c = 4.685 | 50% |
| BIBER (Schweppe) | c = 3.5 | leverage-dependent |
| Modified-M (Koch) | c = 1.5 | leverage-dependent |

Scale parameters: MAD, s₀, a priori, user-defined.

### Additional Features

- Symbolic formula input as strings — parameters and observations are detected automatically
- Compute-on-demand — statistics are only calculated when requested
- Configurable result display with selectable columns and sections
- Outlier detection and removal (`_adj_getOutliers`, `_adj_removeObs`)

→ **[Reference: Full details](docs/reference/feature_overview.md)**

## Installation

### Requirements

- AutoIt v3.3.16+ (x64)
- OpenBLAS DLL (`libopenblas_x64.dll`)

### Setup

1. **Download:** Get the latest release from GitHub
2. **Extract:** Unpack the files into a directory
3. **OpenBLAS:** The `libopenblas_x64.dll` is downloaded automatically on first run, or can be placed manually in the UDF directory
4. **Include:**

```autoit
$__g_hBLAS_DLL = DllOpen("libopenblas_x64.dll")
#include "Adjustment.au3"
```

## Documentation

| Section | Description | Audience |
|---------|------------|----------|
| **[Tutorial](docs/tutorial/01_getting_started.md)** | Step-by-step guide, 9 chapters | Beginners |
| **[Feature Overview](docs/reference/feature_overview.md)** | Full capabilities at a glance | Everyone |
| **[API Reference](docs/reference/api.md)** | All 17 public functions | Developers |
| **[Configuration](docs/reference/configuration.md)** | Solver, display, and robust config | Advanced |
| **[Result Structure](docs/reference/results.md)** | All keys of the result map | Advanced |
| **[Model Types](docs/reference/model_types.md)** | OLS to GGLM with mathematics | Experts |
| **[Solvers](docs/reference/solvers.md)** | GN, LM, QR, SVD in detail | Experts |
| **[Statistics](docs/reference/statistics.md)** | s₀, Qxx, global test, Baarda/Pope | Experts |
| **[Robust Estimation](docs/reference/robust_estimation.md)** | IRLS, 6 estimators, weight functions | Experts |
| **[Error Codes](docs/reference/error_codes.md)** | All $ADJ_ERR_* with solutions | Everyone |

### Tutorial Chapters

1. [Getting Started](docs/tutorial/01_getting_started.md) — What is adjustment?
2. [First Network](docs/tutorial/02_first_network.md) — Trilateration
3. [Weighting](docs/tutorial/03_weighting.md) — Accounting for measurement precision
4. [Mixed Observations](docs/tutorial/04_mixed_observations.md) — Combining distances and angles
5. [Constraints](docs/tutorial/05_constraints.md) — Restrictions and fixed parameters
6. [Regression](docs/tutorial/06_regression.md) — OLS → orthogonal → Deming → York
7. [Covariance Matrix](docs/tutorial/07_covariance_matrix.md) — Correlated measurements
8. [Understanding Results](docs/tutorial/08_results.md) — Configuring and interpreting output
9. [Robust Estimation](docs/tutorial/09_robust_estimation.md) — Handling outliers

## Author

AspirinJunkie

## License

<!-- TODO: Add license -->
