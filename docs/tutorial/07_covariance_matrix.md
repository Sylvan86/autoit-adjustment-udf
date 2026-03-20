# 07 -- Covariance Matrix: Correlations Between Measurements

## Measurements Are Not Always Independent

So far, all observations had only their own standard deviation -- they were **independent** of each other.  In practice this is often not the case:

- Three measurements of the same distance under identical conditions: The errors are related because the same instrument and the same atmosphere affect all three.
- GPS coordinates: The X and Y components are often correlated because the same satellites influence both.
- Transformed measurements: When measured values are converted (e.g. polar points to Cartesian coordinates), the correlation of the original measurements propagates.

These **correlations** are described by the covariance matrix.  It contains the variances (StdDev^2) on the diagonal and the covariances between each pair of measurements in the off-diagonal entries.


## Example: Three Correlated Distance Measurements

Three measurements of the same distance X.  All have a standard deviation of 0.3 m.  But: measurements M1/M2 and M2/M3 are correlated (covariance 0.05) because they were taken close in time with the same instrument.

| Measurement | Value [m] | StdDev [m] |
|-------------|-----------|------------|
| M1          | 10.02     | 0.3        |
| M2          | 10.05     | 0.3        |
| M3          |  9.97     | 0.3        |

The corresponding covariance matrix:

```
         M1      M2      M3
M1    [ 0.09    0.05    0.00 ]
M2    [ 0.05    0.09    0.05 ]
M3    [ 0.00    0.05    0.09 ]
```

Diagonal: 0.09 = 0.3^2 (variance of each measurement).
Off-diagonal: 0.05 (covariance between adjacent measurements).
M1 and M3 are not directly correlated (covariance = 0).


## Without Covariances

For comparison, first the solution **without** covariances -- i.e. as if the measurements were independent:

```autoit
#include "Adjustment.au3"

Local $mSystem = _adj_createSystem()

_adj_addObsFunction($mSystem, "M1", "X", 10.02, 0.3)
_adj_addObsFunction($mSystem, "M2", "X", 10.05, 0.3)
_adj_addObsFunction($mSystem, "M3", "X",  9.97, 0.3)

_adj_setInitialValue($mSystem, "X", 10.0)

_adj_solve($mSystem)
If @error Then Exit MsgBox(48, "error", _adj_getErrorMessage(@error))

ConsoleWrite(_adj_displayResults($mSystem))
```

The result is the simple mean: X = 10.0133 m (= (10.02 + 10.05 + 9.97) / 3).
The standard deviation of X is approximately 0.173 m (= 0.3 / Sqrt(3)).


## With Covariances

Now we add the correlations.  For this we use the function `_adj_addCovariance($mSystem, $sObs1, $sObs2, $fCovariance)`:

```autoit
#include "Adjustment.au3"

Local $mSystem = _adj_createSystem()

_adj_addObsFunction($mSystem, "M1", "X", 10.02, 0.3)
_adj_addObsFunction($mSystem, "M2", "X", 10.05, 0.3)
_adj_addObsFunction($mSystem, "M3", "X",  9.97, 0.3)

; Covariances between the measurements
_adj_addCovariance($mSystem, "M1", "M2", 0.05)
_adj_addCovariance($mSystem, "M2", "M3", 0.05)

_adj_setInitialValue($mSystem, "X", 10.0)

_adj_solve($mSystem)
If @error Then Exit MsgBox(48, "error", _adj_getErrorMessage(@error))

ConsoleWrite(_adj_displayResults($mSystem))
```

### What Changes

- **X shifts slightly**: The correlations change the effective weighting of the measurements.  M2 is correlated with both others and thereby loses some "information content".
- **The standard deviation of X increases**: Correlated measurements together provide less independent information than uncorrelated ones. Three strongly correlated measurements are almost as good as a single one.
- **s0 changes**: The stochastic model (the covariance matrix) affects how well the residuals fit the model.

### Note the Order

The covariances must be added **after** the observations -- both observations must already exist.  The order of the two observation names does not matter: `_adj_addCovariance($mSystem, "M1", "M2", 0.05)` is identical to `_adj_addCovariance($mSystem, "M2", "M1", 0.05)`.


## When Are Covariances Needed?

| Situation | Covariances needed? |
|-----------|---------------------|
| Independent measurements with different instruments | No -- standard deviations suffice |
| Measurements under the same systematic influences | Yes -- model the correlation |
| Transformed measurements (e.g. polar to Cartesian) | Yes -- error propagation creates covariances |
| GPS baselines (X, Y, Z of a baseline) | Yes -- receiver geometry creates correlation |

**Rules of thumb:**

- Positive covariance (> 0): The errors tend in the same direction. If M1 is too large, M2 is probably also too large.
- Negative covariance (< 0): The errors tend in opposite directions.
- Covariance = 0: The measurements are independent.
- The covariance must be smaller in absolute value than the product of the standard deviations: |Cov(M1,M2)| < StdDev(M1) * StdDev(M2).


## Summary

- **`_adj_addCovariance`** registers a covariance between two observations.  Internally, the UDF builds the full covariance matrix from these entries.
- Covariances account for the fact that measurements are **not independent**.
- Without covariances, correlated measurements are overvalued -- the estimated precision is too optimistic.
- In the simplest case (no correlation) the standard deviation per observation suffices -- covariances are optional.


---

Previous chapter: [06 -- Regression: Fitting a Line to Measured Data](06_regression.md)

Next chapter: [08 -- Results: Display and Diagnostics](08_results.md)
