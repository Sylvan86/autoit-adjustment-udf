# 08 -- Results: Display and Diagnostics

## Two Ways to Access the Results

After `_adj_solve` the results are ready.  There are two ways to access them:

| Function | Purpose | Return value |
|----------|---------|--------------|
| `_adj_displayResults($mSystem)` | Formatted text output | String |
| `_adj_getResults($mSystem)` | Programmatic access | Map with all values |

In the previous chapters we always used `_adj_displayResults`. Now we take a closer look: what is contained in the results, and how do you configure the display?


## The Display Configuration

`_adj_displayResults` accepts an optional second parameter -- a display configuration.  Without this parameter, default settings are used.

```autoit
; Create default configuration and customize it
Local $mDisplay = _adj_defaultDisplayConfig()

; What should be displayed?
$mDisplay.showHeader     = True    ; Header line with model, s0, f, vtPv
$mDisplay.showParams     = True    ; Table of estimated parameters
$mDisplay.showObs        = True    ; Table of observations
$mDisplay.showGlobalTest = False   ; Global test (chi-squared)
$mDisplay.showVCE        = False   ; Variance component estimation
$mDisplay.showRobust     = False   ; Robust estimation

; Which columns in the tables?
$mDisplay.paramCols = "name|value|sdx"
$mDisplay.obsCols   = "name|value|v|sdv"

; Number of significant digits
$mDisplay.precision = 6

ConsoleWrite(_adj_displayResults($mSystem, $mDisplay))
```


## Available Columns

### Parameter Table (`paramCols`)

| Column | Meaning |
|--------|---------|
| `name` | Name of the parameter |
| `value` | Estimated value |
| `sdx` | Standard deviation of the parameter |
| `xd` | Change relative to the initial value (dx) |

### Observation Table (`obsCols`)

| Column | Meaning |
|--------|---------|
| `name` | Name of the observation |
| `value` | Measured value |
| `v` | Residual (correction) |
| `sdv` | Standard deviation of the residual |
| `sdyhat` | Standard deviation of the adjusted value |
| `r` | Redundancy number (0..1): How strongly does the network control this measurement? |
| `w` | Baarda test statistic (normalized residual) |
| `T` | Pope test statistic (tau test) |
| `p` | p-value (Baarda): Probability that the residual is due to chance |
| `pPope` | p-value (Pope): like `p`, but using estimated s0 |
| `blunder` | Estimated gross error (nabla) |
| `mdb` | Minimal Detectable Bias |
| `decision` | Test decision: "ok" or "OUTLIER" |
| `robW` | Robust weight (only with robust estimation, see Chapter 09) |

Columns are separated by `|`.  Example for a detailed display:

```autoit
$mDisplay.obsCols = "name|value|v|r|w|decision"
```


## Compute-on-Demand

Many columns require additional computations (e.g. redundancy, test statistics). The UDF computes these **automatically** as soon as you request a corresponding column.  You do not need to call anything extra.

For example: if you include `"r"` in `obsCols`, the UDF internally computes the redundancy matrix -- even if you did not explicitly request it beforehand.


## Example: Complete Diagnostics

```autoit
#include "Adjustment.au3"

Local $mSystem = _adj_createSystem()

; Trilateration from Chapter 02
_adj_addObsFunction($mSystem, "D1", "Sqrt((X - 0)^2 + (Y - 0)^2)",     72.05, 0.05)
_adj_addObsFunction($mSystem, "D2", "Sqrt((X - 100)^2 + (Y - 0)^2)",   63.10, 0.05)
_adj_addObsFunction($mSystem, "D3", "Sqrt((X - 100)^2 + (Y - 100)^2)", 78.20, 0.05)
_adj_addObsFunction($mSystem, "D4", "Sqrt((X - 0)^2 + (Y - 100)^2)",   85.15, 0.05)

_adj_setInitialValue($mSystem, "X", 50.0)
_adj_setInitialValue($mSystem, "Y", 50.0)

_adj_solve($mSystem, _adj_defaultConfig("GN", False))
If @error Then Exit MsgBox(48, "error", _adj_getErrorMessage(@error))

; Detailed display with global test and diagnostic columns
Local $mDisplay = _adj_defaultDisplayConfig()
$mDisplay.showGlobalTest = True
$mDisplay.obsCols = "name|value|v|r|w|decision"

ConsoleWrite(_adj_displayResults($mSystem, $mDisplay))
```

### What the Diagnostics Show

**Global test (chi-squared test):**
Tests whether the overall system is plausible.  "PASSED" means: the measurements are consistent with the assumed model and the specified standard deviations. "NOT PASSED" indicates a model error, incorrect standard deviations, or gross errors.

**Redundancy number r:**
Values close to 1: the measurement is strongly controlled by the network -- an error would be detected.  Values close to 0: the measurement is poorly controlled -- an error may go undetected.

**Test statistic w and decision:**
Large |w| values indicate an outlier.  The `decision` column summarizes the test result: "ok" or "OUTLIER".


## Programmatic Access: `_adj_getResults`

For further processing in your own code there is `_adj_getResults`. The function returns a map with named sub-maps:

```autoit
Local $mRes = _adj_getResults($mSystem)

; --- Metadata ---
ConsoleWrite("Modelltyp: " & $mRes.modelType & @CRLF)
ConsoleWrite("s0:        " & $mRes.s0 & @CRLF)
ConsoleWrite("f:         " & $mRes.f & @CRLF)
ConsoleWrite("vtPv:      " & $mRes.vtpv & @CRLF)
ConsoleWrite("Iterationen: " & $mRes.nIterations & @CRLF)

; --- Parameters ---
ConsoleWrite("X = " & $mRes.x1["X"] & " +/- " & $mRes.sdx["X"] & @CRLF)
ConsoleWrite("Y = " & $mRes.x1["Y"] & " +/- " & $mRes.sdx["Y"] & @CRLF)

; --- Residuals ---
ConsoleWrite("v(D1) = " & $mRes.v["D1"] & @CRLF)
ConsoleWrite("v(D2) = " & $mRes.v["D2"] & @CRLF)
```

### Available Fields in the Result Map

| Field | Type | Content |
|-------|------|---------|
| `.s0` | Number | Standard deviation of unit weight |
| `.f` | Number | Degrees of freedom |
| `.vtpv` | Number | Weighted sum of squared residuals |
| `.nIterations` | Number | Number of iterations until convergence |
| `.modelType` | String | Detected model type (e.g. "WLS", "GLM") |
| `.x1` | Map | Estimated parameters (key = parameter name) |
| `.sdx` | Map | Standard deviations of the parameters |
| `.xd` | Map | Change relative to initial values |
| `.v` | Map | Residuals (key = observation name) |
| `.obsAdj` | Map | Adjusted observation values (measured value + v) |
| `.Qxx` | Matrix | Cofactor matrix of the parameters (for advanced users) |


## Summary

- **`_adj_defaultDisplayConfig()`** creates a configurable display.
- **Columns** can be freely composed -- from simple (`"name|value|v"`) to detailed (`"name|value|v|r|w|decision"`).
- **Compute-on-Demand**: the UDF computes redundancy, test statistics, etc. automatically when you request the corresponding columns.
- **`_adj_getResults()`** returns all values as a map for programmatic processing.
- **Global test and individual tests** detect whether the model fits and whether individual measurements are suspicious.


---

Previous chapter: [07 -- Covariance Matrix: Correlations Between Measurements](07_covariance_matrix.md)

Next chapter: [09 -- Robust Estimation: Automatically Dampening Outliers](09_robust_estimation.md)
