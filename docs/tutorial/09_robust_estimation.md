# 09 -- Robust Estimation: Automatically Dampening Outliers

## The Problem with Outliers

Classical adjustment (least squares) has a weakness: **A single gross error can distort the entire result.**

Imagine a leveling network with 10 precise measurements.  One of them has a reading error of 5 cm.  The adjustment distributes this error across **all** measurements -- the result is systematically shifted.

The global test (Chapter 08) does detect **that** something is wrong. But it does not correct it automatically.  This is exactly what **robust estimation** does: it automatically dampens the influence of suspicious measurements without requiring you to know in advance which measurement is affected.


## How Does It Work?

Robust estimation works iteratively (IRLS -- Iteratively Reweighted Least Squares):

1. **Standard adjustment** as initial solution.
2. **Evaluate residuals**: Which residuals are unusually large?
3. **Adjust weights**: Observations with large residuals receive a smaller weight.
4. **Solve again** with the new weights.
5. Repeat steps 2-4 until the weights no longer change.

At the end, each observation has a **robust weight** between 0 and 1:

- Weight close to 1: The measurement fits well -- full influence.
- Weight close to 0: The measurement is suspicious -- almost no influence.


## Example: Height Measurement with an Outlier

Five height differences in a small network.  One measurement (H3) has a gross error:

| Measurement | Height difference [m] | StdDev [m] | Remark |
|-------------|----------------------|------------|--------|
| H1          | 1.005                | 0.003      | correct   |
| H2          | 1.003                | 0.003      | correct   |
| H3          | 1.050                | 0.003      | **outlier** (+4.5 cm) |
| H4          | 1.002                | 0.003      | correct   |
| H5          | 1.004                | 0.003      | correct   |

Goal: the true height difference **DH**.


### Without Robust Estimation (classical)

```autoit
#include "Adjustment.au3"

Local $mSystem = _adj_createSystem()

_adj_addObsFunction($mSystem, "H1", "DH", 1.005, 0.003)
_adj_addObsFunction($mSystem, "H2", "DH", 1.003, 0.003)
_adj_addObsFunction($mSystem, "H3", "DH", 1.050, 0.003)   ; outlier!
_adj_addObsFunction($mSystem, "H4", "DH", 1.002, 0.003)
_adj_addObsFunction($mSystem, "H5", "DH", 1.004, 0.003)

_adj_setInitialValue($mSystem, "DH", 1.0)

_adj_solve($mSystem, _adj_defaultConfig("GN", False))
If @error Then Exit MsgBox(48, "error", _adj_getErrorMessage(@error))

ConsoleWrite("=== Klassische Loesung ===" & @CRLF)
ConsoleWrite(_adj_displayResults($mSystem))
```

The result: DH = 1.0128 m -- the outlier pulls the mean up by approximately 9 mm.  The true value is around 1.004 m.


### With Robust Estimation (Huber)

```autoit
#include "Adjustment.au3"

Local $mSystem = _adj_createSystem()

_adj_addObsFunction($mSystem, "H1", "DH", 1.005, 0.003)
_adj_addObsFunction($mSystem, "H2", "DH", 1.003, 0.003)
_adj_addObsFunction($mSystem, "H3", "DH", 1.050, 0.003)   ; outlier!
_adj_addObsFunction($mSystem, "H4", "DH", 1.002, 0.003)
_adj_addObsFunction($mSystem, "H5", "DH", 1.004, 0.003)

_adj_setInitialValue($mSystem, "DH", 1.0)

; Configure robust estimation
Local $mConfig = _adj_defaultConfig("GN", False)
$mConfig.robust       = "Huber"
$mConfig.robustParams = _adj_robustDefaults("Huber")

_adj_solve($mSystem, $mConfig)
If @error Then Exit MsgBox(48, "error", _adj_getErrorMessage(@error))

; Display results with robust info and weights
Local $mDisplay = _adj_defaultDisplayConfig()
$mDisplay.showRobust = True
$mDisplay.obsCols    = "name|value|v|robW"

ConsoleWrite("=== Robuste Loesung (Huber) ===" & @CRLF)
ConsoleWrite(_adj_displayResults($mSystem, $mDisplay))
```

### What Changes

- **DH shifts back** toward 1.004 m -- the outlier is downweighted.
- In the **robW** column you can see the robust weights:
  - H1, H2, H4, H5: weight close to 1 (good measurements, full influence)
  - H3: weight significantly below 1 (outlier, dampened)
- The header shows the number of **IRLS iterations** and whether the estimation has **converged**.


## The Estimators

The UDF offers six robust estimators.  They differ in how aggressively they downweight outliers and whether they account for the geometry of the network (leverage).

### L1

Minimizes the **sum of absolute values** instead of squares. For a single unknown this is equivalent to finding the **median** -- the value that is most resistant to outliers.  For multiple unknowns the principle generalizes accordingly.

```autoit
$mConfig.robust       = "L1"
$mConfig.robustParams = _adj_robustDefaults("L1")
```

- Highest breakdown point (50%) -- tolerates up to half the data being outliers.
- Less efficient than other estimators when the data is clean.
- **Caveat:** Even without outliers, the L1 solution can deviate noticeably from the classical least squares (L2) result.  This is inherent to the method -- the median and the mean generally differ unless the data is perfectly symmetric.
- Best for: heavily contaminated data where you have no idea how many outliers exist.

### Huber

A compromise: **small residuals** are treated quadratically (like the classical method).  **Large residuals** (beyond threshold c) are only treated linearly -- dampened but not eliminated.

```autoit
$mConfig.robust       = "Huber"
$mConfig.robustParams = _adj_robustDefaults("Huber")
```

- Default c = 1.345 (95% efficiency at the normal distribution).
- Almost as efficient as the classical method with clean data.
- A common choice when outliers are expected but the data is mostly clean.

### Hampel

A three-zone estimator: small residuals are untouched, medium residuals are dampened, and large residuals are **completely rejected** (weight = 0).

```autoit
$mConfig.robust       = "Hampel"
$mConfig.robustParams = _adj_robustDefaults("Hampel")
; defaults: a=1.7, b=3.4, c=8.5
```

- More aggressive than Huber: observations beyond c get zero weight.
- Three tuning constants (a, b, c) control the transition zones.
- Best for: when you want to fully eliminate gross errors rather than just dampen them.

### Biweight (Tukey)

Like Hampel, a **redescending** estimator -- observations beyond threshold c receive zero weight.  The weight function is smooth (no sharp transitions).

```autoit
$mConfig.robust       = "Biweight"
$mConfig.robustParams = _adj_robustDefaults("Biweight")
; default: c=4.685
```

- High breakdown point (50%) combined with good efficiency (95%).
- Smooth weight function -- well-behaved convergence.
- Best for: a good all-round choice when you want hard rejection of outliers.

### BIBER (Schweppe)

Like Huber, but the standardized residuals are divided by sqrt(r_i) (the square root of the redundancy number).  This makes the estimator **leverage-aware**: observations that are poorly controlled by the network geometry are treated more strictly.

```autoit
$mConfig.robust       = "BIBER"
$mConfig.robustParams = _adj_robustDefaults("BIBER")
; default: c=3.5
```

- Requires an initial computation of redundancy numbers (done automatically).
- Best for: geodetic networks where some observations have low redundancy and could mask outliers in a standard Huber analysis.

### Modified-M (Koch)

Same principle as BIBER (leverage-aware Huber), but with a stricter threshold (c = 1.5 instead of 3.5).  Based on Koch (1999).

```autoit
$mConfig.robust       = "ModifiedM"
$mConfig.robustParams = _adj_robustDefaults("ModifiedM")
; default: c=1.5
```

- More aggressive than BIBER -- flags observations sooner.
- Best for: networks where you suspect multiple small outliers that a standard analysis might miss.

### Which Estimator to Choose?

| Situation | Recommendation |
|-----------|----------------|
| Few outliers expected, mostly clean data | **Huber** (high efficiency) |
| Heavy contamination (>10% outliers) | **L1** or **Biweight** |
| Want to fully reject (not just dampen) outliers | **Hampel** or **Biweight** |
| Geodetic network with varying redundancy | **BIBER** or **Modified-M** |
| Smooth convergence behavior preferred | **Biweight** |
| Simple, no tuning parameters to worry about | **L1** |


## Identifying Outliers

After robust estimation you can query the detected outliers:

```autoit
Local $mRes = _adj_getResults($mSystem)

; Read robust weights directly
For $sKey In MapKeys($mRes.robustWeights)
    Local $fWeight = ($mRes.robustWeights)[$sKey]
    If $fWeight < 0.5 Then
        ConsoleWrite("suspect: " & $sKey & " (weight = " & $fWeight & ")" & @CRLF)
    EndIf
Next
```

Or compactly with `_adj_getOutliers`:

```autoit
Local $aOutliers = _adj_getOutliers($mSystem)
If UBound($aOutliers) > 0 Then
    For $i = 0 To UBound($aOutliers) - 1
        ConsoleWrite("outlier: " & $aOutliers[$i] & @CRLF)
    Next
Else
    ConsoleWrite("no outliers detected." & @CRLF)
EndIf
```


## Configuration Options

| Option | Meaning | Default value |
|--------|---------|---------------|
| `$mConfig.robust` | Estimator name: "L1", "Huber", "Hampel", "Biweight", "BIBER", "ModifiedM" | "" (disabled) |
| `$mConfig.robustParams` | Tuning parameters (from `_adj_robustDefaults`) | Null |
| `$mConfig.robustMaxIter` | Maximum number of IRLS iterations | 30 |
| `$mConfig.robustConvergence` | Convergence threshold (relative weight change) | 0.001 |


## When to Use Robust Estimation?

- When you suspect gross errors in the data but don't know which observations are affected.
- When the global test (Chapter 08) indicates model problems.
- When clean data cannot be guaranteed (field measurements, automated sensors, etc.).
- With clean data and no suspicion of outliers, the classical method is sufficient.

For choosing the right estimator, see the selection table above.


## Summary

- **`$mConfig.robust = "Huber"`** (or any other estimator name) activates robust estimation.
- **`_adj_robustDefaults("...")`** provides the default tuning parameters for the chosen estimator.
- The UDF iterates automatically (IRLS) until the weights are stable.
- Each observation receives a **robust weight** (0..1). Outliers get a weight close to 0.
- Six estimators available: **L1**, **Huber**, **Hampel**, **Biweight**, **BIBER**, **Modified-M** -- from simple to leverage-aware.
- Each estimator has its strengths -- the right choice depends on the specific problem (see the selection table above).
- The **robW** column in the display configuration shows the robust weights.
- **`_adj_getOutliers`** returns a list of detected outliers.


---

Previous chapter: [08 -- Results: Display and Diagnostics](08_results.md)
