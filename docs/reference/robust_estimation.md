# Robust Estimation (IRLS)

The Adjustment module supports robust parameter estimation via Iteratively Reweighted Least Squares (IRLS). Six estimators are available: L1, Huber, Hampel, Biweight (Tukey), BIBER, and ModifiedM. The implementation runs as Phase 1 in the `_adj_solve` workflow -- before the optional VCE (Phase 2).

Related documentation: [Configuration](configuration.md) | [Result Structure](results.md) | [Statistical Measures](statistics.md) | [API Reference](api.md)


## IRLS Algorithm

### Principle

IRLS replaces L2 minimization with iterative reweighting: observations with large residuals receive smaller weights, thereby reducing their influence on the estimate. The procedure converges to the M-estimator solution.

### Workflow

```
Phase 1 (Robust):
  1. Initial L2 adjustment (standard GN/LM)
  2. [BIBER/ModifiedM only] Compute redundancy numbers r_i
  3. Save original weights, upgrade model type to weighted variant
  4. IRLS loop:
     a) Scale estimation sigma_hat (MAD, s0, apriori, or fixed)
     b) Compute standardized residuals u_i
     c) Compute robust weights w(u_i)
     d) Convergence check (max. relative weight change < threshold)
     e) Update whitening vectors
     f) Re-run GN/LM adjustment with updated weights
  5. Save IRLS results

Phase 2 (VCE/Final):
  Reset iteration state, then standard VCE or single-pass solution
```

### Standardized Residuals

The weight functions operate on standardized residuals u_i. Two variants:

**Standard (L1, Huber, Hampel, Biweight):**

    u_i = v_i / (sigma_hat * sigma_i)

**Leverage-based (BIBER, ModifiedM):**

    u_i = v_i / (sigma_hat * sigma_i * sqrt(r_i))

where r_i is the redundancy number from the initial solution. Observations with r_i < 10^(-10) are skipped (weight remains unchanged).

### Convergence

Convergence criterion: maximum relative weight change across all observations:

    max_i |w_i^(k) - w_i^(k-1)| / (1 + |w_i^(k-1)|) < epsilon

Default: epsilon = 10^(-3), maximum iterations = 30.

### Weight Update

The effective weight of an observation is the product:

    p_i^eff = p_i^original * w_i^robust

where p_i^original = 1/sigma_i^2 is the weight from the stochastic model.


## Scale Parameter

The scale parameter sigma_hat normalizes the residuals. Four methods:

| Method | Formula | Description |
|--------|---------|-------------|
| `"MAD"` | sigma_hat = median(\|v_i / sigma_i\|) / 0.6745 | Robust against outliers. Default. |
| `"s0"` | sigma_hat = sqrt(v^T P v / f) | A posteriori standard deviation. Less robust. |
| `"apriori"` | sigma_hat = 1 | Stochastic model is adopted unchanged. For BIBER (Wicki). |
| `"fixed"` | sigma_hat = `.fixedScale` | User-defined fixed value (e.g., Koch: sigma = 0.5). |

### MAD Fallback

When MAD = 0 (e.g., for exact solutions), an automatic fallback to `"s0"` occurs. If s_0 = 0 as well, the IRLS iteration is aborted.

### MAD Constant 0.6745

The factor 0.6745 = Phi^(-1)(3/4) normalizes the MAD so that it is consistent with the standard deviation for normally distributed data:

    MAD = median(|x_i - median(x_i)|) / 0.6745 ≈ sigma


## The 6 Estimators

### L1 (Least Absolute Deviations)

**Weight function:**

    w(u) = 1 / max(|u|, epsilon)     (capped at maxWeight)

| Parameter | Default | Description |
|-----------|---------|-------------|
| -- | -- | No tuning parameters |
| `maxWeight` | 1000 | Upper bound for w (singularity protection) |

**Property:** Minimizes sum |v_i| instead of sum v_i^2. Fully robust, but not smooth -- weights can fluctuate strongly. Good for heavy outliers.

### Huber

**Weight function:**

    w(u) = { 1           when |u| <= c
           { c / |u|     when |u| > c

| Parameter | Default | Description |
|-----------|---------|-------------|
| `c` | 1.345 | Transition point between quadratic and linear |

**Property:** Monotone redescending. Observations within c are treated as L2, beyond linearly. The default c = 1.345 yields 95% asymptotic efficiency under normality. Most widely used.

### Hampel

**Weight function (three zones + null zone):**

    w(u) = { 1                                    when |u| <= a
           { a / |u|                               when a < |u| <= b
           { a * (c - |u|) / ((c - b) * |u|)      when b < |u| <= c
           { 0                                     when |u| > c

| Parameter | Default | Description |
|-----------|---------|-------------|
| `a` | 1.7 | Boundary of the undisturbed region |
| `b` | 3.4 | Start of the linear decay zone |
| `c` | 8.5 | Complete rejection from here |

**Property:** Redescending with smooth transition. Zone [0, a]: full L2 treatment. Zone (a, b]: constant rho derivative (like Huber). Zone (b, c]: linear decrease to zero. Beyond c: complete rejection (weight = 0). More flexible than Huber due to three adjustable zones.

### Biweight (Tukey)

**Weight function:**

    w(u) = { (1 - (u/c)^2)^2     when |u| <= c
           { 0                    when |u| > c

| Parameter | Default | Description |
|-----------|---------|-------------|
| `c` | 4.685 | Rejection threshold |

**Property:** Smoothly redescending. Complete rejection for |u| > c. The default c = 4.685 yields 95% asymptotic efficiency under normality. The smooth weight function avoids discontinuities.

### BIBER (Bounded Influence by Bounded Residual)

**Weight function:** Identical to Huber.

    w(u) = { 1           when |u| <= c
           { c / |u|     when |u| > c

| Parameter | Default | Description |
|-----------|---------|-------------|
| `c` | 3.5 | Transition point |
| `scale` | `"MAD"` | Scale method (often `"apriori"` is used) |

**Difference from Huber:** The standardized residual accounts for leverage:

    u_i = v_i / (sigma_hat * sigma_i * sqrt(r_i))

This penalizes observations with high leverage (small r_i) more strongly. Requires redundancy numbers from the initial solution -- therefore a one-time computation of Q_xx/Q_vv/diag(R) after the first solve.

### ModifiedM

**Weight function:** Identical to Huber.

    w(u) = { 1           when |u| <= c
           { c / |u|     when |u| > c

| Parameter | Default | Description |
|-----------|---------|-------------|
| `c` | 1.5 | Transition point (more conservative than BIBER) |

**Difference from BIBER:** Same leverage consideration, but typically smaller c and different scale method. BIBER and ModifiedM differ only in their tuning parameters -- the weight function and leverage-based standardization are identical.


## Overview: Estimator Comparison

| Estimator | Type | Breakdown | Efficiency | Leverage | Smooth |
|-----------|------|-----------|------------|----------|--------|
| L1 | Convex | 50% | 64% | No | No |
| Huber | Convex, monotone | ~1/n | 95% (c=1.345) | No | No |
| Hampel | Redescending | Adjustable | Adjustable | No | No |
| Biweight | Redescending | ~50% | 95% (c=4.685) | No | Yes |
| BIBER | Convex + leverage | ~1/n | ~95% | Yes | No |
| ModifiedM | Convex + leverage | ~1/n | ~95% | Yes | No |


## Configuration

### Minimal Configuration

```autoit
Local $mConfig = _adj_defaultConfig()
$mConfig.robust = "Huber"
$mConfig.robustParams = _adj_robustDefaults("Huber")
_adj_solve($mSystem, $mConfig)
```

### `_adj_robustDefaults` -- Tuning Parameters

```autoit
Local $mParams = _adj_robustDefaults($sEstimator)
```

Returns the default tuning parameters for each estimator:

| Estimator | Returned Keys |
|-----------|---------------|
| `"L1"` | `.scale = "MAD"`, `.outlierThreshold = 2.5` |
| `"Huber"` | `.c = 1.345`, `.scale = "MAD"`, `.outlierThreshold = 2.5` |
| `"Hampel"` | `.a = 1.7`, `.b = 3.4`, `.c = 8.5`, `.scale = "MAD"`, `.outlierThreshold = 2.5` |
| `"Biweight"` | `.c = 4.685`, `.scale = "MAD"`, `.outlierThreshold = 2.5` |
| `"BIBER"` | `.c = 3.5`, `.scale = "MAD"`, `.outlierThreshold = 2.5` |
| `"ModifiedM"` | `.c = 1.5`, `.scale = "MAD"`, `.outlierThreshold = 2.5` |

All parameter maps contain `.scale = "MAD"` and `.outlierThreshold = 2.5`.

### Customizing Tuning Parameters

```autoit
Local $mConfig = _adj_defaultConfig()
$mConfig.robust = "Hampel"
$mConfig.robustParams = _adj_robustDefaults("Hampel")

; Anpassung der Tuning-Parameter
Local $mRobust = $mConfig.robustParams
$mRobust.a = 2.0         ; breiterer ungestörter Bereich
$mRobust.b = 4.0
$mRobust.c = 10.0        ; spätere vollständige Ablehnung
$mRobust.scale = "s0"    ; s0 statt MAD
$mConfig.robustParams = $mRobust

; IRLS-Konvergenz anpassen
$mConfig.robustMaxIter = 50
$mConfig.robustConvergence = 1e-5

_adj_solve($mSystem, $mConfig)
```

### Using a Fixed Scale

```autoit
Local $mConfig = _adj_defaultConfig()
$mConfig.robust = "Huber"
Local $mRobust = _adj_robustDefaults("Huber")
$mRobust.scale = "fixed"
$mRobust.fixedScale = 0.5    ; Benutzerdefinierter Wert
$mConfig.robustParams = $mRobust
_adj_solve($mSystem, $mConfig)
```


## Accessing Results

### Via `_adj_getResults`

```autoit
Local $mRes = _adj_getResults($mSystem)

; IRLS-Metadaten
ConsoleWrite("IRLS-Iterationen: "   & $mRes.robustIterations & @CRLF)
ConsoleWrite("Konvergiert: "        & $mRes.robustConverged  & @CRLF)
ConsoleWrite("Robuste Skala: "      & $mRes.robustScale      & @CRLF)

; Robuste Gewichte pro Beobachtung
Local $mWeights = $mRes.robustWeights
For $sKey In MapKeys($mWeights)
    ConsoleWrite($sKey & ": w = " & $mWeights[$sKey] & @CRLF)
Next
```

### Via Display

```autoit
Local $mDisplay = _adj_defaultDisplayConfig()
$mDisplay.showRobust = True
$mDisplay.obsCols = "name|value|v|robW"
ConsoleWrite(_adj_displayResults($mSystem, $mDisplay))
```


## Outlier Detection

### Via Robust Weights

```autoit
; Beobachtungen mit w < 0.5 als Ausreißer
Local $aOutliers = _adj_getOutliers($mSystem, 0.5, "robustWeight")

; Beobachtungen mit |u| > 2.5 als Ausreißer
Local $aOutliers = _adj_getOutliers($mSystem, 2.5, "absU")
```

### Iterative Removal

```autoit
; Ausreißer erkennen und entfernen
Local $aOutliers = _adj_getOutliers($mSystem, 0.3, "robustWeight")
If UBound($aOutliers) > 0 Then
    _adj_removeObs($mSystem, $aOutliers)
    ; Erneute Ausgleichung ohne robuste Schätzung
    _adj_solve($mSystem, _adj_defaultConfig())
EndIf
```


## Limitations

- **Generalized models (GLS/GLSE/GGLM/GCLS):** Robust estimation is currently not supported for models with a full covariance matrix (off-diagonal entries). Error: `$ADJ_ERR_INPUT`, `@extended = 7`.
- **Model type upgrade:** For initially unweighted models (OLS/LSE/CLS/GLM), the model type is automatically upgraded to the weighted variant (e.g., OLS --> WLS), so that the whitening transformation accounts for the IRLS weights.
- **Redundancy numbers (BIBER/ModifiedM):** Computed once from the initial solution and kept constant across all IRLS iterations (not updated). This is a common approximation.
- **L1 singularity:** When v_i is near zero, w = 1/|u| is protected by max(|u|, 10^(-10)) and capped at maxWeight (default: 1000).
