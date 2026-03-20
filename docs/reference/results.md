# Result Structure (Result Map)

The results of an adjustment are stored in `$mSystem.results` and are returned via `_adj_getResults($mSystem)` as a flat, user-friendly map. This documentation describes both the internal structure (`.results` sub-map) and the public result map.

Related documentation: [Configuration](configuration.md) | [Statistical Measures](statistics.md) | [Robust Estimation](robust_estimation.md) | [API Reference](api.md)


## Accessing Results

### `_adj_getResults` -- Flat Result Map

```autoit
_adj_solve($mSystem)
Local $mRes = _adj_getResults($mSystem)
If @error Then
    ConsoleWrite("Fehler: " & _adj_getErrorMessage(@error, @extended) & @CRLF)
    Exit
EndIf

; Parameter
ConsoleWrite("X = " & $mRes.x1["X"] & " +/- " & $mRes.sdx["X"] & @CRLF)

; Beobachtungen
ConsoleWrite("v(R1) = " & $mRes.v["R1"] & @CRLF)

; Globalstatistik
ConsoleWrite("s0 = " & $mRes.s0 & ", f = " & $mRes.f & @CRLF)
```

### `_adj_displayResults` -- Formatted Text Output

```autoit
; Standard-Ausgabe
ConsoleWrite(_adj_displayResults($mSystem))

; Erweiterte Ausgabe
Local $mDisplay = _adj_defaultDisplayConfig()
$mDisplay.showGlobalTest = True
$mDisplay.obsCols = "name|value|v|r|decision"
ConsoleWrite(_adj_displayResults($mSystem, $mDisplay))
```


## Complete Result Map

The following table lists all keys of the map returned by `_adj_getResults`.

### Metadata (Always Present)

| Key | Type | Description |
|-----|------|-------------|
| `.modelType` | String | Detected [model type](model_types.md): `"OLS"`, `"WLS"`, `"GLS"`, `"LSE"`, `"WLSE"`, `"GLSE"`, `"GLM"`, `"WGLM"`, `"GGLM"`, `"CLS"`, `"WCLS"`, `"GCLS"` |
| `.algorithm` | String | Algorithm used: `"GN"` or `"LM"` |
| `.s0` | Float | A posteriori standard deviation of unit weight s_0_hat = sqrt(v^T P v / f). See [Statistics](statistics.md#a-posteriori-variance-of-unit-weight). |
| `.f` | Int | Degrees of freedom (redundancy). Computation is model-dependent, see [Statistics](statistics.md#degrees-of-freedom). |
| `.vtpv` | Float | Weighted sum of squared residuals v^T P v |
| `.nIterations` | Int | Number of GN/LM iterations (1 for linear models) |
| `.rankDeficient` | Bool | `True` if rank deficiency was detected |
| `.solver` | String | Solver used: `"QR"` or `"SVD"` |
| `.conditionNumber` | Float/Default | Condition number sigma_max / sigma_min (only with SVD, otherwise `Default`) |

### Parameter Results (Maps: Parameter Name --> Value)

| Key | Type | Description |
|-----|------|-------------|
| `.x1` | Map | Adjusted parameter values. Access: `$mRes.x1["X"]` |
| `.sdx` | Map | Standard deviations of parameters: sd(x_i) = s_0 * sqrt(q_xx_ii). Access: `$mRes.sdx["X"]` |
| `.xd` | Map | Parameter changes Delta_x = x_1 - x_0. Access: `$mRes.xd["X"]`. Only meaningful for nonlinear models. |

For **fixed parameters** (via `_adj_addFixedParam`), no entries exist in these maps -- they were removed before the adjustment.

### Observation Results (Maps: Observation Name --> Value)

| Key | Type | Description |
|-----|------|-------------|
| `.obsValue` | Map | Original observation values l_i |
| `.v` | Map | Residuals v_i = l_hat_i - l_i |
| `.obsAdj` | Map | Adjusted observations l_hat_i = l_i + v_i |

### Conditional Observation Results (Only with Corresponding Compute Flags)

| Key | Type | Compute Flag | Description |
|-----|------|-------------|-------------|
| `.sdv` | Map | `cofactors` | Standard deviation of residual: sd(v_i) = s_0 * sqrt(q_vv_ii) |
| `.sdyhat` | Map | `cofactors` | Standard deviation of adjusted observation: sd(y_hat_i) = s_0 * sqrt(q_yy_ii) |
| `.r` | Map | `redundancy` | Redundancy numbers r_i = p_i * q_vv_ii, values in [0, 1]. See [Statistics](statistics.md#redundancy-numbers). |

### Cofactor Matrices (Only with Corresponding Compute Flags)

| Key | Type | Compute Flag | Description |
|-----|------|-------------|-------------|
| `.Qxx` | Matrix/Null | `qxx` | Cofactor matrix of parameters (u x u). See [Statistics](statistics.md#cofactor-matrix-q_xx). |
| `.Qvv` | Matrix/Null | `cofactors` | Cofactor matrix of residuals (n x n) |
| `.Qyhat` | Matrix/Null | `cofactors` | Cofactor matrix of adjusted observations (n x n) |

### Global Test Results (Only When Computed)

| Key | Type | Description |
|-----|------|-------------|
| `.globalTestPassed` | Bool | `True` if H_0 is not rejected |
| `.globalTestT` | Float | Test statistic T = v^T P v |
| `.globalTestLower` | Float | Lower bound chi2_(alpha/2, f) |
| `.globalTestUpper` | Float | Upper bound chi2_(1-alpha/2, f) |
| `.globalTestAlpha` | Float | Significance level used |

Details: [Statistics -- Global Test](statistics.md#global-test-chi2-test)

### Diagnostics Results (Only When Computed)

All diagnostics maps use observation names as keys. For observations with redundancy number r_i near zero (< 10^-12), the maps contain `Default` instead of numerical values.

| Key | Type | Description |
|-----|------|-------------|
| `.baardaW` | Map | Baarda test statistic \|w_i\| = \|v_i\| / (sigma_i * sqrt(r_i)) |
| `.popeT` | Map | Pope test statistic T_tau_i = \|w_i\| / s_0_hat |
| `.pValue` | Map | p-value from Baarda (normal distribution): 2 * (1 - Phi(\|w_i\|)) |
| `.popePValue` | Map | p-value from Pope (t-distribution): 2 * (1 - F_t(\|T_tau_i\|, f)) |
| `.blunder` | Map | Gross error estimate nabla_hat_i = v_i / r_i |
| `.mdb` | Map | Minimum detectable bias: MDB_i = sigma_i * sqrt(lambda_0 / r_i) |
| `.testDecision` | Map | Decision: `"ok"`, `"suspect"`, or `"outlier"` |
| `.testBasis` | String | Test basis: `"pope"` or `"baarda"` |
| `.baardaWarning` | Bool | `True` if s_0 deviates significantly from 1 (Baarda unreliable) |

Details: [Statistics -- Outlier Diagnostics](statistics.md#outlier-diagnostics)

### VCE Results (Only with Active Variance Component Estimation)

| Key | Type | Description |
|-----|------|-------------|
| `.vceConverged` | Bool | VCE converged |
| `.vceIterations` | Int | Number of VCE iterations |
| `.vceGroups` | Map | Results per group (key = group name) |

The `.vceGroups` map contains a sub-map per variance component group:

| Sub-Key | Type | Description |
|---------|------|-------------|
| `.s0` | Float | Estimated standard deviation of the group: sqrt(sigma2_hat_k) |
| `.sigma2` | Float | Estimated variance component sigma2_hat_k |

### Robust Estimation Results (Only with Active IRLS Estimation)

| Key | Type | Description |
|-----|------|-------------|
| `.robustWeights` | Map | Robust weights w_i per observation (key = observation name). Value range is estimator-specific, see [Robust Estimation](robust_estimation.md). |
| `.robustIterations` | Int | Number of IRLS iterations |
| `.robustScale` | Float | Last robust scale estimate sigma_hat |
| `.robustConverged` | Bool | IRLS converged |


## Outlier Detection (`_adj_getOutliers`)

```autoit
_adj_getOutliers(ByRef $mSystem[, $fThreshold = Default[, $sMethod = "robustWeight"]])
```

Returns an array of observation names identified as outliers.

### Methods and Defaults

| Method | Default Threshold | Criterion | Prerequisite |
|--------|-------------------|-----------|--------------|
| `"robustWeight"` | 0.5 | w_i < threshold | Robust estimation |
| `"absU"` | 2.5 | \|u_i\| > threshold | Robust estimation |
| `"baarda"` | `$mConfig.diagnostics.alpha` | p-value < alpha | Diagnostics computed |
| `"pope"` | `$mConfig.diagnostics.alpha` | p-value (Pope) < alpha | Diagnostics computed |

```autoit
; Robuste Ausreisser
Local $aOutliers = _adj_getOutliers($mSystem, 0.3, "robustWeight")
For $i = 0 To UBound($aOutliers) - 1
    ConsoleWrite("Ausreisser: " & $aOutliers[$i] & @CRLF)
Next

; Baarda-basierte Ausreisser
Local $aOutliers = _adj_getOutliers($mSystem, Default, "baarda")
```


## Removing Observations (`_adj_removeObs`)

```autoit
_adj_removeObs(ByRef $mSystem, $vObs)
```

Removes one or more observations (string or array) from the system. Also removes associated formulas (ObsFunction type) and covariances. Clears the prepared state -- a new `_adj_solve` call is required.

```autoit
; Einzelne Beobachtung entfernen
_adj_removeObs($mSystem, "R3")
_adj_solve($mSystem, $mConfig)

; Mehrere Beobachtungen auf einmal
Local $aOutliers = _adj_getOutliers($mSystem)
_adj_removeObs($mSystem, $aOutliers)
_adj_solve($mSystem, $mConfig)
```


## Internal Result Map (`$mSystem.results`)

The internal `.results` sub-map contains additional keys that are not directly exposed by `_adj_getResults`, but are used internally for compute-on-demand and VCE:

| Key | Type | Description |
|-----|------|-------------|
| `.r2sum` | Float | Copy of `$mState.r2sum` (= v^T P v) |
| `.r2sum0` | Float | Initial sum of squares (before iteration) |
| `.x0` | BLAS vector | Approximate values x_0 (before adjustment) |
| `.x1` | BLAS vector | Adjusted parameters x_1 (BLAS vector, not map!) |
| `.xd` | BLAS vector | Delta_x = x_1 - x_0 (BLAS vector) |
| `.r` | BLAS vector | Residual vector v (BLAS vector) |
| `.sdx` | BLAS vector | Standard deviations (BLAS vector) |
| `.sdv` | BLAS vector | sd(v) (BLAS vector) |
| `.sdy` | BLAS vector | sd(y_hat) (BLAS vector) |
| `.redundancyDiag` | BLAS vector | diag(R) (BLAS vector) |

**Important:** The internal `.results` map contains BLAS vectors (indexed via `$mState.idxObs` / `$mState.idxParams`). The public `_adj_getResults` map converts these into named maps (observation name / parameter name as key). Direct access to `$mSystem.results` is possible for advanced use cases but not recommended.
