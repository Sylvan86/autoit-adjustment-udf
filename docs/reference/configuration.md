# Configuration

Three configuration maps control the behavior of the adjustment: the **Solver Config** (`_adj_defaultConfig`), the embedded **Diagnostics Config** (`.diagnostics`), and the **Display Config** (`_adj_defaultDisplayConfig`). Additionally, there is the **Robust Config** (`_adj_robustDefaults`), which is created separately and attached to the Solver Config.

Related documentation: [API Reference](api.md) | [Statistical Measures](statistics.md) | [Robust Estimation](robust_estimation.md) | [Model Types](model_types.md)


## 1. Solver Config (`_adj_defaultConfig`)

Creates the main configuration for `_adj_solve`. All parameters are optional -- without arguments, a sensible default configuration is used.

### Signature

```autoit
_adj_defaultConfig([$sAlgorithm = "GN"[, $bVCE = True[, $iMaxIter = 100[, $fTolerance = 1e-10[, $fAlpha = 0.05[, $sSolver = "QR"[, $fSolverRCOND = Default]]]]]]])
```

### Parameters

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `.algorithm` | String | `"GN"` | Iteration algorithm: `"GN"` (Gauss-Newton) or `"LM"` (Levenberg-Marquardt). See [Solver Documentation](solvers.md). |
| `.vce` | Bool | `True` | Enable Variance Component Estimation (Helmert). Only effective with multiple variance component groups. |
| `.maxIterations` | Int | `100` | Maximum iterations for GN/LM convergence. |
| `.tolerance` | Float | `1e-10` | Convergence threshold for ||dx||. |
| `.alpha` | Float | `0.05` | Significance level for the [Global Test](statistics.md#global-test-chi2-test). |
| `.solver` | String | `"QR"` | LAPACK solver: `"QR"` (DGELSY) or `"SVD"` (DGELSD). |
| `.solverRCOND` | Float | `1e-5` | Rank threshold for DGELSY/DGELSD. Singular values < RCOND * sigma_max are treated as zero. |
| `.scaling` | Bool | `True` | Jacobi equilibration (column scaling) for improved conditioning. |
| `.deriveMethod` | String | `"Central"` | Numerical differentiation method: `"Central"`, `"Central4"`, `"Forward"`, `"Backward"`, `"Ridder"`, `"Higham"`. |

### LM-Specific Parameters

These parameters are set internally and control the Levenberg-Marquardt damping behavior:

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `.lmTau` | Float | `1e-3` | Initialization factor tau for lambda_0 = tau * max(diag(J^T W J)). |
| `.lmLambdaMin` | Float | `1e-7` | Lower bound for lambda. |
| `.lmLambdaMax` | Float | `1e7` | Upper bound for lambda (termination criterion). |

### VCE-Specific Parameters

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `.vceMaxIter` | Int | `15` | Maximum VCE iterations. |
| `.vceConvergence` | Float | `0.05` | Relative convergence threshold for variance components (5%). |

### Robust Parameters

These keys are set in the config to activate IRLS estimation. Details in [Robust Estimation](robust_estimation.md).

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `.robust` | String | `""` | Estimator name (`""` = no robust estimation). Valid values: `"L1"`, `"Huber"`, `"Hampel"`, `"Biweight"`, `"BIBER"`, `"ModifiedM"`. |
| `.robustParams` | Map/Null | `Null` | Tuning parameter map from `_adj_robustDefaults()`. |
| `.robustMaxIter` | Int | `30` | Maximum IRLS iterations. |
| `.robustConvergence` | Float | `1e-3` | Convergence threshold for robust weight change. |


## 2. Compute Flags (`.compute` Sub-Map)

The `.compute` sub-map controls which statistics are computed after solving. Dependencies are resolved automatically (compute-on-demand via `__adj_ensureComputed`).

| Key | Type | Default | Description | Dependency |
|-----|------|---------|-------------|------------|
| `.qxx` | Bool | `True` | Cofactor matrix Q_xx + standard deviations sd_x | -- |
| `.cofactors` | Bool | `False` | Cofactor matrices Q_vv, Q_yhat + sd(v), sd(y) | qxx |
| `.redundancy` | Bool | `False` | Redundancy numbers diag(R) = diag(P * Q_vv) | cofactors |
| `.globalTest` | Bool | `False` | chi2 global test (H_0: sigma_0^2 = 1) | -- |
| `.diagnostics` | Bool | `False` | Baarda |w|, Pope T_tau, p-values, gross error nabla, MDB | redundancy |

**Automatic activation:** When VCE is active (`$mConfig.vce = True`), `.cofactors` and `.redundancy` are automatically set to `True`, since VCE requires these matrices.

### Dependency Chain

```
diagnostics --> redundancy --> cofactors --> qxx
```

Each request automatically resolves the chain to the left. Example: `.compute.diagnostics = True` automatically computes qxx, cofactors, and redundancy.


## 3. Diagnostics Config (`.diagnostics` Sub-Map)

Controls the parameters for outlier diagnostics. Only evaluated when `compute.diagnostics = True` or the display configuration requests diagnostics columns.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `.alpha` | Float | `0.001` | Significance level for Baarda/Pope test. **Note:** Not identical to `$mConfig.alpha` (0.05, for global test)! |
| `.beta` | Float | `0.20` | Type II error (power = 1 - beta = 0.80) for MDB computation. |
| `.alphaSuspect` | Float/Default | `Default` | Threshold for "suspect" decision. When `Default`: 10 * alpha = 0.01. |
| `.testBasis` | String | `"pope"` | Test basis: `"pope"` (a posteriori s_0_hat, normal approximation) or `"baarda"` (a priori sigma_0 = 1). |

### Baarda vs. Pope

| Property | Baarda | Pope |
|----------|--------|------|
| Variance reference | a priori sigma_0 = 1 | a posteriori s_0_hat |
| Test statistic | \|w_i\| = \|v_i\| / (sigma_i * sqrt(r_i)) | T_tau_i = \|w_i\| / s_0_hat |
| Distribution | Standard normal distribution N(0,1) | Student's t-distribution t(f) |
| Advantage | Classical approach | Robust when sigma_0 deviates |
| Warning | Issued when s_0 > 3.0 or s_0 < 0.3 | -- |

See [Statistical Measures](statistics.md#outlier-diagnostics) for the formulas.


## 4. Display Config (`_adj_defaultDisplayConfig`)

Creates the configuration for `_adj_displayResults`. Independent of the Solver Config -- controls only the text output.

### Signature

```autoit
_adj_defaultDisplayConfig()
```

### Parameters

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `.showHeader` | Bool | `True` | Show header section (model type, s_0, f, v^T P v). |
| `.showParams` | Bool | `True` | Table of adjusted parameters. |
| `.showObs` | Bool | `True` | Table of observations. |
| `.showGlobalTest` | Bool | `False` | Show global test result. |
| `.showVCE` | Bool | `False` | Show VCE result. |
| `.showRobust` | Bool | `False` | Show robust estimation info. |
| `.paramCols` | String | `"name\|value\|sdx"` | Pipe-separated column names for parameter table. |
| `.obsCols` | String | `"name\|value\|v\|sdv"` | Pipe-separated column names for observation table. |
| `.precision` | Int | `6` | Significant digits (StringFormat `%g` format). |

### Available Columns

**Parameter Table (`paramCols`):**

| Column Name | Display Label | Description | Compute Dependency |
|-------------|---------------|-------------|-------------------|
| `name` | Name | Parameter name | -- |
| `value` | Value | Adjusted value x_1 | -- |
| `sdx` | sd | Standard deviation sd(x) = s_0 * sqrt(q_xx_ii) | qxx |
| `xd` | dx | Parameter change Delta_x = x_1 - x_0 | -- |

**Observation Table (`obsCols`):**

| Column Name | Display Label | Description | Compute Dependency |
|-------------|---------------|-------------|-------------------|
| `name` | Name | Observation name | -- |
| `value` | Value | Observation value l | -- |
| `v` | v | Residual v_i | -- |
| `sdv` | sd(v) | Standard deviation of residual | cofactors |
| `sdyhat` | sd(y) | Standard deviation of adjusted observation | cofactors |
| `r` | r | Redundancy number r_i (0..1) | redundancy |
| `w` | \|w\| | Baarda test statistic | diagnostics |
| `T` | Tau | Pope test statistic T_tau | diagnostics |
| `p` | p-value | p-value (Baarda, normal distribution) | diagnostics |
| `pPope` | p(Pope) | p-value (Pope, t-distribution) | diagnostics |
| `blunder` | nabla | Gross error estimate nabla_hat = v_i / r_i | diagnostics |
| `mdb` | MDB | Minimum detectable bias | diagnostics |
| `decision` | Decision | Test decision: `"ok"`, `"suspect"`, `"outlier"` | diagnostics |
| `robW` | rob.w | Robust weight w_i (0..1 for redescending, 0..1000 for L1) | robust |

**Compute-on-Demand:** When a column is requested whose dependency has not yet been computed, the computation is triggered automatically. If the computation fails, the column displays `"---"`.


## 5. Configuration Examples

### Standard Adjustment (Minimal Configuration)

```autoit
; Keine explizite Config nötig -- _adj_solve verwendet Defaults
_adj_solve($mSystem)
; Äquivalent zu:
_adj_solve($mSystem, _adj_defaultConfig())
```

### LM Solver without VCE, with SVD

```autoit
Local $mConfig = _adj_defaultConfig("LM", False)
$mConfig.solver = "SVD"
_adj_solve($mSystem, $mConfig)
```

### Full Statistics with Diagnostics

```autoit
Local $mConfig = _adj_defaultConfig()
Local $mCompute = $mConfig.compute
$mCompute.cofactors   = True
$mCompute.redundancy  = True
$mCompute.globalTest  = True
$mCompute.diagnostics = True
$mConfig.compute = $mCompute
_adj_solve($mSystem, $mConfig)
```

### Requesting Diagnostics via Display Config (Compute-on-Demand)

```autoit
; Solve mit Standard-Config (Diagnostik nicht explizit aktiviert)
_adj_solve($mSystem)

; Display-Config fordert Diagnostik-Spalten an -- Berechnung wird automatisch ausgelöst
Local $mDisplay = _adj_defaultDisplayConfig()
$mDisplay.showGlobalTest = True
$mDisplay.obsCols = "name|value|v|r|w|T|p|decision"
ConsoleWrite(_adj_displayResults($mSystem, $mDisplay))
```

### Configuring Robust Estimation

```autoit
Local $mConfig = _adj_defaultConfig()
$mConfig.robust = "Huber"
$mConfig.robustParams = _adj_robustDefaults("Huber")
$mConfig.robustMaxIter = 50
$mConfig.robustConvergence = 1e-4

; Custom-Tuning:
Local $mRobust = $mConfig.robustParams
$mRobust.c = 2.0         ; statt Default 1.345
$mRobust.scale = "s0"    ; statt Default "MAD"
$mConfig.robustParams = $mRobust

_adj_solve($mSystem, $mConfig)
```

### Customized Display Output

```autoit
Local $mDisplay = _adj_defaultDisplayConfig()
$mDisplay.showRobust     = True
$mDisplay.showVCE        = True
$mDisplay.precision      = 8
$mDisplay.paramCols      = "name|value|sdx|xd"
$mDisplay.obsCols        = "name|value|v|sdv|r|robW"
ConsoleWrite(_adj_displayResults($mSystem, $mDisplay))
```


## Internal Constants

The following constants are defined in `Adjustment.au3` and are used internally:

| Constant | Value | Usage |
|----------|-------|-------|
| `$__ADJ_MAX_ITERATIONS` | 100 | Default maximum for GN/LM |
| `$__ADJ_TOLERANCE` | 1e-10 | Default convergence threshold |
| `$__ADJ_STAGNATION_TOL` | 1e-10 | Stagnation detection |
| `$__ADJ_MARQUARDT_REL_FLOOR` | 1e-10 | LM relative minimum for diag(D) |
| `$__ADJ_MARQUARDT_ABS_FLOOR` | 1e-30 | LM absolute minimum for diag(D) |
| `$__ADJ_LM_TAU_DEFAULT` | 1e-3 | LM initialization factor |
| `$__ADJ_LM_LAMBDA_MIN` | 1e-7 | LM lambda lower bound |
| `$__ADJ_LM_LAMBDA_MAX` | 1e7 | LM lambda upper bound |
| `$__ADJ_VCE_MAX_ITER` | 15 | Maximum VCE iterations |
| `$__ADJ_VCE_CONVERGENCE` | 0.05 | VCE convergence threshold |
| `$__ADJ_VCE_MIN_SIGMA2` | 1e-6 | Minimum variance component (numerical guard) |
