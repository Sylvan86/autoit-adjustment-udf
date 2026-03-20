# API Reference — Adjustment UDF

Compact reference of all 17 public functions. Full description can be found in the function header in the source code.

---

## System Creation

### _adj_createSystem

**Signature:** `_adj_createSystem()`

**Return:** Map — empty adjustment system with `.model` and `.results` sub-maps.

**Note:** The system is populated via `_adj_addObs`, `_adj_addObsFunction`, etc. Also called implicitly when `_adj_addObs`/`_adj_addObsFunction` receive an uninitialized `$mSystem`.

```autoit
Local $mSystem = _adj_createSystem()
```

---

## Model Definition

### _adj_addObs

**Signature:** `_adj_addObs(ByRef $mSystem, $sSymbol = Default, $fValue = 0, $fStdDev = 1, $sVarComp = "s0")`

**Parameters:**
- `$mSystem` — [ByRef] Adjustment system (created if not a map)
- `$sSymbol` — [optional] Observation name. Default: auto-generated `"O1"`, `"O2"`, ...
- `$fValue` — [optional] Observation value. Default: `0`
- `$fStdDev` — [optional] A priori standard deviation. Default: `1`
- `$sVarComp` — [optional] VCE group name. Default: `"s0"`

**Return:** None (modifies `$mSystem` in-place).

**Error:** `$ADJ_ERR_INPUT` with `@extended`:
- `10` — Value is not a finite number
- `11` — Standard deviation must be positive and finite

**Note:** Symbol names are converted to uppercase. When `$fStdDev != 1`, the system is marked as weighted. When `$sVarComp != "s0"`, VCE is activated.

```autoit
_adj_addObs($mSystem, "W1", 60.5, 0.01, "Winkel")
```

---

### _adj_addObsFunction

**Signature:** `_adj_addObsFunction(ByRef $mSystem, $sSymbol = Default, $sFormula = "", $fValue = 0, $fStdDev = 1, $sVarComp = "s0", $fParamInitValue = 0.0, $vDerivatives = Default)`

**Parameters:**
- `$mSystem` — [ByRef] Adjustment system
- `$sSymbol` — [optional] Observation name
- `$sFormula` — [optional] Formula string (e.g., `"sqrt((X-x1)^2+(Y-y1)^2)"`)
- `$fValue` — [optional] Observation value. Default: `0`
- `$fStdDev` — [optional] A priori standard deviation. Default: `1`
- `$sVarComp` — [optional] VCE group. Default: `"s0"`
- `$fParamInitValue` — [optional] Initial value for new parameters. Default: `0.0`
- `$vDerivatives` — [optional] Analytical derivatives: map or pipe string (e.g., `"X=2*X | Y=2*Y"`). Default: numerical differentiation

**Return:** None (modifies `$mSystem` in-place).

**Error:** `$ADJ_ERR_INPUT` (see `_adj_addObs`).

**Note:** Parameters are automatically extracted from the formula via regex. Observation references use the `#` prefix (e.g., `#D1`). The formula is internally converted to uppercase. Creates both an observation and a formula entry.

**Model type impact:**
- 1 observation per formula (only its own symbol) → OLS
- \>1 observation per formula (additional `#` references) → GLM

```autoit
_adj_addObsFunction($mSystem, "S1", "sqrt((X-x1)^2+(Y-y1)^2)", 173.518, 0.01)
```

---

### _adj_addFunction

**Signature:** `_adj_addFunction(ByRef $mSystem, $sFormula, $fValue = 0, $fParamInitValue = 0.0, $vDerivatives = Default, $bRestriction = False)`

**Parameters:**
- `$mSystem` — [ByRef] Adjustment system
- `$sFormula` — Formula string (e.g., `"(#D1*cos(#T1)-XM)^2 + (#D1*sin(#T1)-YM)^2 - R^2"`)
- `$fValue` — [optional] Target value of the equation. Default: `0`
- `$fParamInitValue` — [optional] Initial value for new parameters. Default: `0.0`
- `$vDerivatives` — [optional] Analytical derivatives
- `$bRestriction` — [optional] `True` = treat as restriction instead of formula. Default: `False`

**Return:** None (modifies `$mSystem` in-place).

**Note:** For GLM/CLS models: formulas can mix parameters and observation references (`#` prefix). Called internally by `_adj_addRestriction` with `$bRestriction = True`.

```autoit
_adj_addFunction($mSystem, "#W1 + #W2 + #W3", 180.0)  ; CLS condition
```

---

### _adj_addRestriction

**Signature:** `_adj_addRestriction(ByRef $mSystem, $sFormula, $fValue = 0, $fParamInitValue = 0.0, $vDerivatives = Default)`

**Parameters:**
- `$mSystem` — [ByRef] Adjustment system
- `$sFormula` — Restriction formula (e.g., `"X + Y - 100"`)
- `$fValue` — [optional] Target value. Default: `0`
- `$fParamInitValue` — [optional] Initial value for new parameters. Default: `0.0`
- `$vDerivatives` — [optional] Analytical derivatives

**Return:** None.

**Note:** If the formula consists only of a single parameter name (optionally with sign), `_adj_addFixedParam` is called automatically. Mixed restrictions (parameters + observations) are handled correctly in GLM.

```autoit
_adj_addRestriction($mSystem, "X + Y - 100")     ; Parameter constraint
_adj_addRestriction($mSystem, "HA", 100.0)        ; → _adj_addFixedParam("HA", 100.0)
```

---

### _adj_addFixedParam

**Signature:** `_adj_addFixedParam(ByRef $mSystem, $sParam, $fValue)`

**Parameters:**
- `$mSystem` — [ByRef] Adjustment system
- `$sParam` — Parameter name (case-insensitive, internally uppercase)
- `$fValue` — Fixed value

**Return:** None.

**Note:** The parameter is removed from the adjustable parameter list and stored in `.model.fixed`. The substitution in formulas occurs in `__adj_parseFormulas` during `_adj_solve`.

```autoit
_adj_addFixedParam($mSystem, "HA", 100.0)
```

---

### _adj_addCovariance

**Signature:** `_adj_addCovariance(ByRef $mSystem, $sObs1, $sObs2, $fCovariance)`

**Parameters:**
- `$mSystem` — [ByRef] Adjustment system (observations must already exist)
- `$sObs1` — Name of the first observation
- `$sObs2` — Name of the second observation (= `$sObs1` for variance override)
- `$fCovariance` — Covariance value $\sigma_{12}$ (when `$sObs1 = $sObs2`: variance $\sigma^2$)

**Return:** `True` on success, `False` on error.

**Error:** `$ADJ_ERR_INPUT` with `@extended`:
- `0` — System not initialized
- `2` — `$sObs1` not found
- `3` — `$sObs2` not found
- `4` — Cross-group covariance (incompatible with VCE)

**Note:** When `$sObs1 = $sObs2`: updates `stdDev` to $\sqrt{\sigma^2}$. Off-diagonal covariances switch to a generalized model (GLS/GLSE/GGLM/GCLS).

```autoit
_adj_addCovariance($mSystem, "R1", "R2", 0.0005)  ; Covariance
_adj_addCovariance($mSystem, "R1", "R1", 0.0009)  ; Variance → stdDev = 0.03
```

---

### _adj_setInitialValue

**Signature:** `_adj_setInitialValue(ByRef $mSystem, $sParam, $fValue)`

**Parameters:**
- `$mSystem` — [ByRef] Adjustment system
- `$sParam` — Parameter name (case-insensitive)
- `$fValue` — Initial value / approximate value

**Return:** `True` on success, `False` on error.

**Error:** `$ADJ_ERR_INPUT` — parameter not found.

**Note:** Required for nonlinear models where parameters need good initial values for the iteration. Alternatively, `$fParamInitValue` can be set in `_adj_addObsFunction` / `_adj_addFunction` (applies to all newly discovered parameters of that formula).

```autoit
_adj_setInitialValue($mSystem, "X", 700.0)
_adj_setInitialValue($mSystem, "Y", 300.0)
```

---

## Configuration

### _adj_defaultConfig

**Signature:** `_adj_defaultConfig($sAlgorithm = "GN", $bVCE = True, $iMaxIter = 100, $fTolerance = 1e-10, $fAlpha = 0.05, $sSolver = "QR", $fSolverRCOND = Default)`

**Parameters:**
- `$sAlgorithm` — [optional] `"GN"` (Gauss-Newton) or `"LM"` (Levenberg-Marquardt). Default: `"GN"`
- `$bVCE` — [optional] Enable VCE. Default: `True`
- `$iMaxIter` — [optional] Max. iterations. Default: `100`
- `$fTolerance` — [optional] Convergence tolerance $\|\Delta x\|$. Default: `1e-10`
- `$fAlpha` — [optional] Significance level for global test. Default: `0.05`
- `$sSolver` — [optional] `"QR"` or `"SVD"`. Default: `"QR"`
- `$fSolverRCOND` — [optional] Rank threshold. Default: `1e-5`

**Return:** Map with solver configuration.

**Included sub-maps and fields:**

| Field | Type | Default | Description |
|------|-----|---------|-------------|
| `.compute` | Map | | Compute flags |
| `.compute.qxx` | Bool | `True` | Compute $Q_{xx}$ + $sd_x$ |
| `.compute.cofactors` | Bool | `False` | Compute $Q_{vv}$, $Q_{\hat{y}}$ |
| `.compute.redundancy` | Bool | `False` | Compute redundancy numbers |
| `.compute.globalTest` | Bool | `False` | Compute chi-squared test |
| `.compute.diagnostics` | Bool | `False` | Compute Baarda/Pope/MDB |
| `.diagnostics` | Map | | Diagnostics parameters |
| `.diagnostics.alpha` | Float | `0.001` | Significance level Baarda/Pope |
| `.diagnostics.beta` | Float | `0.20` | Type II error for MDB |
| `.diagnostics.alphaSuspect` | Float | `10*alpha` | Threshold for "suspect" |
| `.diagnostics.testBasis` | String | `"pope"` | `"pope"` or `"baarda"` |
| `.deriveMethod` | String | `"Central"` | Differentiation method |
| `.scaling` | Bool | `True` | Jacobian equilibration |
| `.robust` | String | `""` | Robust estimator (empty = disabled) |
| `.robustParams` | Map/Null | `Null` | Tuning parameters |
| `.robustMaxIter` | Int | `30` | Max. IRLS iterations |
| `.robustConvergence` | Float | `1e-3` | IRLS convergence threshold |

```autoit
Local $mConfig = _adj_defaultConfig("LM", False)
$mConfig.scaling = False
$mConfig.deriveMethod = "Ridder"
```

---

### _adj_defaultDisplayConfig

**Signature:** `_adj_defaultDisplayConfig()`

**Return:** Map with display configuration.

| Field | Type | Default | Description |
|------|-----|---------|-------------|
| `.showHeader` | Bool | `True` | Show header section |
| `.showParams` | Bool | `True` | Show parameter table |
| `.showObs` | Bool | `True` | Show observation table |
| `.showGlobalTest` | Bool | `False` | Show global test results |
| `.showVCE` | Bool | `False` | Show VCE results |
| `.showRobust` | Bool | `False` | Show robust estimation info |
| `.paramCols` | String | `"name\|value\|sdx"` | Parameter table columns |
| `.obsCols` | String | `"name\|value\|v\|sdv"` | Observation table columns |
| `.precision` | Int | `6` | Significant digits |

**Available columns:**

| Table | Column Name | Display | Compute Dependency |
|---------|-----------|---------|---------------------|
| Parameter | `name` | Name | — |
| Parameter | `value` | Value | — |
| Parameter | `sdx` | sd | `qxx` |
| Parameter | `xd` | dx | — |
| Observation | `name` | Name | — |
| Observation | `value` | Value | — |
| Observation | `v` | v | — |
| Observation | `sdv` | sd(v) | `cofactors` |
| Observation | `sdyhat` | sd(y) | `cofactors` |
| Observation | `r` | r | `redundancy` |
| Observation | `w` | \|w\| | `diagnostics` |
| Observation | `T` | Tau | `diagnostics` |
| Observation | `p` | p-value | `diagnostics` |
| Observation | `pPope` | p(Pope) | `diagnostics` |
| Observation | `blunder` | nabla | `diagnostics` |
| Observation | `mdb` | MDB | `diagnostics` |
| Observation | `decision` | Decision | `diagnostics` |
| Observation | `robW` | rob.w | (robust estimation) |

```autoit
Local $mDC = _adj_defaultDisplayConfig()
$mDC.obsCols = "name|value|v|r|w|decision"
$mDC.showGlobalTest = True
$mDC.precision = 8
```

---

### _adj_robustDefaults

**Signature:** `_adj_robustDefaults($sEstimator)`

**Parameters:**
- `$sEstimator` — Estimator name: `"L1"`, `"Huber"`, `"Hampel"`, `"Biweight"`, `"BIBER"`, `"ModifiedM"`

**Return:** Map with tuning parameters:

| Estimator | Fields |
|----------|--------|
| L1 | (no tuning parameters) |
| Huber | `.c = 1.345` |
| Hampel | `.a = 1.7`, `.b = 3.4`, `.c = 8.5` |
| Biweight | `.c = 4.685` |
| BIBER | `.c = 3.5` |
| ModifiedM | `.c = 1.5` |

All additionally contain: `.scale = "MAD"`, `.outlierThreshold = 2.5`.

**Error:** `$ADJ_ERR_INPUT` for unknown estimator.

```autoit
Local $mConfig = _adj_defaultConfig()
$mConfig.robust = "Huber"
$mConfig.robustParams = _adj_robustDefaults("Huber")
$mConfig.robustParams.c = 2.0  ; Adjust tuning
```

---

## Solving

### _adj_solve

**Signature:** `_adj_solve(ByRef $mSystem, $mConfig = Default)`

**Parameters:**
- `$mSystem` — [ByRef] Adjustment system (populated via `_adj_add*` functions)
- `$mConfig` — [optional] Configuration from `_adj_defaultConfig()`. Default: `_adj_defaultConfig()`

**Return:** Implicit (results in `$mSystem.results`). On error: `False`.

**Error:**
- `$ADJ_ERR_INPUT` — Invalid solver type or cross-group covariances
- `$ADJ_ERR_DOF` — Negative degrees of freedom (underdetermined system)
- `$ADJ_ERR_SOLVER` — LAPACK solver error
- `$ADJ_ERR_NO_CONVERGENCE` — No convergence
- `$ADJ_ERR_NOT_POS_DEF` — Covariance matrix $\Sigma_{\ell\ell}$ not positive definite

**Workflow:**
1. Model preparation: validation, formula parsing, index maps, classification, matrix allocation
2. Phase 1 (optional): Robust IRLS estimation
3. Phase 2: VCE loop or single pass

```autoit
_adj_solve($mSystem)                                    ; Default: GN + VCE
_adj_solve($mSystem, _adj_defaultConfig("LM", False))   ; LM without VCE
```

---

## Results and Post-Processing

### _adj_getResults

**Signature:** `_adj_getResults(ByRef $mSystem)`

**Parameters:**
- `$mSystem` — [ByRef] Adjustment system after `_adj_solve()`

**Return:** Map with flat result structure:

| Field | Type | Description |
|------|-----|-------------|
| `.modelType` | String | Model type (e.g., `"WLS"`, `"GGLM"`) |
| `.algorithm` | String | Algorithm used |
| `.s0` | Float | A posteriori $s_0$ |
| `.f` | Int | Degrees of freedom |
| `.vtpv` | Float | $v^T P v$ |
| `.nIterations` | Int | Number of iterations |
| `.rankDeficient` | Bool | Rank deficiency detected |
| `.solver` | String | Solver used |
| `.conditionNumber` | Float/Default | Condition number (SVD only) |
| `.x1` | Map | Adjusted parameters: name $\to$ value |
| `.sdx` | Map | Standard deviations of parameters |
| `.xd` | Map | Parameter corrections $\Delta x$ |
| `.obsValue` | Map | Observation values: name $\to$ value |
| `.v` | Map | Residuals: name $\to$ value |
| `.obsAdj` | Map | Adjusted observations: name $\to$ value |
| `.sdv` | Map | Standard deviations of residuals (if computed) |
| `.sdyhat` | Map | Standard deviations of adjusted observations (if computed) |
| `.r` | Map | Redundancy numbers (if computed) |
| `.Qxx` | Matrix/Null | Cofactor matrix of parameters |
| `.Qvv` | Matrix/Null | Cofactor matrix of residuals (if computed) |
| `.Qyhat` | Matrix/Null | Cofactor matrix of adjusted observations (if computed) |

Conditional fields (only when corresponding compute flags were active):
- `.globalTest*` — Global test results
- `.baardaW`, `.popeT`, `.pValue`, `.popePValue`, `.blunder`, `.mdb`, `.testDecision` — Diagnostics
- `.vceConverged`, `.vceIterations`, `.vceGroups` — VCE
- `.robustWeights`, `.robustIterations`, `.robustScale`, `.robustConverged` — Robust estimation

**Error:** `$ADJ_ERR_INPUT` if no results are available.

```autoit
Local $mRes = _adj_getResults($mSystem)
ConsoleWrite("s0 = " & $mRes.s0 & @CRLF)
ConsoleWrite("X = " & $mRes.x1["X"] & " +/- " & $mRes.sdx["X"] & @CRLF)
```

---

### _adj_displayResults

**Signature:** `_adj_displayResults(ByRef $mSystem, $mDisplayConfig = Default)`

**Parameters:**
- `$mSystem` — [ByRef] Adjustment system after `_adj_solve()`
- `$mDisplayConfig` — [optional] Display configuration. Default: `_adj_defaultDisplayConfig()`

**Return:** Formatted string with result summary. Empty string on error.

**Note:** Automatically triggers compute-on-demand for requested columns. Columns that cannot be computed show `"---"`.

```autoit
ConsoleWrite(_adj_displayResults($mSystem))

; With custom configuration
Local $mDC = _adj_defaultDisplayConfig()
$mDC.obsCols = "name|value|v|r|w|p|decision"
$mDC.showGlobalTest = True
ConsoleWrite(_adj_displayResults($mSystem, $mDC))
```

---

### _adj_getOutliers

**Signature:** `_adj_getOutliers(ByRef $mSystem, $fThreshold = Default, $sMethod = "robustWeight")`

**Parameters:**
- `$mSystem` — [ByRef] Adjustment system after `_adj_solve()`
- `$fThreshold` — [optional] Detection threshold (default: method-dependent)
- `$sMethod` — [optional] Detection method

| Method | Default Threshold | Criterion | Prerequisite |
|---------|-----------------|-----------|---------------|
| `"robustWeight"` | 0.5 | $w_{\text{robust}} < \text{threshold}$ | Robust estimation |
| `"absU"` | 2.5 | $\|u_i\| > \text{threshold}$ | Robust estimation |
| `"baarda"` | `.diagnostics.alpha` | p-value (Baarda) $< \alpha$ | Diagnostics |
| `"pope"` | `.diagnostics.alpha` | p-value (Pope) $< \alpha$ | Diagnostics |

**Return:** Array of observation names (may be empty).

**Error:** `$ADJ_ERR_INPUT` if prerequisites are not met.

```autoit
Local $aOutliers = _adj_getOutliers($mSystem, 0.3, "robustWeight")
For $i = 0 To UBound($aOutliers) - 1
    ConsoleWrite("Outlier: " & $aOutliers[$i] & @CRLF)
Next
```

---

### _adj_removeObs

**Signature:** `_adj_removeObs(ByRef $mSystem, $vObs)`

**Parameters:**
- `$mSystem` — [ByRef] Adjustment system
- `$vObs` — Observation name (string) or array of observation names

**Return:** Number of removed observations (Int).

**Error:** `$ADJ_ERR_INPUT` for invalid input type.

**Note:** Also removes associated formulas (ObsFunction type), covariances, and clears the prepared state. After removal, `_adj_solve()` must be called again.

```autoit
Local $aOutliers = _adj_getOutliers($mSystem)
_adj_removeObs($mSystem, $aOutliers)
_adj_solve($mSystem)  ; Recompute
```

---

### _adj_getErrorMessage

**Signature:** `_adj_getErrorMessage($iErr, $iExt = 0)`

**Parameters:**
- `$iErr` — Error code (`$ADJ_ERR_*` constant or `@error` value)
- `$iExt` — [optional] Extended error code (`@extended` value). Default: `0`

**Return:** String with error description (English).

```autoit
_adj_solve($mSystem)
If @error Then
    ConsoleWrite("Error: " & _adj_getErrorMessage(@error, @extended) & @CRLF)
EndIf
```
