# Error Codes — Adjustment UDF

Complete reference of all error constants with causes, diagnostics, and solutions.

---

## Overview

| Constant | Value | Meaning |
|-----------|------|-----------|
| `$ADJ_ERR_OK` | 0 | No error |
| `$ADJ_ERR_SOLVER` | 1 | LAPACK solver error |
| `$ADJ_ERR_DOF` | 2 | Invalid degrees of freedom |
| `$ADJ_ERR_RANK` | 3 | Rank deficiency |
| `$ADJ_ERR_NO_CONVERGENCE` | 4 | No convergence |
| `$ADJ_ERR_NOT_POS_DEF` | 5 | Matrix not positive definite |
| `$ADJ_ERR_INPUT` | 6 | Invalid input |
| `$ADJ_ERR_RESTR_OVERDETERMINED` | 7 | Restrictions fully determine all parameters |

---

## Error Handling in AutoIt

AutoIt resets `@error` on **every** function call. The error code must therefore be saved in a variable immediately after the call:

```autoit
_adj_solve($mSystem)
Local $iErr = @error, $iExt = @extended   ; SAVE IMMEDIATELY!
If $iErr Then
    ConsoleWrite("Error: " & _adj_getErrorMessage($iErr, $iExt) & @CRLF)
    Return
EndIf
```

---

## $ADJ_ERR_OK (0) — No Error

Successful execution. All results are available.

---

## $ADJ_ERR_SOLVER (1) — LAPACK Solver Error

**Message:** `"LAPACK solver error (LAPACK info = <iExt>)"`

**`@extended`:** LAPACK `INFO` value of the failed routine.

### Typical Causes

| Cause | Diagnostics |
|---------|-----------|
| Singular or rank-deficient system | The whitened system Ã·Δx = w̃ passed to DGELSY/DGELSD is (near-)singular |
| Numerical instability | Very different parameter magnitudes (try enabling scaling or SVD solver) |
| Faulty Jacobian matrix | Zero rows/columns due to incorrect formulas |
| DGELSY/DGELSD internal error | LAPACK convergence failure (rare) |

### Solutions

1. **Use SVD solver:** `_adj_defaultConfig("GN", True, Default, Default, Default, "SVD")` — more tolerant of rank problems
2. **Check equilibration:** `.scaling = True` (default) improves conditioning
3. **Check formulas:** Are all derivatives correct? Zero derivatives indicate missing parameter dependency
4. **Check initial values:** For nonlinear problems, poor initial values can lead to singular Jacobian matrices
5. **Increase RCOND:** `$fSolverRCOND = 1e-3` — more generous rank threshold

---

## $ADJ_ERR_DOF (2) — Invalid Degrees of Freedom

**Message:** `"Invalid degrees of freedom: f = <iExt>"`

**`@extended`:** Computed degree of freedom $f$ (negative).

### Degrees of Freedom Computation

| Model Type | Formula |
|-----------|--------|
| OLS/WLS/GLS | $f = n - u$ |
| LSE/WLSE/GLSE | $f = n - u + r$ |
| CLS/WCLS/GCLS | $f = n_{\text{formulas}}$ |
| GLM/WGLM/GGLM | $f = n_{\text{formulas}} + r - u$ |

with $n$ = observations, $u$ = parameters, $r$ = restrictions.

### Typical Causes

| Cause | Solution |
|---------|--------|
| Too few observations | Add more observations |
| Too many parameters | Fix parameters (`_adj_addFixedParam`) |
| Incorrect model structure | Check model type determination (features in the formula) |
| Unintended parameter detection | Function names are recognized as parameters (e.g., `sin` without parentheses) — correct the formula |

### Note

$f = 0$ (determined system) is accepted — there is no redundancy, but the system is solvable. Only $f < 0$ is reported as an error.

---

## $ADJ_ERR_RANK (3) — Rank Deficiency

**Message:** `"Rank deficiency detected"`

### Typical Causes

| Cause | Diagnostics |
|---------|-----------|
| Linear dependence between parameters | Multiple parameters have identical effects on observations |
| Redundant parameters | Parameters with no influence (derivative = 0 for all observations) |
| Datum defect | Free network without sufficient constraints |
| Collinear observations | Observation geometry does not define all parameters |

### Solutions

1. **Fix parameters:** `_adj_addFixedParam` for datum definition
2. **Add restrictions:** `_adj_addRestriction` for constraint conditions
3. **Use SVD solver:** Provides pseudo-inverse for rank deficiency
4. **Review model:** Identify and remove redundant parameters

### Note

Rank deficiency is also stored in `.rankDeficient` in the results. It is a warning — results may still be available (the solution is then not unique).

---

## $ADJ_ERR_NO_CONVERGENCE (4) — No Convergence

**Message:** `"No convergence after <iExt> iterations"`

**`@extended`:** Number of iterations performed.

### Typical Causes

| Cause | Diagnostics |
|---------|-----------|
| Poor initial values | Linearization at $x_0$ too far from the solution |
| Strong nonlinearity | GN method diverges |
| Ill-conditioned problem | Oscillation around minimum |
| LM stagnation | $\lambda$ reaches bounds ($10^{-7}$ or $10^7$) |

### Solutions

1. **Better initial values:** `_adj_setInitialValue` with more realistic approximations
2. **Use LM method:** `_adj_defaultConfig("LM")` — more robust with poor initial values
3. **Relax tolerance:** `$fTolerance = 1e-6` instead of `1e-10`
4. **Increase max. iterations:** `$iMaxIter = 500`
5. **Change differentiation method:** `"Ridder"` or `"Higham"` instead of `"Central"`
6. **Reformulate problem:** Choose a better parameterization

### Note on LM Stagnation

Levenberg-Marquardt also terminates when $\lambda$ hits its bounds AND the parameter change is below the stagnation tolerance ($10^{-10}$). This can mean:
- $\lambda \to \lambda_{\min}$: System has reached GN convergence, but $\|\Delta x\|$ slightly above tolerance
- $\lambda \to \lambda_{\max}$: No viable descent direction found

---

## $ADJ_ERR_NOT_POS_DEF (5) — Matrix Not Positive Definite

**Message:** `"Matrix not positive definite (Cholesky factorization failed)"`

**`@extended`:** LAPACK `INFO` from `dpotrf` (position of the error in the matrix).

### Where Does the Error Occur?

| Context | Affected Matrix |
|---------|------------------|
| `__adj_allocateMatrices` | Covariance matrix $\Sigma_{\ell\ell}$ (Cholesky factorization) |
| `__adj_solveGLM` | Matrix $M = B P^{-1} B^T$ (GLM normal matrix in equation space) |
| `__adj_estimateVCE` | Scaled covariance matrix after VCE update |

### Typical Causes

| Cause | Solution |
|---------|--------|
| Inconsistent covariances | Check whether $\Sigma_{\ell\ell}$ is symmetric and positive definite |
| Negative variances from VCE | VCE iterations can produce negative variance factors (min. $10^{-6}$) |
| Numerical problems | Very large condition number of the covariance matrix |
| GLM: $B$ matrix deficit | Check observation references in formulas |

### Solutions

1. **Check covariances:** Is the input covariance matrix physically plausible?
2. **Disable VCE:** `_adj_defaultConfig("GN", False)` if VCE is unstable
3. **Simplify model:** Fewer covariances, diagonal instead of generalized model
4. **Convert covariances to correlations:** $\sigma_{ij} < \sigma_i \cdot \sigma_j$ (Cauchy-Schwarz inequality)

---

## $ADJ_ERR_INPUT (6) — Invalid Input

**Message:** Depends on `@extended`.

### @extended Values

| @extended | Message | Triggering Function | Description |
|-----------|---------|--------------------|----|
| 0 | `"Invalid input"` | various | Generally invalid parameter |
| 2 | `"First observation not found in system"` | `_adj_addCovariance` | Observation `$sObs1` does not exist |
| 3 | `"Second observation not found in system"` | `_adj_addCovariance` | Observation `$sObs2` does not exist |
| 4 | `"Cross-group covariance not supported"` | `_adj_addCovariance` | Covariance between different VCE groups |
| 5 | `"Cross-group covariance in system"` | `_adj_solve` | System contains cross-group covariances |
| 7 | (internal) | `__adj_robustIRLS` | Generalized model + covariance matrix + robust estimation not supported |
| 10 | `"Observation value is not a finite number"` | `_adj_addObs`, `_adj_addObsFunction` | Value is NaN, Inf, or not a numeric value |
| 11 | `"Standard deviation must be a positive finite number"` | `_adj_addObs`, `_adj_addObsFunction` | $\sigma \le 0$, NaN, or Inf |

### Detailed Error Causes and Solutions

#### @extended = 0 — General Input Error

**Possible causes:**
- System map not initialized (for functions that expect an existing map)
- Parameter not found in `_adj_setInitialValue`
- Invalid solver type in `_adj_solve` (neither `"QR"` nor `"SVD"`)
- No formulas and no restrictions in the system (in `__adj_validateInput`)
- Unknown estimator name in `_adj_robustDefaults`
- Unknown method in `_adj_getOutliers`

**Solution:** Check input parameters. Add observations and formulas before `_adj_solve`.

---

#### @extended = 2/3 — Observation Not Found

**Causes:**
- Observation name misspelled (case is internally converted to uppercase)
- Observation has not been added yet (order: first `_adj_addObs`, then `_adj_addCovariance`)
- Observation was removed (`_adj_removeObs`)

**Solution:** Check observation names and ensure that the observation is added before the covariance.

---

#### @extended = 4 — Cross-Group Covariance

**Cause:** `_adj_addCovariance` is called with two observations from different VCE groups.

**Background:** The VCE (Helmert method) works with a block-diagonal weight matrix. Off-diagonal covariances between groups violate this assumption.

**Solution:**
- Move observations to the same VCE group (same `$sVarComp`)
- Disable VCE: `_adj_defaultConfig("GN", False)`
- Reconsider the covariance structure

**Note:** The covariance is stored nonetheless and the flag `_hasCrossGroupCovar` is set. `_adj_solve` then rejects the system with `@extended = 5`.

---

#### @extended = 5 — Cross-Group Covariances in System

**Cause:** `_adj_solve` detects that the system contains cross-group covariances (set by `_adj_addCovariance` with `@extended = 4`).

**Solution:** See `@extended = 4`.

---

#### @extended = 7 — Robust Estimation + Generalized Model

**Cause:** Robust estimation (IRLS) is configured, but the system has a full covariance matrix (generalized model GLS/GLSE/GGLM/GCLS).

**Background:** IRLS reweighting modifies the diagonal of the weight matrix. With a full covariance matrix, correct reweighting is non-trivial and currently not implemented.

**Solution:**
- Disable robust estimation (`.robust = ""`)
- Use only weights instead of covariances (diagonal model WLS/WLSE/WGLM/WCLS)

---

#### @extended = 10 — Invalid Observation Value

**Cause:** The observation value `$fValue` is:
- Not a numeric type (`IsNumber` = False)
- NaN ($x \ne x$)
- $+\infty$ or $-\infty$

**Solution:** Ensure that all observation values are finite numbers. Typical error: passing a string instead of a number.

---

#### @extended = 11 — Invalid Standard Deviation

**Cause:** The standard deviation `$fStdDev` is:
- Not a numeric type
- $\le 0$ (must be strictly positive)
- NaN
- $\infty$

**Solution:** Ensure that all standard deviations are positive and finite. A standard deviation of $= 0$ is physically meaningless (a perfect observation does not exist).

---

## $ADJ_ERR_RESTR_OVERDETERMINED (7) — Restrictions Fully Determine Parameters

**Message:** `"Restrictions fully determine all parameters (nRestrictions=<iExt> >= nParams)"`

**`@extended`:** Number of restrictions $r$.

### Description

The number of equality restrictions is greater than or equal to the number of free (adjustable) parameters. In the null-space projection used by LSE and GLM solvers, this leads to $n_{\text{free}} = u - r \le 0$, meaning there are no remaining degrees of freedom in the parameter space. All parameter values are fully determined by the restrictions alone — the observations only contribute residuals.

### Typical Causes

| Cause | Solution |
|---------|--------|
| Too many restrictions for the number of free parameters | Remove restrictions or add more parameters |
| Fixed parameters reduce $u$ below $r$ | Unfix parameters or remove restrictions that become redundant |
| Restriction + fixed parameter together fix the same parameter | Use only one mechanism (restriction OR fixed parameter) |

### Example

```autoit
; Y is fixed at 55, so only X is free (nParams = 1)
_adj_addFixedParam($mSystem, "Y", 55.0)
; Restriction X + Y = 120 fully determines X = 65 (nRestrictions = 1 >= nParams = 1)
_adj_addRestriction($mSystem, "X + Y", 120)
; → Error 7: no free parameters left for adjustment
```

### Solutions

1. **Remove redundant restrictions:** If a fixed parameter already determines the value, the restriction is superfluous
2. **Replace fixed parameter with restriction:** Use `_adj_addRestriction("Y", 55)` instead of `_adj_addFixedParam` — then both X and Y remain adjustable under two restrictions
3. **Add more parameters:** Ensure $u > r$ for the adjustment to have free parameters

---

## Error Code Query

```autoit
_adj_solve($mSystem)
Local $iErr = @error, $iExt = @extended

Switch $iErr
    Case $ADJ_ERR_OK
        ; Success — evaluate results
        ConsoleWrite(_adj_displayResults($mSystem))

    Case $ADJ_ERR_NO_CONVERGENCE
        ; Improve initial values or use LM
        ConsoleWrite("No convergence after " & $iExt & " iterations" & @CRLF)

    Case $ADJ_ERR_INPUT
        ; Evaluate input error in detail
        ConsoleWrite("Input error: " & _adj_getErrorMessage($iErr, $iExt) & @CRLF)

    Case Else
        ConsoleWrite(_adj_getErrorMessage($iErr, $iExt) & @CRLF)
EndSwitch
```
