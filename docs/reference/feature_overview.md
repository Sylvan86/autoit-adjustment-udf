# Feature Overview — Adjustment UDF

Complete overview of all capabilities of the Adjustment UDF for geodetic least-squares adjustment.

---

## Model Types

| Model Type | Description | Weighted | Generalized |
|-----------|-------------|-----------|---------------|
| OLS | Ordinary Least Squares | WLS | GLS |
| LSE | Least Squares with parameter constraints | WLSE | GLSE |
| CLS | Condition adjustment (observation constraints) | WCLS | GCLS |
| GLM | Gauss-Helmert model (general observation equations) | WGLM | GGLM |

### Automatic Model Detection

The UDF determines the model type automatically based on the problem structure:

| Criterion | Result |
|-----------|----------|
| Formulas without `#` observation references, without restrictions | OLS |
| Formulas with exactly 1 observation per formula (`_adj_addObsFunction`) | OLS (observation equations) |
| Formulas with >1 observation per formula or `_adj_addFunction` with `#` references | GLM |
| Formulas with only `#` references, no parameters | CLS |
| Additional restrictions (`_adj_addRestriction`) | LSE variant |
| At least one `stdDev != 1` | W variant (weighted) |
| Covariances via `_adj_addCovariance` (off-diagonal) | G variant (generalized) |

The detection is performed in `__adj_classifyModel` using the internal flags `_hasObsFunctions`, `_hasObsRestrictions`, `_hasParamRestrictions`, `_isWeighted` and `_hasCovar`.

---

## Solvers

| Solver | LAPACK Routine | Description |
|--------|---------------|-------------|
| QR | DGELSY | Rank-revealing QR decomposition with column pivoting (default) |
| SVD | DGELSD | Singular value decomposition, more robust for ill-conditioned systems |

Configuration via `_adj_defaultConfig()`:
- `.solver` = `"QR"` (default) or `"SVD"`
- `.solverRCOND` = rank threshold (default: `1e-5`)

---

## Iteration Methods

| Method | Description |
|-----------|-------------|
| Gauss-Newton (GN) | Linearization + normal equations, quadratic convergence near the solution |
| Levenberg-Marquardt (LM) | Damped GN with Nielsen strategy, more robust with poor initial values |

Configuration: `_adj_defaultConfig("GN")` or `_adj_defaultConfig("LM")`

---

## Statistical Quantities

### Always Computed

| Quantity | Symbol | Description |
|-----------|--------|-------------|
| A posteriori variance factor | $s_0^2 = v^T P v / f$ | Measure of model quality |
| Degrees of freedom | $f$ | Model-dependent |
| Weighted residual sum of squares | $v^T P v$ | Sum of weighted squared residuals |

### Compute-on-Demand (`.compute` flags)

| Flag | Quantities | Default |
|------|-----------|---------|
| `.qxx` | Cofactor matrix $Q_{xx}$, standard deviations $sd_x$ | `True` |
| `.cofactors` | Cofactor matrices $Q_{vv}$, $Q_{\hat{y}}$; standard deviations $sd_v$, $sd_{\hat{y}}$ | `False` |
| `.redundancy` | Redundancy numbers $r_i = \text{diag}(P \cdot Q_{vv})$ | `False` |
| `.globalTest` | Chi-squared global test ($H_0: \sigma_0^2 = 1$) | `False` |
| `.diagnostics` | Baarda $|w|$, Pope $T_\tau$, p-values, gross errors $\hat{\nabla}$, MDB | `False` |

Computation is triggered automatically when the result display (`_adj_displayResults`) requests the corresponding columns.

### Dependency Chain

```
diagnostics → redundancy → cofactors → qxx
```

Each level is computed only once; subsequent calls use cached results.

---

## Diagnostic Tests

| Test | Statistic | Basis | Description |
|------|-----------|-------|-------------|
| Global test | $T = v^T P v$ | $\chi^2(f)$ | Tests whether $\sigma_0^2 = 1$ (null hypothesis) |
| Baarda test | $\|w_i\| = \|v_i\| / (\sigma_i \sqrt{r_i})$ | $N(0,1)$, a priori $\sigma_0 = 1$ | Outlier test for individual observations |
| Pope test | $T_{\tau,i} = \|w_i\| / \hat{s}_0$ | $t(f)$, a posteriori $\hat{s}_0$ | More robust when $\sigma_0$ is unknown |

Both tests provide:
- **p-values** (Baarda: normal distribution, Pope: t-distribution)
- **Gross errors** $\hat{\nabla}_i = v_i / r_i$
- **MDB** (Minimal Detectable Bias): $\text{MDB}_i = \sigma_i \sqrt{\lambda_0 / r_i}$ with $\lambda_0 = (z_{\alpha/2} + z_\beta)^2$
- **Test decision**: `"outlier"`, `"suspect"`, `"ok"`

Configuration in `.diagnostics`:
- `.alpha` = significance level (default: `0.001`)
- `.beta` = type II error for MDB (default: `0.20`)
- `.alphaSuspect` = threshold for `"suspect"` (default: `10 * alpha`)
- `.testBasis` = `"pope"` (default) or `"baarda"`

---

## Variance Component Estimation (VCE)

| Property | Value |
|-------------|------|
| Method | Helmert |
| Groups | Via `$sVarComp` parameter in `_adj_addObs` / `_adj_addObsFunction` |
| Iteration | Until convergence or max. 15 iterations |
| Convergence | Relative change of variance factors < 5% |
| Generalized | Block scaling of the covariance matrix $\Sigma_{\ell\ell}$ |

Activation: `_adj_defaultConfig("GN", True)` (default: active).

---

## Robust Estimation (IRLS)

| Estimator | Weight function $w(u)$ | Tuning parameter | Breakdown point |
|----------|------------------------|------------------|------------|
| L1 | $w = 1 / \max(\|u\|, \varepsilon)$, capped at `maxWeight` | `maxWeight` = 1000 | 50% |
| Huber | $\min(1,\; c/\|u\|)$ | $c = 1.345$ | ~5% |
| Hampel | Three-zone: $1 \to a/\|u\| \to$ linear $\to 0$ | $a=1.7,\; b=3.4,\; c=8.5$ | ~25% |
| Biweight (Tukey) | $(1 - (u/c)^2)^2$ for $\|u\| \le c$, else 0 | $c = 4.685$ | 50% |
| BIBER (Schweppe) | Huber with $\sqrt{r_i}$ correction | $c = 3.5$ | Leverage-dependent |
| ModifiedM (Koch) | Huber with $\sqrt{r_i}$ correction | $c = 1.5$ | Leverage-dependent |

### Scale Parameter

| Method | Description |
|---------|-------------|
| `"MAD"` | $\hat{\sigma} = \text{median}(\|v_i / \sigma_i\|) / 0.6745$ (default) |
| `"s0"` | $\hat{\sigma} = \sqrt{v^T P v / f}$ |
| `"apriori"` | $\hat{\sigma} = 1$ (trusts the stochastic model) |
| `"fixed"` | User-defined fixed value (`.fixedScale`) |

### IRLS Workflow

1. Initial (non-robust) solution
2. Redundancy computation (only for BIBER/ModifiedM)
3. Storage of original weights, upgrade to weighted model variant
4. Iteration: scale estimation $\to$ weight update $\to$ convergence check $\to$ recomputation
5. Convergence: max. relative weight change $< 10^{-3}$ or max. 30 iterations

### Limitations

- Generalized models with full covariance matrix are not yet supported (`@error = $ADJ_ERR_INPUT`, `@extended = 7`)

---

## Numerical Differentiation

| Method | Order | Description |
|---------|---------|-------------|
| `"Central"` / `"Central1"` | $O(h^2)$ | Central difference, 2 function evaluations |
| `"Central2"` | $O(h^4)$ | 4-point central difference |
| `"Central3"` | $O(h^6)$ | 6-point central difference |
| `"Central4"` | $O(h^8)$ | 8-point central difference |
| `"Forward"` | $O(h)$ | Forward difference |
| `"Backward"` | $O(h)$ | Backward difference |
| `"Ridder"` | adaptive | Extrapolation with step-size halving |
| `"Higham"` | adaptive | Iterative refinement to target accuracy |

Configuration: `.deriveMethod` in `_adj_defaultConfig()` (default: `"Central"`).

Step size $h$ is automatically scaled:
- Central differences: $h = \varepsilon^{1/3} \cdot \max(1, |x|)$
- Forward/Backward: $h = \sqrt{\varepsilon} \cdot \max(1, |x|)$

### Analytical Derivatives

Optionally as a map or pipe string:
```autoit
; As map
Local $mDeriv[]
$mDeriv["X"] = "2*X"
$mDeriv["Y"] = "2*Y"

; As pipe string
"X=2*X | Y=2*Y"
```

---

## Additional Features

| Feature | Description |
|---------|-------------|
| Symbolic formula input | Formulas as strings, automatic parameter detection via regex |
| Jacobian equilibration | Automatic column scaling for improved conditioning (`.scaling = True`) |
| Compute-on-demand | Statistics are computed only when requested |
| Configurable result display | Selectable columns for parameter and observation tables |
| Observation removal | `_adj_removeObs` removes observations including formulas and covariances |
| Outlier detection | `_adj_getOutliers` with 4 methods: `robustWeight`, `absU`, `baarda`, `pope` |
