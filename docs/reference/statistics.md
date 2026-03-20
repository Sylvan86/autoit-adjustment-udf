# Statistical Measures

Complete documentation of all statistical measures computed by the Adjustment module: degrees of freedom, variance estimation, cofactor matrices, redundancy numbers, global test, and outlier diagnostics.

Related documentation: [Result Structure](results.md) | [Configuration](configuration.md) | [Model Types](model_types.md) | [Robust Estimation](robust_estimation.md)


## Degrees of Freedom

The degrees of freedom f (redundancy) describe the over-determination of the equation system. Their computation depends on the [model type](model_types.md).

| Model Type | Formula | Variables |
|------------|---------|-----------|
| OLS/WLS/GLS | f = n - u | n = observations, u = parameters |
| LSE/WLSE/GLSE | f = n - u + r | r = restrictions |
| CLS/WCLS/GCLS | f = m | m = condition equations |
| GLM/WGLM/GGLM | f = m + r - u | m = functions, r = restrictions, u = parameters |

**Prerequisite:** f > 0 (otherwise `$ADJ_ERR_DOF`). When f = 0, s_0 cannot be computed.


## A Posteriori Variance of Unit Weight

### Basic Formula

The a posteriori variance of unit weight is the central quality indicator:

    s_0_hat^2 = v^T P v / f

The a posteriori standard deviation of unit weight:

    s_0_hat = sqrt(v^T P v / f)

### Computation of v^T P v

The weighted sum of squares is computed in a model-dependent manner:

**Diagonal weight matrix (WLS/WLSE/WCLS/WGLM):**

    v^T P v = sum_i (v_i / sigma_i)^2 = ||P^(1/2) v||^2

**Full covariance matrix (GLS/GLSE/GCLS/GGLM):**

    v^T P v = ||L^(-1) v||^2

where Sigma_ll = L L^T (Cholesky decomposition) and P = Sigma_ll^(-1).

**Unweighted (OLS):**

    v^T P v = v^T v = sum_i v_i^2

### Interpretation

| s_0_hat | Interpretation |
|---------|----------------|
| close to 1 | Stochastic model fits the data well |
| significantly > 1 | Observation accuracies too optimistic or model errors |
| significantly < 1 | Observation accuracies too pessimistic |


## Cofactor Matrix Q_xx

The cofactor matrix Q_xx of the adjusted parameters is computed via two paths:

### QR Path (Default)

After DGELSY, the upper triangle of the system matrix contains the R factor of the QR decomposition with column pivoting: A P = Q R.

    Q_xx = P R^(-1) R^(-T) P^T

Steps:
1. Extract R (upper n x n triangle)
2. Compute R^(-1) (DTRTRI, in-place)
3. Reverse column permutation: S = P R^(-1)
4. Reverse Jacobi equilibration (if active): S_scaled = diag(1/s_eq) S
5. Q_xx = S S^T (symmetric rank-k update, DSYRK)

### SVD Path

When using DGELSD, the SVD A = U Sigma V^T is computed:

    Q_xx = V diag(1/sigma_i^2) V^T

with rank truncation: singular values sigma_i < RCOND * sigma_max are treated as zero. The condition number sigma_max / sigma_min is stored in `.conditionNumber`.

### LSE/GLM with Restrictions

For systems with restrictions, transformation via the null space Q_2 of the restriction matrix is used:

    Q_xx = Q_2 S S^T Q_2^T = T T^T      with T = Q_2 S

The result is a u x u matrix (full parameter dimension), even though the reduced system has only (u - r) degrees of freedom.

### Standard Deviations of Parameters

    sd(x_i) = s_0_hat * sqrt(q_xx_ii)

Covariance between parameters:

    Sigma_xx = s_0_hat^2 * Q_xx


## Cofactor Matrices Q_vv and Q_yhat

Require compute flag: `$mConfig.compute.cofactors = True`

### Relationship

    P^(-1) = Q_yhat + Q_vv

where P^(-1) = Sigma_ll (covariance matrix of observations) for generalized models, and P^(-1) = diag(sigma_i^2) for diagonal weighting.

### OLS/LSE Computation

    Q_yhat = A Q_xx A^T = (A S)(A S)^T

    Q_vv = P^(-1) - Q_yhat

where A is the (possibly transformed) Jacobian matrix and S is the auxiliary matrix from the Q_xx computation.

### GLM Computation

In the equation space via the hat matrix:

    H_eq = A_tilde Q_xx A_tilde^T

    Q_vv = T^T (I - H_eq) T

where T = L^(-1) B and L is the Cholesky factor of M = B Sigma_ll B^T.

### CLS Computation

Direct computation without Q_xx:

    Q_vv = B^T (B B^T)^(-1) B

via Cholesky decomposition of N = B B^T.

### Back-Transformation

For weighted models, back-transformation from the transformed to the original space is performed:

- **Diagonal (WLS):** Q = diag(sigma) Q_w diag(sigma)
- **Generalized (GLS):** Q = L Q_w L^T

### Standard Deviations

    sd(v_i)     = s_0_hat * sqrt(q_vv_ii)
    sd(y_hat_i) = s_0_hat * sqrt(q_yy_ii)


## Redundancy Numbers

Require compute flag: `$mConfig.compute.redundancy = True`

### Definition

    r_i = (P Q_vv)_ii = p_i * q_vv_ii

For diagonal weighting:

    r_i = q_vv_ii / sigma_i^2

For full covariance matrix:

    r_i = (Sigma_ll^(-1) Q_vv)_ii

via Cholesky: diag(L^(-1) Q_vv L^(-T))

### Properties

| Property | Meaning |
|----------|---------|
| 0 <= r_i <= 1 | Always in the unit interval |
| sum_i r_i = f | Sum equals degrees of freedom |
| r_i close to 1 | Observation well controlled (independently determined) |
| r_i close to 0 | Observation poorly controlled (leverage effect) |

Observations with r_i < 10^(-12) are skipped in diagnostics ("---"), as no meaningful test statistics can be computed.


## Global Test (chi2 Test)

Requires compute flag: `$mConfig.compute.globalTest = True` or display request.

### Hypothesis

    H_0: sigma_0^2 = 1      (stochastic model correct)
    H_1: sigma_0^2 != 1

### Test Statistic

    T = v^T P v

Under H_0: T ~ chi2(f)

### Test Decision

    H_0 not rejected if:   chi2_(alpha/2, f) <= T <= chi2_(1 - alpha/2, f)

The significance level alpha comes from `$mConfig.diagnostics.alpha` (default: 0.001). Quantiles are computed via the inverse regularized incomplete gamma function.

### Configuration

```autoit
Local $mConfig = _adj_defaultConfig()
Local $mDiag = $mConfig.diagnostics
$mDiag.alpha = 0.05    ; Signifikanzniveau für Globaltest und Diagnostik
$mConfig.diagnostics = $mDiag
```


## Outlier Diagnostics

Requires compute flag: `$mConfig.compute.diagnostics = True` or display request via diagnostics columns.

### Baarda Test Statistic

Normalized residual (a priori sigma_0 = 1):

    |w_i| = |v_i| / (sigma_i * sqrt(r_i))

Under H_0: w_i ~ N(0, 1)

p-value: p_i = 2 * (1 - Phi(|w_i|))

### Pope Test Statistic

A posteriori variant with estimated s_0_hat:

    T_tau_i = |w_i| / s_0_hat

Under H_0 approximately: T_tau_i ~ t(f)

p-value: p_i = 2 * (1 - F_t(|T_tau_i|, f))

### Test Decision

Based on the chosen `testBasis` (default: `"pope"`):

| p-value | Decision |
|---------|----------|
| p < alpha | `"outlier"` |
| alpha <= p < alpha_suspect | `"suspect"` |
| p >= alpha_suspect | `"ok"` |

where alpha_suspect = 10 * alpha (default: 0.01) when not explicitly set.

### Gross Error Estimate

    nabla_hat_i = v_i / r_i

Gives the estimated magnitude of a single gross error in observation i.

### Minimum Detectable Bias (MDB)

    MDB_i = sigma_i * sqrt(lambda_0 / r_i)

where lambda_0 = (z_(alpha/2) + z_beta)^2 is the non-centrality parameter.

| Parameter | Default | Description |
|-----------|---------|-------------|
| alpha | 0.001 | Significance level (Type I error) |
| beta | 0.20 | Type II error (power = 1 - beta = 0.80) |

The MDB indicates the minimum size a gross error must have in order to be detected with the chosen power.

### Baarda Warning

When `testBasis = "baarda"` and s_0_hat deviates significantly from 1 (s_0_hat > 3.0 or s_0_hat < 0.3), a warning is issued. In this case, the Baarda test is unreliable, and `testBasis = "pope"` is recommended.

### Configuration Example

```autoit
Local $mConfig = _adj_defaultConfig()
Local $mCompute = $mConfig.compute
$mCompute.diagnostics = True
$mConfig.compute = $mCompute

; Diagnostik-Parameter anpassen
Local $mDiag = $mConfig.diagnostics
$mDiag.alpha        = 0.01     ; 1% statt 0.1%
$mDiag.beta         = 0.10     ; 90% Treffsicherheit
$mDiag.testBasis    = "baarda" ; a priori Test
$mDiag.alphaSuspect = 0.05     ; Suspect-Schwelle
$mConfig.diagnostics = $mDiag
```


## Compute-on-Demand

Statistical measures are not automatically computed after every adjustment, but only when needed. The `__adj_ensureComputed` mechanism resolves dependencies automatically:

### Dependency Chain

```
diagnostics  -->  redundancy  -->  cofactors  -->  qxx
     |                |               |              |
  Baarda/Pope    diag(P*Qvv)     Qvv, Qyhat,    Qxx, sdx
  nabla, MDB                     sdv, sdy
  testDecision
```

### Triggering

1. **Explicit:** Set compute flags in the [Configuration](configuration.md#2-compute-flags-compute-sub-map)
2. **Implicit:** Request columns in the [Display Configuration](configuration.md#4-display-config-adj_defaultdisplayconfig)

Example of implicit triggering:

```autoit
; Solve ohne explizite Compute-Flags
_adj_solve($mSystem)

; Display fordert Redundanz-Spalte an --> berechnet automatisch qxx + cofactors + redundancy
Local $mDisplay = _adj_defaultDisplayConfig()
$mDisplay.obsCols = "name|value|v|r"
ConsoleWrite(_adj_displayResults($mSystem, $mDisplay))
```

### Column Dependencies in the Display Config

| Column | Triggers |
|--------|----------|
| `sdv`, `sdyhat` | cofactors |
| `r` | redundancy |
| `w`, `T`, `p`, `pPope`, `blunder`, `mdb`, `decision` | diagnostics |
| `showGlobalTest = True` | globalTest |

### VCE Special Case

When Variance Component Estimation is active (`$mConfig.vce = True`), `cofactors` and `redundancy` are automatically enforced, since VCE requires these matrices for variance component computation.


## Mathematical Notes

### Weight Matrix P

The weight matrix is not formed explicitly. Instead:

- **Diagonal:** P = diag(1/sigma_i^2), stored as vector `Vector_ObsInvStdDev` (= 1/sigma_i)
- **Full:** Sigma_ll = L L^T (Cholesky), P = L^(-T) L^(-1), whitening via forward/backward substitution

### GLM/CLS Formulas

For the Gauss-Helmert and condition models, the formulas for Q_xx, Q_vv, etc. are more complex than for OLS. The basic structure remains the same, but the system matrices (A_tilde, B) and transformations (M = B Sigma_ll B^T) differ. Details in [Model Types](model_types.md).

### Numerical Stability

- **Jacobi Equilibration:** Column scaling of the Jacobian matrix for improved conditioning. Applied before QR decomposition and reversed during Q_xx computation.
- **SYRK instead of explicit matrix product:** Q_xx = S S^T is computed via BLAS DSYRK (symmetric rank-k update), which is numerically more stable than S * S^T via DGEMM.
- **Rank truncation (SVD):** Singular values below RCOND * sigma_max are treated as zero. Default: RCOND = 10^(-5).
