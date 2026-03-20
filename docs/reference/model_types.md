# Model Types in Detail — Adjustment UDF

The UDF supports four basic model types, each in three weighting variants (unweighted, weighted, generalized). The model type is determined automatically from the problem structure.

---

## Overview

| Basic Type | Unweighted | Weighted ($P$ diagonal) | Generalized ($P = \Sigma^{-1}$, full) |
|----------|-------------|--------------------------|----------------------------------------|
| Observation adjustment | OLS | WLS | GLS |
| With parameter constraints | LSE | WLSE | GLSE |
| Condition adjustment | CLS | WCLS | GCLS |
| Gauss-Helmert model | GLM | WGLM | GGLM |

---

## OLS / WLS / GLS — Observation Adjustment

### Mathematical Model

Functional model (observation equations):

$$\ell + v = f(x)$$

- $\ell$ — observation vector ($n \times 1$)
- $v$ — residual vector ($n \times 1$)
- $x$ — parameter vector ($u \times 1$)
- $f$ — model function (linear or nonlinear)

Linearized (with approximate values $x_0$):

$$A \cdot \Delta x + w = v$$

with Jacobian matrix $A = \partial f / \partial x |_{x_0}$ and misclosure vector $w = f(x_0) - \ell$.

### Stochastic Model

| Variant | Weight matrix $P$ | Covariance matrix |
|----------|-------------------|----------------|
| OLS | $P = I$ | $\Sigma_{\ell\ell} = \sigma_0^2 \cdot I$ |
| WLS | $P = \text{diag}(1/\sigma_i^2)$ | $\Sigma_{\ell\ell} = \sigma_0^2 \cdot \text{diag}(\sigma_i^2)$ |
| GLS | $P = \Sigma_{\ell\ell}^{-1}$ | $\Sigma_{\ell\ell}$ fully populated |

### Normal Equations

$$\hat{x} = (A^T P A)^{-1} \cdot A^T P \ell$$

Solved via DGELSY (QR) or DGELSD (SVD) after whitening transformation.

### Degrees of Freedom

$$f = n - u$$

### When to Use

- **Standard adjustment problems**: Overdetermined observations, parameter estimation
- Examples: Height network adjustment, regression, distance adjustment

### API Mapping

```autoit
; Each observation + formula with exactly 1 observation per formula
_adj_addObsFunction($mSystem, "H1", "HA + dH_AB", 100.5, 0.003)
_adj_addObsFunction($mSystem, "H2", "HA + dH_AC", 101.2, 0.004)
```

Weighted variant: `$fStdDev != 1`. Generalized: additionally `_adj_addCovariance`.

---

## LSE / WLSE / GLSE — With Parameter Constraints

### Mathematical Model

Functional model with equality constraints:

$$\ell + v = f(x) \quad \text{subject to} \quad g(x) = c$$

- $g(x) = c$ — constraint equations ($r$ equations), parameters only (no observations)

Linearized:

$$A \cdot \Delta x + w = v \quad \text{subject to} \quad R_A \cdot \Delta x = w_R$$

### Solution

Extended normal equations with Lagrange multipliers:

$$\begin{pmatrix} A^T P A & R_A^T \\ R_A & 0 \end{pmatrix} \begin{pmatrix} \Delta x \\ \lambda \end{pmatrix} = \begin{pmatrix} A^T P \ell \\ w_R \end{pmatrix}$$

The UDF solves the equivalent null-space projected problem:
1. QR decomposition of $R_A^T$: $R_A^T \cdot P_{col} = Q \cdot \begin{pmatrix} R \\ 0 \end{pmatrix}$
2. Null space $Q_2$ (last $u - r$ columns of $Q$)
3. Reduced system: $(Q_2^T A^T P A Q_2) \cdot z = Q_2^T A^T P \ell$
4. Back-transformation: $\Delta x = Q_2 \cdot z + x_R$

### Degrees of Freedom

$$f = n - u + r$$

### When to Use

- **Known relationships between parameters**: e.g., fixed points, sum conditions
- The constraints are satisfied exactly (no residuals on constraints)

### API Mapping

```autoit
; OLS setup as above, additionally:
_adj_addRestriction($mSystem, "HA + HB - 200")        ; Parameter constraint
_adj_addFixedParam($mSystem, "HA", 100.0)              ; Fix parameter
```

Note: `_adj_addRestriction("X")` is automatically converted to `_adj_addFixedParam("X", 0)`.

---

## CLS / WCLS / GCLS — Condition Adjustment

### Mathematical Model

Condition equations (observations only, no explicit parameters):

$$g(\ell + v) = c$$

Linearized:

$$B \cdot v + w = 0$$

with $B = \partial g / \partial \ell$.

### Solution

$$v = -P^{-1} B^T (B P^{-1} B^T)^{-1} \cdot w$$

### Degrees of Freedom

$$f = n_{\text{formulas}}$$

### When to Use

- **Conditions between observations**: e.g., closing triangle, angle closure condition
- No unknown parameters — only residuals on observations

### API Mapping

```autoit
; Observations without formulas
_adj_addObs($mSystem, "W1", 60.5, 0.01)
_adj_addObs($mSystem, "W2", 59.3, 0.01)
_adj_addObs($mSystem, "W3", 60.2, 0.01)

; Condition equation: observations referenced with #
_adj_addFunction($mSystem, "#W1 + #W2 + #W3", 180.0)
```

---

## GLM / WGLM / GGLM — Gauss-Helmert Model

### Mathematical Model

General implicit observation equations:

$$g(x, \ell + v) = 0$$

This is the most general linear model and contains OLS and CLS as special cases.

Linearized:

$$A \cdot \Delta x + B \cdot v + w = 0$$

with:
- $A = \partial g / \partial x$ (Jacobian with respect to parameters)
- $B = \partial g / \partial \ell$ (Jacobian with respect to observations)

### Solution

Weighted normal matrix in the equation space:

$$M = B \cdot P^{-1} \cdot B^T$$

Cholesky factorization: $M = L \cdot L^T$

Transformed system:

$$\tilde{A} \cdot \Delta x = \tilde{w}$$

with $\tilde{A} = L^{-1} A$ and $\tilde{w} = -L^{-1} w$.

Residuals:

$$v = P^{-1} B^T M^{-1} (w + A \cdot \Delta x)$$

### Degrees of Freedom

$$f = n_{\text{formulas}} + r - u$$

### When to Use

- **Observations in nonlinear relationships**: e.g., circle fitting, plane fitting
- Multiple observations per equation
- Combination of observation and parameter conditions in one formula

### API Mapping

```autoit
; Multiple observations per formula → GLM
_adj_addObsFunction($mSystem, "D1", "(#D1*cos(#T1)-XM)^2+(#D1*sin(#T1)-YM)^2-R^2", 0, 0.01)

; Or: _adj_addFunction with # and parameters mixed
_adj_addFunction($mSystem, "(#D1*cos(#T1)-XM)^2+(#D1*sin(#T1)-YM)^2-R^2")
```

For GLM with constraints, these are handled via null-space projection (as with LSE).

---

## Model Type Determination Logic

The internal classification (`__adj_classifyModel`) uses the following decision logic:

```
Does the system have covariances (_hasCovar)?
├─ Yes → Generalized (G prefix)
│   ├─ Observation references in formulas? → GCLS (without parameters) or GGLM (with parameters)
│   └─ Parameter-only formulas? → GLSE (with constraints) or GLS (without)
├─ Weighted (_isWeighted)?
│   ├─ Observation references? → WCLS / WGLM
│   └─ Parameters only? → WLSE / WLS
└─ Unweighted
    ├─ Observation references? → CLS / GLM
    └─ Parameters only? → LSE / OLS
```

The flags are set automatically when adding observations and formulas:

| Flag | Set by |
|------|--------------|
| `_isWeighted` | `stdDev != 1` in `_adj_addObs` / `_adj_addObsFunction` |
| `_hasCovar` | Off-diagonal covariance in `_adj_addCovariance` |
| `_hasObsFunctions` | `#` references in formulas of `_adj_addObsFunction` / `_adj_addFunction` |
| `_hasObsRestrictions` | `#` references in restrictions |
| `_hasParamFunctions` | Parameter variables in formulas |
| `_hasParamRestrictions` | Parameter variables in restrictions |

---

## Whitening Transformation

For weighted and generalized models, the UDF transforms the equation system into the "white" space (unit covariance):

| Variant | Transformation |
|----------|---------------|
| W variants | Scaling: $\tilde{A} = \text{diag}(1/\sigma_i) \cdot A$, $\tilde{w} = \text{diag}(1/\sigma_i) \cdot w$ |
| G variants | Cholesky: $\Sigma_{\ell\ell} = L L^T$, then $\tilde{A} = L^{-1} A$, $\tilde{w} = L^{-1} w$ |

The back-transformation of the cofactor matrices is performed automatically:
- $Q_{vv} = \text{diag}(\sigma) \cdot \tilde{Q}_{vv} \cdot \text{diag}(\sigma)$ (diagonal)
- $Q_{vv} = L \cdot \tilde{Q}_{vv} \cdot L^T$ (generalized)
