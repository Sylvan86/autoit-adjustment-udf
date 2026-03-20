# Solvers and Iteration Methods — Adjustment UDF

Detailed description of the linear solvers, iteration methods, and convergence control.

---

## Architecture

The solver architecture consists of three levels:

```
Level 1:  _adj_solve()              ─ Entry point, model preparation
            ├── Phase 1: __adj_robustIRLS()        ─ Robust estimation (optional)
            └── Phase 2: __adj_estimateVCE()       ─ VCE loop

Level 2:  __adj_estimateVCE()       ─ Variance component estimation (outer loop)
          __adj_solveNonlinear()     ─ GN/LM iteration loop

Level 3:  __adj_solveIteration()    ─ One linear iteration step
            ├── __adj_computeJacobians()     ─ Jacobian matrices
            ├── __adj_computeContradiction() ─ Misclosure vector
            ├── __adj_applyWhitening()       ─ Weight transformation
            ├── __adj_applyEquilibration()   ─ Column scaling
            ├── __adj_dispatchLinearSolver() ─ LAPACK dispatch
            ├── __adj_computeResiduals()     ─ Residuals
            ├── __adj_updateParameters()     ─ Parameter update
            ├── __adj_updateObservations()   ─ Observation update (GLM)
            └── __adj_updateLMDamping()      ─ LM damping adjustment
```

---

## Linear Solvers

### QR Decomposition (DGELSY) — Default

**LAPACK routine:** `DGELSY` (rank-revealing QR decomposition with column pivoting)

**Principle:**
$$A \cdot P_\pi = Q \cdot R$$

with column pivoting $P_\pi$ for numerical stability and rank determination.

**Properties:**
- Column pivoting reveals rank deficiencies
- Efficient: $O(mn^2)$ for $m \times n$ matrix
- Numerically stable for well-conditioned systems
- `RCOND` parameter controls the rank threshold (default: $10^{-5}$)

**Rank determination:** Diagonal elements of $R$ are tested against `RCOND * max(|R_{ii}|)`. Small diagonal elements indicate rank deficiency.

**Cofactor matrix (Qxx):**
$$Q_{xx} = P_\pi \cdot R^{-1} \cdot R^{-T} \cdot P_\pi^T$$

### Singular Value Decomposition (DGELSD)

**LAPACK routine:** `DGELSD` (SVD-based least-squares solver)

**Principle:**
$$A = U \cdot \Sigma \cdot V^T$$

**Properties:**
- More robust for (nearly) singular systems
- Provides condition number $\kappa = \sigma_{\max} / \sigma_{\min}$
- Rank truncation: $\sigma_i > \text{RCOND} \cdot \sigma_{\max}$
- Slower than QR: $O(mn \cdot \min(m,n))$

**Cofactor matrix (Qxx):**
$$Q_{xx} = V \cdot \text{diag}(1/\sigma_i^2) \cdot V^T$$

where only singular values $\sigma_i > \text{RCOND} \cdot \sigma_{\max}$ are retained.

**When to use SVD:**
- Ill-conditioned systems (high correlation between parameters)
- Suspected rank deficiency
- When the condition number is needed

### Solver Dispatch by Model Type

| Model Type | Solver Function | Equation System |
|-----------|----------------|-----------------|
| OLS/WLS/GLS | `__adj_solveOLS` | $\tilde{A} \cdot \Delta x = \tilde{w}$ |
| LSE/WLSE/GLSE | `__adj_solveLSE` | Null-space projection: $\bar{A} \cdot z = \bar{w}$ |
| CLS/WCLS/GCLS | `__adj_solveCLS` | $\tilde{B}^T \cdot \lambda = 0$, $v = P^{-1} B^T \lambda$ |
| GLM/WGLM/GGLM | `__adj_solveGLM` | Cholesky: $\tilde{A} \cdot \Delta x = \tilde{w}$ |

---

## Gauss-Newton (GN)

### Principle

Iterative linearization of the nonlinear problem:

1. **Linearization** at current approximate values $x_k$:
   $$A_k \cdot \Delta x + w_k = v$$
2. **Solution** of the linear system → $\Delta x_k$
3. **Update**: $x_{k+1} = x_k + \Delta x_k$
4. **Convergence check**: $\|\Delta x_k\| < \text{tolerance}$

### Convergence

- **Quadratic convergence** near the solution (given sufficiently good initial values)
- **May diverge** with poor initial values or ill-conditioned problems

### When to Use

- Good initial values available
- Well-conditioned system
- Fast convergence desired

---

## Levenberg-Marquardt (LM) with Nielsen Damping

### Principle

Damped version of Gauss-Newton — interpolation between GN and gradient descent:

$$(A^T P A + \lambda \cdot D) \cdot \Delta x = A^T P \ell$$

- $\lambda$ — damping parameter
- $D$ — Marquardt scaling matrix: $D = \text{diag}(A^T P A)$, dynamically updated

The augmentation is realized by adding $\sqrt{\lambda} \cdot \sqrt{D}$ as additional rows (least-squares form).

### Nielsen Strategy

The adjustment of $\lambda$ is based on the gain ratio $\rho$:

$$\rho = \frac{F(x_k) - F(x_k + \Delta x)}{L(0) - L(\Delta x)}$$

with:
- $F(x) = v^T P v$ — actual cost function
- $L(\Delta x)$ — linearized model (predicted reduction)

**Adjustment rules:**

| Condition | Action | Meaning |
|-----------|--------|-----------|
| $\rho > 0$ | $\lambda \leftarrow \lambda \cdot \max(1/3, 1 - (2\rho - 1)^3)$, $\nu = 2$ | Good step → less damping (→ GN) |
| $\rho \le 0$ | $\lambda \leftarrow \lambda \cdot \nu$, $\nu \leftarrow 2 \nu$, rollback to $x_k$ | Bad step → more damping (→ gradient) |

### Initialization of $\lambda_0$

$$\lambda_0 = \tau \cdot \max(\text{diag}(A^T P A))$$

with $\tau = 10^{-3}$ (default).

### Marquardt Scaling Matrix $D$

- Dynamic: $D_{ii} = \max(D_{ii}^{\text{old}}, (A^T P A)_{ii})$ — only increases
- Floor: $D_{ii} \ge \varepsilon_{\text{rel}} \cdot \max(D) + \varepsilon_{\text{abs}}$ with $\varepsilon_{\text{rel}} = 10^{-10}$, $\varepsilon_{\text{abs}} = 10^{-30}$

### Stagnation Detection

Termination when $\lambda$ hits its bounds:
- $\lambda < \lambda_{\min} = 10^{-7}$ (practically reached GN → stagnating in GN mode)
- $\lambda > \lambda_{\max} = 10^7$ (pure gradient descent → no meaningful direction)

Checked via stagnation tolerance: $|\Delta x| < 10^{-10}$.

### LM Rollback

When $\rho \le 0$, the state is rolled back to the previous iteration step:
- **Important:** $\lambda$ is saved BEFORE the rollback and applied AFTER the rollback
- Otherwise the $\lambda$ increase is lost through the state copy

### When to Use

- Poor or unknown initial values
- Ill-conditioned system
- Robustness more important than speed

---

## Convergence Criteria

| Parameter | Default | Description |
|-----------|---------|-------------|
| `.tolerance` | $10^{-10}$ | Termination when $\|\Delta x\| < \text{tolerance}$ |
| `.maxIterations` | 100 | Maximum number of GN/LM iterations |

### LM-Specific Termination Conditions

| Condition | Meaning |
|-----------|-----------|
| $\lambda < \lambda_{\min}$ and stagnation | GN convergence reached |
| $\lambda > \lambda_{\max}$ and stagnation | No meaningful direction |
| $\|\Delta x\| < \text{tolerance}$ | Parameter change below threshold |

### Linear Systems

For fully linear problems (all formulas linear), GN converges in **one iteration**. The UDF detects linearity automatically via `__adj_isLinear` and sets `maxIterations = 1`.

---

## Jacobian Equilibration (Column Scaling)

### Problem

When parameters have different orders of magnitude (e.g., coordinates ~1000 m and angles ~0.001 rad), the normal equation matrix $A^T P A$ is ill-conditioned.

### Solution

Automatic column scaling of the Jacobian matrix:

$$\tilde{A}_{:,j} = A_{:,j} / s_j \quad \text{with} \quad s_j = \|A_{:,j}\|_2$$

This improves the conditioning of the normal equation matrix.

### Back-Transformation

Parameter corrections and cofactor matrices are back-transformed after the solution:

$$\Delta x_j = \Delta \tilde{x}_j / s_j$$
$$(Q_{xx})_{ij} = (\tilde{Q}_{xx})_{ij} / (s_i \cdot s_j)$$

### Configuration

```autoit
Local $mConfig = _adj_defaultConfig()
$mConfig.scaling = False  ; Disable equilibration
```

Default: enabled (`True`).

---

## Differentiation Methods

The Jacobian matrices are computed numerically by default. The method is configurable via `.deriveMethod`.

### Available Methods

| Method | Formula | Function Evaluations | Accuracy |
|---------|--------|----------------------|-------------|
| `"Central"` | $(f(x+h) - f(x-h)) / 2h$ | 2 | $O(h^2)$ |
| `"Central2"` | 4-point stencil | 4 | $O(h^4)$ |
| `"Central3"` | 6-point stencil | 6 | $O(h^6)$ |
| `"Central4"` | 8-point stencil | 8 | $O(h^8)$ |
| `"Forward"` | $(f(x+h) - f(x)) / h$ | 2 | $O(h)$ |
| `"Backward"` | $(f(x) - f(x-h)) / h$ | 2 | $O(h)$ |
| `"Ridder"` | Richardson extrapolation | adaptive | adaptive |
| `"Higham"` | Iterative refinement | adaptive | $\sim \varepsilon^{6/7}$ |

### Step Size

Automatically scaled to the parameter magnitude:

$$h = \begin{cases} \varepsilon^{1/3} \cdot \max(1, |x|) & \text{Central} \\ \sqrt{\varepsilon} \cdot \max(1, |x|) & \text{Forward/Backward} \end{cases}$$

with machine epsilon $\varepsilon \approx 2.2 \times 10^{-16}$.

### Analytical Derivatives

Optionally per formula or map:

```autoit
; Pipe string (per formula)
_adj_addObsFunction($mSystem, "S1", "sqrt(X^2+Y^2)", 100.0, 0.01, "s0", 0, _
    "X = X/sqrt(X^2+Y^2) | Y = Y/sqrt(X^2+Y^2)")

; Map (per formula)
Local $mDeriv[]
$mDeriv["X"] = "X/sqrt(X^2+Y^2)"
$mDeriv["Y"] = "Y/sqrt(X^2+Y^2)"
_adj_addObsFunction($mSystem, "S1", "sqrt(X^2+Y^2)", 100.0, 0.01, "s0", 0, $mDeriv)
```

Analytical derivatives undergo the same formula transformations as the main formulas (parameter names → map accesses, fixed parameters → values).

---

## Variance Component Estimation (VCE)

### Helmert Method

Outer iteration loop around the solver:

1. Solve the adjustment problem with current weights
2. Compute statistics (cofactor matrices, redundancy numbers)
3. Estimate variance factors $\hat{\sigma}_k^2$ per group $k$
4. Update weights: $p_i \leftarrow p_i / \hat{\sigma}_k^2$
5. Repeat until convergence or max. 15 iterations

### Convergence

Relative change of all variance factors $< 5\%$.

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `.vce` | `True` | Enable VCE |
| `.vceMaxIter` | 15 | Max. VCE iterations |
| `.vceConvergence` | 0.05 | Convergence threshold (relative change) |

### Limitation

- Cross-group covariances ($\sigma_{ij}$ with observations from different VCE groups) are rejected (`$ADJ_ERR_INPUT`, `@extended = 4/5`).
- For generalized models with only one group: VCE is skipped (equivalent to $s_0^2$ estimation).
