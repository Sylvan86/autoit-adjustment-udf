# 06 -- Regression: From Simple to Errors-in-Variables

## Adjustment Is Regression

In the previous chapters we estimated coordinates or physical quantities. But adjustment can also do what statisticians call **regression**: fitting a curve (or line) through data points such that the deviations are minimized.

This is not a special case -- it is **the same principle** as before.  Instead of "distance as a function of coordinates" we write "y as a function of x".  Everything else stays the same.

This chapter walks through a hierarchy of regression models, each generalizing the previous one.  All variants use the **same 6 data points** so you can directly compare the results.


## The Data

Six measurements of a physical quantity y at different values of x:

| Point | x    | y     |
|-------|------|-------|
| P1    | 1.0  |  2.1  |
| P2    | 2.0  |  3.9  |
| P3    | 3.0  |  6.2  |
| P4    | 4.0  |  7.8  |
| P5    | 5.0  | 10.3  |
| P6    | 6.0  | 11.7  |

Goal: the best-fit line **y = A + B * x**.


## The Hierarchy

| # | Method | Errors in x? | Weights | UDF model type |
|---|--------|-------------|---------|----------------|
| 1 | OLS (Ordinary Least Squares) | no | equal | OLS |
| 2 | WLS (Weighted Least Squares) | no | individual σ_y | WLS |
| 3 | TLS / Orthogonal Regression | yes, σ_x = σ_y | equal | GLM |
| 4 | Deming Regression | yes, σ_x/σ_y = const | constant ratio | WGLM |
| 5 | York Regression | yes, individual σ_x, σ_y | per point | WGLM |
| 6 | York with correlations | yes, individual + correlated | per point + r_xy | GGLM |

Each step generalizes the previous one.  The UDF handles them all -- you just change how you set up the observations.


## Variant 1: Simple Regression (OLS)

All measurements have the same standard deviation -- an unweighted regression.  Only y is subject to error; x is treated as exact.

We use `_adj_addObsFunction`: each observation has its own formula relating it to the parameters A and B.

```autoit
#include "Adjustment.au3"

Local $mSystem = _adj_createSystem()

; Each observation: y_i = A + B * x_i
_adj_addObsFunction($mSystem, "P1", "A + B * 1.0",  2.1, 0.5)
_adj_addObsFunction($mSystem, "P2", "A + B * 2.0",  3.9, 0.5)
_adj_addObsFunction($mSystem, "P3", "A + B * 3.0",  6.2, 0.5)
_adj_addObsFunction($mSystem, "P4", "A + B * 4.0",  7.8, 0.5)
_adj_addObsFunction($mSystem, "P5", "A + B * 5.0", 10.3, 0.5)
_adj_addObsFunction($mSystem, "P6", "A + B * 6.0", 11.7, 0.5)

_adj_solve($mSystem)
If @error Then Exit MsgBox(48, "error", _adj_getErrorMessage(@error))

ConsoleWrite(_adj_displayResults($mSystem))
```

Result: A (y-intercept), B (slope), residuals v(P1)..v(P6), f = 6 - 2 = 4 degrees of freedom.


## Variant 2: Weighted Regression (WLS)

Same model, but the measurements have **different precisions**. More precise points pull the line more strongly toward themselves.

```autoit
Local $mSystem = _adj_createSystem()

; Same data, but P5 and P6 measured more precisely (smaller stdDev)
_adj_addObsFunction($mSystem, "P1", "A + B * 1.0",  2.1, 1.0)
_adj_addObsFunction($mSystem, "P2", "A + B * 2.0",  3.9, 1.0)
_adj_addObsFunction($mSystem, "P3", "A + B * 3.0",  6.2, 1.0)
_adj_addObsFunction($mSystem, "P4", "A + B * 4.0",  7.8, 1.0)
_adj_addObsFunction($mSystem, "P5", "A + B * 5.0", 10.3, 0.2)
_adj_addObsFunction($mSystem, "P6", "A + B * 6.0", 11.7, 0.2)

_adj_solve($mSystem)
If @error Then Exit MsgBox(48, "error", _adj_getErrorMessage(@error))

ConsoleWrite(_adj_displayResults($mSystem))
```

The line shifts compared to Variant 1: P5 and P6 carry 25x more weight (stdDev 0.2 vs. 1.0 → weight 1/0.04 vs. 1/1.0).


## Variant 3: Orthogonal / Total Least Squares (TLS)

Up to now only y was subject to error.  But what if **both x and y are measured** -- and therefore both subject to error?

This is an **errors-in-variables** model.  In the simplest case (TLS / orthogonal regression), x and y have the **same standard deviation**.

The key API difference: we now use `_adj_addObs` for the observations (both x and y) and `_adj_addFunction` for the condition equations that link them.  This is a **Gauss-Helmert model - or in english "Generalized Linear Model" (GLM)** -- observations appear inside the functions, referenced by the `#` prefix.

```autoit
Local $mSystem = _adj_createSystem()

; All observations — both x and y are measured (same stdDev = 0.5)
_adj_addObs($mSystem, "X1", 1.0, 0.5)
_adj_addObs($mSystem, "X2", 2.0, 0.5)
_adj_addObs($mSystem, "X3", 3.0, 0.5)
_adj_addObs($mSystem, "X4", 4.0, 0.5)
_adj_addObs($mSystem, "X5", 5.0, 0.5)
_adj_addObs($mSystem, "X6", 6.0, 0.5)
_adj_addObs($mSystem, "Y1", 2.1, 0.5)
_adj_addObs($mSystem, "Y2", 3.9, 0.5)
_adj_addObs($mSystem, "Y3", 6.2, 0.5)
_adj_addObs($mSystem, "Y4", 7.8, 0.5)
_adj_addObs($mSystem, "Y5", 10.3, 0.5)
_adj_addObs($mSystem, "Y6", 11.7, 0.5)

; Condition equations: y_i - A - B * x_i = 0
; #X1, #Y1 etc. reference the observations above
_adj_addFunction($mSystem, "#Y1 - A - B * #X1", 0)
_adj_addFunction($mSystem, "#Y2 - A - B * #X2", 0)
_adj_addFunction($mSystem, "#Y3 - A - B * #X3", 0)
_adj_addFunction($mSystem, "#Y4 - A - B * #X4", 0)
_adj_addFunction($mSystem, "#Y5 - A - B * #X5", 0)
_adj_addFunction($mSystem, "#Y6 - A - B * #X6", 0)

_adj_solve($mSystem)
If @error Then Exit MsgBox(48, "error", _adj_getErrorMessage(@error))

ConsoleWrite(_adj_displayResults($mSystem))
```

### What changes

- There are now **12 observations** (6 x-values + 6 y-values) but still only 2 parameters.
- Both x and y receive residuals -- the line adjusts in **both** directions.
- Because σ_x = σ_y, the line is fitted orthogonally (minimizing perpendicular distances rather than vertical ones).

### Why `_adj_addObs` + `_adj_addFunction`?

In Variants 1-2, each formula depended on exactly **one** observation -- so `_adj_addObsFunction` was appropriate.

Here, each condition equation depends on **two** observations (#X_i and #Y_i).  In such cases, observations and functions must be added separately:

- `_adj_addObs` — registers an observation with its value and standard deviation
- `_adj_addFunction` — defines a condition equation referencing observations via `#`

Basically, all models can be described using _adj_addObs/_adj_addFunction. _adj_addObsFunction is simply a shorthand version of this for cases where the function contains only one observation.

## Variant 4: Deming Regression

Like TLS, but the **ratio** of standard deviations between x and y is constant (and known), not necessarily 1:1.

For example, if y-measurements are less precise than x-measurements by a factor of about 1.7 (σ_y/σ_x ≈ 1.7):

```autoit
Local $mSystem = _adj_createSystem()

; x-values with stdDev = 0.3, y-values with stdDev = 0.5
; The ratio sigma_y/sigma_x = 0.5/0.3 ≈ 1.67 is constant across all points
_adj_addObs($mSystem, "X1", 1.0, 0.3)
_adj_addObs($mSystem, "X2", 2.0, 0.3)
_adj_addObs($mSystem, "X3", 3.0, 0.3)
_adj_addObs($mSystem, "X4", 4.0, 0.3)
_adj_addObs($mSystem, "X5", 5.0, 0.3)
_adj_addObs($mSystem, "X6", 6.0, 0.3)
_adj_addObs($mSystem, "Y1",  2.1, 0.5)
_adj_addObs($mSystem, "Y2",  3.9, 0.5)
_adj_addObs($mSystem, "Y3",  6.2, 0.5)
_adj_addObs($mSystem, "Y4",  7.8, 0.5)
_adj_addObs($mSystem, "Y5", 10.3, 0.5)
_adj_addObs($mSystem, "Y6", 11.7, 0.5)

; Same condition equations as TLS
_adj_addFunction($mSystem, "#Y1 - A - B * #X1", 0)
_adj_addFunction($mSystem, "#Y2 - A - B * #X2", 0)
_adj_addFunction($mSystem, "#Y3 - A - B * #X3", 0)
_adj_addFunction($mSystem, "#Y4 - A - B * #X4", 0)
_adj_addFunction($mSystem, "#Y5 - A - B * #X5", 0)
_adj_addFunction($mSystem, "#Y6 - A - B * #X6", 0)

_adj_solve($mSystem)
If @error Then Exit MsgBox(48, "error", _adj_getErrorMessage(@error))

ConsoleWrite(_adj_displayResults($mSystem))
```

The only difference to TLS: the standard deviations of x and y differ. This shifts the line because the solver now weights the x-corrections differently than the y-corrections.

Special cases:
- σ_x = σ_y → Orthogonal / TLS (Variant 3)
- σ_x → 0 (x exact) → OLS (Variant 1)


## Variant 5: York Regression

The most general case without correlations: each point has its own **individual** standard deviations for both x and y.  This is known as York regression (York, 1966).

```autoit
Local $mSystem = _adj_createSystem()

; Each point has individual uncertainties in x and y
_adj_addObs($mSystem, "X1", 1.0, 0.10)
_adj_addObs($mSystem, "X2", 2.0, 0.20)
_adj_addObs($mSystem, "X3", 3.0, 0.15)
_adj_addObs($mSystem, "X4", 4.0, 0.30)
_adj_addObs($mSystem, "X5", 5.0, 0.25)
_adj_addObs($mSystem, "X6", 6.0, 0.40)
_adj_addObs($mSystem, "Y1",  2.1, 0.5)
_adj_addObs($mSystem, "Y2",  3.9, 0.3)
_adj_addObs($mSystem, "Y3",  6.2, 0.4)
_adj_addObs($mSystem, "Y4",  7.8, 0.6)
_adj_addObs($mSystem, "Y5", 10.3, 0.2)
_adj_addObs($mSystem, "Y6", 11.7, 0.3)

_adj_addFunction($mSystem, "#Y1 - A - B * #X1", 0)
_adj_addFunction($mSystem, "#Y2 - A - B * #X2", 0)
_adj_addFunction($mSystem, "#Y3 - A - B * #X3", 0)
_adj_addFunction($mSystem, "#Y4 - A - B * #X4", 0)
_adj_addFunction($mSystem, "#Y5 - A - B * #X5", 0)
_adj_addFunction($mSystem, "#Y6 - A - B * #X6", 0)

_adj_solve($mSystem)
If @error Then Exit MsgBox(48, "error", _adj_getErrorMessage(@error))

ConsoleWrite(_adj_displayResults($mSystem))
```

Each point now contributes to the solution according to its own precision in both dimensions.  Points with small σ_x and σ_y constrain the line more than imprecise ones.


## Variant 6: York with Correlations (GGLM)

Sometimes x and y of the same point are **correlated** -- for example because they originate from the same measurement process.  This is the full York regression (York, 1969).

Add covariances between the x and y observations of each point:

```autoit
; ... (same observations and functions as Variant 5) ...

; Correlation between x and y of each point
_adj_addCovariance($mSystem, "X1", "Y1", 0.02)
_adj_addCovariance($mSystem, "X2", "Y2", 0.03)
_adj_addCovariance($mSystem, "X3", "Y3", 0.02)
_adj_addCovariance($mSystem, "X4", "Y4", 0.05)
_adj_addCovariance($mSystem, "X5", "Y5", 0.01)
_adj_addCovariance($mSystem, "X6", "Y6", 0.04)

_adj_solve($mSystem)
```

The covariance information shifts the result because the solver now accounts for the statistical dependency between x and y of each point. Chapter 07 explains covariances in detail.


## Summary

| Variant | Errors in | Weights | API pattern |
|---------|----------|---------|-------------|
| 1 - OLS | y only | equal | `_adj_addObsFunction` |
| 2 - WLS | y only | individual σ_y | `_adj_addObsFunction` |
| 3 - TLS | x and y | equal σ | `_adj_addObs` + `_adj_addFunction` |
| 4 - Deming | x and y | constant ratio σ_y/σ_x | `_adj_addObs` + `_adj_addFunction` |
| 5 - York | x and y | individual σ_x, σ_y | `_adj_addObs` + `_adj_addFunction` |
| 6 - York+corr. | x and y | individual + correlated | + `_adj_addCovariance` |

Key takeaway: when only **one** observation feeds into a formula, use `_adj_addObsFunction`.  When **multiple** observations appear in a condition equation, use `_adj_addObs` + `_adj_addFunction` separately.


---

Previous chapter: [05 -- Constraints and Fixed Parameters](05_constraints.md)

Next chapter: [07 -- Covariance Matrix: Correlations Between Measurements](07_covariance_matrix.md)