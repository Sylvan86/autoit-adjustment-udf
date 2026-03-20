# 02 -- First Network: Determining Position from Distances (Trilateration)

## Multiple Unknowns

In Chapter 01, all measurements shared the same unknown (X).
Now it gets more interesting: **several different measurements** jointly determine **multiple unknowns**.

Imagine you are standing at an unknown point and measure the distance to four known stations.  Each individual distance describes a circle on which you could be located.  Only the combination of all four circles yields your position - and the adjustment finds the best possible intersection point.


## The Example: Trilateration

Four reference points with known coordinates:

| Point | X [m] | Y [m] |
|-------|-------|-------|
| P1    |   0   |   0   |
| P2    | 100   |   0   |
| P3    | 100   | 100   |
| P4    |   0   | 100   |

Measured distances to the unknown point (X, Y):

| Measurement | Distance [m] | StdDev [m] |
|-------------|--------------|------------|
| D1          | 72.05        | 0.05       |
| D2          | 63.10        | 0.05       |
| D3          | 78.20        | 0.05       |
| D4          | 85.15        | 0.05       |

And we estimate that, over these distances, we can measure a distance with our tape measure to within about 5 cm. (these are the 0.05 StdDev values).

## The Functional Relationship

The distance from the unknown point (X, Y) to a known point (x0, y0) is:

    Distance = Sqrt((X - x0)^2 + (Y - y0)^2)

This formula is used for each measurement with the respective reference point coordinates.  Rather than hardcoding the coordinates directly into the formula strings, you can declare them as **fixed parameters** via `_adj_addFixedParam()`.  The solver substitutes their values automatically — keeping the formulas readable and the coordinates easy to change.


## The Code

```autoit
#include "Adjustment.au3"

Local $mSystem = _adj_createSystem()

; Known reference point coordinates — declared as fixed parameters
; Alternatively, you can enter the numerical value directly into the formula - but this way is clearer.
_adj_addFixedParam($mSystem, "X1",   0)
_adj_addFixedParam($mSystem, "Y1",   0)
_adj_addFixedParam($mSystem, "X2", 100)
_adj_addFixedParam($mSystem, "Y2",   0)
_adj_addFixedParam($mSystem, "X3", 100)
_adj_addFixedParam($mSystem, "Y3", 100)
_adj_addFixedParam($mSystem, "X4",   0)
_adj_addFixedParam($mSystem, "Y4", 100)

; Distance measurements with formulas
; Each formula describes the relationship measured value <-> unknowns (X, Y)
; Fixed parameters (X1, Y1, ...) are substituted automatically.
_adj_addObsFunction($mSystem, "D1", "Sqrt((X - X1)^2 + (Y - Y1)^2)", 72.05, 0.05)
_adj_addObsFunction($mSystem, "D2", "Sqrt((X - X2)^2 + (Y - Y2)^2)", 63.10, 0.05)
_adj_addObsFunction($mSystem, "D3", "Sqrt((X - X3)^2 + (Y - Y3)^2)", 78.20, 0.05)
_adj_addObsFunction($mSystem, "D4", "Sqrt((X - X4)^2 + (Y - Y4)^2)", 85.15, 0.05)

; Initial values -- rough estimate of the position
_adj_setInitialValue($mSystem, "X", 50.0)
_adj_setInitialValue($mSystem, "Y", 50.0)

; Solve
_adj_solve($mSystem)
If @error Then Exit MsgBox(48, "error", _adj_getErrorMessage(@error))

; Results
ConsoleWrite(_adj_displayResults($mSystem))

; Programmatic access to individual values (see `\docs\reference\results.md`)
Local $mRes = _adj_getResults($mSystem)
ConsoleWrite("X = " & $mRes.x1["X"] & " m" & @CRLF)
ConsoleWrite("Y = " & $mRes.x1["Y"] & " m" & @CRLF)
ConsoleWrite("s0 = " & $mRes.s0 & @CRLF)
```


## Why Are Initial Values Needed?

The formula `Sqrt((X - x0)^2 + (Y - y0)^2)` is **nonlinear** - it contains squares and a square root.  The solver (Gauss-Newton) works iteratively: it starts at the initial values and improves the solution step by step.  Without reasonable initial values, it may search in the wrong region or fail to converge entirely.

**Rule of thumb:** The initial values do not need to be exact, but they should be in the right order of magnitude.  In the example, we know that the point lies somewhere within the 100x100 square -- so (50, 50) is a reasonable estimate.


## Understanding the Results

| Quantity | Meaning |
|----------|---------|
| **X, Y** | Adjusted coordinates of the unknown point |
| **sdx(X), sdx(Y)** | Standard deviations -- how precisely are X and Y determined? |
| **s0** | Quality measure of the overall solution (s0 close to 1: model and measurement precision are consistent) |
| **f** | Degrees of freedom = 4 measurements - 2 unknowns = 2 |
| **v(D1)..v(D4)** | Residuals -- show how much each measurement was corrected |
| **nIterations** | Number of iterations until convergence |


## Summary

- Each observation gets its **own formula**, but all formulas share the **same unknowns** (X, Y).
- For nonlinear formulas, **initial values** are required.
- The adjustment finds the combination of X and Y that best fits **all** measurements simultaneously.


---

Next chapter: [03 -- Weighting: Not All Measurements Are Equally Precise](03_weighting.md)
