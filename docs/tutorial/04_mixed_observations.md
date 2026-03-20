# 04 -- Mixed Observations: Combining Distances and Angles

## Different Measurement Types in One System

In practice, rarely only one measurement method is available.  A surveyor measures distances and directions simultaneously.  A robot combines GPS positions with wheel odometry.

The Adjustment UDF can **mix arbitrary functional relationships** - as long as each measurement can be expressed as a formula of the unknowns and parameters. Distances, angles, height differences, coordinates: all in the same system.


## The Example: Position from Distances and Directions

Two known points, one unknown point (X, Y):

| Fixed Point | X [m] | Y [m] |
|-------------|-------|-------|
| P1          |   0   |   0   |
| P2          | 100   |   0   |

Measured distances:

| Measurement | Value [m] | StdDev [m] | Formula |
|-------------|-----------|------------|---------|
| S1          | 72.00     | 0.03       | Sqrt((X-0)^2 + (Y-0)^2) |
| S2          | 63.00     | 0.03       | Sqrt((X-100)^2 + (Y-0)^2) |

Measured direction angles (from fixed point to unknown point, in radians):

| Measurement | Value [rad] | StdDev [rad] | Formula |
|-------------|-------------|--------------|---------|
| R1          | 0.9828      | 0.002        | ATan((Y-0) / (X-0)) |
| R2          | 1.9720      | 0.002        | _atan2(Y-0, X-100) |

**Note on ATan and _atan2:** AutoIt only provides `ATan` (arctangent with one argument), which returns values in the range -pi/2 to +pi/2.  For angles in all four quadrants, a helper function `_atan2` is needed.  In the example, P1 lies in the first quadrant (X > 0, Y > 0), so `ATan(Y/X)` is sufficient. For P2 (where X-100 can be negative), we use the helper function (you can use them also directly in the formula).


## The Code

```autoit
#include "Adjustment.au3"

Local $mSystem = _adj_createSystem()

; --- Distance measurements ---
_adj_addObsFunction($mSystem, "S1", "Sqrt(X^2 + Y^2)",           72.00, 0.03)
_adj_addObsFunction($mSystem, "S2", "Sqrt((X - 100)^2 + Y^2)",   63.00, 0.03)

; --- Direction measurements [rad] ---
; R1: Angle from P1(0,0) to the point -- X and Y are positive, ATan suffices
_adj_addObsFunction($mSystem, "R1", "ATan(Y / X)", 0.9828, 0.002)

; R2: Angle from P2(100,0) to the point -- _atan2 for correct quadrant
_adj_addObsFunction($mSystem, "R2", "_atan2(Y, X - 100)", 1.9720, 0.002)

; Initial values
_adj_setInitialValue($mSystem, "X", 50.0)
_adj_setInitialValue($mSystem, "Y", 50.0)

; Solve
_adj_solve($mSystem)
If @error Then Exit MsgBox(48, "error", _adj_getErrorMessage(@error))

; show results
ConsoleWrite(_adj_displayResults($mSystem))

; --- Helper function: atan2 for all quadrants ---
Func _atan2($fY, $fX)
    If $fX > 0 Then Return ATan($fY / $fX)
    If $fX < 0 And $fY >= 0 Then Return ATan($fY / $fX) + ACos(-1)
    If $fX < 0 And $fY < 0 Then Return ATan($fY / $fX) - ACos(-1)
    If $fX = 0 And $fY > 0 Then Return ACos(-1) / 2
    If $fX = 0 And $fY < 0 Then Return -ACos(-1) / 2
    Return 0
EndFunc
```


## What Happens Here

- **Distances** and **angles** have completely different units (meters vs. radians) and different standard deviations.
- The adjustment handles this automatically: through the standard deviations, the measurement types are correctly weighted relative to each other.
- The formulas may contain any AutoIt expressions, as long as they use the unknowns (here X and Y) as variables and are evaluable by AutoIt's `Execute()`.
- Custom functions (such as `_atan2`) can be used in formulas, provided they are defined globally.


## What You Will See in the Results

| Quantity | Meaning |
|----------|---------|
| **f = 2** | 4 measurements - 2 unknowns = 2 degrees of freedom |
| **v(S1), v(S2)** | Residuals of the distance measurements [m] |
| **v(R1), v(R2)** | Residuals of the direction measurements [rad] |
| **sdx(X), sdx(Y)** | The position is determined more precisely than with only distances or only angles alone |

Combining different measurement types generally improves both the accuracy and the reliability of the result.


## Summary

- The UDF mixes **arbitrary observation types** in the same system.
- Each observation only needs: **name**, **formula**, **measured value**, **standard deviation**.
- The weighting between different measurement types is derived automatically from the standard deviations.
- **Custom AutoIt functions** can be used in formulas.


---

Next chapter: [05 -- Constraints and Fixed Parameters](05_constraints.md)