# 05 -- Constraints and Fixed Parameters

## Incorporating Additional Knowledge

Sometimes you know more than just the measurements.  For example:

- "The sum of two distances must equal 120 m" (geometric condition)
- "The Y coordinate is exactly known" (fixed parameter)
- "Two points lie on a straight line" (linear constraint)

Such **additional knowledge** can be incorporated into the system as a constraint (restriction) or fixed parameter.  The adjustment then considers both simultaneously: the measurements **and** the conditions.


## Example 1: Constraint -- "X + Y = 120"

We take the trilateration network from Chapter 02 and add a constraint: from an independent source, we know that the sum of the coordinates X + Y must equal exactly 120.

```autoit
#include "Adjustment.au3"

Local $mSystem = _adj_createSystem()

; Distance measurements (same as Chapter 02)
_adj_addObsFunction($mSystem, "D1", "Sqrt((X - 0)^2 + (Y - 0)^2)",     72.05, 0.05)
_adj_addObsFunction($mSystem, "D2", "Sqrt((X - 100)^2 + (Y - 0)^2)",   63.10, 0.05)
_adj_addObsFunction($mSystem, "D3", "Sqrt((X - 100)^2 + (Y - 100)^2)", 78.20, 0.05)
_adj_addObsFunction($mSystem, "D4", "Sqrt((X - 0)^2 + (Y - 100)^2)",   85.15, 0.05)

; Constraint: X + Y = 120
_adj_addRestriction($mSystem, "X + Y", 120)

_adj_setInitialValue($mSystem, "X", 50.0)
_adj_setInitialValue($mSystem, "Y", 50.0)

_adj_solve($mSystem, _adj_defaultConfig("GN", False))
If @error Then Exit MsgBox(48, "error", _adj_getErrorMessage(@error))

ConsoleWrite(_adj_displayResults($mSystem))
```

### What Changes

- The solution satisfies the condition **exactly**: X + Y = 120.000
- The result shifts compared to the unconstrained solution (Chapter 02), because the constraint forces the solution onto a line in the X-Y space.
- The degrees of freedom increase: f = 4 measurements - 2 unknowns + 1 constraint = 3.


## Example 2: Fixing a Parameter -- "Y is Known"

Sometimes a parameter is not unknown but **exactly specified**. Example: The Y coordinate was precisely determined by another method and should be held fixed at 55.0 m.

```autoit
#include "Adjustment.au3"

Local $mSystem = _adj_createSystem()

; Distance measurements (same as Chapter 02)
_adj_addObsFunction($mSystem, "D1", "Sqrt((X - 0)^2 + (Y - 0)^2)",     72.05, 0.05)
_adj_addObsFunction($mSystem, "D2", "Sqrt((X - 100)^2 + (Y - 0)^2)",   63.10, 0.05)
_adj_addObsFunction($mSystem, "D3", "Sqrt((X - 100)^2 + (Y - 100)^2)", 78.20, 0.05)
_adj_addObsFunction($mSystem, "D4", "Sqrt((X - 0)^2 + (Y - 100)^2)",   85.15, 0.05)

; Fix Y -- no longer an unknown, but a fixed value
_adj_addFixedParam($mSystem, "Y", 55.0)

_adj_setInitialValue($mSystem, "X", 50.0)

_adj_solve($mSystem, _adj_defaultConfig("GN", False))
If @error Then Exit MsgBox(48, "error", _adj_getErrorMessage(@error))

ConsoleWrite(_adj_displayResults($mSystem))
```

### What Changes

- **Y no longer appears as an unknown** in the results -- it is substituted everywhere as the fixed value 55.0.
- Only X is estimated, with f = 4 - 1 = 3 degrees of freedom.
- The standard deviation of X becomes smaller, because fewer unknowns need to be determined.


## When to Use Which Approach?

| Situation | Method | Example |
|-----------|--------|---------|
| Value is known with certainty to be exact | `_adj_addFixedParam` | Fixed point coordinates |
| A relationship between unknowns must hold | `_adj_addRestriction` | Geometric condition |
| Value is known but with uncertainty | Additional observation | Pseudo-observation with StdDev |

**Tip:** If you are not sure whether a value is truly exact, you can incorporate it as an **additional observation** with a very small standard deviation instead of as a fixed parameter´(this is called a "*pseudo observation*"):

```autoit
; Instead of: _adj_addFixedParam($mSystem, "Y", 55.0)
; Better, if Y has an uncertainty of 0.01 m:
_adj_addObsFunction($mSystem, "Yfix", "Y", 55.0, 0.001)
```


## Summary

- **`_adj_addRestriction`**: Enforces an equation between unknowns. The condition is satisfied exactly.  Increases the degrees of freedom.
- **`_adj_addFixedParam`**: Removes a parameter from the estimation and sets it to a fixed value.  Reduces the number of unknowns.
- Both methods can be combined with any number of observations.
- For "nearly fixed" values with residual uncertainty: use a pseudo-observation instead of a fixed parameter.


---

Next chapter: [06 -- Regression: Fitting a Line to Measurement Data](06_regression.md)