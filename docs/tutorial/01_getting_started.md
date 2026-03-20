# 01 -- Getting Started: What is Adjustment?

## The Basic Idea

Imagine you measure a distance five times with a tape measure.
Each measurement yields a slightly different value -- because every measurement is subject to random measurement errors.  So which value is the "correct" one?

The answer: **none of them**.  But from all five measurements together, a **best possible estimate** can be computed -- one that keeps the measurement errors as small as possible.  That is exactly what an adjustment does.

**Key principle:** More measurements than necessary --> best possible estimate + quality measure for the result.


## The Workflow in Four Steps

```
1. Create a system            _adj_createSystem()
2. Add observations           _adj_addObsFunction(...)
3. Solve                      _adj_solve(...)
4. Evaluate results           _adj_displayResults(...) / _adj_getResults(...)
```

This pattern stays the same in every chapter -- only the observations and the functional relationship become more complex.


## Example: Five Measurements, One Unknown

Five distance measurements (in meters):

| Measurement | Value [m] |
|-------------|-----------|
| M1          | 10.02     |
| M2          |  9.98     |
| M3          | 10.01     |
| M4          | 10.03     |
| M5          |  9.99     |

Goal: find the best possible value **X** for the distance.

Each measurement "observes" the same quantity X.  The formula is therefore simply `"X"`.

```autoit
#include "Adjustment.au3"

; --- 1. Create system ---
Local $mSystem = _adj_createSystem()

; --- 2. Add observations ---
;                          Symbol Formula Measured value
_adj_addObsFunction($mSystem, "M1", "X", 10.02)
_adj_addObsFunction($mSystem, "M2", "X",  9.98)
_adj_addObsFunction($mSystem, "M3", "X", 10.01)
_adj_addObsFunction($mSystem, "M4", "X", 10.03)
_adj_addObsFunction($mSystem, "M5", "X",  9.99)

; --- 3. Solve ---
_adj_solve($mSystem)
If @error Then Exit MsgBox(48, "error", _adj_getErrorMessage(@error))

; --- 4. Display results ---
ConsoleWrite(_adj_displayResults($mSystem))
```


## What Do the Results Mean?

| Quantity | Meaning | Expected Value |
|----------|---------|----------------|
| **X**    | Estimated distance value | 10.006 m (= mean of the 5 measurements) |
| **s0**   | Standard deviation of the unit weight - the compatibility between the actual measurements, the mathematical functional model, and the a priori stochastic weights | approx. 0.021 |
| **f**    | Degrees of freedom = number of measurements minus number of unknowns = 5 - 1 = 4 | 4 |
| **v**    | Residuals -- the difference between adjusted and measured value for each measurement | e.g. v(M1) = -0.014 |

**s0 close to 1** means: the specified standard deviations match the actual scatter of the measurements well.  If s0 is significantly larger than 1, the measurements were worse than assumed (or the model is incorrect).


## Setup Notes

For the code to work, you need:

1. AutoIt
2. The file `Adjustment.au3` (and the internal files it includes) in the same directory or reachable via a path specification.
3. The OpenBLAS DLL (`libopenblas.dll` for 64-bit, `libopenblas_x86.dll` for 32-bit) in the same directory.  The DLL is searched for automatically on the first `#include`.  Alternatively, it can be loaded beforehand:

```autoit
$__g_hBLAS_DLL = DllOpen("libopenblas.dll")
#include "Adjustment.au3"
```


---

Next chapter: [02 -- First Network: Trilateration](02_first_network.md)
