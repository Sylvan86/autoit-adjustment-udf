# 03 -- Weighting: Not All Measurements Are Equally Precise

## Why Weighting?

In reality, measurements have different levels of precision. A high-precision distance measuring device delivers better values than a simple tape measure. A measurement taken in clear visibility is more reliable than one taken in rain.

The adjustment should take this into account: **More precise measurements should have a stronger influence on the result than imprecise ones.** This is exactly what the standard deviation (`$fStdDev`) in `_adj_addObsFunction` controls.

**Small standard deviation = high precision = strong influence.**


## The Example: Trilateration with Different Precision Levels

The same network as in Chapter 02, but with more realistic assumptions:

| Measurement | Distance [m] | StdDev [m] | Remark |
|-------------|--------------|------------|--------|
| D1          | 72.05        | 0.02       | Precision instrument |
| D2          | 63.10        | 0.05       | Standard |
| D3          | 78.20        | 0.05       | Standard |
| D4          | 85.15        | 0.10       | Poor visibility conditions |

Measurement D1 is 5x more precise than D4. Its weight (internally: 1/StdDev^2) is therefore 25x higher.


## The Code

```autoit
#include "Adjustment.au3"

Local $mSystem = _adj_createSystem()

; Same formulas as Chapter 02 -- only the standard deviations change
_adj_addObsFunction($mSystem, "D1", "Sqrt((X - 0)^2 + (Y - 0)^2)",     72.05, 0.02)
_adj_addObsFunction($mSystem, "D2", "Sqrt((X - 100)^2 + (Y - 0)^2)",   63.10, 0.05)
_adj_addObsFunction($mSystem, "D3", "Sqrt((X - 100)^2 + (Y - 100)^2)", 78.20, 0.05)
_adj_addObsFunction($mSystem, "D4", "Sqrt((X - 0)^2 + (Y - 100)^2)",   85.15, 0.10)

_adj_setInitialValue($mSystem, "X", 50.0)
_adj_setInitialValue($mSystem, "Y", 50.0)

_adj_solve($mSystem, _adj_defaultConfig("GN", False))
If @error Then
    ConsoleWrite("Fehler: @error=" & @error & @CRLF)
    Exit 1
EndIf

ConsoleWrite(_adj_displayResults($mSystem))
```


## Comparison: Equal vs. Different Weighting

| Aspect | Equal StdDev (Ch. 02) | Different StdDev |
|--------|----------------------|------------------|
| Influence | All measurements count equally | Precise measurement D1 dominates |
| Position | Result lies "in the middle" | Result is pulled toward D1 |
| Residuals | Residuals v are evenly distributed | D4 receives a larger residual |

You will notice that the adjusted coordinates shift slightly. The more precise measurement D1 "pulls" the result toward itself, while the imprecise measurement D4 can only contribute weakly.


## Reading the Residuals

```autoit
Local $mRes = _adj_getResults($mSystem)
Local $mV   = $mRes.v

ConsoleWrite("v(D1) = " & $mV["D1"] & " m" & @CRLF)
ConsoleWrite("v(D2) = " & $mV["D2"] & " m" & @CRLF)
ConsoleWrite("v(D3) = " & $mV["D3"] & " m" & @CRLF)
ConsoleWrite("v(D4) = " & $mV["D4"] & " m" & @CRLF)
```

With equal weighting, the residuals are similar in magnitude. With different weighting, the imprecise measurement D4 is corrected more strongly than the precise measurement D1 -- this is intentional.


## Summary

- The **standard deviation** controls how strongly a measurement influences the result.
- You do not need to compute weights manually -- the UDF handles this internally.
- **Realistic standard deviations are crucial** for a good result. When in doubt, it is better to estimate slightly too large rather than too optimistically.
- When all measurements have the same standard deviation, the adjustment corresponds to simple ordinary least squares (OLS).


---

Next chapter: [04 -- Mixed Observations: Combining Distances and Angles](04_mixed_observations.md)
