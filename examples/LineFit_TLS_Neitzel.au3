#AutoIt3Wrapper_Change2CUI=y

; ===============================================================================
;   Line Fit (TLS) -- Neitzel & Petrovic
; ===============================================================================
;
; Source:
;   Neitzel, F., Petrovic, S. (2008): Total Least Squares (TLS) im Kontext
;   der Ausgleichung nach kleinsten Quadraten am Beispiel der ausgleichenden
;   Geraden, zfv 133(3), pp. 141-148.
;   Table 3 (strict Gauss-Helmert solution).
;
; Problem:
;   Best-fit line y = a*x + b through 4 points on a parabola (y = x^2).
;   Both x and y are error-prone observations -> Total Least Squares.
;   Equal unit weights for all observations (P = I, sigma = 1).
;
;   Points:
;     (0, 0)   (1, 1)   (2, 4)   (3, 9)
;
; Classification:
;   Model:       Gauss-Helmert (implicit functional model)
;   Weighting:   Unweighted (P = I) — orthogonal regression
;   Linearity:   Linear condition equations (y - a*x - b = 0)
;   Algorithm:   Gauss-Newton
;   VCE:         No
;   Redundancy:  f = c - u = 4 - 2 = 2
;
; ===============================================================================


; -- Initialization -----------------------------------------------------------
#include "../Adjustment.au3"


; -- Model setup --------------------------------------------------------------
Global $mSystem

; Observations: x and y coordinates (all with sigma = 1 -> unit weights)
;
;   Obs.     Value     sigma
;   x₁        0        1
;   y₁        0        1
;   x₂        1        1
;   y₂        1        1
;   x₃        2        1
;   y₃        4        1
;   x₄        3        1
;   y₄        9        1

_adj_addObs($mSystem, "x1", 0, 1)
_adj_addObs($mSystem, "y1", 0, 1)
_adj_addObs($mSystem, "x2", 1, 1)
_adj_addObs($mSystem, "y2", 1, 1)
_adj_addObs($mSystem, "x3", 2, 1)
_adj_addObs($mSystem, "y3", 4, 1)
_adj_addObs($mSystem, "x4", 3, 1)
_adj_addObs($mSystem, "y4", 9, 1)

; Condition equations:  yᵢ - a·xᵢ - b = 0
; Parameters (A, B) have NO # prefix; observations have # prefix.
_adj_addFunction($mSystem, "#y1 - A * #x1 - B", 0)
_adj_addFunction($mSystem, "#y2 - A * #x2 - B", 0)
_adj_addFunction($mSystem, "#y3 - A * #x3 - B", 0)
_adj_addFunction($mSystem, "#y4 - A * #x4 - B", 0)

; Initial values for parameters (from OLS pre-solution)
_adj_setInitialValue($mSystem, "A", 3.0)
_adj_setInitialValue($mSystem, "B", -1.0)


; -- Adjustment ---------------------------------------------------------------
_adj_solve($mSystem, _adj_defaultConfig("GN", False))
If @error Then
	ConsoleWrite("!> ERROR in _adj_solve: @error=" & @error & " @extended=" & @extended & @CRLF)
	Exit 1
EndIf


; -- Results ------------------------------------------------------------------
ConsoleWrite(_adj_displayResults($mSystem))


; -- Validation vs. reference values (Neitzel & Petrovic 2008, Table 3) -------
Global $mRes = _adj_getResults($mSystem)
Global $mX1  = $mRes.x1
Global $mSdx = $mRes.sdx

ConsoleWrite(@CRLF)
ConsoleWrite("===================================================" & @CRLF)
ConsoleWrite("  Validation vs. Reference (Neitzel 2008)          " & @CRLF)
ConsoleWrite("===================================================" & @CRLF & @CRLF)

; --- Parameters ---
ConsoleWrite("Adjusted parameters:" & @CRLF)
_check("a",  $mX1["A"],   3.241804, 6)
_check("b",  $mX1["B"],  -1.362705, 6)

; --- Standard deviations ---
ConsoleWrite(@CRLF & "Standard deviations:" & @CRLF)
_check("sigma(a)", $mSdx["A"], 0.678679, 6)
_check("sigma(b)", $mSdx["B"], 1.254155, 6)

; --- Statistical summary ---
ConsoleWrite(@CRLF & "Statistical summary:" & @CRLF)
_check("vtPv", $mRes.vtpv, 0.372946, 6)
_check("f",    $mRes.f,    2, 0)

; --- Method comparison (informational) ---
ConsoleWrite(@CRLF & "Method comparison (for reference):" & @CRLF)
ConsoleWrite("  OLS (error-free x):  a = 3.000000,  b = -1.000000,  vtPv = 2.000000" & @CRLF)
ConsoleWrite("  OLS (error-free y):  a = 0.285714,  b =  0.333333,  vtPv = 0.761905" & @CRLF)
ConsoleWrite(StringFormat("  TLS / GHM (strict):  a = %f,  b = %f,  vtPv = %f", _
	$mX1["A"], $mX1["B"], $mRes.vtpv) & @CRLF)


; -- Helper -------------------------------------------------------------------
Func _check($sName, $fActual, $fExpected, $iDez)
	Local $sFmt  = "%." & $iDez & "f"
	Local $fTol  = 0.5 * (10 ^ (-$iDez))
	Local $fDiff = Abs($fActual - $fExpected)
	Local $sOK   = ($fDiff <= $fTol) ? "✓" : "✗ MISMATCH"
	ConsoleWrite(StringFormat("  %-12s  expect: " & $sFmt & "  got: " & $sFmt & "  [%s]", _
		$sName, $fExpected, $fActual, $sOK) & @CRLF)
EndFunc
