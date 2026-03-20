#AutoIt3Wrapper_Change2CUI=y

; ===============================================================================
;   2D Trilateration Network -- Krumm / Benning
; ===============================================================================
;
; Source:
;   Krumm, F. (2020): Geodetic Network Adjustment Examples, Rev. 3.5,
;   Geodätisches Institut, Universität Stuttgart.
;   Example 3.1.1 (pp. 48–50), based on Benning (2011), Ex. 8-2.
;
; Problem:
;   2D trilateration network with 4 points.
;   Points 1 and 2 are fixed (datum).
;   Points 3 and 4 are new points (4 unknowns: x3, y3, x4, y4).
;   5 horizontal distance observations.
;   Nonlinear model solved by Gauss–Newton iteration.
;
;          1 ─────────── 2
;          │╲           ╱│
;          │  ╲       ╱  │
;          │    ╲   ╱    │
;          │      ╳      │
;          │    ╱   ╲    │
;          │  ╱       ╲  │
;          │╱           ╲│
;          3 ─────────── 4
;
; Classification:
;   Model:       Observation equations (Gauss–Markov)
;   Weighting:   Weighted (WLS) — uniform σ = 10 mm
;   Linearity:   Nonlinear  (sᵢⱼ = √((xⱼ−xᵢ)² + (yⱼ−yᵢ)²))
;   Algorithm:   Gauss–Newton (multiple iterations)
;   VCE:         No
;   Datum:       Points 1 and 2 fixed
;   Redundancy:  f = n − u = 5 − 4 = 1
;
; ===============================================================================


; -- Initialization -----------------------------------------------------------------
#include "../Adjustment.au3"


; -- Fixed points -------------------------------------------------------------------
;   Point 1:  x = 0.0000 m,  y = 1000.0000 m
;   Point 2:  x = 1000.0000 m,  y = 1000.0000 m
Local Const $fX1 = 0.0, $fY1 = 1000.0
Local Const $fX2 = 1000.0, $fY2 = 1000.0


; -- Model setup --------------------------------------------------------------------
Local $mSystem

; Observations: horizontal distances [m] with σ = 0.010 m
; Functional model:  sᵢⱼ = Sqrt((xⱼ − xᵢ)² + (yⱼ − yᵢ)²)
;
;   Obs.    from  to     s [m]      σ [m]
;   s_13     1    3    1000.0200    0.010
;   s_14     1    4    1414.2000    0.010
;   s_23     2    3    1414.2400    0.010
;   s_24     2    4     999.9800    0.010
;   s_34     3    4    1000.0000    0.010
Local $sFunc

$sFunc = StringFormat("Sqrt((x3 - %.4f)^2 + (y3 - %.4f)^2)", $fX1, $fY1)
_adj_addObsFunction($mSystem, "s_13", $sFunc, 1000.0200, 0.010)

$sFunc = StringFormat("Sqrt((x4 - %.4f)^2 + (y4 - %.4f)^2)", $fX1, $fY1)
_adj_addObsFunction($mSystem, "s_14", $sFunc, 1414.2000, 0.010)

$sFunc = StringFormat("Sqrt((x3 - %.4f)^2 + (y3 - %.4f)^2)", $fX2, $fY2)
_adj_addObsFunction($mSystem, "s_23", $sFunc, 1414.2400, 0.010)

$sFunc = StringFormat("Sqrt((x4 - %.4f)^2 + (y4 - %.4f)^2)", $fX2, $fY2)
_adj_addObsFunction($mSystem, "s_24", $sFunc, 999.9800, 0.010)

$sFunc = "Sqrt((x4 - x3)^2 + (y4 - y3)^2)"
_adj_addObsFunction($mSystem, "s_34", $sFunc, 1000.0000, 0.010)

; Approximate values for new points
_adj_setInitialValue($mSystem, "x3",    0.0)
_adj_setInitialValue($mSystem, "y3",    0.0)
_adj_setInitialValue($mSystem, "x4", 1000.0)
_adj_setInitialValue($mSystem, "y4",    0.0)


; -- Adjustment ---------------------------------------------------------------------
_adj_solve($mSystem, _adj_defaultConfig("GN", False))
If @error Then
	ConsoleWrite("!> ERROR in _adj_solve: @error=" & @error & " @extended=" & @extended & @CRLF)
	Exit 1
EndIf


; -- Results ------------------------------------------------------------------------
ConsoleWrite(_adj_displayResults($mSystem))


; -- Validation vs. reference values (Krumm 2020, pp. 48-50) ------------------------
;
; Residual sign convention:
;   UDF:   v = l̂ − l   (adjusted minus observed)
;   Krumm: ê = l − l̂   (observed minus adjusted)  →  v_UDF = −ê_Krumm

Local $mRes = _adj_getResults($mSystem)
Local $mX1  = $mRes.x1
Local $mSdx = $mRes.sdx
Local $mV   = $mRes.v

ConsoleWrite(@CRLF)
ConsoleWrite("===================================================" & @CRLF)
ConsoleWrite("  Validation vs. Reference (Krumm 2020)            " & @CRLF)
ConsoleWrite("===================================================" & @CRLF & @CRLF)

; --- Statistical summary ---
ConsoleWrite("Statistical summary:" & @CRLF)
_check("s₀ [mm]",  $mRes.s0 * 1000, 6.882, 1)
_check("f",        $mRes.f,         1,     0)

; --- Adjusted coordinates [m] ---
ConsoleWrite(@CRLF & "Adjusted coordinates [m]:" & @CRLF)
_check("x₃", $mX1["X3"],  -0.0096, 4)
_check("y₃", $mX1["Y3"],  -0.0226, 4)
_check("x₄", $mX1["X4"], 999.9930, 4)
_check("y₄", $mX1["Y4"],   0.0174, 4)

; --- Standard deviations [cm] ---
ConsoleWrite(@CRLF & "Standard deviations σ̂ [cm]:" & @CRLF)
_check("σ̂(x₃)", $mSdx["X3"] * 100, 0.901, 2)
_check("σ̂(y₃)", $mSdx["Y3"] * 100, 0.637, 2)
_check("σ̂(x₄)", $mSdx["X4"] * 100, 0.901, 2)
_check("σ̂(y₄)", $mSdx["Y4"] * 100, 0.637, 2)

; --- Residuals v [cm]  (v = l^ - l; signs reversed vs. Krumm) ---
ConsoleWrite(@CRLF & "Residuals v [cm]  (v = l̂ − l; Krumm: ê = −v):" & @CRLF)
_check("v(s_13)", $mV["S_13"] * 100,  0.260, 3)
_check("v(s_14)", $mV["S_14"] * 100, -0.368, 3)
_check("v(s_23)", $mV["S_23"] * 100, -0.368, 3)
_check("v(s_24)", $mV["S_24"] * 100,  0.260, 3)
_check("v(s_34)", $mV["S_34"] * 100,  0.260, 3)


; -- Helper -------------------------------------------------------------------------
; Compares a computed value against a reference value.
; $iDez controls display precision and tolerance (±0.5 units of last displayed digit).
Func _check($sName, $fActual, $fExpected, $iDez)
	Local $sFmt  = "%." & $iDez & "f"
	Local $fTol  = 0.5 * (10 ^ (-$iDez))
	Local $fDiff = Abs($fActual - $fExpected)
	Local $sOK   = ($fDiff <= $fTol) ? "✓" : "✗ MISMATCH"
	ConsoleWrite(StringFormat("  %-12s  expect: " & $sFmt & "  got: " & $sFmt & "  [%s]", _
		$sName, $fExpected, $fActual, $sOK) & @CRLF)
EndFunc
