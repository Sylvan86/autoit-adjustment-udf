#AutoIt3Wrapper_Change2CUI=y

; ===============================================================================
;   Direction Network with Distance Constraint -- Wolf / Stuttgart
; ===============================================================================
;
; Source:
;   Sneeuw, N., Krumm, F., Roth, M. (2022): Lecture Notes Adjustment
;   Theory, Rev. 4.48, Geodaetisches Institut, Univ. Stuttgart.
;   Example 5.1.6 (pp. 65-77), based on Wolf (1979, pp. 66-78).
;
; Problem:
;   10-point network with 6 fixed points (A-F) and 3 new points (G,H,I).
;   36 direction observations + 1 distance constraint (G-I = 2121.90 m).
;   9 orientation unknowns + 6 coordinate unknowns = 15 unknowns total.
;   Nonlinear model, solved by Gauss-Newton iteration.
;
; Adaptation:
;   Original source uses 35 directions + 1 distance obs + 1 angle obs
;   (37 obs, DOF = 22).  Here: 36 directions as observations, distance
;   G-I as equality constraint (LSE).  DOF = 36 - 15 + 1 = 22.
;
; Direction computation:
;   The UDF computes Jacobians by numerical differentiation.  Using
;   Mod(bearing, 400) to wrap directions into [0, 400) creates a
;   discontinuity at the 0/400 boundary that produces wrong derivatives
;   when the computed direction is near 0.  Solution: return the raw
;   (bearing - R0) value without Mod, and "unwrap" observed directions
;   to match the range of the formula output.
;
; Classification:
;   Model:       Observation equations with equality constraint (LSE)
;   Weighting:   Weighted (WLSE) -- sigma_r = 2.5 mgon
;   Linearity:   Nonlinear (direction formulas)
;   Algorithm:   Gauss-Newton (iterative)
;   VCE:         No
;   Datum:       6 fixed points (A-F)
;   Constraint:  Distance G-I = 2121.90 m
;   Redundancy:  f = n - u + c = 36 - 15 + 1 = 22
;
; ===============================================================================


; -- Initialization --------------------------------------------------------------
#include "../Adjustment.au3"


; -- Point coordinates [m] (Table 5.42) ------------------------------------------
;    Fixed (A-F): known coordinates.  New (G-I): approximate coordinates.
Local $mPt[]
$mPt["A"] = _pt(184423.28, 726419.33)
$mPt["B"] = _pt(186444.18, 726476.66)
$mPt["C"] = _pt(183257.84, 725490.35)
$mPt["D"] = _pt(184292.00, 723313.00)
$mPt["E"] = _pt(185487.00, 721829.00)
$mPt["F"] = _pt(186708.72, 722104.58)
$mPt["G"] = _pt(184868.20, 725139.70)
$mPt["H"] = _pt(186579.30, 725336.60)
$mPt["I"] = _pt(185963.07, 723322.02)

; New points: coordinates are estimated parameters (XG, YG, XH, YH, XI, YI)
Local $mNewPt[]
$mNewPt["G"] = True
$mNewPt["H"] = True
$mNewPt["I"] = True

; -- Approximate orientation unknowns [gon] (Table 5.43, omega_0) ----------------
Local $mR0[]
$mR0["A"] = 98.196
$mR0["B"] = 192.486
$mR0["C"] = 57.164
$mR0["D"] = 19.448
$mR0["E"] = 19.636
$mR0["F"] = 285.864
$mR0["G"] = 55.215
$mR0["H"] = 197.451
$mR0["I"] = 18.903

Local Const $fSigmaR = 0.0025  ; 2.5 mgon


; -- Direction observations (Table 5.43) -----------------------------------------
;    Each row: [station, target, observed direction in gon]
Local $aDirs[] = [ _
	"A","B",  0.0000,  "A","G", 80.5000,  "A","C",158.9610, _
	"B","H",  0.0000,  "B","G", 62.7260,  "B","A",105.7120, _
	"C","A",  0.0000,  "C","G", 56.4960,  "C","I", 85.8450,  "C","D",114.5950, _
	"D","G",  0.0000,  "D","H", 34.4500,  "D","I", 80.2110,  "D","E",137.4020,  "D","C",352.3090, _
	"E","I",  0.0000,  "E","F", 66.2450,  "E","D",337.2160, _
	"F","E",  0.0000,  "F","I", 79.1690,  "F","H",111.5820, _
	"G","B",  0.0000,  "G","H", 37.4980,  "G","I",110.2580,  "G","D",164.2320,  "G","C",258.4410,  "G","A",323.4860, _
	"H","F",  0.0000,  "H","I", 21.4450,  "H","D", 56.4420, _
	"I","H",  0.0000,  "I","F",146.1430,  "I","E",200.7330,  "I","D",280.7560,  "I","C",324.1050,  "I","G",346.5690 _
]


; -- Model setup -----------------------------------------------------------------
Local $mSystem

; Add all direction observations
For $i = 0 To UBound($aDirs) - 1 Step 3
	Local $sFrom = $aDirs[$i], $sTo = $aDirs[$i + 1], $fDir = $aDirs[$i + 2]
	Local $mF = $mPt[$sFrom], $mT = $mPt[$sTo]

	; Formula: fixed-point coords are literal numbers, new-point coords are variable names
	Local $sFc = MapExists($mNewPt, $sFrom) ? ("X" & $sFrom & ", Y" & $sFrom) : StringFormat("%.2f, %.2f", $mF.x, $mF.y)
	Local $sTc = MapExists($mNewPt, $sTo)   ? ("X" & $sTo   & ", Y" & $sTo)   : StringFormat("%.2f, %.2f", $mT.x, $mT.y)
	Local $sFunc = "_calcDirection_gon(" & $sFc & ", " & $sTc & ", R0" & $sFrom & ")"

	; Unwrap observed direction to match the Mod-free formula range
	Local $fRaw = _calcBearing($mF.x, $mF.y, $mT.x, $mT.y) - $mR0[$sFrom]
	_adj_addObsFunction($mSystem, "r_" & $sFrom & "_" & $sTo, $sFunc, _unwrap($fDir, $fRaw), $fSigmaR)
Next

; Equality constraint: distance G-I = 2121.90 m
_adj_addRestriction($mSystem, "sqrt((XI - XG)^2 + (YI - YG)^2)", 2121.90)

; Initial values for coordinates and orientations
For $s In MapKeys($mNewPt)
	_adj_setInitialValue($mSystem, "X" & $s, ($mPt[$s]).x)
	_adj_setInitialValue($mSystem, "Y" & $s, ($mPt[$s]).y)
Next
For $s In MapKeys($mR0)
	_adj_setInitialValue($mSystem, "R0" & $s, $mR0[$s])
Next


; -- Adjustment ------------------------------------------------------------------
_adj_solve($mSystem, _adj_defaultConfig("GN", False))
If @error Then
	ConsoleWrite("!> ERROR in _adj_solve: @error=" & @error & " @extended=" & @extended & @CRLF)
	Exit 1
EndIf


; -- Results ---------------------------------------------------------------------
ConsoleWrite(_adj_displayResults($mSystem))


; -- Validation vs. reference values (Stuttgart LN, Ex. 5.1.6, pp. 68-69) -------
;
; Note: Reference values are from the full model (35 dirs + 1 dist + 1 angle,
;       DOF = 22).  This adapted model uses 36 dirs + 1 distance constraint
;       (LSE, DOF = 22).  Small deviations from reference values are expected.

Local $mRes = _adj_getResults($mSystem)
Local $mX1  = $mRes.x1

ConsoleWrite(@CRLF)
ConsoleWrite("===================================================" & @CRLF)
ConsoleWrite("  Validation vs. Reference (Stuttgart LN, approx.)  " & @CRLF)
ConsoleWrite("===================================================" & @CRLF & @CRLF)

ConsoleWrite("Adjusted coordinates [m]:" & @CRLF)
_check("X_G", $mX1["XG"], 184868.038, 1)
_check("Y_G", $mX1["YG"], 725139.657, 1)
_check("X_H", $mX1["XH"], 186579.4, 1)
_check("Y_H", $mX1["YH"], 725336.414, 1)
_check("X_I", $mX1["XI"], 185963.215, 1)
_check("Y_I", $mX1["YI"], 723322.303, 1)

ConsoleWrite(@CRLF & "Orientation unknowns [gon]:" & @CRLF)
_check("R0_A", $mX1["R0A"], 98.1987, 2)
_check("R0_B", $mX1["R0B"], 192.4866, 2)
_check("R0_C", $mX1["R0C"], 57.1634, 2)
_check("R0_D", $mX1["R0D"], 19.4452, 2)
_check("R0_E", $mX1["R0E"], 19.6364, 2)
_check("R0_F", $mX1["R0F"], 285.8684, 2)
_check("R0_G", $mX1["R0G"], 55.2150, 2)
_check("R0_H", $mX1["R0H"], 197.4525, 2)
_check("R0_I", $mX1["R0I"], 18.9001, 2)

ConsoleWrite(@CRLF & "Statistical summary:" & @CRLF)
_check("f", $mRes.f, 22, 0)

Local $fDistGI = Sqrt(($mX1["XI"] - $mX1["XG"])^2 + ($mX1["YI"] - $mX1["YG"])^2)
ConsoleWrite(@CRLF & "Constraint check (dist G-I = 2121.90 m):" & @CRLF)
_check("d_GI", $fDistGI, 2121.90, 2)


; -- Helpers ---------------------------------------------------------------------

Func _check($sName, $fActual, $fExpected, $iDez)
	Local $sFmt  = "%." & $iDez & "f"
	Local $fTol  = 0.5 * (10 ^ (-$iDez))
	Local $fDiff = Abs($fActual - $fExpected)
	Local $sOK   = ($fDiff <= $fTol) ? "pass" : "MISMATCH"
	ConsoleWrite(StringFormat("  %-12s  expect: " & $sFmt & "  got: " & $sFmt & "  [%s]", _
		$sName, $fExpected, $fActual, $sOK) & @CRLF)
EndFunc

; Creates a point map with .x and .y fields.
Func _pt($fX, $fY)
	Local $m[]
	$m.x = $fX
	$m.y = $fY
	Return $m
EndFunc

; Geodetic bearing from (Xfrom,Yfrom) to (Xto,Yto) in gon.
; Half-angle ATan formula, returns values in (-200, 200] without discontinuity.
Func _calcBearing($fXfrom, $fYfrom, $fXto, $fYto)
	Local Static $fRho = 200 / ACos(-1)
	Local $dN = $fXto - $fXfrom
	Local $dE = $fYto - $fYfrom
	If Abs($dN) < 1e-15 And Abs($dE) < 1e-15 Then Return 0
	Return 2 * $fRho * ATan($dN / ($dE + Sqrt($dN * $dN + $dE * $dE)))
EndFunc

; Direction = bearing - orientation (no Mod, avoids 0/400 discontinuity).
Func _calcDirection_gon($fXfrom, $fYfrom, $fXto, $fYto, $fR0)
	Return _calcBearing($fXfrom, $fYfrom, $fXto, $fYto) - $fR0
EndFunc

; Unwraps observed direction to match the Mod-free formula range.
Func _unwrap($fObs, $fComputed)
	Return $fObs + 400 * Round(($fComputed - $fObs) / 400)
EndFunc
