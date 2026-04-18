#AutoIt3Wrapper_Change2CUI=y

; ===============================================================================
;   Thurber — NIST Nonlinear Regression, Higher Difficulty
; ===============================================================================
;
; Source:
;   NIST/ITL Statistical Reference Datasets (StRD) for Nonlinear Least Squares.
;   Dataset: Thurber.
;   URL: https://www.itl.nist.gov/div898/strd/nls/data/thurber.shtml
;   Original: Thurber, R., NIST — Semiconductor electron mobility modeling.
;
; Problem:
;   Electron mobility in semiconductors as function of log density.
;   Rational model (cubic/cubic) with 7 parameters:
;     y = (b1 + b2*x + b3*x^2 + b4*x^3) / (1 + b5*x + b6*x^2 + b7*x^3)
;   37 observations, 7 unknowns.
;
;   Classified as "Higher Difficulty" by NIST — numerically challenging
;   due to high parameter count and rational model structure.
;
; Classification:
;   Model:       Observation equations (Gauss-Markov), nonlinear
;   Weighting:   Unweighted (OLS)
;   Linearity:   Nonlinear
;   Algorithm:   Levenberg-Marquardt
;   VCE:         No
;   Redundancy:  f = n - u = 37 - 7 = 30
;
; ===============================================================================


; -- Initialization --------------------------------------------------------------
#include "../Adjustment.au3"


; -- Model setup -----------------------------------------------------------------
Global $mSystem

; Observations: y = f(x) + e
; Functional model:  y = (b1 + b2*x + b3*x^2 + b4*x^3) / (1 + b5*x + b6*x^2 + b7*x^3)
;
; 37 data points (electron mobility vs. log density)
Global $aData[][2] = [ _
	[ 80.574,  -3.067], _
	[ 84.248,  -2.981], _
	[ 87.264,  -2.921], _
	[ 87.195,  -2.912], _
	[ 89.076,  -2.840], _
	[ 89.608,  -2.797], _
	[ 89.868,  -2.702], _
	[ 90.101,  -2.699], _
	[ 92.405,  -2.633], _
	[ 95.854,  -2.481], _
	[100.696,  -2.363], _
	[101.060,  -2.322], _
	[401.672,  -1.501], _
	[390.724,  -1.460], _
	[567.534,  -1.274], _
	[635.316,  -1.212], _
	[733.054,  -1.100], _
	[759.087,  -1.046], _
	[894.206,  -0.915], _
	[990.785,  -0.714], _
	[1090.109, -0.566], _
	[1080.914, -0.545], _
	[1122.643, -0.400], _
	[1178.351, -0.309], _
	[1260.531, -0.109], _
	[1273.514, -0.103], _
	[1288.339,  0.010], _
	[1327.543,  0.119], _
	[1353.863,  0.377], _
	[1414.509,  0.790], _
	[1425.208,  0.963], _
	[1421.384,  1.006], _
	[1442.962,  1.115], _
	[1464.350,  1.572], _
	[1468.705,  1.841], _
	[1447.894,  2.047], _
	[1457.628,  2.200]  _
]

; Build observation equations with analytical derivatives
; Formula: y = p/q  where  p = B1+B2*x+B3*x²+B4*x³,  q = 1+B5*x+B6*x²+B7*x³
; Derivatives (quotient rule):
;   dF/dBi (i=1..4) = x^(i-1) / q               (numerator params)
;   dF/dBi (i=5..7) = -x^(i-4) * p / q²         (denominator params)
For $i = 0 To UBound($aData) - 1
	Global $fX = $aData[$i][1]
	Global $fY = $aData[$i][0]
	Global $sX = StringFormat("%.3f", $fX)

	; Numerator p and denominator q as subexpressions
	Global $sP = "(B1 + B2*(" & $sX & ") + B3*(" & $sX & ")^2 + B4*(" & $sX & ")^3)"
	Global $sQ = "(1 + B5*(" & $sX & ") + B6*(" & $sX & ")^2 + B7*(" & $sX & ")^3)"

	Global $sFormula = $sP & " / " & $sQ

	; Analytical derivatives (note: use -1* instead of - prefix to avoid precedence issues with ^)
	Global $sDeriv = "B1 = 1 / " & $sQ _
	    & " | B2 = (" & $sX & ") / " & $sQ _
	    & " | B3 = (" & $sX & ")^2 / " & $sQ _
	    & " | B4 = (" & $sX & ")^3 / " & $sQ _
	    & " | B5 = -1 * (" & $sX & ") * " & $sP & " / " & $sQ & "^2" _
	    & " | B6 = -1 * (" & $sX & ")^2 * " & $sP & " / " & $sQ & "^2" _
	    & " | B7 = -1 * (" & $sX & ")^3 * " & $sP & " / " & $sQ & "^2"
	_adj_addObsFunction($mSystem, "y" & ($i + 1), $sFormula, $fY, 1, "s0", 0.0, $sDeriv)
Next

; Starting values (NIST Start 2 — "far" / higher difficulty)
_adj_setInitialValue($mSystem, "B1", 1300)
_adj_setInitialValue($mSystem, "B2", 1500)
_adj_setInitialValue($mSystem, "B3", 500)
_adj_setInitialValue($mSystem, "B4", 75)
_adj_setInitialValue($mSystem, "B5", 1)
_adj_setInitialValue($mSystem, "B6", 0.4)
_adj_setInitialValue($mSystem, "B7", 0.05)


; -- Adjustment ------------------------------------------------------------------
Global $mConfig = _adj_defaultConfig("LM", False)
$mConfig.solver = "SVD"
Global $mDiagCfg = $mConfig.diagnostics
$mDiagCfg.testBasis = "pope"
$mConfig.diagnostics = $mDiagCfg
_adj_solve($mSystem, $mConfig)
If @error Then
	ConsoleWrite("!> ERROR in _adj_solve: @error=" & @error & " @extended=" & @extended & @CRLF)
	Exit 1
EndIf


; -- Results ---------------------------------------------------------------------
Global $mDisplay = _adj_defaultDisplayConfig()
$mDisplay.obsCols = "name|value|v|r|w|T|pPope|decision|mdb"
$mDisplay.showGlobalTest = True
ConsoleWrite(_adj_displayResults($mSystem, $mDisplay))


; -- Validation vs. NIST certified values ----------------------------------------
; NIST certified values (11 significant digits):
;   b1 = 1.2881396800E+03    sd(b1) = 4.6647963344E+00
;   b2 = 1.4910792535E+03    sd(b2) = 3.9571156086E+01
;   b3 = 5.8323836877E+02    sd(b3) = 2.8698696102E+01
;   b4 = 7.5416644291E+01    sd(b4) = 5.5675370270E+00
;   b5 = 9.6629502864E-01    sd(b5) = 3.1333340687E-02
;   b6 = 3.9797285797E-01    sd(b6) = 1.4984928198E-02
;   b7 = 4.9727297349E-02    sd(b7) = 6.5842344623E-03
;   Residual SS = 5.6427082397E+03

Global $mRes = _adj_getResults($mSystem)
Global $mX1  = $mRes.x1
Global $mSdx = $mRes.sdx

ConsoleWrite(@CRLF)
ConsoleWrite("+===================================================+" & @CRLF)
ConsoleWrite("|  Validation vs. NIST Certified Values              |" & @CRLF)
ConsoleWrite("+===================================================+" & @CRLF & @CRLF)

; --- Adjusted parameters ---
ConsoleWrite("Adjusted parameters:" & @CRLF)
_check("b1", $mX1["B1"], 1.2881396800E+03, 3)
_check("b2", $mX1["B2"], 1.4910792535E+03, 3)
_check("b3", $mX1["B3"],  5.8323836877E+02, 3)
_check("b4", $mX1["B4"],   7.5416644291E+01, 3)
_check("b5", $mX1["B5"],    9.6629502864E-01, 6)
_check("b6", $mX1["B6"],    3.9797285797E-01, 6)
_check("b7", $mX1["B7"],    4.9727297349E-02, 6)

; --- Standard deviations ---
ConsoleWrite(@CRLF & "Standard deviations:" & @CRLF)
_check("sd(b1)", $mSdx["B1"],  4.6647963344E+00, 2)
_check("sd(b2)", $mSdx["B2"], 3.9571156086E+01, 1)
_check("sd(b3)", $mSdx["B3"], 2.8698696102E+01, 1)
_check("sd(b4)", $mSdx["B4"],  5.5675370270E+00, 2)
_check("sd(b5)", $mSdx["B5"],  3.1333340687E-02, 3)
_check("sd(b6)", $mSdx["B6"],  1.4984928198E-02, 3)
_check("sd(b7)", $mSdx["B7"],  6.5842344623E-03, 3)

; --- Statistical summary ---
ConsoleWrite(@CRLF & "Statistical summary:" & @CRLF)
_check("s0",   $mRes.s0,   1.3714600784E+01, 4)
_check("vtPv", $mRes.vtpv, 5.6427082397E+03, 2)
_check("f",    $mRes.f,    30, 0)
_check("iter", $mRes.nIterations, -1, 0)   ; show iteration count (no reference)


; -- Helper ----------------------------------------------------------------------
; Compares a computed value against a reference value.
; $iDez controls display precision and tolerance (+-0.5 units of last displayed digit).
; Also shows the absolute deviation.
Func _check($sName, $fActual, $fExpected, $iDez)
	Local $sFmt  = "%." & $iDez & "f"
	Local $fTol  = 0.5 * (10 ^ (-$iDez))
	Local $fDiff = $fActual - $fExpected
	Local $sOK   = (Abs($fDiff) <= $fTol) ? "OK" : ($fExpected = -1 ? "" : "MISMATCH")
	If $fExpected = -1 Then
		ConsoleWrite(StringFormat("  %-12s  got: " & $sFmt, $sName, $fActual) & @CRLF)
	Else
		ConsoleWrite(StringFormat("  %-12s  expect: " & $sFmt & "  got: " & $sFmt & "  dev: %+.2e  [%s]", _
			$sName, $fExpected, $fActual, $fDiff, $sOK) & @CRLF)
	EndIf
EndFunc
