#AutoIt3Wrapper_Change2CUI=y

; ===============================================================================
;   MGH10 — NIST Nonlinear Regression, Higher Difficulty
; ===============================================================================
;
; Source:
;   NIST/ITL Statistical Reference Datasets (StRD) for Nonlinear Least Squares.
;   Dataset: MGH10.
;   URL: https://www.itl.nist.gov/div898/strd/nls/data/mgh10.shtml
;   Original: More, Garbow & Hillstrom (1981), Testing unconstrained optimization software, ACM TOMS 7(1), pp. 17-41.
;
; Problem:
;   Meyer function (exponential).
;   Exponential model with 3 parameters:
;     y = B1*Exp(B2/(x + B3))
;   16 observations, 3 unknowns.
;
;   Classified as "Higher Difficulty" by NIST — exponential model with
;   large dynamic range causes convergence difficulties.
;
; Classification:
;   Model:       Observation equations (Gauss-Markov), nonlinear
;   Weighting:   Unweighted (OLS)
;   Linearity:   Nonlinear
;   Algorithm:   Levenberg-Marquardt
;   VCE:         No
;   Redundancy:  f = n - u = 16 - 3 = 13
;
; ===============================================================================


; -- Initialization --------------------------------------------------------------
#include "../Adjustment.au3"


; -- Model setup -----------------------------------------------------------------
Global $mSystem

; 16 data points (y, x)
Global $aData[][2] = [ _
	[         34780,             50], _
	[         28610,             55], _
	[         23650,             60], _
	[         19630,             65], _
	[         16370,             70], _
	[         13720,             75], _
	[         11540,             80], _
	[          9744,             85], _
	[          8261,             90], _
	[          7030,             95], _
	[          6005,            100], _
	[          5147,            105], _
	[          4427,            110], _
	[          3820,            115], _
	[          3307,            120], _
	[          2872,            125]  _
]

; Build observation equations
For $i = 0 To UBound($aData) - 1
	Global $sX = StringFormat("%.15g", $aData[$i][1])
	Global $sFormula = "B1*Exp(B2/(" & $sX & " + B3))"
	_adj_addObsFunction($mSystem, "y" & ($i + 1), $sFormula, $aData[$i][0])
Next

; Starting values (NIST Start 2)
_adj_setInitialValue($mSystem, "B1", 0.02)
_adj_setInitialValue($mSystem, "B2", 4000)
_adj_setInitialValue($mSystem, "B3", 250)


; -- Adjustment ------------------------------------------------------------------
Global $mConfig = _adj_defaultConfig("LM", False)
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
; NIST certified values (3 parameters):
;   b1 =     5.6096364710E-03    sd(b1) = 1.5687892471E-04
;   b2 =     6.1813463463E+03    sd(b2) = 2.3309021107E+01
;   b3 =     3.4522363462E+02    sd(b3) = 7.8486103508E-01
;   Residual SS = 8.7945855171E+01

Global $mRes = _adj_getResults($mSystem)
Global $mX1  = $mRes.x1
Global $mSdx = $mRes.sdx

ConsoleWrite(@CRLF)
ConsoleWrite("+===================================================+" & @CRLF)
ConsoleWrite("|  Validation vs. NIST Certified Values              |" & @CRLF)
ConsoleWrite("+===================================================+" & @CRLF & @CRLF)

; --- Adjusted parameters ---
ConsoleWrite("Adjusted parameters:" & @CRLF)
_check("b1", $mX1["B1"], 5.6096364710E-03, 8)
_check("b2", $mX1["B2"], 6.1813463463E+03, 2)
_check("b3", $mX1["B3"], 3.4522363462E+02, 3)

; --- Standard deviations ---
ConsoleWrite(@CRLF & "Standard deviations:" & @CRLF)
_check("sd(b1)", $mSdx["B1"], 1.5687892471E-04, 8)
_check("sd(b2)", $mSdx["B2"], 2.3309021107E+01, 3)
_check("sd(b3)", $mSdx["B3"], 7.8486103508E-01, 5)

; --- Statistical summary ---
ConsoleWrite(@CRLF & "Statistical summary:" & @CRLF)
_check("s0",   $mRes.s0,   2.6009740065E+00, 4)
_check("vtPv", $mRes.vtpv, 8.7945855171E+01, 3)
_check("f",    $mRes.f,    13, 0)


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
