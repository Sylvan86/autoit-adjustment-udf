#AutoIt3Wrapper_Change2CUI=y

; ===============================================================================
;   MGH09 — NIST Nonlinear Regression, Higher Difficulty
; ===============================================================================
;
; Source:
;   NIST/ITL Statistical Reference Datasets (StRD) for Nonlinear Least Squares.
;   Dataset: MGH09.
;   URL: https://www.itl.nist.gov/div898/strd/nls/data/mgh09.shtml
;   Original: More, Garbow & Hillstrom (1981), Testing unconstrained optimization software, ACM TOMS 7(1), pp. 17-41.
;
; Problem:
;   Kowalik and Osborne rational function.
;   Rational model with 4 parameters:
;     y = B1*(x^2 + x*B2)/(x^2 + x*B3 + B4)
;   11 observations, 4 unknowns.
;
;   Classified as "Higher Difficulty" by NIST — small dataset with
;   rational model structure makes optimization sensitive to starting values.
;
; Classification:
;   Model:       Observation equations (Gauss-Markov), nonlinear
;   Weighting:   Unweighted (OLS)
;   Linearity:   Nonlinear
;   Algorithm:   Levenberg-Marquardt
;   VCE:         No
;   Redundancy:  f = n - u = 11 - 4 = 7
;
; ===============================================================================


; -- Initialization --------------------------------------------------------------
#include "../Adjustment.au3"


; -- Model setup -----------------------------------------------------------------
Local $mSystem

; 11 data points (y, x)
Local $aData[][2] = [ _
	[        0.1957,              4], _
	[        0.1947,              2], _
	[        0.1735,              1], _
	[          0.16,            0.5], _
	[        0.0844,           0.25], _
	[        0.0627,          0.167], _
	[        0.0456,          0.125], _
	[        0.0342,            0.1], _
	[        0.0323,         0.0833], _
	[        0.0235,         0.0714], _
	[        0.0246,         0.0625]  _
]

; Build observation equations
For $i = 0 To UBound($aData) - 1
	Local $sX = StringFormat("%.15g", $aData[$i][1])
	Local $sFormula = "B1*(" & $sX & "^2 + " & $sX & "*B2)/(" & $sX & "^2 + " & $sX & "*B3 + B4)"
	_adj_addObsFunction($mSystem, "y" & ($i + 1), $sFormula, $aData[$i][0])
Next

; Starting values (NIST Start 2)
_adj_setInitialValue($mSystem, "B1", 0.25)
_adj_setInitialValue($mSystem, "B2", 0.39)
_adj_setInitialValue($mSystem, "B3", 0.415)
_adj_setInitialValue($mSystem, "B4", 0.39)


; -- Adjustment ------------------------------------------------------------------
Local $mConfig = _adj_defaultConfig("LM", False)
Local $mDiagCfg = $mConfig.diagnostics
$mDiagCfg.testBasis = "pope"
$mConfig.diagnostics = $mDiagCfg
_adj_solve($mSystem, $mConfig)
If @error Then
	ConsoleWrite("!> ERROR in _adj_solve: @error=" & @error & " @extended=" & @extended & @CRLF)
	Exit 1
EndIf


; -- Results ---------------------------------------------------------------------
Local $mDisplay = _adj_defaultDisplayConfig()
$mDisplay.obsCols = "name|value|v|r|w|T|pPope|decision|mdb"
$mDisplay.showGlobalTest = True
ConsoleWrite(_adj_displayResults($mSystem, $mDisplay))


; -- Validation vs. NIST certified values ----------------------------------------
; NIST certified values (4 parameters):
;   b1 =     1.9280693458E-01    sd(b1) = 1.1435312227E-02
;   b2 =     1.9128232873E-01    sd(b2) = 1.9633220911E-01
;   b3 =     1.2305650693E-01    sd(b3) = 8.0842031232E-02
;   b4 =     1.3606233068E-01    sd(b4) = 9.0025542308E-02
;   Residual SS = 3.0750560385E-04

Local $mRes = _adj_getResults($mSystem)
Local $mX1  = $mRes.x1
Local $mSdx = $mRes.sdx

ConsoleWrite(@CRLF)
ConsoleWrite("+===================================================+" & @CRLF)
ConsoleWrite("|  Validation vs. NIST Certified Values              |" & @CRLF)
ConsoleWrite("+===================================================+" & @CRLF & @CRLF)

; --- Adjusted parameters ---
ConsoleWrite("Adjusted parameters:" & @CRLF)
_check("b1", $mX1["B1"], 1.9280693458E-01, 6)
_check("b2", $mX1["B2"], 1.9128232873E-01, 6)
_check("b3", $mX1["B3"], 1.2305650693E-01, 6)
_check("b4", $mX1["B4"], 1.3606233068E-01, 6)

; --- Standard deviations ---
ConsoleWrite(@CRLF & "Standard deviations:" & @CRLF)
_check("sd(b1)", $mSdx["B1"], 1.1435312227E-02, 6)
_check("sd(b2)", $mSdx["B2"], 1.9633220911E-01, 5)
_check("sd(b3)", $mSdx["B3"], 8.0842031232E-02, 6)
_check("sd(b4)", $mSdx["B4"], 9.0025542308E-02, 6)

; --- Statistical summary ---
ConsoleWrite(@CRLF & "Statistical summary:" & @CRLF)
_check("s0",   $mRes.s0,   6.6279236551E-03, 7)
_check("vtPv", $mRes.vtpv, 3.0750560385E-04, 8)
_check("f",    $mRes.f,    7, 0)


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
