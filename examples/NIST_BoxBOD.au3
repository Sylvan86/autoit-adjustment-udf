#AutoIt3Wrapper_Change2CUI=y

; ===============================================================================
;   BoxBOD — NIST Nonlinear Regression, Higher Difficulty
; ===============================================================================
;
; Source:
;   NIST/ITL Statistical Reference Datasets (StRD) for Nonlinear Least Squares.
;   Dataset: BoxBOD.
;   URL: https://www.itl.nist.gov/div898/strd/nls/data/boxbod.shtml
;   Original: Box, Hunter & Hunter (1978), Statistics for Experimenters, pp. 483-487.
;
; Problem:
;   Biochemical oxygen demand (BOD) as function of incubation time.
;   Exponential saturation model with 2 parameters:
;     y = b1 * (1 - exp(-b2 * x))
;   6 observations, 2 unknowns (b1, b2).
;
;   Classified as "Higher Difficulty" by NIST — distant starting values
;   cause naive Gauss-Newton to diverge. Requires Levenberg-Marquardt.
;
; Classification:
;   Model:       Observation equations (Gauss-Markov), nonlinear
;   Weighting:   Unweighted (OLS)
;   Linearity:   Nonlinear
;   Algorithm:   Levenberg-Marquardt
;   VCE:         No
;   Redundancy:  f = n - u = 6 - 2 = 4
;
; ===============================================================================


; -- Initialization --------------------------------------------------------------
#include "../Adjustment.au3"


; -- Model setup -----------------------------------------------------------------
Global $mSystem

; Observations: y = f(x) + e
; Functional model:  y = b1 * (1 - exp(-b2 * x))
;
;   Obs.     x (days)    y (BOD, mg/l)     sigma
;   y1       1           109               1
;   y2       2           149               1
;   y3       3           149               1
;   y4       5           191               1
;   y5       7           213               1
;   y6       10          224               1

_adj_addObsFunction($mSystem, "y1", "B1 * (1 - Exp(-B2 * 1))",   109)
_adj_addObsFunction($mSystem, "y2", "B1 * (1 - Exp(-B2 * 2))",   149)
_adj_addObsFunction($mSystem, "y3", "B1 * (1 - Exp(-B2 * 3))",   149)
_adj_addObsFunction($mSystem, "y4", "B1 * (1 - Exp(-B2 * 5))",   191)
_adj_addObsFunction($mSystem, "y5", "B1 * (1 - Exp(-B2 * 7))",   213)
_adj_addObsFunction($mSystem, "y6", "B1 * (1 - Exp(-B2 * 10))",  224)

; Starting values (NIST Start 2 — "far" / higher difficulty)
_adj_setInitialValue($mSystem, "B1", 100)
_adj_setInitialValue($mSystem, "B2", 0.75)


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
Global $mRes = _adj_getResults($mSystem)
Global $mX1  = $mRes.x1
Global $mSdx = $mRes.sdx

ConsoleWrite(@CRLF)
ConsoleWrite("+===================================================+" & @CRLF)
ConsoleWrite("|  Validation vs. NIST Certified Values              |" & @CRLF)
ConsoleWrite("+===================================================+" & @CRLF & @CRLF)

; --- Adjusted parameters ---
ConsoleWrite("Adjusted parameters:" & @CRLF)
_check("b1", $mX1["B1"], 213.80940889, 4)
_check("b2", $mX1["B2"],   0.54723749, 6)

; --- Standard deviations ---
ConsoleWrite(@CRLF & "Standard deviations:" & @CRLF)
_check("sd(b1)", $mSdx["B1"], 12.3545, 2)
_check("sd(b2)", $mSdx["B2"],  0.1046, 3)

; --- Statistical summary ---
ConsoleWrite(@CRLF & "Statistical summary:" & @CRLF)
_check("s0",   $mRes.s0,   17.0881, 2)
_check("vtPv", $mRes.vtpv, 1168.0089, 2)
_check("f",    $mRes.f,    4, 0)


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
