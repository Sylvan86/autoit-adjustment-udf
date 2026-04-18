#AutoIt3Wrapper_Change2CUI=y

; ===============================================================================
;   Eckerle4 — NIST Nonlinear Regression, Higher Difficulty
; ===============================================================================
;
; Source:
;   NIST/ITL Statistical Reference Datasets (StRD) for Nonlinear Least Squares.
;   Dataset: Eckerle4.
;   URL: https://www.itl.nist.gov/div898/strd/nls/data/eckerle4.shtml
;   Original: Eckerle, K., NIST (1979). Circular Interference Transmittance Study.
;
; Problem:
;   Circular interference transmittance study.
;   Gaussian peak model with 3 parameters:
;     y = (B1/B2)*Exp(-0.5*((x - B3)/B2)^2)
;   35 observations, 3 unknowns.
;
;   Classified as "Higher Difficulty" by NIST — numerically challenging
;   due to narrow peak and sensitive parameterization.
;
; Classification:
;   Model:       Observation equations (Gauss-Markov), nonlinear
;   Weighting:   Unweighted (OLS)
;   Linearity:   Nonlinear
;   Algorithm:   Levenberg-Marquardt
;   VCE:         No
;   Redundancy:  f = n - u = 35 - 3 = 32
;
; ===============================================================================


; -- Initialization --------------------------------------------------------------

#include "../Adjustment.au3"


; -- Model setup -----------------------------------------------------------------
Global $mSystem

; 35 data points (y, x)
Global $aData[][2] = [ _
	[     0.0001575,            400], _
	[     0.0001699,            405], _
	[      0.000235,            410], _
	[     0.0003102,            415], _
	[     0.0004917,            420], _
	[      0.000871,            425], _
	[     0.0017418,            430], _
	[       0.00464,            435], _
	[     0.0065895,          436.5], _
	[     0.0097302,            438], _
	[     0.0149002,          439.5], _
	[      0.023731,            441], _
	[     0.0401683,          442.5], _
	[     0.0712559,            444], _
	[     0.1264458,          445.5], _
	[     0.2073413,            447], _
	[     0.2902366,          448.5], _
	[     0.3445623,            450], _
	[     0.3698049,          451.5], _
	[     0.3668534,            453], _
	[     0.3106727,          454.5], _
	[     0.2078154,            456], _
	[     0.1164354,          457.5], _
	[     0.0616764,            459], _
	[       0.03372,          460.5], _
	[     0.0194023,            462], _
	[     0.0117831,          463.5], _
	[     0.0074357,            465], _
	[     0.0022732,            470], _
	[       0.00088,            475], _
	[     0.0004579,            480], _
	[     0.0002345,            485], _
	[     0.0001586,            490], _
	[     0.0001143,            495], _
	[       7.1e-05,            500]  _
]

; Build observation equations
For $i = 0 To UBound($aData) - 1
	Global $sX = StringFormat("%.15g", $aData[$i][1])
	Global $sFormula = "(B1/B2)*Exp(-0.5*((" & $sX & " - B3)/B2)^2)"
	_adj_addObsFunction($mSystem, "y" & ($i + 1), $sFormula, $aData[$i][0])
Next

; Starting values (NIST Start 2)
_adj_setInitialValue($mSystem, "B1", 1.5)
_adj_setInitialValue($mSystem, "B2", 5)
_adj_setInitialValue($mSystem, "B3", 450)


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
;   b1 =     1.5543827178E+00    sd(b1) = 1.5408051163E-02
;   b2 =     4.0888321754E+00    sd(b2) = 4.6803020753E-02
;   b3 =     4.5154121844E+02    sd(b3) = 4.6800518816E-02
;   Residual SS = 1.4635887487E-03

Global $mRes = _adj_getResults($mSystem)
Global $mX1  = $mRes.x1
Global $mSdx = $mRes.sdx

ConsoleWrite(@CRLF)
ConsoleWrite("+===================================================+" & @CRLF)
ConsoleWrite("|  Validation vs. NIST Certified Values              |" & @CRLF)
ConsoleWrite("+===================================================+" & @CRLF & @CRLF)

; --- Adjusted parameters ---
ConsoleWrite("Adjusted parameters:" & @CRLF)
_check("b1", $mX1["B1"], 1.5543827178E+00, 5)
_check("b2", $mX1["B2"], 4.0888321754E+00, 5)
_check("b3", $mX1["B3"], 4.5154121844E+02, 3)

; --- Standard deviations ---
ConsoleWrite(@CRLF & "Standard deviations:" & @CRLF)
_check("sd(b1)", $mSdx["B1"], 1.5408051163E-02, 6)
_check("sd(b2)", $mSdx["B2"], 4.6803020753E-02, 6)
_check("sd(b3)", $mSdx["B3"], 4.6800518816E-02, 6)

; --- Statistical summary ---
ConsoleWrite(@CRLF & "Statistical summary:" & @CRLF)
_check("s0",   $mRes.s0,   6.7629245447E-03, 7)
_check("vtPv", $mRes.vtpv, 1.4635887487E-03, 8)
_check("f",    $mRes.f,    32, 0)


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
