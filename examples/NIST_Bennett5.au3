#AutoIt3Wrapper_Change2CUI=y

; ===============================================================================
;   Bennett5 — NIST Nonlinear Regression, Higher Difficulty
; ===============================================================================
;
; Source:
;   NIST/ITL Statistical Reference Datasets (StRD) for Nonlinear Least Squares.
;   Dataset: Bennett5.
;   URL: https://www.itl.nist.gov/div898/strd/nls/data/bennett5.shtml
;   Original: Bennett, Swartzendruber & Brown, NIST (1994). Superconductivity Magnetization Modeling.
;
; Problem:
;   Magnetization modeling of superconductors.
;   Power-law model with 3 parameters:
;     y = B1*(B2 + x)^(-1/B3)
;   154 observations, 3 unknowns.
;
;   Classified as "Higher Difficulty" by NIST — power-law model with
;   large parameter ranges makes optimization challenging.
;
; Classification:
;   Model:       Observation equations (Gauss-Markov), nonlinear
;   Weighting:   Unweighted (OLS)
;   Linearity:   Nonlinear
;   Algorithm:   Levenberg-Marquardt
;   VCE:         No
;   Redundancy:  f = n - u = 154 - 3 = 151
;
; ===============================================================================


; -- Initialization --------------------------------------------------------------
#include "../Adjustment.au3"


; -- Model setup -----------------------------------------------------------------
Global $mSystem

; 154 data points (y, x)
Global $aData[][2] = [ _
	[    -34.834702,       7.447168], _
	[      -34.3932,       8.102586], _
	[    -34.152901,       8.452547], _
	[    -33.979099,       8.711278], _
	[    -33.845901,       8.916774], _
	[    -33.732899,       9.087155], _
	[    -33.640301,        9.23259], _
	[      -33.5592,       9.359535], _
	[    -33.486801,       9.472166], _
	[      -33.4231,       9.573384], _
	[    -33.365101,       9.665293], _
	[       -33.313,       9.749461], _
	[    -33.260899,       9.827092], _
	[      -33.2174,       9.899128], _
	[    -33.176899,       9.966321], _
	[    -33.139198,       10.02928], _
	[    -33.101601,       10.08851], _
	[    -33.066799,       10.14443], _
	[       -33.035,       10.19738], _
	[    -33.003101,       10.24767], _
	[    -32.971298,       10.29556], _
	[    -32.942299,       10.34125], _
	[    -32.916302,       10.38495], _
	[    -32.890202,       10.42682], _
	[    -32.864101,         10.467], _
	[       -32.841,       10.50564], _
	[    -32.817799,       10.54283], _
	[    -32.797501,       10.57869], _
	[      -32.7743,       10.61331], _
	[       -32.757,       10.64678], _
	[    -32.733799,       10.67915], _
	[      -32.7164,       10.71052], _
	[      -32.6991,       10.74092], _
	[    -32.678799,       10.77044], _
	[      -32.6614,        10.7991], _
	[    -32.644001,       10.82697], _
	[    -32.626701,       10.85408], _
	[    -32.612202,       10.88047], _
	[    -32.597698,       10.90619], _
	[    -32.583199,       10.93126], _
	[    -32.568699,       10.95572], _
	[    -32.554298,       10.97959], _
	[    -32.539799,       11.00291], _
	[    -32.525299,        11.0257], _
	[    -32.510799,       11.04798], _
	[    -32.499199,       11.06977], _
	[    -32.487598,        11.0911], _
	[    -32.473202,       11.11198], _
	[    -32.461601,       11.13244], _
	[    -32.435501,       11.15248], _
	[    -32.435501,       11.17213], _
	[      -32.4268,       11.19141], _
	[      -32.4123,       11.21031], _
	[    -32.400799,       11.22887], _
	[    -32.392101,       11.24709], _
	[    -32.380501,       11.26498], _
	[    -32.366001,       11.28256], _
	[      -32.3573,       11.29984], _
	[    -32.348598,       11.31682], _
	[    -32.339901,       11.33352], _
	[      -32.3284,       11.34994], _
	[    -32.319698,        11.3661], _
	[    -32.311001,         11.382], _
	[      -32.2994,       11.39766], _
	[    -32.290699,       11.41307], _
	[    -32.282001,       11.42824], _
	[      -32.2733,        11.4432], _
	[    -32.264599,       11.45793], _
	[    -32.256001,       11.47244], _
	[    -32.247299,       11.48675], _
	[    -32.238602,       11.50086], _
	[      -32.2299,       11.51477], _
	[    -32.224098,       11.52849], _
	[    -32.215401,       11.54202], _
	[      -32.2038,       11.55538], _
	[    -32.198002,       11.56855], _
	[      -32.1894,       11.58156], _
	[    -32.183601,       11.59442], _
	[      -32.1749,      11.607121], _
	[    -32.169102,       11.61964], _
	[      -32.1633,         11.632], _
	[    -32.154598,       11.64421], _
	[    -32.145901,       11.65628], _
	[    -32.140099,        11.6682], _
	[    -32.131401,       11.67998], _
	[    -32.125599,       11.69162], _
	[    -32.119801,       11.70313], _
	[    -32.111198,       11.71451], _
	[      -32.1054,       11.72576], _
	[    -32.096699,       11.73688], _
	[      -32.0909,       11.74789], _
	[    -32.088001,       11.75878], _
	[      -32.0793,       11.76955], _
	[    -32.073502,        11.7802], _
	[    -32.067699,       11.79073], _
	[    -32.061901,       11.80116], _
	[    -32.056099,       11.81148], _
	[    -32.050301,        11.8217], _
	[    -32.044498,       11.83181], _
	[    -32.038799,       11.84182], _
	[    -32.033001,       11.85173], _
	[    -32.027199,       11.86155], _
	[      -32.0243,       11.87127], _
	[    -32.018501,       11.88089], _
	[    -32.012699,       11.89042], _
	[    -32.004002,       11.89987], _
	[    -32.001099,       11.90922], _
	[      -31.9953,       11.91849], _
	[      -31.9895,       11.92768], _
	[      -31.9837,       11.93678], _
	[      -31.9779,       11.94579], _
	[    -31.972099,       11.95473], _
	[    -31.969299,       11.96359], _
	[    -31.963501,       11.97237], _
	[    -31.957701,       11.98107], _
	[      -31.9519,        11.9897], _
	[      -31.9461,       11.99826], _
	[      -31.9403,       12.00674], _
	[    -31.937401,       12.01515], _
	[    -31.931601,       12.02349], _
	[      -31.9258,       12.03176], _
	[    -31.922899,       12.03997], _
	[    -31.917101,        12.0481], _
	[    -31.911301,       12.05617], _
	[      -31.9084,       12.06418], _
	[    -31.902599,       12.07212], _
	[      -31.8969,       12.08001], _
	[    -31.893999,       12.08782], _
	[    -31.888201,       12.09558], _
	[      -31.8853,       12.10328], _
	[    -31.882401,       12.11092], _
	[      -31.8766,        12.1185], _
	[    -31.873699,       12.12603], _
	[    -31.867901,        12.1335], _
	[    -31.862101,       12.14091], _
	[      -31.8592,       12.14827], _
	[      -31.8563,       12.15557], _
	[      -31.8505,       12.16283], _
	[      -31.8447,       12.17003], _
	[    -31.841801,       12.17717], _
	[      -31.8389,       12.18427], _
	[    -31.833099,       12.19132], _
	[      -31.8302,       12.19832], _
	[    -31.827299,       12.20527], _
	[      -31.8216,       12.21217], _
	[    -31.818701,       12.21903], _
	[    -31.812901,       12.22584], _
	[    -31.809999,        12.2326], _
	[      -31.8071,       12.23932], _
	[      -31.8013,       12.24599], _
	[    -31.798401,       12.25262], _
	[      -31.7955,        12.2592], _
	[      -31.7897,       12.26575], _
	[      -31.7868,       12.27224]  _
]

; Build observation equations
For $i = 0 To UBound($aData) - 1
	Global $sX = StringFormat("%.15g", $aData[$i][1])
	Global $sFormula = "B1*(B2 + " & $sX & ")^(-1/B3)"
	_adj_addObsFunction($mSystem, "y" & ($i + 1), $sFormula, $aData[$i][0])
Next

; Starting values (NIST Start 2)
_adj_setInitialValue($mSystem, "B1", -1500)
_adj_setInitialValue($mSystem, "B2", 45)
_adj_setInitialValue($mSystem, "B3", 0.85)


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
; NIST certified values (3 parameters):
;   b1 =    -2.5235058043E+03    sd(b1) = 2.9715175411E+02
;   b2 =     4.6736564644E+01    sd(b2) = 1.2448871856E+00
;   b3 =     9.3218483193E-01    sd(b3) = 2.0272299378E-02
;   Residual SS = 5.2404744073E-04

Global $mRes = _adj_getResults($mSystem)
Global $mX1  = $mRes.x1
Global $mSdx = $mRes.sdx

ConsoleWrite(@CRLF)
ConsoleWrite("+===================================================+" & @CRLF)
ConsoleWrite("|  Validation vs. NIST Certified Values              |" & @CRLF)
ConsoleWrite("+===================================================+" & @CRLF & @CRLF)

; --- Adjusted parameters ---
ConsoleWrite("Adjusted parameters:" & @CRLF)
_check("b1", $mX1["B1"], -2.5235058043E+03, 2)
_check("b2", $mX1["B2"], 4.6736564644E+01, 4)
_check("b3", $mX1["B3"], 9.3218483193E-01, 6)

; --- Standard deviations ---
ConsoleWrite(@CRLF & "Standard deviations:" & @CRLF)
_check("sd(b1)", $mSdx["B1"], 2.9715175411E+02, 2)
_check("sd(b2)", $mSdx["B2"], 1.2448871856E+00, 4)
_check("sd(b3)", $mSdx["B3"], 2.0272299378E-02, 6)

; --- Statistical summary ---
ConsoleWrite(@CRLF & "Statistical summary:" & @CRLF)
_check("s0",   $mRes.s0,   1.8629312528E-03, 7)
_check("vtPv", $mRes.vtpv, 5.2404744073E-04, 8)
_check("f",    $mRes.f,    151, 0)


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
