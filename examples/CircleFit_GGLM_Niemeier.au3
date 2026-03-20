#AutoIt3Wrapper_Change2CUI=y

; ===============================================================================
;   Circle Fitting from Polar Coordinates -- Niemeier (WGLM + GGLM)
; ===============================================================================
;
; Source:
;   Niemeier, W. (2008): Ausgleichungsrechnung, 2nd ed., de Gruyter.
;   Section 5.2.3, pp. 178-181.
;
; Problem:
;   Determine circle parameters (center XM, YM and radius R) from 12 polar
;   coordinate measurements (distance + direction). Observations are error-
;   prone and enter the condition equations via #-references -> GLM.
;
;   Two variants:
;     WGLM: Diagonal weight matrix (no covariances)
;     GGLM: Additional covariances between adjacent distances (cyclic)
;
;   Condition equation (for each measurement i):
;     (Dᵢ·cos(Tᵢ·pi/200) - XM)² + (Dᵢ·sin(Tᵢ·pi/200) - YM)² - R² = 0
;
; Classification:
;   Model:       Gauss-Helmert (implicit functional model)
;   Weighting:   WGLM (diagonal) / GGLM (with covariances)
;   Linearity:   Nonlinear (quadratic condition equations)
;   Algorithm:   Gauss-Newton
;   VCE:         No
;   Redundancy:  f = c - u = 12 - 3 = 9
;
; ===============================================================================


; -- Initialization -----------------------------------------------------------
#include "../Adjustment.au3"

; -- Observation data ---------------------------------------------------------
; Distances [m]
Local $aD[] = [ _
	156.5221, 161.9524, 159.2147, 148.6576, 131.7163, 111.2273, _
	 92.1938,  82.2948,  87.4644, 104.4039, 125.1007, 143.6252]

; Directions [gon]
Local $aT[] = [ _
	70.48282, 62.48844, 54.34595, 46.97062, 41.55401, 39.98572, _
	45.11189, 58.42577, 73.50564, 81.44536, 81.75780, 77.36336]

; Standard deviations
Local $fSigmaD = 0.002    ; [m] for distances
Local $fSigmaT = 0.0003   ; [gon] for directions

; Covariance between adjacent distances (GGLM only)
Local $fCovDD = 0.000001   ; [m²]  (correlation rho ~ 0.25)


; ===============================================================================
;   WGLM: Diagonal weight matrix (no covariances)
; ===============================================================================

ConsoleWrite("===============================================================" & @CRLF)
ConsoleWrite("  WGLM: Diagonal covariance matrix" & @CRLF)
ConsoleWrite("===============================================================" & @CRLF & @CRLF)

Local $mWGLM
Local $i

; Add distance observations  D₁ .. D₁₂
For $i = 0 To 11
	_adj_addObs($mWGLM, "D" & ($i + 1), $aD[$i], $fSigmaD)
Next

; Add direction observations  T₁ .. T₁₂
For $i = 0 To 11
	_adj_addObs($mWGLM, "T" & ($i + 1), $aT[$i], $fSigmaT)
Next

; Condition equations:
;   (Dᵢ·cos(Tᵢ·pi/200) - XM)² + (Dᵢ·sin(Tᵢ·pi/200) - YM)² - R² = 0
;
; Note: gon -> rad conversion:  rad = gon * pi / 200
; In AutoIt/UDF formulas, use  ACOS(-1)/200  for pi/200.
Local Const $sPiOver200 = "acos(-1)/200"
Local $sD, $sT, $sFormula
For $i = 1 To 12
	$sD = "#D" & $i
	$sT = "#T" & $i
	$sFormula = "(" & $sD & "*cos(" & $sT & "*" & $sPiOver200 & ") - XM)^2 + " & _
	            "(" & $sD & "*sin(" & $sT & "*" & $sPiOver200 & ") - YM)^2 - R^2"
	_adj_addFunction($mWGLM, $sFormula, 0)
Next

; Initial values for parameters
_adj_setInitialValue($mWGLM, "XM", 70.050)
_adj_setInitialValue($mWGLM, "YM", 99.970)
_adj_setInitialValue($mWGLM, "R",  40.030)

; Solve
_adj_solve($mWGLM, _adj_defaultConfig("GN", False))
If @error Then
	ConsoleWrite("!> ERROR WGLM: @error=" & @error & " @extended=" & @extended & @CRLF)
	Exit 1
EndIf

Local $mDisplay = _adj_defaultDisplayConfig()
$mDisplay.obsCols = "name|value|v|r|w|p|decision|mdb"
$mDisplay.showGlobalTest = True
ConsoleWrite(_adj_displayResults($mWGLM, $mDisplay))

; Validation (Niemeier 2008, p. 181)
Local $mResW = _adj_getResults($mWGLM)
Local $mXW   = $mResW.x1

ConsoleWrite(@CRLF & "Validation vs. Reference (Niemeier 2008, WGLM):" & @CRLF)
_check("XM",  $mXW["XM"], 70.000, 3)
_check("YM",  $mXW["YM"], 100.000, 3)
_check("R",   $mXW["R"],  40.000, 3)
_check("s₀",  $mResW.s0,  1.13, 2)
_check("f",   $mResW.f,   9, 0)


; ===============================================================================
;   GGLM: With covariances between adjacent distances
; ===============================================================================

ConsoleWrite(@CRLF & @CRLF)
ConsoleWrite("===============================================================" & @CRLF)
ConsoleWrite("  GGLM: With off-diagonal covariances" & @CRLF)
ConsoleWrite("===============================================================" & @CRLF & @CRLF)

Local $mGGLM

; Add same observations
For $i = 0 To 11
	_adj_addObs($mGGLM, "D" & ($i + 1), $aD[$i], $fSigmaD)
Next
For $i = 0 To 11
	_adj_addObs($mGGLM, "T" & ($i + 1), $aT[$i], $fSigmaT)
Next

; Add covariances between adjacent distances (cyclic)
; cov(Dᵢ, Dᵢ₊₁) = 0.000001 m²  for i = 1..11
; cov(D₁₂, D₁) = 0.000001 m²   (cyclic closure)
For $i = 1 To 11
	_adj_addCovariance($mGGLM, "D" & $i, "D" & ($i + 1), $fCovDD)
Next
_adj_addCovariance($mGGLM, "D12", "D1", $fCovDD)

; Same condition equations
For $i = 1 To 12
	$sD = "#D" & $i
	$sT = "#T" & $i
	$sFormula = "(" & $sD & "*cos(" & $sT & "*" & $sPiOver200 & ") - XM)^2 + " & _
	            "(" & $sD & "*sin(" & $sT & "*" & $sPiOver200 & ") - YM)^2 - R^2"
	_adj_addFunction($mGGLM, $sFormula, 0)
Next

; Same initial values
_adj_setInitialValue($mGGLM, "XM", 70.050)
_adj_setInitialValue($mGGLM, "YM", 99.970)
_adj_setInitialValue($mGGLM, "R",  40.030)

; Solve
_adj_solve($mGGLM, _adj_defaultConfig("GN", False))
If @error Then
	ConsoleWrite("!> ERROR GGLM: @error=" & @error & " @extended=" & @extended & @CRLF)
	Exit 1
EndIf

Local $mDisplay = _adj_defaultDisplayConfig()
$mDisplay.obsCols = "name|value|v|r|w|p|decision|mdb"
$mDisplay.showGlobalTest = True
ConsoleWrite(_adj_displayResults($mGGLM, $mDisplay))

; Validation — Niemeier provides the same parameter values for both variants
; (simulated data -> round values). The Q_xx matrix differs.
Local $mResG = _adj_getResults($mGGLM)
Local $mXG   = $mResG.x1

ConsoleWrite(@CRLF & "Validation vs. Reference (Niemeier 2008, GGLM):" & @CRLF)
_check("XM",  $mXG["XM"], 70.000, 3)
_check("YM",  $mXG["YM"], 100.000, 3)
_check("R",   $mXG["R"],  40.000, 3)
_check("f",   $mResG.f,   9, 0)

; Compare Q_xx between WGLM and GGLM (informational)
ConsoleWrite(@CRLF & "Note: Both variants yield the same parameters (simulated data)." & @CRLF)
ConsoleWrite("The difference lies in the cofactor matrix Q_xx, reflecting" & @CRLF)
ConsoleWrite("the impact of observation covariances on parameter precision." & @CRLF)


; -- Helper -------------------------------------------------------------------
Func _check($sName, $fActual, $fExpected, $iDez)
	Local $sFmt  = "%." & $iDez & "f"
	Local $fTol  = 0.5 * (10 ^ (-$iDez))
	Local $fDiff = Abs($fActual - $fExpected)
	Local $sOK   = ($fDiff <= $fTol) ? "✓" : "✗ MISMATCH"
	ConsoleWrite(StringFormat("  %-12s  expect: " & $sFmt & "  got: " & $sFmt & "  [%s]", _
		$sName, $fExpected, $fActual, $sOK) & @CRLF)
EndFunc
