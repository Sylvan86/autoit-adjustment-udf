#AutoIt3Wrapper_Change2CUI=y

; ===============================================================================
;   Variance Component Estimation -- Koch (Helmert Method)
; ===============================================================================
;
; Source:
;   Koch, K.R. (1999): Parameter Estimation and Hypothesis Testing
;   in Linear Models, 2nd ed., Springer.
;   Example I (Section 3.6.3, pp. 235–237).
;
; Problem:
;   Estimation of a single mean (β₁) from 5 observations.
;   Two variance component groups:
;     Group 1: observations 1–3  (σ₁² · α₁², α₁² = 1×10⁻⁴ → σ = 0.01)
;     Group 2: observations 4–5  (σ₂² · α₂², α₂² = 4×10⁻⁴ → σ = 0.02)
;   VCE estimates σ₁² and σ₂² iteratively (Helmert method).
;
;   Obs.   Value    σ         Group
;   y₁     1.25     0.01      grp1
;   y₂     1.26     0.01      grp1
;   y₃     1.24     0.01      grp1
;   y₄     1.22     0.02      grp2
;   y₅     1.27     0.02      grp2
;
; Classification:
;   Model:       Observation equations (Gauss–Markov)
;   Weighting:   Weighted, with variance component estimation
;   Linearity:   Linear  (yᵢ = β₁)
;   Algorithm:   Gauss–Newton (single iteration for linear model)
;   VCE:         Yes (Helmert method)
;   Groups:      grp1 (obs 1–3), grp2 (obs 4–5)
;   Redundancy:  f = n − u = 5 − 1 = 4
;
; ===============================================================================


; -- Initialization -----------------------------------------------------------------
#include "../Adjustment.au3"


; -- Model setup --------------------------------------------------------------------
Global $mSystem

; Observations: 5 measurements of a single mean β₁
; Functional model:  yᵢ = β₁
;
;   Obs.   Formula    Value    σ        Group
;   y₁     beta1      1.25    0.01     grp1
;   y₂     beta1      1.26    0.01     grp1
;   y₃     beta1      1.24    0.01     grp1
;   y₄     beta1      1.22    0.02     grp2
;   y₅     beta1      1.27    0.02     grp2

_adj_addObsFunction($mSystem, "y1", "beta1", 1.25, 0.01, "grp1")
_adj_addObsFunction($mSystem, "y2", "beta1", 1.26, 0.01, "grp1")
_adj_addObsFunction($mSystem, "y3", "beta1", 1.24, 0.01, "grp1")
_adj_addObsFunction($mSystem, "y4", "beta1", 1.22, 0.02, "grp2")
_adj_addObsFunction($mSystem, "y5", "beta1", 1.27, 0.02, "grp2")


; -- Adjustment with VCE --------------------------------------------------------
_adj_solve($mSystem, _adj_defaultConfig("GN", True))
If @error Then
	ConsoleWrite("!> ERROR in _adj_solve: @error=" & @error & " @extended=" & @extended & @CRLF)
	Exit 1
EndIf


; -- Results ------------------------------------------------------------------------
ConsoleWrite(_adj_displayResults($mSystem))


; -- Validation vs. reference values (Koch 1999, pp. 235-237) -----------------------
;
; Koch provides iterative VCE convergence:
;   Iteration 1: σ₁² = 0.8889, σ₂² = 1.7917
;   Iteration 2: σ₁² = 1.0914, σ₂² = 0.9318
;   Iteration 3: σ₁² = 0.9892, σ₂² = 1.0127
; After convergence, both variance factors should be near 1.0.
;
; The adjusted mean β₁ can be verified against the weighted mean.

Global $mRes = _adj_getResults($mSystem)
Global $mX1  = $mRes.x1

ConsoleWrite(@CRLF)
ConsoleWrite("===================================================" & @CRLF)
ConsoleWrite("  Validation vs. Reference (Koch 1999)             " & @CRLF)
ConsoleWrite("===================================================" & @CRLF & @CRLF)

; --- Statistical summary ---
ConsoleWrite("Statistical summary:" & @CRLF)
_check("f", $mRes.f, 4, 0)

; --- Adjusted parameter ---
ConsoleWrite(@CRLF & "Adjusted parameter:" & @CRLF)
; The weighted mean of {1.25, 1.26, 1.24, 1.22, 1.27} with VCE-updated weights
; Koch doesn't give explicit final beta1; validate VCE convergence instead
_check("β₁", $mX1["BETA1"], 1.249, 3)

; --- VCE convergence ---
ConsoleWrite(@CRLF & "VCE convergence:" & @CRLF)
_check("VCE converged", ($mRes.vceConverged ? 1 : 0), 1, 0)

; After VCE, both group variance factors should be near 1.0
; Koch final values (iteration 3): σ₁² ≈ 0.99, σ₂² ≈ 1.01
Global $mGroups = $mRes.vceGroups
Global $mGrp1 = $mGroups["grp1"]
Global $mGrp2 = $mGroups["grp2"]

ConsoleWrite(@CRLF & "Per-group variance factors (should be near 1.0):" & @CRLF)
_check("σ²(grp1)", $mGrp1.sigma2, 1.0, 0)
_check("σ²(grp2)", $mGrp2.sigma2, 1.0, 0)

; The group s0 values (= sqrt of sigma2)
ConsoleWrite(@CRLF & "Per-group s₀:" & @CRLF)
_check("s₀(grp1)", $mGrp1.s0, 1.0, 0)
_check("s₀(grp2)", $mGrp2.s0, 1.0, 0)


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
