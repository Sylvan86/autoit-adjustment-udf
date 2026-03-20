#AutoIt3Wrapper_Change2CUI=y

; ===============================================================================
;   Leveling Network -- Ghilani
; ===============================================================================
;
; Source:
;   Krumm, F. (2020): Geodetic Network Adjustment Examples, Rev. 3.5,
;   Geodätisches Institut, Universität Stuttgart.
;   Example 2.1 (pp. 6–9), based on Ghilani (2010), Ex. 12.6.
;
; Problem:
;   Leveling network with 4 height points.
;   Point 1 serves as datum and is held fixed (H₁ = 437.596 m).
;   6 leveled height differences → 3 unknowns (H₂, H₃, H₄).
;   (Points labeled A–D in the source, renamed to 1–4 for subscript notation.)
;
;            1(A) ───── 2(B)
;              │╲       ╱│
;              │  ╲   ╱  │
;              │    ╳    │
;              │  ╱   ╲  │
;              │╱       ╲│
;            4(D) ───── 3(C)
;
; Classification:
;   Model:       Observation equations (Gauss–Markov)
;   Weighting:   Weighted (WLS) — individual σ per observation
;   Linearity:   Linear  (Δhᵢⱼ = Hⱼ − Hᵢ)
;   Algorithm:   Gauss–Newton (single iteration for linear model)
;   VCE:         No
;   Datum:       Point 1 fixed (H₁ = 437.596 m)
;   Redundancy:  f = n − u = 6 − 3 = 3
;
; ===============================================================================


; -- Initialization -----------------------------------------------------------------
#include "../Adjustment.au3"


; -- Model setup --------------------------------------------------------------------
Local $mSystem

; Observations: leveled height differences Δh [m] with standard deviation σ [m]
; Functional model:  Δhᵢⱼ = Hⱼ − Hᵢ
;
;   Obs.     Formula      Δh [m]       σ [m]
;   Δh₁₂    H₂ − H₁      10.5090     0.006
;   Δh₁₃    H₃ − H₁      15.8810     0.012
;   Δh₂₃    H₃ − H₂       5.3600     0.004
;   Δh₂₄    H₄ − H₂      −3.1670     0.004
;   Δh₃₄    H₄ − H₃      −8.5230     0.005
;   Δh₄₁    H₁ − H₄      −7.3480     0.003

_adj_addObsFunction($mSystem, "Δh₁₂", "H2 - H1",  10.5090, 0.006)
_adj_addObsFunction($mSystem, "Δh₁₃", "H3 - H1",  15.8810, 0.012)
_adj_addObsFunction($mSystem, "Δh₂₃", "H3 - H2",   5.3600, 0.004)
_adj_addObsFunction($mSystem, "Δh₂₄", "H4 - H2",  -3.1670, 0.004)
_adj_addObsFunction($mSystem, "Δh₃₄", "H4 - H3",  -8.5230, 0.005)
_adj_addObsFunction($mSystem, "Δh₄₁", "H1 - H4",  -7.3480, 0.003)

; Datum definition: point 1 is held fixed
_adj_addFixedParam($mSystem, "H1", 437.596)

; Approximate values (optional for linear models — included for documentation)
_adj_setInitialValue($mSystem, "H2", 448.105) ; 448,10871
_adj_setInitialValue($mSystem, "H3", 453.465) ; 453,46847
_adj_setInitialValue($mSystem, "H4", 444.942) ; 444,94361


; -- Adjustment ---------------------------------------------------------------------
_adj_solve($mSystem, _adj_defaultConfig("GN", False))
If @error Then
	ConsoleWrite("!> ERROR in _adj_solve: @error=" & @error & " @extended=" & @extended & @CRLF)
	Exit 1
EndIf


; -- Results ------------------------------------------------------------------------
ConsoleWrite(_adj_displayResults($mSystem))


; -- Validation vs. reference values (Krumm 2020, pp. 8-9) -------------------------
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
_check("s₀",    $mRes.s0,   0.6512, 4)
_check("vᵀPv",  $mRes.vtpv, 1.27212, 5)
_check("f",     $mRes.f,    3, 0)

; --- Adjusted heights [m] ---
ConsoleWrite(@CRLF & "Adjusted heights [m]:" & @CRLF)
_check("H₂", $mX1["H2"], 448.10871, 4)
_check("H₃", $mX1["H3"], 453.46847, 4)
_check("H₄", $mX1["H4"], 444.94361, 4)

; --- Standard deviations [mm] ---
ConsoleWrite(@CRLF & "Standard deviations σ̂ [mm]:" & @CRLF)
_check("σ̂(H₂)", $mSdx["H2"] * 1000, 2.30, 2)
_check("σ̂(H₃)", $mSdx["H3"] * 1000, 2.64, 2)
_check("σ̂(H₄)", $mSdx["H4"] * 1000, 1.76, 2)

; --- Residuals v [mm]  (v = l^ - l; signs reversed vs. Krumm) ---
ConsoleWrite(@CRLF & "Residuals v [mm]  (v = l̂ − l):" & @CRLF)
_check("v(Δh₁₂)", $mV["ΔH₁₂"] * 1000,   3.71, 2)
_check("v(Δh₁₃)", $mV["ΔH₁₃"] * 1000,  -8.53, 2)
_check("v(Δh₂₃)", $mV["ΔH₂₃"] * 1000,  -0.24, 2)
_check("v(Δh₂₄)", $mV["ΔH₂₄"] * 1000,   1.89, 2)
_check("v(Δh₃₄)", $mV["ΔH₃₄"] * 1000,  -1.86, 2)
_check("v(Δh₄₁)", $mV["ΔH₄₁"] * 1000,   0.39, 2)


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
