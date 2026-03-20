#AutoIt3Wrapper_Change2CUI=y

; ===============================================================================
;   Dynamic Height Network -- Krumm
; ===============================================================================
;
; Source:
;   Krumm, F. (2020): Geodetic Network Adjustment Examples, Rev. 3.5,
;   Geodaetisches Institut, Universitaet Stuttgart.
;   Example 2.3 (pp. 13–16).
;
; Problem:
;   Leveling network with 5 height points (2, 3, 6, 7, 8).
;   Points 2 and 3 have known heights with a FULL covariance matrix
;   (dynamic datum — not fixed, but stochastically constrained).
;   Points 6, 7, 8 are new points to be determined.
;   5 leveled height differences + 2 pseudo-observations (datum heights).
;   Total: n = 7 observations, u = 5 unknowns, f = 2.
;
;          2 ──────── 8
;          │           │
;          │           │
;          3 ──── 6 ── 7
;
; GLS Feature:
;   The two datum heights H₂ and H₃ share a non-diagonal 2×2
;   covariance matrix Σ_datum, creating off-diagonal entries in the
;   full observation covariance matrix Σ_ll.
;   This makes the weight matrix P = Σ_ll⁻¹ non-diagonal → GLS.
;
; Classification:
;   Model:       Observation equations (Gauss–Markov)
;   Weighting:   Generalized (GLS) — full covariance matrix Σ_ll
;   Linearity:   Linear  (Δhᵢⱼ = Hⱼ − Hᵢ,  H_pseudo = Hᵢ)
;   Algorithm:   Gauss–Newton (single iteration for linear model)
;   VCE:         No
;   Datum:       Dynamic (H₂, H₃ as pseudo-observations with covariance)
;   Redundancy:  f = n − u = 7 − 5 = 2
;
; ===============================================================================


; -- Initialization -----------------------------------------------------------------
#include "../Adjustment.au3"


; -- Model setup --------------------------------------------------------------------
Local $mSystem

; Datum pseudo-observations: known heights with full covariance matrix
; Stochastic model for datum points (dynamic datum):
;
;   Σ_datum = │ 0.0025  −0.0015 │  [m²]
;             │−0.0015   0.0036 │
;
;   σ(H₂) = √0.0025 = 0.0500 m
;   σ(H₃) = √0.0036 = 0.0600 m
;   Cov(H₂, H₃)     = −0.0015 m²

_adj_addObsFunction($mSystem, "H2obs", "H2", 107.7541, 0.05)
_adj_addObsFunction($mSystem, "H3obs", "H3", 103.4535, 0.06)

; Off-diagonal covariance between the two datum observations
_adj_addCovariance($mSystem, "H2obs", "H3obs", -0.0015)

; Observations: leveled height differences Δh [m]
; Stochastic model: σᵢ = √(sᵢ / 1000) [m]  where sᵢ is the distance in meters
; (standard leveling variance model: σ² ∝ distance)
;
;   Obs.     Formula        Δh [m]    s [m]    σ [m]
;   Δh₂₈    H8 − H2        5.1280     700     0.83666
;   Δh₃₆    H6 − H3        2.1830     500     0.70711
;   Δh₃₇    H7 − H3       12.2540     500     0.70711
;   Δh₆₇    H7 − H6       10.0710     800     0.89443
;   Δh₈₇    H7 − H8        2.8240     800     0.89443

_adj_addObsFunction($mSystem, "Dh28", "H8 - H2",  5.1280, Sqrt(700.0 / 1000.0))
_adj_addObsFunction($mSystem, "Dh36", "H6 - H3",  2.1830, Sqrt(500.0 / 1000.0))
_adj_addObsFunction($mSystem, "Dh37", "H7 - H3", 12.2540, Sqrt(500.0 / 1000.0))
_adj_addObsFunction($mSystem, "Dh67", "H7 - H6", 10.0710, Sqrt(800.0 / 1000.0))
_adj_addObsFunction($mSystem, "Dh87", "H7 - H8",  2.8240, Sqrt(800.0 / 1000.0))

; Approximate values (optional for linear models — included for documentation)
_adj_setInitialValue($mSystem, "H2", 107.7541)
_adj_setInitialValue($mSystem, "H3", 103.4535)
_adj_setInitialValue($mSystem, "H6", 105.6400)
_adj_setInitialValue($mSystem, "H7", 115.7110)
_adj_setInitialValue($mSystem, "H8", 112.8850)


; -- Adjustment ---------------------------------------------------------------------
_adj_solve($mSystem, _adj_defaultConfig("GN", False))
If @error Then
	ConsoleWrite("!> ERROR in _adj_solve: @error=" & @error & " @extended=" & @extended & @CRLF)
	Exit 1
EndIf


; -- Results ------------------------------------------------------------------------
ConsoleWrite(_adj_displayResults($mSystem))


; -- Validation vs. reference values (Krumm 2020, pp. 14-16) ------------------------
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
_check("f",     $mRes.f,    2, 0)

; --- Adjusted heights [m] ---
ConsoleWrite(@CRLF & "Adjusted heights [m]:" & @CRLF)
_check("H₂", $mX1["H2"], 107.7541, 4)
_check("H₃", $mX1["H3"], 103.4535, 4)
_check("H₆", $mX1["H6"], 105.6364, 4)
_check("H₇", $mX1["H7"], 115.7072, 4)
_check("H₈", $mX1["H8"], 112.8826, 4)

; --- Standard deviations [mm] ---
ConsoleWrite(@CRLF & "Standard deviations σ̂ [mm]:" & @CRLF)
_check("σ̂(H₂)", $mSdx["H2"] * 1000, 0.04, 2)
_check("σ̂(H₃)", $mSdx["H3"] * 1000, 0.04, 2)
_check("σ̂(H₆)", $mSdx["H6"] * 1000, 0.43, 2)
_check("σ̂(H₇)", $mSdx["H7"] * 1000, 0.39, 2)
_check("σ̂(H₈)", $mSdx["H8"] * 1000, 0.48, 2)

; --- Residuals v [mm]  (v = l^ - l; signs reversed vs. Krumm) ---
; Krumm gives ê = l − l̂ for height differences:
;   ê(Δh₂₈) = −0.52 mm  →  v_UDF = +0.52 mm
;   ê(Δh₃₆) =  0.10 mm  →  v_UDF = −0.10 mm
;   ê(Δh₃₇) =  0.27 mm  →  v_UDF = −0.27 mm
;   ê(Δh₆₇) =  0.17 mm  →  v_UDF = −0.17 mm
;   ê(Δh₈₇) = −0.60 mm  →  v_UDF = +0.60 mm
ConsoleWrite(@CRLF & "Residuals v [mm]  (v = l̂ − l):" & @CRLF)
_check("v(Δh₂₈)", $mV["DH28"] * 1000,  0.52, 2)
_check("v(Δh₃₆)", $mV["DH36"] * 1000, -0.10, 2)
_check("v(Δh₃₇)", $mV["DH37"] * 1000, -0.27, 2)
_check("v(Δh₆₇)", $mV["DH67"] * 1000, -0.17, 2)
_check("v(Δh₈₇)", $mV["DH87"] * 1000,  0.60, 2)


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
