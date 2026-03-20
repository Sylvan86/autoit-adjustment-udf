#AutoIt3Wrapper_Change2CUI=y

; =============================================================================
;   Robust Estimation (BIBER) — Wicki Leveling Network with 2 Gross Errors
; =============================================================================
;
; Source:
;   Wicki, F. (1998): Robuste Schätzverfahren für die Parameterschätzung
;   in geodätischen Netzen, Kap. 18.2, pp. 151–157.
;
; Problem:
;   Leveling network with 4 unknowns (H6, H8, H10, H11).
;   Datum: point 9 held fixed (H = 0).
;   9 height difference observations with individual weights.
;   Two gross errors injected: Obs 1 (6→8) += 100 mm, Obs 7 (9→10) -= 100 mm.
;
;   BIBER estimator (Schweppe-type M-estimator) with tuning constant c = 3.5.
;
; Expected results:
;   H6 ≈ −27.816, H8 ≈ 4.246, H10 ≈ −2.315, H11 ≈ 30.415  [mm]
;   Obs 1 and Obs 7 detected as outliers (|w_ROB| = 3.50 = c)
;   Max. coordinate deviation from outlier-free solution: ~5 mm
;
; Classification:
;   Model:       Observation equations (Gauss–Markov)
;   Weighting:   Weighted (WLS) — σᵢ = 1/√pᵢ
;   Linearity:   Linear  (Δhᵢⱼ = Hⱼ − Hᵢ)
;   Algorithm:   Gauss–Newton with IRLS (robust)
;   VCE:         No
;   Robust:      BIBER (Schweppe-type), c = 3.5, a-priori scale
;   Datum:       Point 9 fixed (H = 0, implicit)
;   Redundancy:  f = n − u = 9 − 4 = 5
;
; =============================================================================


; ── Initialization ────────────────────────────────────────────────────────────
#include "../Adjustment.au3"


; ── Model setup ───────────────────────────────────────────────────────────────
Local $mSystem

; Observations: leveled height differences Δh [mm] with σ = 1/√pᵢ [mm]
; Functional model:  Δhᵢⱼ = Hⱼ − Hᵢ  (datum: H9 = 0, so H9 drops out)
;
; Wicki weights pᵢ imply σᵢ = 1/√pᵢ [mm].
; Observations given in mm so that weights and units are consistent.
;
;   Obs.    from→to     Δh [mm]    pᵢ       Formula
;   Obs1    6→8          32159     0.1276    H8 − H6      (contaminated: +100 mm)
;   Obs2    8→10         −6556     0.1372    H10 − H8
;   Obs3    8→11         26170     0.0772    H11 − H8
;   Obs4    10→11        32726     0.1041    H11 − H10
;   Obs5    9→6         −27809     0.0693    H6
;   Obs6    9→11         30419     0.1041    H11
;   Obs7    9→10         −2417     0.0918    H10           (contaminated: −100 mm)
;   Obs8    9→8           4246     0.1372    H8
;   Obs9    6→10         25496     0.0865    H10 − H6

_adj_addObsFunction($mSystem, "Obs1",  "H8 - H6",      32159, 1.0 / Sqrt(0.1276))
_adj_addObsFunction($mSystem, "Obs2",  "H10 - H8",     -6556, 1.0 / Sqrt(0.1372))
_adj_addObsFunction($mSystem, "Obs3",  "H11 - H8",     26170, 1.0 / Sqrt(0.0772))
_adj_addObsFunction($mSystem, "Obs4",  "H11 - H10",    32726, 1.0 / Sqrt(0.1041))
_adj_addObsFunction($mSystem, "Obs5",  "H6",          -27809, 1.0 / Sqrt(0.0693))
_adj_addObsFunction($mSystem, "Obs6",  "H11",          30419, 1.0 / Sqrt(0.1041))
_adj_addObsFunction($mSystem, "Obs7",  "H10",          -2417, 1.0 / Sqrt(0.0918))
_adj_addObsFunction($mSystem, "Obs8",  "H8",            4246, 1.0 / Sqrt(0.1372))
_adj_addObsFunction($mSystem, "Obs9",  "H10 - H6",     25496, 1.0 / Sqrt(0.0865))

; Initial values (zeros — linear system, converges in 1 iteration)
_adj_setInitialValue($mSystem, "H6",  0.0)
_adj_setInitialValue($mSystem, "H8",  0.0)
_adj_setInitialValue($mSystem, "H10", 0.0)
_adj_setInitialValue($mSystem, "H11", 0.0)


; ── Adjustment ────────────────────────────────────────────────────────────────
Local $mParams = _adj_robustDefaults("BIBER")
$mParams.c = 3.5
$mParams.scale = "apriori"   ; Wicki uses σ₀=1 (a-priori model), no scale estimation
$mParams.outlierThreshold = 3.49   ; just below c so that outliers are flagged
Local $mConfig = _adj_defaultConfig("GN", False)
$mConfig.robust = "BIBER"
$mConfig.robustParams = $mParams

_adj_solve($mSystem, $mConfig)
If @error Then
	ConsoleWrite("!> ERROR in _adj_solve: @error=" & @error & " @extended=" & @extended & @CRLF)
	Exit 1
EndIf


; ── Results ───────────────────────────────────────────────────────────────────
Local $mDisplay = _adj_defaultDisplayConfig()
$mDisplay.showRobust = True
$mDisplay.showGlobalTest = True
$mDisplay.obsCols = "name|value|v|sdv|r|w|robW"
ConsoleWrite(_adj_displayResults($mSystem, $mDisplay) & @CRLF)


; ── Validation vs. reference values (Wicki 1998, pp. 151–157) ─────────────────
Local $mRes = _adj_getResults($mSystem)
Local $mX1  = $mRes.x1

ConsoleWrite(@CRLF)
ConsoleWrite("╔═══════════════════════════════════════════════════╗" & @CRLF)
ConsoleWrite("║  Validation vs. Reference (Wicki 1998)            ║" & @CRLF)
ConsoleWrite("╚═══════════════════════════════════════════════════╝" & @CRLF & @CRLF)

; ─── Adjusted heights [mm] — Wicki reference is approximate (≈) ───
ConsoleWrite("Adjusted heights [mm]:" & @CRLF)
_check("H6",  $mX1["H6"],  -27819, 0)
_check("H8",  $mX1["H8"],    4247, 0)
_check("H10", $mX1["H10"],  -2317, 0)
_check("H11", $mX1["H11"],  30415, 0)

; ─── Deviation from outlier-free L2 solution [mm] ───
; Outlier-free L2 solution: [-27811, 4246, -2313, 30416]
ConsoleWrite(@CRLF & "Deviation from outlier-free L2 solution [mm]:" & @CRLF)
ConsoleWrite(StringFormat("  dH6:  %.1f mm", $mX1["H6"]  - (-27811)) & @CRLF)
ConsoleWrite(StringFormat("  dH8:  %.1f mm", $mX1["H8"]  - 4246)     & @CRLF)
ConsoleWrite(StringFormat("  dH10: %.1f mm", $mX1["H10"] - (-2313))  & @CRLF)
ConsoleWrite(StringFormat("  dH11: %.1f mm", $mX1["H11"] - 30416)    & @CRLF)

; ─── Robust estimation summary ───
ConsoleWrite(@CRLF & "Robust estimation summary:" & @CRLF)
ConsoleWrite("  IRLS iterations: " & $mRes.robustIterations & @CRLF)
ConsoleWrite("  Converged:       " & $mRes.robustConverged & @CRLF)
ConsoleWrite("  Robust scale:    " & StringFormat("%.4f", $mRes.robustScale) & @CRLF)

; ─── Robust weights — outlier detection ───
ConsoleWrite(@CRLF & "Robust weights (outliers marked with *):" & @CRLF)
Local $mWeights = $mRes.robustWeights
For $i = 1 To 9
	Local $sName = "OBS" & $i
	Local $sFlag = ""
	If $mWeights[$sName] < 1.0 Then $sFlag = " *"
	ConsoleWrite(StringFormat("  %5s: rob.w = %.4f%s", $sName, $mWeights[$sName], $sFlag) & @CRLF)
Next
ConsoleWrite(@CRLF)
ConsoleWrite("Expected: Obs1 and Obs7 detected as outliers (rob.w < 1)" & @CRLF)


; ── Helper ────────────────────────────────────────────────────────────────────
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
