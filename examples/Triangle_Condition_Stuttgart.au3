#AutoIt3Wrapper_Change2CUI=y

; =====================================================================
;  Triangle --- Nonlinear Condition Equation (Sine Rule)
; =====================================================================
;
; Source:
;   Sneeuw, N., Krumm, F., Roth, M. (2022): Lecture Notes Adjustment
;   Theory, Rev. 4.48, Geodaetisches Institut, Univ. Stuttgart.
;   Example 5.2.2 (pp. 79-80).
;
; Problem:
;   Triangle with 2 observed angles (alpha, beta) and 2 sides (a, b).
;   Nonlinear condition equation from the sine rule:
;
;     a * sin(beta) - b * sin(alpha) = 0
;
;   This is a B-model (pure condition equations) with a nonlinear
;   condition and mixed observation types (lengths and angles).
;
; Classification:
;   Model:       Condition equations (CLS, B-model)
;   Weighting:   Weighted --- P = diag(1, 0.1, 1, 0.2)
;   Linearity:   Nonlinear (sine rule)
;   Algorithm:   Gauss-Newton (iterative)
;   VCE:         No
;   Redundancy:  f = 1 (1 condition, 0 parameters)
;
; =====================================================================


; -- Initialization ---------------------------------------------------
#include "../Adjustment.au3"

Global Const $fDeg2Rad = 3.14159265358979323846 / 180.0


; -- Model setup ------------------------------------------------------
Local $mSystem

; Observations with weights P = diag(1, 0.1, 1, 0.2)
; sigma_i = 1/sqrt(P_i) (assuming sigma_0 = 1)
; Angles are stored in radians for compatibility with Sin().
;
;   Obs.     Value          P       sigma = 1/sqrt(P)
;   a        10 m           1       1.0
;   b         5 m           0.1     sqrt(10) = 3.16228
;   alpha    60 deg          1       1.0       [rad]
;   beta     23.7 deg        0.2     sqrt(5) = 2.23607  [rad]

_adj_addObs($mSystem, "A",     10.0,               1.0)
_adj_addObs($mSystem, "B",      5.0,               1.0 / Sqrt(0.1))
_adj_addObs($mSystem, "ALPHA", 60.0 * $fDeg2Rad,   1.0)
_adj_addObs($mSystem, "BETA",  23.7 * $fDeg2Rad,   1.0 / Sqrt(0.2))

; Nonlinear condition equation (sine rule):
;   a * sin(beta) - b * sin(alpha) = 0
_adj_addFunction($mSystem, "#A * Sin(#BETA) - #B * Sin(#ALPHA)")


; -- Adjustment -------------------------------------------------------
_adj_solve($mSystem, _adj_defaultConfig("GN", False))
If @error Then
	ConsoleWrite("!> ERROR in _adj_solve: @error=" & @error & " @extended=" & @extended & @CRLF)
	Exit 1
EndIf


; -- Results ----------------------------------------------------------
ConsoleWrite(_adj_displayResults($mSystem))


; -- Validation vs. reference values (Stuttgart, Ex. 5.2.2) -----------
Local $mRes = _adj_getResults($mSystem)
Local $mV   = $mRes.v

ConsoleWrite(@CRLF)
ConsoleWrite("===================================================" & @CRLF)
ConsoleWrite("  Validation vs. Reference (Stuttgart, Ex. 5.2.2)   " & @CRLF)
ConsoleWrite("===================================================" & @CRLF & @CRLF)

; --- Statistical summary ---
ConsoleWrite("Statistical summary:" & @CRLF)
_check("f", $mRes.f, 1, 0)

; Reference: e^T P e = 4.017e-6
; In the UDF, vtPv is computed from the internal weight matrix (P = sigma0^2/sigma_i^2).
; The reference vtPv = 4.017e-6 assumes P = diag(1, 0.1, 1, 0.2).
ConsoleWrite(StringFormat("  vtPv (UDF):  %.6e", $mRes.vtpv) & @CRLF)

; --- Residuals ---
ConsoleWrite(@CRLF & "Residuals:" & @CRLF)
ConsoleWrite(StringFormat("  v(a)     = %+.6e m", $mV["A"]) & @CRLF)
ConsoleWrite(StringFormat("  v(b)     = %+.6e m", $mV["B"]) & @CRLF)

; Convert angle residuals from radians to degrees/minutes/seconds
Local $fVAlpha_deg = $mV["ALPHA"] / $fDeg2Rad
Local $fVBeta_deg  = $mV["BETA"]  / $fDeg2Rad

ConsoleWrite(StringFormat("  v(alpha) = %+.6e rad = ", $mV["ALPHA"]) & _fmtDMS($fVAlpha_deg) & @CRLF)
ConsoleWrite(StringFormat("  v(beta)  = %+.6e rad = ", $mV["BETA"])  & _fmtDMS($fVBeta_deg) & @CRLF)

; Reference residuals:
;   e_a     = -5.63e-6
;   e_b     = +1.13e-4
;   e_alpha = +6' 26.241"
;   e_beta  = -1d 55' 42.492"
ConsoleWrite(@CRLF & "Reference residuals (Stuttgart):" & @CRLF)
ConsoleWrite("  e_a     = -5.63e-6 m" & @CRLF)
ConsoleWrite("  e_b     = +1.13e-4 m" & @CRLF)
ConsoleWrite("  e_alpha = +6' 26.241""" & @CRLF)
ConsoleWrite("  e_beta  = -1d 55' 42.492""" & @CRLF)

; Verification: sine rule must hold for adjusted observations
Local $fAadj     = $mRes.obsAdj["A"]
Local $fBadj     = $mRes.obsAdj["B"]
Local $fAlphaAdj = $mRes.obsAdj["ALPHA"]
Local $fBetaAdj  = $mRes.obsAdj["BETA"]
Local $fCheck    = $fAadj * Sin($fBetaAdj) - $fBadj * Sin($fAlphaAdj)
ConsoleWrite(@CRLF & StringFormat("Sine rule check (should be 0): %.6e", $fCheck) & @CRLF)


; -- Helpers ----------------------------------------------------------
Func _check($sName, $fActual, $fExpected, $iDez)
	Local $sFmt  = "%." & $iDez & "f"
	Local $fTol  = 0.5 * (10 ^ (-$iDez))
	Local $fDiff = Abs($fActual - $fExpected)
	Local $sOK   = ($fDiff <= $fTol) ? "OK" : "MISMATCH"
	ConsoleWrite(StringFormat("  %-12s  expect: " & $sFmt & "  got: " & $sFmt & "  [%s]", _
		$sName, $fExpected, $fActual, $sOK) & @CRLF)
EndFunc

; Formats a decimal-degree value as  [+/-]D d M' S.sss"
Func _fmtDMS($fDeg)
	Local $sSign = ($fDeg < 0) ? "-" : "+"
	Local $fAbs  = Abs($fDeg)
	Local $iD    = Int($fAbs)
	Local $fRem  = ($fAbs - $iD) * 60
	Local $iM    = Int($fRem)
	Local $fS    = ($fRem - $iM) * 60
	If $iD > 0 Then
		Return StringFormat("%s%dd %d' %.3f""", $sSign, $iD, $iM, $fS)
	Else
		Return StringFormat("%s%d' %.3f""", $sSign, $iM, $fS)
	EndIf
EndFunc
