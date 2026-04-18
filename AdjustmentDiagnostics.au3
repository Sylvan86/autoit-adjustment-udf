#include-once

; #INDEX# =======================================================================================================================
; Title .........: AdjustmentDiagnostics
; AutoIt Version : 3.3.16.1
; Description ...: Compute-on-demand dispatcher and observation diagnostics for geodetic adjustment:
;                  Baarda test, Pope test, p-values, blunder estimation, MDB, and normal distribution.
; Author(s) .....: AspirinJunkie
; Dll ...........: ucrtbase.dll
; ===============================================================================================================================


#Region Compute-on-demand dispatcher

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_ensureComputed
; Description ...: Ensures a specific computation has been performed; resolves dependencies and computes
;                  missing steps. Results are cached in $mSystem.results to avoid recomputation.
; Syntax.........: __adj_ensureComputed(ByRef $mSystem, $sWhat)
; Parameters ....: $mSystem     - [ByRef] The adjustment system
;                  $sWhat       - Computation to ensure: "qxx", "cofactors", "redundancy",
;                                 "globalTest", or "diagnostics"
; Return values .: Success      - True
;                  Failure      - False, @error set:
;                               |$ADJ_ERR_SOLVER - Solver failure during Qxx computation
;                               |$ADJ_ERR_INPUT  - Unknown $sWhat value
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: Dependency chain: diagnostics → redundancy → cofactors → qxx
;                  Each step is only computed once; subsequent calls return cached results.
; Related .......: __adj_computeQxx, __adj_computeCofactors, __adj_computeRedundancy,
;                  __adj_computeGlobalTest, __adj_computeDiagnostics
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __adj_ensureComputed(ByRef $mSystem, $sWhat)
	Local $mResults = $mSystem.results
	Local $mState = $mSystem.state

	Switch $sWhat
		Case "qxx"
			If MapExists($mResults, "Qxx") Then Return True
			If StringRegExp($mState.modelType, 'CLS$') Then Return True ; CLS has no Qxx
			Local $mS = __adj_computeQxx($mSystem, $mState)
			If @error Then Return SetError(@error, @extended, False)
			$mState.S = $mS
			$mSystem.state = $mState
			; sdx (only when s0 already computed — __adj_computeStatistics handles this separately)
			$mResults = $mSystem.results
			If MapExists($mResults, "Qxx") And MapExists($mResults, "s0") And $mResults.s0 > 0 Then
				Local $mDiagQxx = _la_getDiag($mResults.Qxx)
				Local $mSdx = _la_sqrtElements($mDiagQxx)
				$mResults.sdx = _la_scale($mSdx, $mResults.s0)
				$mSystem.results = $mResults
			EndIf
			Return True

		Case "cofactors"
			If MapExists($mResults, "Qvv") Then Return True
			; dependency: qxx (for S-matrix)
			__adj_ensureComputed($mSystem, "qxx")
			If @error Then Return SetError(@error, @extended, False)
			$mState = $mSystem.state ; re-read after qxx
			Local $mSCof = MapExists($mState, "S") ? $mState.S : Null
			__adj_computeCofactors($mSystem, $mState, $mSCof)
			If @error Then Return SetError(@error, @extended, False)
			Return True

		Case "redundancy"
			If MapExists($mResults, "redundancyDiag") Then Return True
			; dependency: cofactors
			__adj_ensureComputed($mSystem, "cofactors")
			If @error Then Return SetError(@error, @extended, False)
			$mState = $mSystem.state ; re-read
			__adj_computeRedundancy($mSystem, $mState)
			If @error Then Return SetError(@error, @extended, False)
			Return True

		Case "correlation"
			If MapExists($mResults, "correlation") Then Return True
			; dependency: qxx
			__adj_ensureComputed($mSystem, "qxx")
			If @error Then Return SetError(@error, @extended, False)
			__adj_computeCorrelation($mSystem)
			If @error Then Return SetError(@error, @extended, False)
			Return True

		Case "globalTest"
			If MapExists($mResults, "globalTestPassed") Then Return True
			__adj_computeGlobalTest($mSystem)
			Return True

		Case "diagnostics"
			If MapExists($mResults, "baardaW") Then Return True
			; dependency: redundancy
			__adj_ensureComputed($mSystem, "redundancy")
			If @error Then Return SetError(@error, @extended, False)
			__adj_computeDiagnostics($mSystem)
			If @error Then Return SetError(@error, @extended, False)
			Return True

		Case Else
			Return SetError($ADJ_ERR_INPUT, 0, False)
	EndSwitch
EndFunc

#EndRegion ; Compute-on-demand dispatcher

#Region Observation diagnostics

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_computeDiagnostics
; Description ...: Computes per-observation diagnostics: Baarda |w|, Pope Tτ, p-values, blunder
;                  estimate, MDB (Marginal Detectable Blunder), and test decision.
; Syntax.........: __adj_computeDiagnostics(ByRef $mSystem)
; Parameters ....: $mSystem     - [ByRef] The adjustment system (redundancyDiag must exist in results)
; Return values .: Success      - Implicit (results stored in $mSystem.results)
;                  Failure      - SetError($ADJ_ERR_INPUT) if redundancyDiag not available
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: Precondition: call __adj_ensureComputed("redundancy") first.
;                  Baarda |w| = |v_i| / (σ_i · √r_i) uses a priori σ₀ = 1.
;                  Pope Tτ = |w| / ŝ₀ uses a posteriori standard deviation.
;                  testDecision based on config: diagnostics.testBasis ("baarda" or "pope"),
;                  diagnostics.alpha, diagnostics.alphaSuspect.
; Related .......: __adj_ensureComputed, __adj_normCdf, __adj_normQuantile
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __adj_computeDiagnostics(ByRef $mSystem)
	Local $mResults = $mSystem.results
	Local $mState   = $mSystem.state
	Local $mModel   = $mSystem.model
	Local $mConfig  = $mSystem.config

	If Not MapExists($mResults, "redundancyDiag") Then Return SetError($ADJ_ERR_INPUT, 0, False)

	Local $fS0       = $mResults.s0
	Local $iF        = $mResults.f
	Local $mVecR     = $mState.r                  ; BLAS solver residual vector (lives in state, synced to results only after _adj_solve returns)
	Local $mVecRed   = $mResults.redundancyDiag   ; redundancy diagonal (BLAS vector, indexed by solver position)
	Local $mIdxObs   = $mState.idxObs
	Local $mObs      = $mModel.obs

	; diagnostics config
	Local $mDiag      = $mConfig.diagnostics
	Local $fAlpha     = $mDiag.alpha
	Local $fBeta      = $mDiag.beta
	Local $fAlphaSusp = IsKeyword($mDiag.alphaSuspect) ? (10 * $fAlpha) : $mDiag.alphaSuspect
	Local $sTestBasis = MapExists($mDiag, "testBasis") ? $mDiag.testBasis : "baarda"

	; non-centrality parameter λ₀ = (z_{α/2} + z_{β})²
	Local $fZAlpha = __adj_normQuantile(1 - $fAlpha / 2)
	Local $fZBeta  = __adj_normQuantile(1 - $fBeta)
	Local $fLambda0 = ($fZAlpha + $fZBeta) ^ 2

	; result maps
	Local $mBaardaW[], $mPopeT[], $mPValue[], $mPopePValue[], $mBlunder[], $mMDB[], $mTestDecision[]

	; warning flag: Baarda unreliable when s₀ ≫ 1 or s₀ ≪ 1
	Local $bBaardaWarning = ($sTestBasis = "baarda" And ($fS0 > 3.0 Or ($fS0 > 0 And $fS0 < 0.3)))

	; VCE group scaling factors for generalized models (covariance path):
	; Helmert VCE updates .weight cumulatively → 1/√weight already correct.
	; Generalized VCE updates Σₗₗ but NOT .weight → need explicit group scaling.
	Local $mVCEGroups = ($mState.hasCovariances And MapExists($mResults, "vceGroups")) ? $mResults.vceGroups : Null

	For $sKey In MapKeys($mIdxObs)
		Local $iIdx   = $mIdxObs[$sKey]
		Local $fV     = __adj_vecGet($mVecR, $iIdx)
		Local $fRi    = __adj_vecGet($mVecRed, $iIdx)

		; effective σ for test statistics:
		; - IRLS stores pre-robust weight in ._weightOriginal (includes VCE corrections but NOT IRLS downweighting)
		; - Helmert VCE updates .weight cumulatively → 1/√weight reflects VCE corrections
		; - Generalized VCE updates Σₗₗ but NOT .weight → need explicit group scaling
		; Priority: _weightOriginal > weight > stdDev
		Local $mObsEntry = $mObs[$sKey]
		Local $fSigma
		If MapExists($mObsEntry, "_weightOriginal") Then
			$fSigma = 1 / Sqrt($mObsEntry._weightOriginal)
		ElseIf MapExists($mObsEntry, "weight") Then
			$fSigma = 1 / Sqrt($mObsEntry.weight)
		Else
			$fSigma = $mObsEntry.stdDev
		EndIf
		; generalized VCE: weight not updated by VCE → scale by group's estimated variance factor
		If IsMap($mVCEGroups) Then
			Local $sVCGroup = $mObsEntry.varComp
			If MapExists($mVCEGroups, $sVCGroup) Then $fSigma *= Sqrt(($mVCEGroups[$sVCGroup]).sigma2)
		EndIf

		; guard against zero/negative redundancy
		If $fRi < 1e-12 Or $fSigma < 1e-15 Then
			$mBaardaW[$sKey]      = Default ; not computable
			$mPopeT[$sKey]        = Default
			$mPValue[$sKey]       = Default
			$mPopePValue[$sKey]   = Default
			$mBlunder[$sKey]      = Default
			$mMDB[$sKey]          = Default
			$mTestDecision[$sKey] = "---"
			ContinueLoop
		EndIf

		; Baarda |w| = |v_i| / (σ_i · √r_i)  (a priori σ₀ = 1)
		Local $fW = Abs($fV) / ($fSigma * Sqrt($fRi))
		$mBaardaW[$sKey] = $fW

		; Pope Tτ = |w| / ŝ₀
		Local $fTau = ($fS0 > 1e-15) ? ($fW / $fS0) : Default
		$mPopeT[$sKey] = $fTau

		; p-value from Baarda (normal distribution)
		; use 2·Φ(-|w|) instead of 2·(1-Φ(|w|)) to avoid catastrophic cancellation for large |w|
		Local $fPBaarda = 2 * __adj_normCdf(-$fW)
		$mPValue[$sKey] = $fPBaarda

		; p-value from Pope (τ-distribution, NOT plain Student-t — see __adj_popeCdf header)
		; cancellation-safe form: 2·F_τ(-|τ|), evaluated via the small-tail incomplete beta.
		Local $fPPope = (Not IsKeyword($fTau)) ? 2 * __adj_popeCdf(-$fTau, $iF) : Default
		$mPopePValue[$sKey] = $fPPope

		; blunder estimate ∇̂ = v_i / r_i
		$mBlunder[$sKey] = $fV / $fRi

		; MDB = σ_i · √(λ₀ / r_i)
		$mMDB[$sKey] = $fSigma * Sqrt($fLambda0 / $fRi)

		; testDecision based on selected test basis
		Local $fPDecision = ($sTestBasis = "pope" And Not IsKeyword($fPPope)) ? $fPPope : $fPBaarda
		If $fPDecision < $fAlpha Then
			$mTestDecision[$sKey] = "outlier"
		ElseIf $fPDecision < $fAlphaSusp Then
			$mTestDecision[$sKey] = "suspect"
		Else
			$mTestDecision[$sKey] = "ok"
		EndIf
	Next

	$mResults.baardaW        = $mBaardaW
	$mResults.popeT          = $mPopeT
	$mResults.pValue         = $mPValue
	$mResults.popePValue     = $mPopePValue
	$mResults.blunder        = $mBlunder
	$mResults.mdb            = $mMDB
	$mResults.testDecision   = $mTestDecision
	$mResults.testBasis      = $sTestBasis
	$mResults.baardaWarning  = $bBaardaWarning
	$mSystem.results = $mResults
EndFunc

#EndRegion ; Observation diagnostics

#Region Normal distribution functions

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_normCdf
; Description ...: Standard normal CDF Φ(x) = P(Z ≤ x)
; Syntax.........: __adj_normCdf($x)
; Parameters ....: $x           - Input value
; Return values .: Φ(x) in [0, 1]
; Author ........: AspirinJunkie
; ===============================================================================================================================
Func __adj_normCdf($x)
	Return 0.5 * __adj_erfc(-0.707106781186547524 * $x) ; 1/√2
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_normQuantile
; Description ...: Inverse standard normal CDF: returns z such that Φ(z) = p
; Syntax.........: __adj_normQuantile($p)
; Parameters ....: $p           - Probability in (0, 1)
; Return values .: Success      - z-quantile (Float)
;                  Failure      - SetError(1) if p out of range
; Author ........: AspirinJunkie
; ===============================================================================================================================
Func __adj_normQuantile($p)
	If $p <= 0 Or $p >= 1 Then Return SetError(1, 0, 0)
	Return -1.41421356237309505 * __adj_inverfc(2 * $p) ; √2
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_erfc
; Description ...: Complementary error function erfc(x) via C-Runtime
; Syntax.........: __adj_erfc($x)
; Parameters ....: $x           - Input value
; Return values .: erfc(x) as Double
; Author ........: AspirinJunkie
; ===============================================================================================================================
Func __adj_erfc($x)
	Return DllCall("ucrtbase.dll", "double:cdecl", "erfc", "double", $x)[0]
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_inverfc
; Description ...: Inverse complementary error function using polynomial approximation + Newton-Raphson
; Syntax.........: __adj_inverfc($p)
; Parameters ....: $p           - Input in (0, 2)
; Return values .: x such that erfc(x) = p
; Author ........: AspirinJunkie
; ===============================================================================================================================
Func __adj_inverfc($p)
	If $p >= 2 Then Return -100
	If $p <= 0 Then Return 100

	Local $pp = ($p < 1) ? $p : 2 - $p
	Local $t = Sqrt(-2 * Log($pp / 2))
	Local $x = -0.70711 * ((2.30753 + $t * 0.27061) / (1 + $t * (0.99229 + $t * 0.04481)) - $t)

	; 2x Newton-Raphson refinement
	Local $fErr
	For $j = 0 To 1
		$fErr = __adj_erfc($x) - $pp
		$x += $fErr / (1.12837916709551257 * Exp(-($x * $x)) - $x * $fErr) ; 2/√π
	Next

	Return ($p < 1) ? $x : -$x
EndFunc

#EndRegion ; Normal distribution functions

#Region Student's t-distribution functions

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_tCdf
; Description ...: Student's t-distribution CDF: P(T ≤ t) for T ~ t(f)
; Syntax.........: __adj_tCdf($t, $f)
; Parameters ....: $t           - Test statistic value
;                  $f           - Degrees of freedom (> 0)
; Return values .: CDF value in [0, 1]
; Author ........: AspirinJunkie
; Remarks .......: Uses regularized incomplete beta function. Ported from Stat_Distributions.au3.
; Related .......: __adj_betaRI, __adj_computeDiagnostics
; ===============================================================================================================================
Func __adj_tCdf($t, $f)
	If $f <= 0 Then Return 0.5
	Local $p = 0.5 * __adj_betaRI($f / ($f + $t * $t), 0.5 * $f, 0.5)
	Return ($t >= 0) ? (1 - $p) : $p
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_popeCdf
; Description ...: Pope-τ-distribution CDF: P(T ≤ t) for T ~ tau(f)
; Syntax.........: __adj_popeCdf($t, $f)
; Parameters ....: $t           - Test statistic value
;                  $f           - Degrees of freedom (> 1)
; Return values .: CDF value in [0, 1]
; Author ........: AspirinJunkie
; Remarks .......: Pope (1976): τ²·(f-1)/(f-τ²) ~ F(1, f-1) — i.e. the τ-distribution
;                  is supported on (-√f, √f); outside the interval the CDF is 0 or 1.
;                  Used for outlier diagnostics with a-posteriori ŝ₀ (numerator and
;                  denominator of τ_i = w_i/ŝ₀ are NOT independent, so plain Student-t
;                  with f DoF is wrong — Pope showed the correct distribution has f-1 DoF
;                  inside this transform).
;                  Implementation uses the symmetry I_x(a,b) = 1 − I_{1-x}(b,a) of the
;                  incomplete beta to evaluate the small-tail directly, avoiding
;                  catastrophic cancellation for large |t|.
; Related .......: __adj_betaRI, __adj_tCdf, __adj_computeDiagnostics
; ===============================================================================================================================
Func __adj_popeCdf($t, $f)
	If $f <= 1 Then Return ($t >= 0) ? 1 : 0   ; degenerate (no a-posteriori scale info)
	Local $fT2 = $t * $t
	If $fT2 >= $f Then Return ($t > 0) ? 1 : 0  ; outside support → tail is exactly 0
	; one-sided tail: P(|τ| ≥ |t|) = I_{(f-t²)/f}((f-1)/2, 1/2)
	; (form chosen to keep the small quantity small — no 1−x with cancellation)
	Local $fHalfTail = 0.5 * __adj_betaRI(($f - $fT2) / $f, 0.5 * ($f - 1), 0.5)
	Return ($t >= 0) ? (1 - $fHalfTail) : $fHalfTail
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_betaRI
; Description ...: Regularized incomplete beta function I_x(a,b) = B(x;a,b) / B(a,b)
; Syntax.........: __adj_betaRI($x, $a, $b)
; Parameters ....: $x           - Upper limit of integration in [0, 1]
;                  $a           - Shape parameter a > 0
;                  $b           - Shape parameter b > 0
; Return values .: I_x(a,b) in [0, 1]
; Author ........: AspirinJunkie
; Remarks .......: Uses continued fraction representation. Ported from Stat_Distributions.au3 / Stat_Basics.au3.
; Related .......: __adj_betacf, __adj_tCdf
; ===============================================================================================================================
Func __adj_betaRI($x, $a, $b)
	If $a <= 0 Or $b <= 0 Then Return SetError(1, 0, -1)
	If $x <= 0 Then Return 0
	If $x >= 1 Then Return 1

	Local $bt = Exp(__adj_gammaLn($a + $b) - __adj_gammaLn($a) - __adj_gammaLn($b) + $a * Log($x) + $b * Log(1 - $x))

	If $x < ($a + 1) / ($a + $b + 2) Then Return $bt * __adj_betacf($x, $a, $b) / $a
	Return 1 - $bt * __adj_betacf(1 - $x, $b, $a) / $b
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_betacf
; Description ...: Continued fraction representation for regularized incomplete beta function
; Syntax.........: __adj_betacf($x, $a, $b)
; Parameters ....: $x           - Upper limit in (0, 1)
;                  $a           - Shape parameter a
;                  $b           - Shape parameter b
; Return values .: Continued fraction value
; Author ........: AspirinJunkie
; Remarks .......: Lentz's method. Ported from Stat_Basics.au3.
; Related .......: __adj_betaRI
; ===============================================================================================================================
Func __adj_betacf($x, $a, $b)
	Local Const $EPS = 2.2204460492503131e-16  ; machine epsilon
	Local Const $FPMIN = 1e-30

	Local $m, $m2, $aa, $del
	Local $qab = $a + $b, $qap = $a + 1.0, $qam = $a - 1.0, $c = 1, $d = 1 - $qab * $x / $qap
	If Abs($d) < $FPMIN Then $d = $FPMIN
	$d = 1 / $d
	Local $h = $d

	For $m = 1 To 9999
		$m2 = 2 * $m
		$aa = $m * ($b - $m) * $x / (($qam + $m2) * ($a + $m2))
		$d = 1 + $aa * $d
		If Abs($d) < $FPMIN Then $d = $FPMIN
		$c = 1 + $aa / $c
		If Abs($c) < $FPMIN Then $c = $FPMIN
		$d = 1 / $d
		$h *= $d * $c
		$aa = -($a + $m) * ($qab + $m) * $x / (($a + $m2) * ($qap + $m2))
		$d = 1 + $aa * $d
		If Abs($d) < $FPMIN Then $d = $FPMIN
		$c = 1 + $aa / $c
		If Abs($c) < $FPMIN Then $c = $FPMIN
		$d = 1 / $d
		$del = $d * $c
		$h *= $del
		If Abs($del - 1) <= $EPS Then ExitLoop
	Next
	Return $h
EndFunc

#EndRegion ; Student's t-distribution functions
