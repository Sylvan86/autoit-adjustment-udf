#include-once

; #INDEX# =======================================================================================================================
; Title .........: AdjustmentStats
; AutoIt Version : 3.3.16.1
; Description ...: Statistical inference for geodetic adjustment: cofactor matrices (Qxx, Qvv, Qyhat),
;                  standard deviations, redundancy, chi-squared global test, and chi-squared distribution.
; Author(s) .....: AspirinJunkie
; Dll ...........: ucrtbase.dll
; ===============================================================================================================================

#Region Statistics orchestration

; ══════════════════════════════════════════════════════════════════
; Phase 5: Statistical inference — Qxx, s₀, sdx, Qvv, redundancy
; ══════════════════════════════════════════════════════════════════

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_computeStatistics
; Description ...: Orchestrates all statistical computations: DOF, Qxx, s₀, sdx, cofactors, redundancy, global test, diagnostics
; Syntax.........: __adj_computeStatistics($mSystem, $mState)
; Parameters ....: $mSystem     - [ByRef] The adjustment system
;                  $mState      - [ByRef] Solver state
; Return values .: Success      - None (modifies $mSystem.results and $mState in-place)
;                  Failure      - None
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: Forces cofactors+redundancy on for VCE. Uses compute-on-demand via __adj_ensureComputed.
; Related .......: __adj_ensureComputed, __adj_computeDOF, __adj_computeQxx
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __adj_computeStatistics(ByRef $mSystem, ByRef $mState)
	; NOTE: AutoIt Maps are value types. Each sub-function that modifies $mSystem.results
	; requires re-reading the local copy afterwards, otherwise the local copy is stale.
	; IMPORTANT: __adj_ensureComputed reads $mSystem.state, so we must write back
	; the ByRef $mState before calling it, and re-read $mState afterwards.

	Local $mConfig = $mSystem.config
	Local $mCompute = $mConfig.compute

	; VCE needs cofactors + redundancy — force them on regardless of compute flags
	If $mConfig.vce Then
		$mCompute.cofactors  = True
		$mCompute.redundancy = True
		$mConfig.compute = $mCompute
		$mSystem.config = $mConfig
	EndIf

	; Step 1: degrees of freedom — always
	Local $iF = __adj_computeDOF($mState)
	Local $mResults = $mSystem.results
	$mResults.f = $iF
	$mSystem.results = $mResults

	; Step 2: Qxx — when compute.qxx = True (default)
	If $mCompute.qxx Then
		$mSystem.state = $mState ; write-back for ensureComputed
		__adj_ensureComputed($mSystem, "qxx")
		$mState = $mSystem.state ; re-read (ensureComputed may store S-matrix)
	EndIf

	; Step 3: s₀ = √(vᵀPv / f) — always (re-read results after Qxx)
	$mResults = $mSystem.results
	$mResults.vtpv = $mState.r2sum
	If $iF > 0 Then
		$mResults.s0 = Sqrt($mState.r2sum / $iF)
	Else
		$mResults.s0 = 0
	EndIf

	; Step 4: sdx = s₀ · √(diag(Qxx)) — when Qxx was computed
	If MapExists($mResults, "Qxx") And $mResults.s0 > 0 Then
		Local $mDiagQxx = _la_getDiag($mResults.Qxx)
		Local $mSdx = _la_sqrtElements($mDiagQxx)
		$mResults.sdx = _la_scale($mSdx, $mResults.s0)
	EndIf
	$mSystem.results = $mResults

	; Step 5-7: delegate to resolver based on compute flags
	If $mCompute.cofactors Or $mCompute.redundancy Or $mCompute.diagnostics Then
		$mSystem.state = $mState ; write-back for ensureComputed
	EndIf
	If $mCompute.cofactors   Then __adj_ensureComputed($mSystem, "cofactors")
	If $mCompute.redundancy  Then __adj_ensureComputed($mSystem, "redundancy")
	If $mCompute.globalTest  Then __adj_ensureComputed($mSystem, "globalTest")
	If $mCompute.diagnostics Then __adj_ensureComputed($mSystem, "diagnostics")
EndFunc

#EndRegion ; Statistics orchestration

#Region Cofactor matrices (Qxx, Qvv, Qyhat)

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_computeDOF
; Description ...: Computes degrees of freedom f based on model type
; Syntax.........: __adj_computeDOF($mState)
; Parameters ....: $mState      - [ByRef] Solver state
; Return values .: Int — degrees of freedom (OLS: n-u, LSE: n-u+r, CLS: nFormulas, GLM: nFormulas+r-u)
; Author ........: AspirinJunkie
; ===============================================================================================================================
Func __adj_computeDOF(ByRef $mState)
	Local $sModel = $mState.modelType

	If StringRegExp($sModel, '^[OWG]LS$') Then
		Return $mState.nObs - $mState.nParams

	ElseIf StringRegExp($sModel, '^[WG]?LSE$') Then
		Return $mState.nObs - $mState.nParams + $mState.nRestrictions

	ElseIf StringRegExp($sModel, '^[WG]?CLS$') Then
		Return $mState.nFormulas

	ElseIf StringRegExp($sModel, '^[WG]?GLM$') Then
		Return $mState.nFormulas + $mState.nRestrictions - $mState.nParams

	Else
		Return 0
	EndIf
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_computeQxx
; Description ...: Computes cofactor matrix Qxx of adjusted parameters
; Syntax.........: __adj_computeQxx($mSystem, $mState)
; Parameters ....: $mSystem     - [ByRef] The adjustment system
;                  $mState      - [ByRef] Solver state
; Return values .: Success      - S matrix (n×n or n×nFree) for reuse in cofactor computation
;                  Failure      - Null if rank deficient
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: QR branch: Qxx = P·R⁻¹·R⁻ᵀ·Pᵀ. SVD branch: Qxx = V·diag(1/σᵢ²)·Vᵀ.
;                  Handles LSE/GLM+restrictions via Q₂ back-transformation.
;                  Includes Jacobi-Equilibration back-transform for QR branch.
; Related .......: __adj_computeQxxFromSVD, __adj_computeStatistics
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __adj_computeQxx(ByRef $mSystem, ByRef $mState)
	Local $mResults = $mSystem.results
	Local $sModel = $mState.modelType
	Local $iN, $mSource, $iLDA

	; determine source matrix (DGELSY output with R in upper triangle) and dimension
	If MapExists($mState, "Q2") Then
		; LSE or GLM+Restrictions: reduced system in nullspace of restrictions
		$iN = $mState.nParams - $mState.nRestrictions
		$mSource = $mState.Abar
		$iLDA = $mSource.rows
	ElseIf StringRegExp($sModel, 'GLM$') Then
		; GLM without restrictions: Cholesky-transformed system
		$iN = $mState.nParams
		$mSource = $mState.Atilde
		$iLDA = $mSource.rows
	Else ; OLS
		$iN = $mState.nParams
		$mSource = $mState.Matrix_A
		$iLDA = $mSource.rows
	EndIf

	; ══════ SVD branch: Qxx = V · diag(1/σᵢ²) · Vᵀ ══════
	If $mState.solver = $ADJ_SOLVER_SVD Then
		; determine source matrix: original (pre-solver) version
		Local $mSvdSource
		If MapExists($mState, "Q2") Then
			$mSvdSource = _la_duplicate($mState.Abar_orig)
		ElseIf StringRegExp($sModel, 'GLM$') Then
			$mSvdSource = _la_duplicate($mState.Atilde_orig)
		Else
			$mSvdSource = _la_duplicate($mState.A_orig)
		EndIf

		; Apply Jacobi-Equilibration to source matrix (same as QR branch uses equilibrated R).
		; Without this, the RCOND-based rank truncation operates on un-equilibrated singular
		; values and may falsely truncate small-but-valid singular values in poorly scaled systems.
		Local $bHasEquilibration = MapExists($mState, "EquilibrationScale")
		If $bHasEquilibration Then
			Local $mEqScale = $mState.EquilibrationScale
			Local $tEqS = $mEqScale.struct
			Local $iSvdRows = $mSvdSource.rows
			; scale column j by 1/S_eq_j (same scaling that was applied to the solver's working matrix)
			For $__j = 0 To $iN - 1
				_blas_scal($mSvdSource, 1.0 / DllStructGetData($tEqS, 1, $__j + 1), $__j * $iSvdRows, 1, $iSvdRows)
			Next
		EndIf

		Local $mSVDResult = __adj_computeQxxFromSVD($mSvdSource, $mState, $iN, $mSvdSource.rows, $mState.solverRCOND)
		Local $iSvdErr = @error
		If $iSvdErr Then Return SetError($iSvdErr, @extended, Null)

		Local $mQxx = $mSVDResult.Qxx
		Local $mS = $mSVDResult.V

		; Equilibration back-transform: scale V rows by 1/S_eq
		; so that Qxx = S_scaled · S_scaledᵀ = diag(S⁻¹) · Qxx_eq · diag(S⁻¹)
		; (same as QR branch lines below)
		If $bHasEquilibration Then
			For $__i = 0 To $iN - 1
				_blas_scal($mS, 1.0 / DllStructGetData($tEqS, 1, $__i + 1), $__i, $iN, $iN)
			Next
			; recompute Qxx with back-transformed V
			$mQxx = _blas_createMatrix($iN, $iN)
			_blas_syrk($mS, $mQxx, 1, 0, "U", "N", $iN, $iN)
		EndIf

		If MapExists($mState, "Q2") Then
			; LSE / GLM+Restrictions back-transformation: Qxx = Q₂ · S · Sᵀ · Q₂ᵀ = T · Tᵀ  where T = Q₂ · S
			Local $iNpar = $mState.nParams
			Local $mQ2 = $mState.Q2
			Local $mT = _blas_createMatrix($iNpar, $iN)
			_blas_gemm($mQ2, $mS, $mT, 1, 0, "N", "N", $iNpar, $iN, $iN)

			Local $mQxxFull = _blas_createMatrix($iNpar, $iNpar)
			_blas_syrk($mT, $mQxxFull, 1, 0, "U", "N", $iNpar, $iN)
			__adj_fillLowerFromUpper($mQxxFull, $iNpar)

			$mResults.Qxx = $mQxxFull
			$mSystem.results = $mResults
			Return $mT
		Else
			__adj_fillLowerFromUpper($mQxx, $iN)
			$mResults.Qxx = $mQxx
			$mSystem.results = $mResults
			Return $mS
		EndIf
	EndIf

	; ══════ QR branch: Qxx = P · R⁻¹ · R⁻ᵀ · Pᵀ ══════
	; 1. Extract R (upper n×n triangle)
	Local $mR = _blas_createMatrix($iN, $iN)
	_lp_lacpy($mSource, $mR, "U", $iN, $iN, $iLDA, $iN)

	; 2. R⁻¹ (in-place, upper triangle)
	_lp_trtri($mR, "U", "N", $iN)
	If @error Then
		$mState.rankDeficient = True
		Return Null
	EndIf

	; 3. R⁻ᵀ (transpose: upper → lower, but _la_transpose creates full n×n)
	Local $mRinvT = _la_transpose($mR)

	; 4. Copy JPVT (lapmt modifies it!)
	Local $tJPVT = $mState.JPVT
	Local $iBytes = DllStructGetSize($tJPVT)
	Local $tJPVT_copy = DllStructCreate("INT[" & $iN & "]")
	DllCall("kernel32.dll", "NONE", "RtlCopyMemory", "PTR", DllStructGetPtr($tJPVT_copy), "PTR", DllStructGetPtr($tJPVT), "ULONG_PTR", $iBytes)

	; 5. Column permutation: R⁻ᵀ · Pᵀ (backward permutation)
	_lp_lapmt($mRinvT, $tJPVT_copy, False, $iN, $iN, $iN)

	; 6. Transpose → S = P · R⁻¹
	$mS = _la_transpose($mRinvT)

	; Jacobi-Equilibration back-transform: scale S rows by 1/S_eq
	; so that Qxx = S_scaled · S_scaledᵀ = diag(S⁻¹) · Qxx_eq · diag(S⁻¹)
	; Only for QR branch — SVD branch uses pre-equilibration A_orig (no scaling needed)
	If MapExists($mState, "EquilibrationScale") Then
		$mEqScale = $mState.EquilibrationScale
		$tEqS = $mEqScale.struct
		; S is n×n (OLS/GLM) or nFree×nFree (LSE) — scale each row i by 1/S_eq_i
		For $__i = 0 To $iN - 1
			_blas_scal($mS, 1.0 / DllStructGetData($tEqS, 1, $__i + 1), $__i, $iN, $iN)
		Next
	EndIf

	; 7. Qxx = S · Sᵀ (symmetric rank-k update)
	If MapExists($mState, "Q2") Then
		; LSE / GLM+Restrictions back-transformation: Qxx = Q₂ · S · Sᵀ · Q₂ᵀ = T · Tᵀ  where T = Q₂ · S
		$iNpar = $mState.nParams
		$mQ2 = $mState.Q2
		$mT = _blas_createMatrix($iNpar, $iN)
		_blas_gemm($mQ2, $mS, $mT, 1, 0, "N", "N", $iNpar, $iN, $iN)

		$mQxx = _blas_createMatrix($iNpar, $iNpar)
		_blas_syrk($mT, $mQxx, 1, 0, "U", "N", $iNpar, $iN)

		; fill lower triangle from upper (for symmetric display/access)
		__adj_fillLowerFromUpper($mQxx, $iNpar)

		$mResults.Qxx = $mQxx
		$mSystem.results = $mResults
		Return $mT  ; return T for cofactor matrix computation
	Else
		; OLS/GLM: Qxx = S · Sᵀ directly
		$mQxx = _blas_createMatrix($iN, $iN)
		_blas_syrk($mS, $mQxx, 1, 0, "U", "N", $iN, $iN)

		__adj_fillLowerFromUpper($mQxx, $iN)

		$mResults.Qxx = $mQxx
		$mSystem.results = $mResults
		Return $mS
	EndIf
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_computeQxxFromSVD
; Description ...: Computes Qxx via SVD: A = UΣVᵀ with rank truncation
; Syntax.........: __adj_computeQxxFromSVD($mSource, $mState, $iN, $iM, $fRCOND)
; Parameters ....: $mSource     - Source matrix
;                  $mState      - [ByRef] Solver state
;                  $iN          - Number of columns
;                  $iM          - Number of rows
;                  $fRCOND      - Rank threshold for singular value truncation
; Return values .: Success      - Map with .Qxx and .V
;                  Failure      - SetError($ADJ_ERR_SOLVER)
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: Also computes condition number stored in $mState.fConditionNumber.
; Related .......:
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __adj_computeQxxFromSVD($mSource, ByRef $mState, $iN, $iM, $fRCOND)
	; 1. Full SVD: A = U · Σ · Vᵀ
	Local $mSVD = _la_SVD(_la_duplicate($mSource))
	Local $iSvdErr = @error
	If $iSvdErr Or IsKeyword($mSVD) Then Return SetError($ADJ_ERR_SOLVER, 2, Null)

	; 2. V = Transpose(VT)  (VT is n×n, V is n×n)
	Local $mV = _la_transpose($mSVD.VT)

	; 3. Rank truncation: σᵢ > RCOND · σ_max
	;    S has min(M,N) entries
	Local $iMinMN = ($iM < $iN) ? $iM : $iN
	Local $fSigmaMax = DllStructGetData($mSVD.S.struct, 1, 1)
	Local $fThreshold = $fRCOND * $fSigmaMax
	Local $iEffectiveRank = 0

	; 4. Scale columns of V: vᵢ *= 1/σᵢ (only for σᵢ > threshold)
	For $i = 1 To $iMinMN
		Local $fSigma = DllStructGetData($mSVD.S.struct, 1, $i)
		If $fSigma > $fThreshold Then
			_blas_scal($mV, 1.0 / $fSigma, ($i - 1) * $iN, 1, $iN)
			$iEffectiveRank += 1
		Else
			_blas_scal($mV, 0, ($i - 1) * $iN, 1, $iN)
		EndIf
	Next
	; Zero columns beyond min(M,N) (underdetermined case)
	For $i = $iMinMN + 1 To $iN
		_blas_scal($mV, 0, ($i - 1) * $iN, 1, $iN)
	Next

	; 5. Qxx = V_scaled · V_scaledᵀ (symmetric rank-k update, upper triangle)
	Local $mQxx = _blas_createMatrix($iN, $iN)
	_blas_syrk($mV, $mQxx, 1, 0, "U", "N", $iN, $iN)

	; 6. Condition number (guard against rank 0)
	If $iEffectiveRank > 0 Then
		Local $fSigmaMin = DllStructGetData($mSVD.S.struct, 1, $iEffectiveRank)
		$mState.fConditionNumber = $fSigmaMax / $fSigmaMin
	Else
		$mState.fConditionNumber = Default
	EndIf

	; 7. Return Qxx + V_scaled (V is equivalent to S = P·R⁻¹ from QR path)
	Local $mResult[]
	$mResult.Qxx = $mQxx
	$mResult.V = $mV
	Return $mResult
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_computeCofactors
; Description ...: Computes Qvv (residuals) and Qŷ (adjusted observations) cofactor matrices
; Syntax.........: __adj_computeCofactors($mSystem, $mState, $mS)
; Parameters ....: $mSystem     - [ByRef] The adjustment system
;                  $mState      - [ByRef] Solver state
;                  $mS          - S matrix from __adj_computeQxx (or Null for CLS)
; Return values .: Success      - None (stores Qvv, Qyhat, sdv, sdy in $mSystem.results)
;                  Failure      - None
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: Three branches: CLS (direct via B), GLM (via hat matrix), OLS/LSE (via A·S·Sᵀ·Aᵀ).
;                  Back-transforms from whitened to original space.
; Related .......:
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __adj_computeCofactors(ByRef $mSystem, ByRef $mState, $mS)
	Local $mResults = $mSystem.results
	Local $sModel = $mState.modelType
	Local $iPobs = $mState.nObs  ; number of observations
	Local $fS0 = $mResults.s0

	If StringRegExp($sModel, 'CLS$') Then
		; ══════ CLS: direct Qvv computation without Qxx ══════
		; Qvv_w = Bᵀ · (B·Bᵀ)⁻¹ · B  via Cholesky
		Local $mB = $mState.B_orig
		Local $iMeq = $mB.rows

		; N = B · Bᵀ (meq × meq, symmetric)
		Local $mN = _blas_createMatrix($iMeq, $iMeq)
		_blas_syrk($mB, $mN, 1, 0, "L", "N", $iMeq, $iPobs)

		; Cholesky: N = L · Lᵀ
		_lp_potrf($mN, "L")
		If @error Then Return

		; S_cls = L⁻¹ · B (meq × pobs)
		Local $mScls = _la_duplicate($mB)
		_blas_trsm($mN, $mScls, 1, "L", "L", "N", "N", $iMeq, $iPobs)

		; Qvv_w = S_clsᵀ · S_cls (pobs × pobs, symmetric)
		Local $mQvv = _blas_createMatrix($iPobs, $iPobs)
		_blas_syrk($mScls, $mQvv, 1, 0, "U", "T", $iPobs, $iMeq)

		; back-transformation: Qvv = diag(σ)·Qvv_w·diag(σ) or L·Qvv_w·Lᵀ
		If StringRegExp($sModel, "^G(?!LM)") Then
			__adj_fillLowerFromUpper($mQvv, $iPobs)
			__adj_scaleSymmetricMatrixCholesky($mQvv, $mState.CovCholeskyL, $iPobs)
		ElseIf StringRegExp($sModel, "^W") Then
			__adj_scaleSymmetricMatrix($mQvv, $mState.Vector_ObsStdDev, $iPobs)
		EndIf

		__adj_fillLowerFromUpper($mQvv, $iPobs)
		$mResults.Qvv = $mQvv

		; Qŷ = P⁻¹ - Qvv  (P⁻¹ = Σₗₗ for generalized, diag(σ²) for diagonal)
		Local $mQyhat = _la_duplicate($mQvv)
		_la_scale($mQyhat, -1, True)
		If StringRegExp($sModel, "^G(?!LM)") Then
			_blas_axpy($mState.Matrix_Sigma, $mQyhat, 1, 0, 0, 1, 1, $iPobs * $iPobs)
		Else
			__adj_addDiagVariance($mQyhat, $mState.Vector_ObsStdDev, $iPobs)
		EndIf
		$mResults.Qyhat = $mQyhat

	ElseIf StringRegExp($sModel, 'GLM$') Then
		; ══════ GLM: Qvv via hat matrix in equation space ══════
		Local $mBglm = $mState.B_orig
		Local $mCholeskyM = $mState.CholeskyM
		$iMeq = $mState.nFormulas
		Local $iNpar = $mState.nParams
		Local $iScols = $mS.cols  ; nParams (no restrictions) or nFree (with restrictions, S = T = Q2·S)

		; Use original Ã saved before DGELSY (not the QR-factored Atilde)
		Local $mAtilde_orig = $mState.Atilde_orig

		; H_eq = Ã·Qxx·Ãᵀ = (Ã·S)·(Ã·S)ᵀ  where S from Qxx computation
		; U_eq = Ã · S (meq × iScols)
		Local $mUeq = _blas_createMatrix($iMeq, $iScols)
		_blas_gemm($mAtilde_orig, $mS, $mUeq, 1, 0, "N", "N", $iMeq, $iScols, $iNpar)

		; H_eq = U_eq · U_eqᵀ (meq × meq)
		Local $mHeq = _blas_createMatrix($iMeq, $iMeq)
		_blas_syrk($mUeq, $mHeq, 1, 0, "U", "N", $iMeq, $iScols)
		__adj_fillLowerFromUpper($mHeq, $iMeq)

		; IH = I - H_eq
		Local $mIH = _blas_createMatrix($iMeq, $iMeq)
		_lp_laset($mIH, "X", 0, 1, $iMeq, $iMeq)  ; IH = I
		_blas_axpy($mHeq, $mIH, -1, 0, 0, 1, 1, $iMeq * $iMeq)  ; IH = I - H_eq

		; Qvv_w = Tᵀ · (I-H̃) · T  where T = L⁻¹ · B̃
		; This correctly applies M⁻¹ = L⁻ᵀL⁻¹ from both sides of (I-H̃)
		; T = L⁻¹ · B_orig (meq × pobs)
		Local $mT_glm = _la_duplicate($mBglm)
		_blas_trsm($mCholeskyM, $mT_glm, 1, "L", "L", "N", "N", $iMeq, $iPobs)

		; Temp = (I-H̃) · T (meq × pobs)
		Local $mTemp = _blas_createMatrix($iMeq, $iPobs)
		_blas_gemm($mIH, $mT_glm, $mTemp, 1, 0, "N", "N", $iMeq, $iPobs, $iMeq)

		; Qvv_w = Tᵀ · Temp (pobs × pobs)
		$mQvv = _blas_createMatrix($iPobs, $iPobs)
		_blas_gemm($mT_glm, $mTemp, $mQvv, 1, 0, "T", "N", $iPobs, $iPobs, $iMeq)

		; back-transformation: Qvv = diag(σ)·Qvv_w·diag(σ) or L·Qvv_w·Lᵀ
		If StringRegExp($sModel, "^G(?!LM)") Then
			__adj_scaleSymmetricMatrixCholesky($mQvv, $mState.CovCholeskyL, $iPobs)
		ElseIf StringRegExp($sModel, "^W") Then
			__adj_scaleSymmetricMatrix($mQvv, $mState.Vector_ObsStdDev, $iPobs)
		EndIf

		$mResults.Qvv = $mQvv

		; Qŷ = P⁻¹ - Qvv
		$mQyhat = _la_duplicate($mQvv)
		_la_scale($mQyhat, -1, True)
		If StringRegExp($sModel, "^G(?!LM)") Then
			_blas_axpy($mState.Matrix_Sigma, $mQyhat, 1, 0, 0, 1, 1, $iPobs * $iPobs)
		Else
			__adj_addDiagVariance($mQyhat, $mState.Vector_ObsStdDev, $iPobs)
		EndIf
		$mResults.Qyhat = $mQyhat

	Else ; OLS/LSE
		; ══════ OLS/LSE: Qŷ = A·Qxx·Aᵀ = (A·S)·(A·S)ᵀ ══════
		If $mS = Null Then Return

		Local $mA = $mState.A_orig  ; NOTE: for WLS/WLSE this is the whitened A
		Local $iM = $iPobs  ; use nObs, not A.rows (may be augmented with LM rows)
		Local $iNcols = $mS.cols  ; n for OLS, nFree for LSE (but S/T already back-transformed)

		; U = A_orig · S (m × ncols)
		Local $mU = _blas_createMatrix($iM, $iNcols)
		_blas_gemm($mA, $mS, $mU, 1, 0, "N", "N", $iM, $iNcols, $mS.rows, $mA.rows)

		; Qŷ_w = U · Uᵀ (m × m, symmetric) — in whitened space for WLS/WLSE
		$mQyhat = _blas_createMatrix($iM, $iM)
		_blas_syrk($mU, $mQyhat, 1, 0, "U", "N", $iM, $iNcols)
		__adj_fillLowerFromUpper($mQyhat, $iM)

		; back-transformation: Qŷ = diag(σ)·Qŷ_w·diag(σ) or L·Qŷ_w·Lᵀ
		If StringRegExp($sModel, "^G(?!LM)") Then
			__adj_fillLowerFromUpper($mQyhat, $iM)
			__adj_scaleSymmetricMatrixCholesky($mQyhat, $mState.CovCholeskyL, $iM)
		ElseIf StringRegExp($sModel, "^W") Then
			__adj_scaleSymmetricMatrix($mQyhat, $mState.Vector_ObsStdDev, $iM)
		EndIf
		$mResults.Qyhat = $mQyhat

		; Qvv = P⁻¹ - Qŷ
		$mQvv = _la_duplicate($mQyhat)
		_la_scale($mQvv, -1, True)
		If StringRegExp($sModel, "^G(?!LM)") Then
			_blas_axpy($mState.Matrix_Sigma, $mQvv, 1, 0, 0, 1, 1, $iM * $iM)
		Else
			__adj_addDiagVariance($mQvv, $mState.Vector_ObsStdDev, $iM)
		EndIf
		$mResults.Qvv = $mQvv
	EndIf

	; standard deviations of adjusted observations and residuals
	If MapExists($mResults, "Qyhat") And $fS0 > 0 Then
		$mResults.sdy = _la_scale(_la_sqrtElements(_la_getDiag($mResults.Qyhat)), $fS0)
	EndIf
	If MapExists($mResults, "Qvv") And $fS0 > 0 Then
		$mResults.sdv = _la_scale(_la_sqrtElements(_la_getDiag($mResults.Qvv)), $fS0)
	EndIf

	$mSystem.results = $mResults
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_computeRedundancy
; Description ...: Computes redundancy diagonal diag(R) = diag(P·Qvv)
; Syntax.........: __adj_computeRedundancy($mSystem, $mState)
; Parameters ....: $mSystem     - [ByRef] The adjustment system
;                  $mState      - [ByRef] Solver state
; Return values .: None (stores redundancyDiag in $mSystem.results)
; Author ........: AspirinJunkie
; ===============================================================================================================================
Func __adj_computeRedundancy(ByRef $mSystem, ByRef $mState)
	Local $mResults = $mSystem.results
	If Not MapExists($mResults, "Qvv") Then Return

	If $mState.hasCovariances Then
		; generalized: diag(R) = diag(P · Qvv) where P = Σ⁻¹
		Local $iPobs = $mState.nObs
		Local $mPQvv = _la_duplicate($mResults.Qvv)
		_blas_trsm($mState.CovCholeskyL, $mPQvv.ptr, 1.0, "L", "L", "N", "N", $iPobs, $iPobs, $iPobs, $iPobs)
		_blas_trsm($mState.CovCholeskyL, $mPQvv.ptr, 1.0, "L", "L", "N", "T", $iPobs, $iPobs, $iPobs, $iPobs)
		$mResults.redundancyDiag = _la_getDiag($mPQvv)
	Else
		; diagonal: diag(R) = diag(Qvv) · p = diag(Qvv) / σ²
		Local $mDiagQvv = _la_getDiag($mResults.Qvv)
		$mResults.redundancyDiag = _la_mul($mDiagQvv, _la_squareElements($mState.Vector_ObsInvStdDev))
	EndIf

	$mSystem.results = $mResults
EndFunc

#EndRegion ; Cofactor matrices

#Region Matrix helpers

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_fillLowerFromUpper
; Description ...: Copies upper triangle to lower triangle of symmetric matrix (column-major)
; Syntax.........: __adj_fillLowerFromUpper($mMatrix, $iN)
; Parameters ....: $mMatrix     - [ByRef] Matrix map
;                  $iN          - Matrix dimension
; Return values .: None (in-place modification)
; Author ........: AspirinJunkie
; ===============================================================================================================================
Func __adj_fillLowerFromUpper(ByRef $mMatrix, $iN)
	; transpose → lower triangle of Mᵀ contains original upper triangle values
	; then lacpy "L" copies that lower triangle back to M
	Local $mT = _la_transpose($mMatrix)
	_lp_lacpy($mT, $mMatrix, "L", $iN, $iN, $iN, $iN)
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_scaleSymmetricMatrix
; Description ...: Scales symmetric matrix: M[i,j] *= σ[i]·σ[j] (= diag(σ)·M·diag(σ))
; Syntax.........: __adj_scaleSymmetricMatrix($mMatrix, $mSigmaVec, $iN)
; Parameters ....: $mMatrix     - [ByRef] Matrix map
;                  $mSigmaVec   - Vector of σ values
;                  $iN          - Matrix dimension
; Return values .: None (in-place modification)
; Author ........: AspirinJunkie
; ===============================================================================================================================
Func __adj_scaleSymmetricMatrix(ByRef $mMatrix, $mSigmaVec, $iN)
	; M := diag(σ) · M · diag(σ) via trmm (2 BLAS calls instead of 2n scal calls)
	Local $mD = _la_VectorToDiag($mSigmaVec)
	_blas_trmm($mD, $mMatrix.ptr, 1.0, "L", "L", "N", "N", $iN, $iN, $iN, $iN)   ; M := D · M
	_blas_trmm($mD, $mMatrix.ptr, 1.0, "R", "L", "N", "N", $iN, $iN, $iN, $iN)   ; M := M · D
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_scaleSymmetricMatrixCholesky
; Description ...: Back-transforms M := L·M·Lᵀ for cofactor recovery from whitened space
; Syntax.........: __adj_scaleSymmetricMatrixCholesky($mMatrix, $mCholeskyL, $iN)
; Parameters ....: $mMatrix     - [ByRef] Matrix map
;                  $mCholeskyL  - Lower Cholesky factor L
;                  $iN          - Matrix dimension
; Return values .: None (in-place modification)
; Author ........: AspirinJunkie
; ===============================================================================================================================
Func __adj_scaleSymmetricMatrixCholesky(ByRef $mMatrix, $mCholeskyL, $iN)
	_blas_trmm($mCholeskyL, $mMatrix.ptr, 1.0, "L", "L", "N", "N", $iN, $iN, $iN, $iN)   ; M := L · M
	_blas_trmm($mCholeskyL, $mMatrix.ptr, 1.0, "R", "L", "N", "T", $iN, $iN, $iN, $iN)   ; M := M · Lᵀ
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_addDiagVariance
; Description ...: Adds P⁻¹ = diag(σ²) to matrix diagonal: M[i,i] += σ[i]²
; Syntax.........: __adj_addDiagVariance($mMatrix, $mStdDevVec, $iN)
; Parameters ....: $mMatrix     - [ByRef] Matrix map
;                  $mStdDevVec  - Vector of standard deviations
;                  $iN          - Matrix dimension
; Return values .: None (in-place modification)
; Author ........: AspirinJunkie
; ===============================================================================================================================
Func __adj_addDiagVariance(ByRef $mMatrix, $mStdDevVec, $iN)
	Local $tSigma = $mStdDevVec.struct
	Local $tM = $mMatrix.struct
	For $i = 0 To $iN - 1
		Local $fSigma = DllStructGetData($tSigma, 1, $i + 1)
		Local $iPos = $i * $iN + $i + 1  ; diagonal position (1-indexed for DllStructGetData)
		DllStructSetData($tM, 1, DllStructGetData($tM, 1, $iPos) + $fSigma * $fSigma, $iPos)
	Next
EndFunc

#EndRegion ; Matrix helpers

#Region Global test

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_computeGlobalTest
; Description ...: Computes χ² global test (H₀: σ₀² = 1)
; Syntax.........: __adj_computeGlobalTest($mSystem)
; Parameters ....: $mSystem     - [ByRef] The adjustment system
; Return values .: Success      - None (stores globalTestPassed, globalTestT, globalTestLower, globalTestUpper, globalTestAlpha in $mSystem.results)
;                  Failure      - None
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: Test statistic T = vᵀPv. Bounds from chi² quantiles at alpha/2 and 1-alpha/2. Alpha from $mConfig.diagnostics.alpha.
; Related .......:
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __adj_computeGlobalTest(ByRef $mSystem)
	Local $mResults = $mSystem.results
	If $mResults.f <= 0 Then Return

	Local $fAlpha = ($mSystem.config).diagnostics.alpha
	Local $fT     = $mResults.vtpv
	Local $fLower = __adj_chi2Quantile($fAlpha / 2,     $mResults.f)
	Local $fUpper = __adj_chi2Quantile(1 - $fAlpha / 2, $mResults.f)

	$mResults.globalTestPassed = ($fT >= $fLower And $fT <= $fUpper)
	$mResults.globalTestT      = $fT
	$mResults.globalTestLower  = $fLower
	$mResults.globalTestUpper  = $fUpper
	$mResults.globalTestAlpha  = $fAlpha
	$mSystem.results = $mResults
EndFunc

#EndRegion ; Global test

#Region Chi2-distribution functions

; ═══════════════════════════════════════════════════════════════
; Chi²-quantile function (for global test)
; Ported from Stat_Distributions.au3 — dependency chain:
;   __adj_chi2Quantile → __adj_invgammp → __adj_gammp
;   → __adj_gser/__adj_gcf/__adj_gammpApprox → __adj_gammaLn
; ═══════════════════════════════════════════════════════════════

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_chi2Quantile
; Description ...: Chi² quantile: P(χ² ≤ x) = p → x
; Syntax.........: __adj_chi2Quantile($p, $iDF)
; Parameters ....: $p           - Probability in [0,1)
;                  $iDF         - Degrees of freedom
; Return values .: Quantile value, or SetError(1) if p out of range
; Author ........: AspirinJunkie
; ===============================================================================================================================
Func __adj_chi2Quantile($p, $iDF)
	If $p < 0 Or $p >= 1 Then Return SetError(1, 0, 0)
	Return 2 * __adj_invgammp($p, 0.5 * $iDF)
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_gammaLn
; Description ...: Log-Gamma function via ucrtbase.dll lgamma
; Syntax.........: __adj_gammaLn($z)
; Parameters ....: $z           - Input value
; Return values .: ln(Γ(z)) as Double
; Author ........: AspirinJunkie
; ===============================================================================================================================
Func __adj_gammaLn($z)
	Return DllCall("ucrtbase.dll", "double:cdecl", "lgamma", "double", $z)[0]
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_gammp
; Description ...: Regularized incomplete gamma function P(a,x)
; Syntax.........: __adj_gammp($a, $x)
; Parameters ....: $a           - Shape parameter (> 0)
;                  $x           - Upper limit (≥ 0)
; Return values .: P(a,x) in [0,1], or SetError(1)
; Author ........: AspirinJunkie
; ===============================================================================================================================
Func __adj_gammp($a, $x)
	If $x < 0 Or $a <= 0 Then Return SetError(1, 0, 0)
	If $x = 0 Then Return 0
	If Int($a) >= 100 Then Return __adj_gammpApprox($a, $x, 1)
	If $x < $a + 1 Then
		Return __adj_gser($a, $x)
	Else
		Return 1 - __adj_gcf($a, $x)
	EndIf
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_gser
; Description ...: Series representation of P(a,x) for x < a+1
; Syntax.........: __adj_gser($a, $x)
; Parameters ....: $a           - Shape parameter
;                  $x           - Upper limit
; Return values .: P(a,x) as Float
; Author ........: AspirinJunkie
; ===============================================================================================================================
Func __adj_gser($a, $x)
	Local Const $fEPS = 2.2204460492503131e-16
	Local $gln = __adj_gammaLn($a)
	Local $ap = $a, $sum = 1 / $a, $del = $sum
	Do
		$ap += 1
		$del *= $x / $ap
		$sum += $del
		If Abs($del) < Abs($sum) * $fEPS Then Return $sum * Exp(-$x + $a * Log($x) - $gln)
	Until 0
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_gcf
; Description ...: Continued fraction for Q(a,x) = 1-P(a,x) for x ≥ a+1
; Syntax.........: __adj_gcf($a, $x)
; Parameters ....: $a           - Shape parameter
;                  $x           - Upper limit
; Return values .: Q(a,x) as Float
; Author ........: AspirinJunkie
; ===============================================================================================================================
Func __adj_gcf($a, $x)
	Local Const $fEPS = 2.2204460492503131e-16
	Local Const $fFPMIN = 1.0020841800044864e-292
	Local $an, $del
	Local $gln = __adj_gammaLn($a), $b = $x + 1 - $a, $c = 1 / $fFPMIN, $d = 1 / $b, $h = $d

	For $i = 1 To 100000
		$an = -$i * ($i - $a)
		$b += 2
		$d = $an * $d + $b
		If Abs($d) < $fFPMIN Then $d = $fFPMIN
		$c = $b + $an / $c
		If Abs($c) < $fFPMIN Then $c = $fFPMIN
		$d = 1 / $d
		$del = $d * $c
		$h *= $del
		If Abs($del - 1) <= $fEPS Then ExitLoop
	Next
	Return Exp(-$x + $a * Log($x) - $gln) * $h
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_gammpApprox
; Description ...: Gauss-Legendre quadrature approximation for P(a,x) when a ≥ 100
; Syntax.........: __adj_gammpApprox($a, $x, $iPsig)
; Parameters ....: $a           - Shape parameter
;                  $x           - Upper limit
;                  $iPsig       - 1 for P, other for Q
; Return values .: Approximation of P(a,x) or Q(a,x)
; Author ........: AspirinJunkie
; ===============================================================================================================================
Func __adj_gammpApprox($a, $x, $iPsig = 1)
	Local Static $aY[18] = [0.0021695375159141994, 0.011413521097787704, 0.027972308950302116, _
		0.051727015600492421, 0.082502225484340941, 0.12007019910960293, 0.16415283300752470, _
		0.21442376986779355, 0.27051082840644336, 0.33199876341447887, 0.39843234186401943, _
		0.46931971407375483, 0.54413605556657973, 0.62232745288031077, 0.70331500465597174, _
		0.78649910768313447, 0.87126389619061517, 0.95698180152629142]
	Local Static $aW[18] = [0.0055657196642445571, 0.012915947284065419, 0.020181515297735382, _
		0.027298621498568734, 0.034213810770299537, 0.040875750923643261, 0.047235083490265582, _
		0.053244713977759692, 0.058860144245324798, 0.064039797355015485, 0.068745323835736408, _
		0.072941885005653087, 0.076598410645870640, 0.079687828912071670, 0.082187266704339706, _
		0.084078218979661945, 0.085346685739338721, 0.085983275670394821]

	Local $a1 = $a - 1
	Local $lna1 = $a1 >= 0 ? Log($a1) : 0
	Local $sqrta1 = Sqrt($a1)
	Local $gln = __adj_gammaLn($a)
	Local $xu = $x > $a1 ? _
		($a1 + 11.5 * $sqrta1 > $x + 6 * $sqrta1 ? $a1 + 11.5 * $sqrta1 : $x + 6 * $sqrta1) : _
		(0 > $a1 - 7.5 * $sqrta1 ? (0 > $x - 5 * $sqrta1 ? 0 : $x - 5 * $sqrta1) : _
		($a1 - 7.5 * $sqrta1 > $x - 5 * $sqrta1 ? $x - 5 * $sqrta1 : $a1 - 7.5 * $sqrta1))

	Local $sum = 0, $t
	For $j = 0 To 17
		$t = $x + ($xu - $x) * $aY[$j]
		$sum += $aW[$j] * Exp(-($t - $a1) + $a1 * (Log($t) - $lna1))
	Next
	Local $ans = $sum * ($xu - $x) * Exp($a1 * ($lna1 - 1) - $gln)
	Return $iPsig = 1 ? ($ans > 0 ? 1 - $ans : -$ans) : ($ans >= 0 ? $ans : 1 + $ans)
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_invgammp
; Description ...: Inverse regularized incomplete gamma function: P(a,x) = p → x
; Syntax.........: __adj_invgammp($p, $a)
; Parameters ....: $p           - Target probability
;                  $a           - Shape parameter
; Return values .: Success      - x such that P(a,x) = p
;                  Failure      - SetError(1) if a ≤ 0
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: Initial approximation followed by up to 12 Newton-Raphson iterations.
; Related .......:
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __adj_invgammp($p, $a)
	Local $x, $err, $t, $u, $pp, $lna1, $afac, $a1 = $a - 1
	Local $gln = __adj_gammaLn($a)

	If $a <= 0 Then Return SetError(1, 0, 0)
	If $p >= 1 Then Return $a + 100 * Sqrt($a) > 100 ? $a + 100 * Sqrt($a) : 100
	If $p <= 0 Then Return 0

	If $a > 1 Then
		$lna1 = Log($a1)
		$afac = Exp($a1 * ($lna1 - 1) - $gln)
		$pp = $p < 0.5 ? $p : 1 - $p
		$t = Sqrt(-2 * Log($pp))
		$x = (2.30753 + $t * 0.27061) / (1 + $t * (0.99229 + $t * 0.04481)) - $t
		If $p < 0.5 Then $x = -$x
		$x = $a * ((1 - 1 / (9 * $a) - $x / (3 * Sqrt($a))) ^ 3)
		If $x < 1e-3 Then $x = 1e-3
	Else
		$t = 1 - $a * (0.253 + $a * 0.12)
		$x = $p < $t ? ($p / $t) ^ (1 / $a) : 1 - Log(1 - ($p - $t) / (1 - $t))
	EndIf

	For $j = 0 To 11
		If $x <= 0 Then Return 0
		$err = __adj_gammp($a, $x) - $p
		$t = ($a > 1) ? $afac * Exp(-($x - $a1) + $a1 * (Log($x) - $lna1)) : Exp(-$x + $a1 * Log($x) - $gln)
		$u = $err / $t
		Local $fMin = $u * (($a - 1) / $x - 1)
		$t = $u / (1 - 0.5 * ($fMin < 1 ? $fMin : 1))
		$x -= $t
		If $x <= 0 Then $x = 0.5 * ($x + $t)
		If Abs($t) < 1e-8 * $x Then ExitLoop
	Next
	Return $x
EndFunc

#EndRegion ; Chi2-distribution functions
