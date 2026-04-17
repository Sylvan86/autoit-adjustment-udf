#include-once

; #INDEX# =======================================================================================================================
; Title .........: AdjustmentSolver
; AutoIt Version : 3.3.16.1
; Description ...: Three-tier solver architecture for geodetic adjustment: VCE outer loop,
;                  Gauss-Newton/Levenberg-Marquardt iteration, and LAPACK linear solver dispatch.
;                  Supports OLS/WLS/GLS, LSE, CLS, GLM and their generalized variants.
; Author(s) .....: AspirinJunkie
; Dll ...........: libopenblas.dll
; ===============================================================================================================================


#Region VCE outer loop (Tier 1)

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_computeVtPV
; Description ...: Computes vᵀPv = ‖L⁻¹·v‖² (generalized) or ‖v/σ‖² (diagonal)
; Syntax.........: __adj_computeVtPV($mState)
; Parameters ....: $mState      - [ByRef] Solver state map
; Return values .: Float — weighted sum of squared residuals
; Author ........: AspirinJunkie
; ===============================================================================================================================
Func __adj_computeVtPV(ByRef $mState)
	Local $mResiduals = $mState.r
	If $mState.hasCovariances Then
		; generalized: vᵀPv = ‖L⁻¹·v‖² where Σₗₗ = L·Lᵀ
		Local $mTmp = _la_duplicate($mResiduals)
		_blas_trsv($mState.CovCholeskyL, $mTmp, "L", "N", "N", $mState.nObs, 1, $mState.nObs)
		Return _blas_nrm2($mTmp) ^ 2
	Else
		; diagonal: vᵀPv = Σ(vᵢ/σᵢ)²
		Return _blas_nrm2(_la_mul($mState.Vector_ObsInvStdDev, $mResiduals)) ^ 2
	EndIf
EndFunc


; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_estimateVCE
; Description ...: Helmert variance component estimation outer loop
; Syntax.........: __adj_estimateVCE($mSystem, $mState)
; Parameters ....: $mSystem     - [ByRef] Adjustment system map
;                  $mState      - [ByRef] Solver state map
; Return values .: Success      - True (VCE converged or single-pass without VCE)
;                  Failure      - SetError($ADJ_ERR_*) on solver or factorization error
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: Runs GN/LM solve + statistics iteratively. Updates group weights until
;                  convergence or max 15 iterations. For generalized models: block-scales Σₗₗ.
;                  For diagonal: weight update pᵢ/σ̂²_k.
; Related .......: __adj_solveNonlinear, __adj_computeStatistics, __adj_updateWhiteningVectors
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __adj_estimateVCE(ByRef $mSystem, ByRef $mState)
	Local $mConfig = $mSystem.config
	Local $bDoVCE  = $mConfig.vce

	; VKS with covariances: skip if only 1 group (equivalent to s₀² estimation)
	If $mState.hasCovariances And $bDoVCE Then
		If Not MapExists($mState, "vceGroupOffsets") Or UBound(MapKeys($mState.vceGroupOffsets)) <= 1 Then
			$bDoVCE = False
		EndIf
	EndIf

	; guard: cross-group covariances detected (should not reach here, safety net)
	If $bDoVCE And ($mSystem.model)._hasCrossGroupCovar Then
		Return SetError($ADJ_ERR_INPUT, 5, False)
	EndIf

	; build group index map: varComp group → array of observation indices
	; and reverse index map: index → observation name (for direct access)
	Local $mGroupIndices[], $aIdxToName[0]
	Local $mModel = $mSystem.model
	If $bDoVCE Then
		Local $mObservations = $mModel.obs, $mIdxObs = $mState.idxObs
		ReDim $aIdxToName[$mState.nObs]
		For $sObsName In MapKeys($mIdxObs)
			Local $iOIdx = $mIdxObs[$sObsName]
			$aIdxToName[$iOIdx] = $sObsName
			Local $sGroup = ($mObservations[$sObsName]).varComp
			If Not MapExists($mGroupIndices, $sGroup) Then
				Local $aTmp[1] = [$iOIdx]
				$mGroupIndices[$sGroup] = $aTmp
			Else
				Local $aIdx = $mGroupIndices[$sGroup]
				ReDim $aIdx[UBound($aIdx) + 1]
				$aIdx[UBound($aIdx) - 1] = $iOIdx
				$mGroupIndices[$sGroup] = $aIdx
			EndIf
		Next
	EndIf

	Local $iVCEIter = 0, $mGroupResults[]
	For $iVCE = 1 To 15
		$iVCEIter = $iVCE
		; reset iteration state for fresh GN/LM convergence (VKS iterations > 1)
		If $iVCE > 1 Then
			If MapExists($mState, "nIterations") Then MapRemove($mState, "nIterations")
			If MapExists($mState, "r_accumulated") Then MapRemove($mState, "r_accumulated")

			; reset Marquardt scaling state (D depends on weights which VKS changes)
			If MapExists($mState, "LM_D") Then MapRemove($mState, "LM_D")
			If MapExists($mState, "LM_gradient") Then MapRemove($mState, "LM_gradient")
			If MapExists($mState, "LM_step") Then MapRemove($mState, "LM_step")
			If MapExists($mState, "EquilibrationScale") Then MapRemove($mState, "EquilibrationScale")
		EndIf

		__adj_solveNonlinear($mSystem, $mState)
		If @error Then Return SetError(@error, @extended, False)

		; compute the residuals vector r
		__adj_computeResiduals($mSystem, $mState)

		; compute statistics: Qxx, s₀, sdx, Qvv, redundancy (Phase 5)
		__adj_computeStatistics($mSystem, $mState)

		; don't do a variance component estimation if not chosen by the user
		If Not $bDoVCE Then ExitLoop

		; shared variables for both VKS branches (covariance + diagonal)
		Local $mResults = $mSystem.results
		Local $mVecV = $mState.r
		Local $mVecR = $mResults.redundancyDiag
		Local $bConverged = True
		Local $fVtPV, $fRedundancySum, $fSigmaSq, $bClipped
		Local $mGrp[]

		; === Generalized VKS: block-diagonal Σₗₗ ===
		If $mState.hasCovariances Then
			Local $mGroupOffsets = $mState.vceGroupOffsets
			Local $iNObs = $mState.nObs

			; Step 1: z = L_original⁻¹ · v (TRSV with original Cholesky)
			Local $mZ = _la_duplicate($mState.r)
			_blas_trsv($mState.CovCholeskyL_Original, $mZ, "L", "N", "N", $iNObs, 1, $iNObs)

			; Step 2+3: per-group vᵀPv_k and σ̂²_k
			$bConverged = True
			Local $mVCELastSigma2 = MapExists($mState, "vceLastSigma2") ? $mState.vceLastSigma2 : Null

			For $sGroup In MapKeys($mGroupOffsets)
				Local $mOff = $mGroupOffsets[$sGroup]
				Local $iStart = $mOff.start, $iCount = $mOff.count

				; group vtPv via nrm2 with offset (no temporary vector needed)
				$fVtPV = _blas_nrm2($mZ, $iStart, 1, $iCount) ^ 2

				; group redundancy sum via asum (contiguous, non-negative r_i)
				$fRedundancySum = _blas_asum($mVecR, $iStart, 1, $iCount)

				; σ̂²_k = vₖᵀQₖ⁻¹vₖ / rₖ  (absolute estimate against original Qₖ)
				$fSigmaSq = 1.0
				If $fRedundancySum > 1e-10 Then
					$fSigmaSq = $fVtPV / $fRedundancySum
				EndIf

				; protection against very small / negative variance estimates:
				; floor at $__ADJ_VCE_MIN_SIGMA2 instead of resetting to 1.0 (which would
				; throw away the partial information that the group's variance is in fact small)
				$bClipped = ($fSigmaSq < $__ADJ_VCE_MIN_SIGMA2)
				If $bClipped Then $fSigmaSq = $__ADJ_VCE_MIN_SIGMA2

				; store per-group results
				$mGrp = Null
				Local $mGrpCov[]
				$mGrpCov.s0      = Sqrt($fSigmaSq)
				$mGrpCov.sigma2  = $fSigmaSq
				$mGrpCov.vtpv    = $fVtPV
				$mGrpCov.r       = $fRedundancySum
				$mGrpCov.clipped = $bClipped
				$mGroupResults[$sGroup] = $mGrpCov

				; convergence test: |σ̂²_k - σ̂²_k_old| / σ̂²_k_old < 0.05
				Local $fOld = 1.0
				If IsMap($mVCELastSigma2) And MapExists($mVCELastSigma2, $sGroup) Then $fOld = $mVCELastSigma2[$sGroup]
				If Abs($fSigmaSq - $fOld) / $fOld >= 0.05 Then $bConverged = False

				; store for next iteration
				If Not IsMap($mVCELastSigma2) Then
					Local $mTmpLast[]
					$mVCELastSigma2 = $mTmpLast
				EndIf
				$mVCELastSigma2[$sGroup] = $fSigmaSq
			Next
			$mState.vceLastSigma2 = $mVCELastSigma2

			; Step 4: rebuild Σₗₗ from original with block-scaling
			Local $mSigmaNew = _la_duplicate($mState.Matrix_Sigma_Original)
			For $sGroup In MapKeys($mGroupOffsets)
				$mOff = $mGroupOffsets[$sGroup]
				$iStart = $mOff.start
				$iCount = $mOff.count
				$fSigmaSq = ($mGroupResults[$sGroup]).sigma2

				; extract block, scale, write back via _lp_lacpy
				Local $mBlock = _la_extractBlock($mSigmaNew, $iStart + 1, $iStart + 1, $iStart + $iCount, $iStart + $iCount)
				_blas_scal($mBlock, $fSigmaSq)
				Local $iOffset = ($iStart * $iNObs + $iStart) * ($mSigmaNew.datatype = "DOUBLE" ? 8 : 4)
				_lp_lacpy($mBlock.ptr, $mSigmaNew.ptr + $iOffset, "X", $iCount, $iCount, $iCount, $iNObs)
			Next

			; Step 5: Cholesky factorization of new Σₗₗ
			Local $mCholeskyNew = _la_duplicate($mSigmaNew)
			_lp_potrf($mCholeskyNew, "L")
			If @error Then Return SetError($ADJ_ERR_NOT_POS_DEF, @error, False)

			; Step 6: update state
			$mState.Matrix_Sigma = $mSigmaNew
			$mState.CovCholeskyL = $mCholeskyNew

			If $bConverged Then
				; update vtpv and s₀ with final Σₗₗ
				; uses CovCholeskyL = L_neu (correct: vtpv = vᵀΣₗₗ_neu⁻¹v)
				$mResults = $mSystem.results
				$mResults.vtpv = __adj_computeVtPV($mState)
				$mResults.r2sum = $mResults.vtpv
				If $mResults.f > 0 Then $mResults.s0 = Sqrt($mResults.vtpv / $mResults.f)
				$mResults.vceGroups = $mGroupResults
				$mSystem.results = $mResults
				ExitLoop
			EndIf

			; next VKS iteration: no weight/whitening-vector update needed
			; (generalized models use CovCholeskyL directly in whitening)
			ContinueLoop
		EndIf

		; --- Helmert VCE: estimate variance components per group ---
		$mResults = $mSystem.results
		$mVecV = $mState.r
		$mVecR = $mResults.redundancyDiag
		$bConverged = True

		$mModel = $mSystem.model
		$mObservations = $mModel.obs
		$mIdxObs = $mState.idxObs

		For $sGroup In MapKeys($mGroupIndices)
			$aIdx = $mGroupIndices[$sGroup]

			; compute group redundancy sum rₖ and vₖᵀPₖvₖ
			$fRedundancySum = 0
			$fVtPV = 0
			For $j = 0 To UBound($aIdx) - 1
				Local $iIdx = $aIdx[$j]
				$fRedundancySum += __adj_vecGet($mVecR, $iIdx)

				Local $fVi = __adj_vecGet($mVecV, $iIdx)
				Local $fPi = ($mObservations[$aIdxToName[$iIdx]]).weight
				$fVtPV += $fVi^2 * $fPi
			Next

			; variance component estimate σ̂²_k = vₖᵀPₖvₖ / rₖ
			$fSigmaSq = 1.0
			If $fRedundancySum > 1e-10 Then
				$fSigmaSq = $fVtPV / $fRedundancySum
			EndIf

			; protection against negative/very small variance estimates (Step 6.3):
			; floor at $__ADJ_VCE_MIN_SIGMA2 instead of resetting to 1.0
			$bClipped = ($fSigmaSq < $__ADJ_VCE_MIN_SIGMA2)
			If $bClipped Then
				$fSigmaSq = $__ADJ_VCE_MIN_SIGMA2
			EndIf

			; store per-group VCE results (overwritten each iteration, final values kept)
			$mGrp.s0      = Sqrt($fSigmaSq)
			$mGrp.sigma2  = $fSigmaSq
			$mGrp.vtpv    = $fVtPV
			$mGrp.r       = $fRedundancySum
			$mGrp.clipped = $bClipped
			$mGroupResults[$sGroup] = $mGrp

			; convergence test: |1 - σ̂²_k| < 0.05
			If Abs(1 - $fSigmaSq) >= 0.05 Then $bConverged = False

			; update weights: pᵢ,neu = pᵢ / σ̂²_k
			For $j = 0 To UBound($aIdx) - 1
				$iIdx = $aIdx[$j]
				Local $sName = $aIdxToName[$iIdx]
				Local $mObs = $mObservations[$sName]
				$mObs.weight = $mObs.weight / $fSigmaSq
				$mObservations[$sName] = $mObs
			Next
		Next

		$mModel.obs = $mObservations
		$mSystem.model = $mModel ;! MAP WRITE-BACK

		; update whitening vectors from new weights (Step 6.4)
		; σᵢ,eff = 1/√pᵢ, 1/σᵢ,eff = √pᵢ
		__adj_updateWhiteningVectors($mSystem, $mState)

		If $bConverged Then
			; recompute vtpv and s₀ with final weights so results are consistent
			Local $fVtPvFinal = 0
			For $sGroup In MapKeys($mGroupIndices)
				Local $aIdx2 = $mGroupIndices[$sGroup]
				For $j = 0 To UBound($aIdx2) - 1
					Local $iIdx2 = $aIdx2[$j]
					Local $fVi2 = __adj_vecGet($mVecV, $iIdx2)
					Local $fPi2 = ($mObservations[$aIdxToName[$iIdx2]]).weight
					$fVtPvFinal += $fVi2^2 * $fPi2
				Next
			Next
			$mResults["vtpv"] = $fVtPvFinal
			$mResults["s0"] = Sqrt($fVtPvFinal / $mResults.f)
			$mSystem.results = $mResults
			ExitLoop
		EndIf

	Next

	; --- Store VCE results in $mSystem.results ---
	If $bDoVCE Then
		Local $mVCEResults = $mSystem.results
		$mVCEResults.vceConverged  = $bConverged
		$mVCEResults.vceIterations = $iVCEIter
		$mVCEResults.vceGroups     = $mGroupResults
		$mSystem.results = $mVCEResults
	EndIf
EndFunc

#EndRegion ; VCE outer loop (Tier 1)

#Region Nonlinear iteration (Tier 2)

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_updateWhiteningVectors
; Description ...: Updates Vector_ObsStdDev and Vector_ObsInvStdDev from current observation weights
; Syntax.........: __adj_updateWhiteningVectors($mSystem, $mState)
; Parameters ....: $mSystem     - [ByRef] Adjustment system map
;                  $mState      - [ByRef] Solver state map
; Return values .: None
; Author ........: AspirinJunkie
; ===============================================================================================================================
Func __adj_updateWhiteningVectors(ByRef $mSystem, ByRef $mState)
	Local $mObservations = ($mSystem.model).obs, $mIdxObs = $mState.idxObs
	Local $tS  = $mState.Vector_ObsStdDev.struct
	Local $tSI = $mState.Vector_ObsInvStdDev.struct

	For $sObsName In MapKeys($mIdxObs)
		Local $iIdx = $mIdxObs[$sObsName]
		Local $fWeight = ($mObservations[$sObsName]).weight
		Local $fSigma = 1 / Sqrt($fWeight)
		DllStructSetData($tS,  1, $fSigma,      $iIdx + 1)
		DllStructSetData($tSI, 1, Sqrt($fWeight), $iIdx + 1)
	Next
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_solveNonlinear
; Description ...: Gauss-Newton or Levenberg-Marquardt iteration loop
; Syntax.........: __adj_solveNonlinear($mSystem, $mState)
; Parameters ....: $mSystem     - [ByRef] Adjustment system map
;                  $mState      - [ByRef] Solver state map
; Return values .: Success      - True (converged within maxIterations)
;                  Failure      - SetError($ADJ_ERR_*) on solver error
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: LM: gain ratio evaluation, step acceptance/rejection, rollback, deferred GLM
;                  gain ratio, final clean GN iterations. GN: convergence via parameter correction
;                  (OLS/LSE) or relative residual change (CLS/GLM). Adaptive stagnation detection.
; Related .......: __adj_solveIteration, __adj_updateParameters, __adj_computeResiduals
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __adj_solveNonlinear(ByRef $mSystem, ByRef $mState)
	Local $mConfig = $mSystem.config
	Local $iMaxIterations = $mConfig.maxIterations
	Local $fTolerance     = $mConfig.tolerance

	; CLS still unsupported with LM (no parameters to dampen) — fallback to GN
	If $mConfig.algorithm = "LM" And StringRegExp($mState.modelType, 'CLS$') Then
		$mConfig.algorithm = "GN"
		$mSystem.config = $mConfig
	EndIf

	If $mConfig.algorithm = "LM" Then
		Local $mParamsOld, $mStateOld
		Local $bIsLSE = StringRegExp($mState.modelType, '^W?G?LSE$')
		Local $bIsGLM = StringRegExp($mState.modelType, '(CLS|GLM)$')
		Local $fR2sumPrev_GLM = 0     ; previous r2sum for GLM/CLS convergence test

		; GLM deferred gain ratio state: F(x+h) is only available after the NEXT
		; iteration's __adj_solveIteration (which computes r2sum0 = ‖w̃(x)‖² at the new point).
		; For OLS/WLS/GLS: re-evaluate functions at x+h and use vᵀPv immediately.
		Local $bGLM_pending = False    ; deferred check waiting for next r2sum0
		Local $fGLM_prevCost = 0      ; F(x) from the step being evaluated
		Local $fGLM_prevFPred = 0     ; predicted decrease from the step
		Local $mGLM_rollbackParams    ; params before the tentative step
		Local $mGLM_rollbackState     ; state before the tentative step

		; use the Levenberg-Marquardt algorithm for iterative solving the nonlinear model
		For $i = 1 To $iMaxIterations
			; backup current state for potential rollback (before solve modifies state)
			$mParamsOld = ($mSystem.model).params
			$mStateOld  = $mState

			; solve system with current values for the parameters and observations
			__adj_solveIteration($mSystem, $mState)
			If @error Then Return SetError(@error, @extended, False)

			; --- GLM deferred gain ratio check ---
			; After solve, r2sum0 = F(x_current) = actual cost at the current point.
			; If a previous step is pending evaluation, compare costs now.
			If $bIsGLM And $bGLM_pending Then
				Local $fGLM_actualCost = $mState.r2sum0
				Local $fGLM_fActual = $fGLM_prevCost - $fGLM_actualCost
				Local $fGLM_rho = ($fGLM_prevFPred > 0) ? ($fGLM_fActual / $fGLM_prevFPred) : 0

				; λ update
				If $fGLM_rho > 0.75 Then
					$mState.fLMLambda = _Max(1e-7, $mState.fLMLambda / 3)
				ElseIf $fGLM_rho < 0.25 Then
					$mState.fLMLambda = _Min(1e7, $mState.fLMLambda * 2)
				EndIf

				If $fGLM_rho <= 0 Then
					; Reject step: rollback to pre-step state, but keep increased lambda
					Local $fRejLambda = $mState.fLMLambda
					Local $mModel_rb1 = $mSystem.model
					$mModel_rb1.params = $mGLM_rollbackParams
					$mSystem.model = $mModel_rb1
					$mState         = $mGLM_rollbackState
					$mState.fLMLambda = $fRejLambda  ; preserve increased lambda after rollback
					$bGLM_pending = False
					ContinueLoop ; re-solve with increased λ at the old point
				EndIf
				$bGLM_pending = False
			EndIf

			; solution converged - stop iteration
			; GLM/CLS: |xd|→0 alone is unreliable because relinearisation at l+v
			; can shift the system; require r2sum stagnation instead (analogous to GN path,
			; cf. Neitzel/Petrovic 2008).  All other model types: classic |xd|<tol.
			If $bIsGLM Then
				If $fR2sumPrev_GLM > 1e-30 _
					And Abs($mState.r2sum - $fR2sumPrev_GLM) / $fR2sumPrev_GLM < $fTolerance Then ExitLoop
				$fR2sumPrev_GLM = $mState.r2sum
			Else
				If Abs(_blas_amax($mState.xd)) < $fTolerance Then ExitLoop
			EndIf

			; tentatively apply parameter update
			__adj_updateParameters($mSystem, $mState)

			If Not $bIsLSE And Not $bIsGLM Then
				; For unconstrained OLS/WLS/GLS: evaluate step via observation residuals
				__adj_computeResiduals($mSystem, $mState)
				$mState.r2sum = __adj_computeVtPV($mState)
			EndIf

			; --- Gain ratio for non-GLM models (immediate) ---
			If Not $bIsGLM Then
				; termination: check for stagnation
				If $mState.r2sum0 > 0 Then
					If ((Abs($mState.r2sum0 - $mState.r2sum) / $mState.r2sum0)) < $__ADJ_STAGNATION_TOL Then ExitLoop
				EndIf

				; evaluate step quality via gain ratio
				If $mState.r2sum0 > 0 Then
					Local $fNewLambda = __adj_updateLMDamping($mState)
					Local $bAccepted  = @extended
					$mState.fLMLambda = $fNewLambda
					If Not $bAccepted Then
						; reject step: rollback parameters and state, but keep updated lambda
						Local $mModel_rb2 = $mSystem.model
						$mModel_rb2.params = $mParamsOld
						$mSystem.model = $mModel_rb2
						$mState         = $mStateOld
						$mState.fLMLambda = $fNewLambda  ; preserve increased lambda after rollback
					EndIf
				EndIf
			EndIf

			; --- GLM: save state for deferred gain ratio ---
			If $bIsGLM And MapExists($mState, "fLMLambda") Then
				$fGLM_prevCost = $mState.r2sum0  ; F(x) before step

				; compute fPred = hᵀg + λ‖Dh‖² from current gradient/step
				Local $mStep_g = MapExists($mState, "LM_step") ? $mState.LM_step : $mState.xd
				Local $iStepN_g = $mStep_g.size
				$fGLM_prevFPred = 0
				If MapExists($mState, "LM_gradient") And MapExists($mState, "LM_D") Then
					Local $fHtG_g = _blas_dot($mStep_g, $mState.LM_gradient, 0, 0, 1, 1, $iStepN_g)
					; ‖Dh‖² via sbmv: treat D as diagonal band matrix (K=0), compute Dh = D.*h, then ‖Dh‖²
					Local $mDh_g = _blas_createVector($iStepN_g)
					_blas_sbmv($mState.LM_D, $mStep_g, $mDh_g, 1.0, 0.0, 0, "U", $iStepN_g, 1, 1, 1)
					Local $fDh2_g = _blas_nrm2($mDh_g, 0, 1, $iStepN_g) ^ 2
					$fGLM_prevFPred = $fHtG_g + $mState.fLMLambda * $fDh2_g
				EndIf
				If $fGLM_prevFPred <= 0 Then
					$fGLM_prevFPred = $mState.fLMLambda * _blas_nrm2($mStep_g, 0, 1, $iStepN_g)^2
				EndIf

				$mGLM_rollbackParams = $mParamsOld   ; params BEFORE the step
				$mGLM_rollbackState  = $mStateOld    ; state BEFORE the step
				$bGLM_pending = True
			EndIf

		Next

		; final clean GN iterations without LM damping for correct statistics (vᵀPv, s₀, Qxx)
		If MapExists($mState, "fLMLambda") Then
			; Switch to GN so __adj_solveIteration does not re-initialize LM
			$mConfig.algorithm = "GN"
			$mSystem.config = $mConfig
			MapRemove($mState, "fLMLambda")
			If MapExists($mState, "LM_D") Then MapRemove($mState, "LM_D")
			If MapExists($mState, "LM_gradient") Then MapRemove($mState, "LM_gradient")
			If MapExists($mState, "LM_step") Then MapRemove($mState, "LM_step")
			If MapExists($mState, "EquilibrationScale") Then MapRemove($mState, "EquilibrationScale")

			; For GLM/CLS: reset accumulated observation residuals for fresh GN convergence.
			; The LM path accumulates v incrementally via damped steps, which for nonlinear
			; problems leads to suboptimal v (path-dependent linearization). Resetting lets
			; the clean GN phase converge to the correct (x*, v*) from scratch.
			If StringRegExp($mState.modelType, '(CLS|GLM)$') Then
				If MapExists($mState, "nIterations") Then MapRemove($mState, "nIterations")
				If MapExists($mState, "r_accumulated") Then MapRemove($mState, "r_accumulated")
			EndIf

			Local $fR2sumPrevLM = 0
			Local $fMaxDxPrevLM = 1e308, $iStagnLM = 0
			For $i = 1 To $iMaxIterations
				__adj_solveIteration($mSystem, $mState)
				If @error Then Return SetError(@error, @extended, False)
				If StringRight($mState.modelType, 3) <> "CLS" Then __adj_updateParameters($mSystem, $mState)
				; GLM/CLS: r2sum-based convergence (same reasoning as GN path)
				If StringRegExp($mState.modelType, '(CLS|GLM)$') Then
					Local $fDeltaR2LM = Abs($mState.r2sum - $fR2sumPrevLM)
					If $fR2sumPrevLM > 1e-30 And $fDeltaR2LM / $fR2sumPrevLM < $fTolerance Then ExitLoop
					$fR2sumPrevLM = $mState.r2sum
				Else
					Local $fMaxDxLM = Abs(_blas_amax($mState.xd))
					If $fMaxDxLM < $fTolerance Then ExitLoop
					; adaptive stagnation: corrections no longer decreasing → numerical accuracy floor
					If $fMaxDxLM > 0.5 * $fMaxDxPrevLM Then
						$iStagnLM += 1
						If $iStagnLM >= 3 Then ExitLoop
					Else
						$iStagnLM = 0
					EndIf
					$fMaxDxPrevLM = $fMaxDxLM
				EndIf
			Next

			; restore original algorithm name for display/results
			$mConfig.algorithm = "LM"
			$mSystem.config = $mConfig
		EndIf

	Else
		; use gauss-newton algorithm for iterative solving the nonlinear model

		Local $fR2sumPrev = 0  ; for CLS/GLM convergence test (relative residual change)
		Local $bIsCLS = StringRegExp($mState.modelType, 'CLS$')
		Local $bIsGLM = StringRegExp($mState.modelType, 'GLM$')
		Local $fMaxDxPrev = 1e308, $iStagnation = 0

		For $i = 1 To $iMaxIterations
			; solve the current linear state of the system --> xd and/or r
			__adj_solveIteration($mSystem, $mState)
			If @error Then Return SetError(@error, @extended, False)

			; improve the parameter values
			If Not $bIsCLS Then __adj_updateParameters($mSystem, $mState)    ; x = x + xd

			; in linear case no further processing needed
			; Exception: GLM models are always effectively nonlinear (parameter × observation interaction)
			If Not $mState.isNonlinear And Not $bIsGLM Then ExitLoop

			; convergence test (model-dependent)
			If $bIsCLS Or $bIsGLM Then
				; CLS/GLM: convergence via relative change in residual sum vᵀPv
				; For GLM, xd=0 does NOT mean convergence — the relinearization at l+v
				; can change the system even when parameter corrections are zero.
				; See: Neitzel/Petrovic (2008), rigorous vs. approximate GH evaluation.
				Local $fDeltaR2 = Abs($mState.r2sum - $fR2sumPrev)
				If $fR2sumPrev > 1e-30 And $fDeltaR2 / $fR2sumPrev < $fTolerance Then ExitLoop
				$fR2sumPrev = $mState.r2sum
			Else
				; OLS/WLS/LSE: convergence via parameter correction
				Local $fMaxDx = Abs(_blas_amax($mState.xd))
				If $fMaxDx < $fTolerance Then ExitLoop
				; adaptive stagnation: if corrections no longer halving → numerical accuracy floor reached
				If $fMaxDx > 0.5 * $fMaxDxPrev Then
					$iStagnation += 1
					If $iStagnation >= 3 Then ExitLoop
				Else
					$iStagnation = 0
				EndIf
				$fMaxDxPrev = $fMaxDx
			EndIf
		Next

	EndIf

EndFunc

#EndRegion ; Nonlinear iteration (Tier 2)

#Region Single iteration step (Tier 3)

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_solveIteration
; Description ...: Single linear adjustment step
; Syntax.........: __adj_solveIteration($mSystem, $mState)
; Parameters ....: $mSystem     - [ByRef] Adjustment system map
;                  $mState      - [ByRef] Solver state map
; Return values .: Success      - True (linear solve completed)
;                  Failure      - SetError($ADJ_ERR_*) on solver error
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: Computes Jacobians, contradiction vector, applies whitening, equilibration,
;                  LM damping, dispatches linear solver. Handles GLM/CLS residual accumulation.
; Related .......: __adj_computeJacobians, __adj_computeContradiction, __adj_applyWhitening, __adj_dispatchLinearSolver
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __adj_solveIteration(ByRef $mSystem, ByRef $mState)

	; calculate the jacobian matrices
	__adj_computeJacobians($mSystem, $mState) ; A, B Matrizen

	; Calculate the contradiction vector w (or -w)
	__adj_computeContradiction($mSystem, $mState, True)

	; apply the "Whitening" to convert a weighted problem into an unweighted one
	__adj_applyWhitening($mState)

	; Jacobi-Equilibration: normalize column scales of whitened A
	; Applied BEFORE Marquardt-D computation so D̃ captures true curvature, not scale artifacts
	; Only for OLS/WLS/GLS — GLM and LSE handle equilibration in their own solvers
	If $mState.scaling And StringRegExp($mState.modelType, '(OLS|WLS|GLS)$') Then
		; save pre-equilibration A for residuals and SVD Qxx
		$mState.A_orig = _la_duplicate($mState.Matrix_A)
		; equilibrate: A_*,j /= S_j (only observation rows, not LM augmentation)
		; nFormulas = observation equation rows (correct for both GN and LM:
		; GN: Matrix_A.rows = nFormulas; LM: Matrix_A.rows = nFormulas + nParams)
		__adj_applyEquilibration($mState, $mState.Matrix_A, $mState.nFormulas, $mState.nParams)
	EndIf

	; initialise Levenberg-Marquardt damping factor
	; if LM is chosen:
	Local $mConfig = $mSystem.config
	If $mConfig.algorithm = "LM" And (Not MapExists($mState, "nIterations")) Then
		__adj_initLMDamping($mState)
	ElseIf $mConfig.algorithm = "LM" And MapExists($mState, "fLMLambda") _
		And StringRegExp($mState.modelType, '(OLS|WLS|GLS)$') Then
		; Marquardt D-update for subsequent iterations (OLS/WLS/GLS only — GLM/LSE handle this in solveLinearSystem)
		Local $mA_lm = $mState.Matrix_A
		Local $nP_lm = $mState.nParams
		Local $iNobs_lm = $mA_lm.rows - $nP_lm  ; real observation rows (excl. augmentation)

		; compute D_new = √diag(Aᵀ·A) from whitened Jacobian (only real observation rows)
		Local $mC_lm = _blas_createMatrix($nP_lm, $nP_lm)
		_blas_syrk($mA_lm, $mC_lm, 1, 0, "U", "T", $nP_lm, $iNobs_lm, $mA_lm.rows, $nP_lm)
		Local $mDnew = _blas_copy($mC_lm, 0, $mC_lm.rows + 1, 0, 1, $nP_lm)
		__adj_applyMarquardtFloorAndSqrt($mDnew, $nP_lm)

		; update D = max(D_old, D_new) element-wise
		__adj_updateMarquardtD($mState.LM_D, $mDnew, $nP_lm)

		; overwrite augmentation rows with updated √λ · D
		Local $iOff_lm = $mState.nFormulas + $mState.nRestrictions
		__adj_fillLMaugmentation($mA_lm.ptr, $mState.LM_D, $mState.fLMLambda, $iOff_lm, $mA_lm.rows + 1, $nP_lm)

		; update gradient g = Aᵀ · w (whitened, only first nObs rows)
		Local $mGrad_lm = _blas_createVector($nP_lm)
		_blas_gemv($mA_lm, $mState.Vector_W, $mGrad_lm, 1, 0, "T", 1, 1, $iNobs_lm, $nP_lm, $mA_lm.rows)
		$mState.LM_gradient = $mGrad_lm
	EndIf

	; solve the current linear iteration
	__adj_dispatchLinearSolver($mState)
	If @error Then Return SetError(@error, @extended, False)

	; For GLM: w at l₀ → solver yields TOTAL v. Store directly.
	; For CLS: w at l₀+v (from iter 2) → solver yields INCREMENTAL Δv. Accumulate.
	If StringRegExp($mState.modelType, '(CLS|GLM)$') Then
		If StringRegExp($mState.modelType, 'CLS$') And MapExists($mState, "nIterations") Then
			; CLS iteration ≥ 2: $mState.r contains Δv → accumulate v = v_prev + Δv
			_blas_axpy($mState.r, $mState.r_accumulated) ; r_accumulated += Δv
			$mState.r = _la_duplicate($mState.r_accumulated)
		Else
			; First iteration (CLS/GLM) or GLM: r is total v → store directly
			$mState.r_accumulated = _la_duplicate($mState.r)
		EndIf
		; Recompute r2sum from observation residuals vᵀPv — but ONLY when NOT in LM phase.
		; During LM, r2sum must stay as the DGELSY residual norm in Cholesky space for
		; consistent gain ratio computation (predicted decrease also uses Cholesky space).
		; After LM convergence, the clean GN phase will use vᵀPv for correct statistics.
		If Not MapExists($mState, "fLMLambda") Then
			$mState.r2sum = __adj_computeVtPV($mState)
		EndIf
	EndIf

	$mState.nIterations = MapExists($mState, "nIterations") ? $mState.nIterations + 1 : 1
EndFunc



; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_computeJacobians
; Description ...: Builds Jacobian matrices A and B via numerical/analytical differentiation
; Syntax.........: __adj_computeJacobians($mSystem, $mState [, $sDeriveMethod])
; Parameters ....: $mSystem        - [ByRef] Adjustment system map
;                  $mState         - [ByRef] Solver state map
;                  $sDeriveMethod  - [optional] Differentiation method (Default from config)
; Return values .: Success      - True (Jacobians filled in $mState)
;                  Failure      - None (derivatives always computable)
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: Handles all model types. For CLS/GLM iter>=2: evaluates at l₀+v.
;                  Fills LM augmentation diagonal.
; Related .......: __adj_evalDerivative, __adj_solveIteration
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __adj_computeJacobians(ByRef $mSystem, ByRef $mState, $sDeriveMethod = Default)
	; beachte Modelltyp, LM

	Local $iRFormulas       = $mState.nFormulas, _
	      $iRRestrictions   = $mState.nRestrictions

	If IsKeyword($sDeriveMethod) = 1 Then
		Local $mCfg_jac = $mSystem.config
		If MapExists($mCfg_jac, "deriveMethod") Then
			$sDeriveMethod = $mCfg_jac.deriveMethod
		Else
			$sDeriveMethod = $mState.isNonlinear ? "Forward" : "Central"
		EndIf
	EndIf

	; prepare lists of parameter and observations maps
	Local $mModel   = $mSystem.model
	Local $mObsFull = $mModel.obs, $mObs[], $mPar = $mModel.params
	Local $mIdxObs  = $mState.idxObs

	; for CLS/GLM cases the jacobian matrix B for the observations must be developed at the point o_i + r_i
	If StringRegExp($mState.modelType, '(CLS|GLM)$') And MapExists($mState, "nIterations") Then ; further iterations require improvement of the observed values with the residuals
		; in CLS/GLM cases __adj_computeResiduals() has already been implicitly implemented in the solution
		Local $mR = $mState.r
		For $sObsName In MapKeys($mObsFull)
			$mObs[$sObsName] = $mObsFull[$sObsName].value + __adj_vecGet($mR, $mIdxObs[$sObsName])
		Next
	Else
		For $sObsName In MapKeys($mObsFull)
			$mObs[$sObsName] = $mObsFull[$sObsName].value
		Next
	EndIf

	Local $mJacobiA, $tJacobiA, $mJacobiB, $tJacobiB, _
		$iSkip, $iOffset, $mIdxParam, $mFormula, $idxF, $idxP, $idxEntry, $sFunc, $sSymbol, $mDerivExec
	Switch $mState.modelType
	  	Case  "LSE", "WLSE", "GLSE" ; Linear Equality-Constrained Least Squares: f(x) = y under restrictions g(x) = z
			Local $mRA        = $mState.Matrix_RA, _
				$tRA_struct = $mRA.struct
				$mIdxParam  = $mState.idxParams
				$idxF       = 1

			; clear RA
			_la_clear($mRA)

		  	For $mFormula In $mModel.restrictions
				  $sFunc     = $mFormula.executable
				  $mDerivExec = MapExists($mFormula, "derivativesExec") ? $mFormula.derivativesExec : False

				For $sSymbol in $mFormula.params
					$idxP = $mIdxParam[$sSymbol]
					$idxEntry = $iRRestrictions * $idxP + $idxF
					DllStructSetData($tRA_struct, 1, __adj_evalDerivative($mDerivExec, $sSymbol, $sFunc, $mPar, $mObs, False, $sDeriveMethod), $idxEntry)
				Next

				$idxF += 1
			Next
			ContinueCase ; to handle the parameter jacobian for the formulas

		Case "OLS", "WLS", "GLS"
			$mJacobiA  = $mState.Matrix_A
			$tJacobiA  = $mJacobiA.struct
			$iSkip     = $mJacobiA.rows - $iRFormulas
			$mIdxParam = $mState.idxParams
			$idxF      = 1

			; clear A
			_la_clear($mJacobiA)

			For $mFormula In $mModel.formulas
				$sFunc = $mFormula.executable
				$mDerivExec = MapExists($mFormula, "derivativesExec") ? $mFormula.derivativesExec : False

				For $sSymbol in $mFormula.params
					$idxP     = $mIdxParam[$sSymbol]
					$idxEntry = $iRFormulas * $idxP + $idxF + $idxP * $iSkip
					DllStructSetData($tJacobiA, 1, __adj_evalDerivative($mDerivExec, $sSymbol, $sFunc, $mPar, $mObs, False, $sDeriveMethod), $idxEntry)
				Next

				$idxF += 1
			Next

		Case "GLM", "WGLM", "GGLM" ; General Gauss-Markov Linear Model: f(x,y) = 0
			$mJacobiA = $mState.Matrix_A
			$tJacobiA = $mJacobiA.struct
			$mJacobiB = $mState.Matrix_B
			$tJacobiB = $mJacobiB.struct
			$mIdxParam = $mState.idxParams

			; clear A & B
			_la_clear($mJacobiA)
			_la_clear($mJacobiB)

			; process the jacobians for the functions
			$iSkip = $mJacobiA.rows - $iRFormulas
			$idxF  = 1
			For $mFormula In $mModel.formulas
		  		$sFunc = $mFormula.executable
				$mDerivExec = MapExists($mFormula, "derivativesExec") ? $mFormula.derivativesExec : False

				; write parameter jacobian
		  		For $sSymbol in $mFormula.params
					$idxP = $mIdxParam[$sSymbol]
					$idxEntry = $iRFormulas * $idxP + $idxF + $idxP * $iSkip
					DllStructSetData($tJacobiA, 1, __adj_evalDerivative($mDerivExec, $sSymbol, $sFunc, $mPar, $mObs, False, $sDeriveMethod), $idxEntry)
				Next

				; write observation jacobian
				For $sSymbol in $mFormula.obs
					$idxP = $mIdxObs[$sSymbol]
					$idxEntry = $iRFormulas * $idxP + $idxF + $idxP * $iSkip
					DllStructSetData($tJacobiB, 1, __adj_evalDerivative($mDerivExec, $sSymbol, $sFunc, $mPar, $mObs, True, $sDeriveMethod), $idxEntry)
				Next

				$idxF += 1
			Next

			; process the parameter jacobian for the restrictions into Matrix_RA (null-space projection)
			If $mState.nRestrictions > 0 Then
				$mRA       = $mState.Matrix_RA
				$tRA_struct = $mRA.struct
				_la_clear($mRA)

				$idxF = 1
				For $mFormula In $mModel.restrictions
					$sFunc      = $mFormula.executable
					$mDerivExec = MapExists($mFormula, "derivativesExec") ? $mFormula.derivativesExec : False

					For $sSymbol In $mFormula.params
						$idxP = $mIdxParam[$sSymbol]
						$idxEntry = $iRRestrictions * $idxP + $idxF
						DllStructSetData($tRA_struct, 1, __adj_evalDerivative($mDerivExec, $sSymbol, $sFunc, $mPar, $mObs, False, $sDeriveMethod), $idxEntry)
					Next

					$idxF += 1
				Next
			EndIf

		Case "CLS", "WCLS", "GCLS" ; Conditional Least Squares: g(y) = 0
			$mJacobiB = $mState.Matrix_B
			$tJacobiB = $mJacobiB.struct

			; clear B
			_la_clear($mJacobiB)

			; process the jacobians for the functions
			$iSkip = $mJacobiB.rows - $iRFormulas
			$idxF  = 1
			For $mFormula In $mModel.formulas
				$sFunc = $mFormula.executable
				$mDerivExec = MapExists($mFormula, "derivativesExec") ? $mFormula.derivativesExec : False

				; write observation jacobian
				For $sSymbol in $mFormula.obs
					$idxP = $mIdxObs[$sSymbol]
					$idxEntry = $iRFormulas * $idxP + $idxF + $idxP * $iSkip
					DllStructSetData($tJacobiB, 1, __adj_evalDerivative($mDerivExec, $sSymbol, $sFunc, $mPar, $mObs, True, $sDeriveMethod), $idxEntry)
				Next

				$idxF += 1
	  		Next

			; process the parameter jacobian for the restrictions
			If $mState.nRestrictions > 0 Then
				Local $mRB        = $mState.Matrix_RB, _
				      $tRB_struct = $mRB.struct
				      $idxF       = 1
				For $mFormula In $mModel.restrictions
					$sFunc = $mFormula.executable
					$mDerivExec = MapExists($mFormula, "derivativesExec") ? $mFormula.derivativesExec : False

					For $sSymbol in $mFormula.obs
						$idxP = $mIdxObs[$sSymbol]
						$idxEntry = $iRRestrictions * $idxP + $idxF
						DllStructSetData($tRB_struct, 1, __adj_evalDerivative($mDerivExec, $sSymbol, $sFunc, $mPar, $mObs, True, $sDeriveMethod), $idxEntry)
					Next
					$idxF += 1
				Next
			EndIf
	EndSwitch

	; add Marquardt-scaled diagonal √λ·D for Levenberg-Marquardt augmentation
	; Note: LSE augmentation happens in __adj_dispatchLinearSolver (on reduced Ā, not on A).
	; The rows filled here for LSE are unused (gemm only reads first nObs rows of A).
	If MapExists($mState, "fLMLambda") And StringRegExp($mState.modelType, '(OLS|WLS|GLS)$') Then
		$iOffset = $iRFormulas + $iRRestrictions
		If StringRegExp($mState.modelType, "^[WG]?CLS$") Then
			__blas_fillWithScalar($tJacobiB, Sqrt($mState.fLMLambda), $iOffset, $mJacobiB.rows + 1, $mState.nObs)
		ElseIf MapExists($mState, "LM_D") Then
			; Marquardt: √λ · Dᵢ (D will be updated after whitening in solve_iter)
			__adj_fillLMaugmentation($tJacobiA, $mState.LM_D, $mState.fLMLambda, $iOffset, $mJacobiA.rows + 1, $mState.nParams)
		Else
			; Fallback (should not happen after initializeLMdampingfactor)
			__blas_fillWithScalar($tJacobiA, Sqrt($mState.fLMLambda), $iOffset, $mJacobiA.rows + 1, $mState.nParams)
		EndIf
	EndIf

EndFunc


; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_computeContradiction
; Description ...: Computes contradiction vector w = f(x,l)
; Syntax.........: __adj_computeContradiction($mSystem, $mState [, $bInvert [, $bReturn]])
; Parameters ....: $mSystem     - [ByRef] Adjustment system map
;                  $mState      - [ByRef] Solver state map
;                  $bInvert     - If True, compute -w (default for solver)
;                  $bReturn     - If True, return copy instead of modifying in-place
; Return values .: Success      - If $bReturn: BLAS vector copy, else None
;                  Failure      - None
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: GLM/first: evaluates at l₀. CLS iter>=2: evaluates at l₀+v.
; Related .......: __adj_solveIteration
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __adj_computeContradiction(ByRef $mSystem, ByRef $mState, $bInvert = False, $bReturn = False)
	; prepare lists of parameter and observations maps
	Local $mModel   = $mSystem.model
	Local $mObsFull = $mModel.obs, $mObs[]
	Local $mIdxObs  = $mState.idxObs

	; Evaluation point for the contradiction vector w = -g(…):
	;
	; GLM (has parameters): w evaluated at original observations l₀ but current
	;   parameters x_k.  The parameter updates provide convergence feedback.
	;   Corresponds to total-v formulation (Lenzmann/Lenzmann 2004, Neitzel/Petrovic 2008).
	;
	; CLS (no parameters): w MUST be evaluated at l₀+v (current adjusted obs)
	;   after the first iteration, because without parameter updates there is no
	;   other convergence feedback.  The solver then yields the incremental Δv,
	;   which is accumulated in __adj_solveIteration.
	If StringRegExp($mState.modelType, 'CLS$') And MapExists($mState, "nIterations") Then
		; CLS iterations ≥ 2: evaluate at improved observations l₀ + v
		Local $mR = $mState.r
		For $sObsName In MapKeys($mObsFull)
			$mObs[$sObsName] = $mObsFull[$sObsName].value + __adj_vecGet($mR, $mIdxObs[$sObsName])
		Next
	Else
		; First iteration (any model) and GLM iterations: evaluate at l₀
		For $sObsName In MapKeys($mObsFull)
			$mObs[$sObsName] = $mObsFull[$sObsName].value
		Next
	EndIf
	Local $mPar = $mModel.params
	#forceref $mObs, $mPar

	Local $mW    = $mState.Vector_W,  _
	      $tW    = $mW.struct,         _
	      $iIdxV = 1, $mFormula, $sFunc, $mWR, $tWR

	; set all values to 0 (especially for Levenberg-Marquardt)
	_la_clear($mW)

	If $bReturn Then
		$mW = _blas_duplicate($mW)
		$tW    = $mW.struct
	EndIf

	Switch $mState.modelType
		Case "LSE", "WLSE", "GLSE", "CLS", "WCLS", "GCLS"
			If $mState.nRestrictions > 0 Then
				$mWR = $mState.Vector_WR
		 		$tWR = $mWR.struct

				For $mFormula In $mModel.restrictions
					$sFunc = $mFormula.executable
					DllStructSetData($tWR, 1, $bInvert ? - Execute($sFunc) : Execute($sFunc), $iIdxV)
					$iIdxV += 1
				Next
			EndIf

			ContinueCase ; continue with calculating the values for the formulas

		Case "OLS", "WLS", "GLS", "GLM", "WGLM", "GGLM"
			$iIdxV = 1
			For $mFormula In $mModel.formulas
				$sFunc = $mFormula.executable
				DllStructSetData($tW, 1, $bInvert ? - Execute($sFunc) : Execute($sFunc), $iIdxV)
				$iIdxV += 1
			Next

			; fill restriction contradiction vector (Vector_WR) for GLM with restrictions
			If StringRegExp($mState.modelType, "^[WG]?GLM$") And $mState.nRestrictions > 0 Then
				Local $mWR_glm = $mState.Vector_WR
				Local $tWR_glm = $mWR_glm.struct
				Local $iIdxWR = 1
				For $mFormula In $mModel.restrictions
					$sFunc = $mFormula.executable
					DllStructSetData($tWR_glm, 1, $bInvert ? - Execute($sFunc) : Execute($sFunc), $iIdxWR)
					$iIdxWR += 1
				Next
			EndIf
	EndSwitch

	If $bReturn Then Return $mW
EndFunc


; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_applyWhitening
; Description ...: Transforms weighted/generalized problem into unweighted form
; Syntax.........: __adj_applyWhitening($mState)
; Parameters ....: $mState      - [ByRef] Solver state map
; Return values .: Success      - True (whitening applied in-place)
;                  Failure      - SetError on BLAS error
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: WLS: diag(1/σ)·A, diag(1/σ)·w. GLS: L⁻¹·A, L⁻¹·w.
;                  WGLM/WCLS: B·diag(σ). GGLM/GCLS: B·L.
; Related .......: __adj_solveIteration
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __adj_applyWhitening(ByRef $mState)
	If Not StringRegExp($mState.modelType, '^[WG]') Then Return


	; dimensions
	Local $iM = $mState.nObs, _
	      $mL, $mTmp, $mA, $mB, $mStdDevs
;~ 	Local $mS, $mD, $mV, $mC

	Switch $mState.modelType
		Case "OLS", "LSE", "CLS", "GLM"
			; do nothing - all observations have similar weight

		; weighted case - observations have different a-priori standard deviations but there are no covariances between them
		Case "WLS", "WLSE"
			; - create matrix L out of the standard deviations: L = diag(1/σ₁, 1/σ₂, 1/σ₃, ..., 1/σₙ) (std deviatons - NOT variances!)
			; - A' = L ⋅ A, b' = L ⋅ b
			; - do the normal OLS/LSE but with A' and b' instead of A and b

			$mA = $mState.Matrix_A

			; P · b --> b
			$mL = $mState.Vector_ObsInvStdDev
			$mTmp = _blas_duplicate($mState.Vector_W)
			_blas_sbmv($mL.ptr, $mState.Vector_W, $mTmp.ptr, 1, 0, 0, "L", $iM, 1, 1, 1, "DOUBLE")
			If @error Then Return SetError(@error + 70, @extended, Null)
			$mState.Vector_W = $mTmp

			; L · A --> A
			_blas_trmm(_la_VectorToDiag($mL), $mA.ptr, 1, "L", "U", "N", "N", $iM, $mA.cols, $iM, $mA.rows, "DOUBLE")
			If @error Then Return SetError(@error + 80, @extended, Null)
			; alternative approach for L · A: Loop over L-elements and scale corresponding row of A with this value (with _blas_scal()). But this needs a AutoIt loop so i choose the other way.

		Case "WGLM", "WCLS" ; CLS/WCLS/GCLS is only handled as a special case of GLM/WGLM/GGLM
			; Approach: "Whitening": The Cholesky decomposition P=L⋅Lᵀ is calculated from the P matrix. Then B'=B⋅L⁻ᵀ is calculated.
			;            dggglm is then calculated further with this B'. This corresponds to the weighted form.
			;            for weighted case where Σₗₗ is a diagonal matrix this simply leads to:  L⁻¹ = diag(σ₁, σ₂, σ₃, ..., σₙ) = vector of apriori standard deviations for the observations

			; P=Σₗₗ⁻¹ --> P=L⋅Lᵀ --> L⁻¹ leads to:
			$mStdDevs = $mState.Vector_ObsStdDev
			$mB = $mState.Matrix_B

			; B ⋅ L⁻ᵀ
			$mL = _la_VectorToDiag($mStdDevs)
			_blas_trmm($mL.ptr, $mB.ptr, 1, "R", "U", "N", "T", $mB.rows, $mB.cols, $mL.rows, $mB.rows)
			If @error Then Return SetError(@error + 80, @extended, Null)

			$mState.Matrix_B = $mB

		;~ Case "WCLS"
			;~ ; - Because in the CLS model a function can be dependent on several observations (f(x,y) = 0), 
			;~ ;   we must first determine aggregated standard deviations for the functions using the error propagation law.
			;~ ; 	- the formula for this looks like: σ²ₕ = (∂ₕ/∂l₁)² ⋅ σ₁² + (∂ₕ/∂l₂)² ⋅ σ₂² + ... + (∂ₕ/∂lₙ)² ⋅ σₙ²
			;~ ; 	- we can use existing parts for this:
			;~ ;     P = (B ⋅ Σₗₗ ⋅ Bᵀ)⁻¹
			;~ ;     for weighted case this means: P = (B ⋅ S_B2 ⋅ Bᵀ)⁻¹ with S_B2 = diag(σ₁², σ₂², σ₃², ..., σₙ²)
			;~ ;     hint: the inverse for weighted case is nothing as the 
			;~ ; - do a cholesky 

			;~ $iF = $mSystem.nFormulas
			;~ $iO = $mSystem.nObs
			;~ $iR = $mSystem.nRestrictions
			;~ $mB = $mSystem.Matrix_B

			;~ ; D = diag(σ₁, σ₂, σ₃, ..., σₙ)
			;~ $mS = $mSystem.Vector_ObsStdDev
			;~ $mD = _la_VectorToDiag($mS, False)
			
			;~ ; V = B ⋅ D
			;~ $mV = _la_mul($mB, $mD)

			;~ ; C = V ⋅ Vᵀ
			;~ $mC = _blas_createMatrix($mV.rows, $mV.rows)
			;~ _blas_syrk($mV, $mC)
			;~ $mC = _la_getDiag($mC)

			;~ ; P = diag(1/C₁₁, 1/C₂₂, 1/C₃₃, ..., 1/Cₙₙ)
			;~ $mP = _la_invElements($mC)

			;~ ; workaround to prevent weighting of restrictions (Todo: find a way to directly scale formulas only)
			;~ If $iR > 0 Then __blas_fillWithScalar($mP, 1.0, $iF, 1, $iR)

			;~ ; P * d --> d
			;~ $mTmp = _blas_duplicate($mSystem.Vector_W)
			;~ _blas_sbmv($mP.ptr, $mSystem.Vector_W, $mTmp.ptr, 1, 0, 0, "L", $iF, 1, 1, 1, "DOUBLE")
			;~ If @error Then Return SetError(@error + 70, @extended, Null)
			;~ ;~ _la_scale($mTmp, -1, True) 
			;~ $mSystem.Vector_W = $mTmp

			;~ ; alternative approach for P · A and P · B: Loop over L-elements and scale corresponding row of A with this value (with _blas_scal()).
			;~ ; But this needs a AutoIt loop so i choose the following way.

			;~ ; P-vector to P- diagonal matrix
			;~ _la_VectorToDiag($mP, True)

			;~ ; P * B --> B
			;~ _blas_trmm($mP.ptr, $mB.ptr, 1, "L", "U", "N", "N", $iF + $iR, $iO, $iF + $iR, $iF + $iR)
			;~ If @error Then Return SetError(@error + 80, @extended, Null)
			;~ $mSystem.Matrix_B = $mB

		; general case - observations have different a-priori standard deviations and there exist covariances between them
		Case "GLS", "GLSE"
			; Whitening via Cholesky: Σₗₗ = L·Lᵀ (cached in $mState.CovCholeskyL)
			; Solve L·A' = A  →  A' = L⁻¹·A   (DTRSM, SIDE=L)
			; Solve L·w' = w  →  w' = L⁻¹·w   (DTRSV)
			$mA = $mState.Matrix_A
			Local $mCholL = $mState.CovCholeskyL

			; w' = L⁻¹ · w  (only first nObs elements, LM-augmentation rows stay untouched)
			$mTmp = _blas_duplicate($mState.Vector_W)
			_blas_trsv($mCholL, $mTmp, "L", "N", "N", $iM, 1, $iM)
			If @error Then Return SetError(@error + 70, @extended, Null)
			$mState.Vector_W = $mTmp

			; A' = L⁻¹ · A  (only first nObs rows; LDA = A.rows to skip LM rows)
			_blas_trsm($mCholL, $mA.ptr, 1.0, "L", "L", "N", "N", $iM, $mA.cols, $iM, $mA.rows)
			If @error Then Return SetError(@error + 80, @extended, Null)

		Case "GGLM", "GCLS"
			; Whitening: Cholesky Σₗₗ = L·Lᵀ (L cached in $mState.CovCholeskyL)
			; B' = B · L  (TRMM, SIDE=R) — analogous to WGLM: B' = B·diag(σ)
			; This transforms min ‖ṽ‖₂² into min vᵀPv where ṽ = L⁻¹·v
			$mB = $mState.Matrix_B
			Local $mCholL = $mState.CovCholeskyL

			; B' = B · L: TRMM with SIDE=R, UPLO=L, TRANS=N
			; B has dimension (nFormulas+nRestrictions) × nObs, L is nObs × nObs
			_blas_trmm($mCholL, $mB.ptr, 1.0, "R", "L", "N", "N", $mB.rows, $iM, $iM, $mB.rows)
			If @error Then Return SetError(@error + 80, @extended, Null)
			$mState.Matrix_B = $mB
	EndSwitch

EndFunc

#EndRegion ; Single iteration step (Tier 3)

#Region Linear solvers

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_dispatchLinearSolver
; Description ...: Dispatches to model-specific solver (OLS/LSE/CLS/GLM)
; Syntax.........: __adj_dispatchLinearSolver($mState)
; Parameters ....: $mState      - [ByRef] Solver state map
; Return values .: Success      - True
;                  Failure      - SetError on solver failure
; Author ........: AspirinJunkie
; ===============================================================================================================================
Func __adj_dispatchLinearSolver(ByRef $mState)
	; save previous vᵀ·P·v
	$mState.r2sum0 = $mState.r2sum > 0 ? $mState.r2sum : _blas_nrm2($mState.Vector_W)^2

	Switch $mState.modelType
		Case "OLS", "WLS", "GLS"
			__adj_solveOLS($mState)
		Case "LSE", "WLSE", "GLSE"
			__adj_solveLSE($mState)
		Case "CLS", "WCLS", "GCLS"
			__adj_solveCLS($mState)
		Case "GLM", "WGLM", "GGLM"
			__adj_solveGLM($mState)
	EndSwitch
	If @error Then Return SetError(@error, @extended, False)
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_solveOLS
; Description ...: OLS/WLS/GLS solver via DGELSY or DGELSD
; Syntax.........: __adj_solveOLS($mState)
; Parameters ....: $mState      - [ByRef] Solver state map
; Return values .: Success      - True (xd and r2sum stored in $mState)
;                  Failure      - SetError($ADJ_ERR_SOLVER) on LAPACK error
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: Extracts xd and r2sum. Handles equilibration back-transform and LM step caching.
; Related .......: __adj_dispatchLinearSolver, __adj_reverseEquilibration
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __adj_solveOLS(ByRef $mState)
	Local $iN = $mState.nParams, $mXd, $mW

	; Ordinary Least Squares: f(x) = y
	; Adjustment by indirect observations
	; DGELSY: QR with column pivoting and rank detection

	; A_orig is saved in __adj_solveIteration before equilibration (pre-equilibration copy)
	; For non-scaling mode, save here as fallback
	If Not MapExists($mState, "A_orig") Then $mState.A_orig = _la_duplicate($mState.Matrix_A)

	; solve system
	Local $iLpErr, $iRank
	If $mState.solver = $ADJ_SOLVER_SVD Then
		Local $fRCOND = $mState.solverRCOND
		; Use Default dimensions — _lp_gelsd reads actual rows/cols from matrix map
		; (critical for LM: Matrix_A includes augmentation rows beyond nObs)
		Local $tSV = _lp_gelsd($mState.Matrix_A, $mState.Vector_W, "N", Default, Default, Default, Default, Default, $fRCOND)
		$iLpErr = @error
		$iRank = @extended
		If $iLpErr Then Return SetError($ADJ_ERR_SOLVER, $iLpErr, False)
		$mState.tSingularValues = $tSV
	Else
		Local $tJPVT = _lp_gelsy($mState.Matrix_A, $mState.Vector_W)
		$iLpErr = @error
		$iRank = @extended
		If $iLpErr Then Return SetError($ADJ_ERR_SOLVER, $iLpErr, False)
		$mState.JPVT = $tJPVT
	EndIf

	; rank deficiency warning (continue but flag for statistics)
	If $iRank < $iN Then $mState.rankDeficient = True

	$mW = $mState.Vector_W

	; ausgeglichene ParameterDIFFERENZEN(!) xd
	$mXd = _blas_createVector($iN)
	_lp_lacpy($mW, $mXd, "X", $iN, 1)
	$mState.xd = $mXd

	; cache dx̃ in equilibrated space for LM gain ratio BEFORE back-transform
	; (LM_gradient and LM_D are in equilibrated space — LM_step must match)
	If MapExists($mState, "fLMLambda") And MapExists($mState, "EquilibrationScale") Then
		$mState.LM_step = _la_duplicate($mXd)
	EndIf

	; reverse equilibration: dx = S⁻¹ · dx̃ (back to original parameter space)
	__adj_reverseEquilibration($mState, $mState.xd, $iN)

	; vᵀ·P·v (residuals in positions N+1 to M of the solution vector)
	$mState.r2sum = _blas_nrm2($mW, $iN, 1, $mW.size - $iN)^2
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_solveLSE
; Description ...: LSE/WLSE/GLSE solver via null-space projection + DGELSY/DGELSD
; Syntax.........: __adj_solveLSE($mState)
; Parameters ....: $mState      - [ByRef] Solver state map
; Return values .: Success      - True (xd and r2sum stored in $mState)
;                  Failure      - SetError($ADJ_ERR_SOLVER) on LAPACK error
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: QR of RAᵀ -> Q₂ projection -> reduced solve -> back-transform
;                  Δx = α·x_p + Q₂·z. Includes LM augmentation in reduced space.
; Related .......: __adj_dispatchLinearSolver
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __adj_solveLSE(ByRef $mState)
	Local $iM = $mState.nObs, $iN = $mState.nParams, $iR = $mState.nRestrictions, $mXd

	; Linear Equality-Constrained Least Squares: f(x) = y under restrictions g(x) = z
	; save original matrices before solver overwrites them (needed for Qxx, redundancy matrix)
	$mState.A_orig  = _la_duplicate($mState.Matrix_A)
	$mState.RA_orig = _la_duplicate($mState.Matrix_RA)

	; ══════════════════════════════════════════════════════════════
	; Nullraum-Projektion + DGELSY
	; ══════════════════════════════════════════════════════════════
	; min ‖Ax - w‖₂  subject to RA·x = WR
	;
	; 1. QR of RAᵀ → Q = [Q₁|Q₂], R_c
	; 2. Particular solution: x_p = Q₁ · R_c⁻ᵀ · WR
	; 3. Project: Ā = A·Q₂, w̄ = w - A·x_p
	; 4. DGELSY(Ā, w̄) → z
	; 5. Δx = x_p + Q₂·z

	; guard: restrictions must not fully determine all parameters
	If $iR >= $iN Then Return SetError($ADJ_ERR_RESTR_OVERDETERMINED, $iR, False)

	; Step 1: QR factorization of RAᵀ (n×p matrix)
	Local $mRAt = _la_transpose($mState.Matrix_RA)  ; RAᵀ: n×p
	Local $tTau = _lp_geqrf($mRAt, $iN, $iR)
	If @error Then Return SetError($ADJ_ERR_SOLVER, @error, False)

	; Extract R_c (upper p×p triangular) from the QR result before orgqr overwrites it
	Local $mRc = _blas_createMatrix($iR, $iR)
	_lp_lacpy($mRAt, $mRc, "U", $iR, $iR)

	; Reconstruct Q from QR factorization (n×n orthogonal matrix)
	; orgqr needs n×n storage, but mRAt is only n×p → copy into larger buffer
	Local $mQ = _blas_createMatrix($iN, $iN)
	_lp_lacpy($mRAt, $mQ, "X", $iN, $iR, $iN, $iN)  ; copy n×p data into n×n matrix
	_lp_orgqr($mQ, $tTau, $iN, $iN, $iR)
	If @error Then Return SetError($ADJ_ERR_SOLVER, @error, False)
	; $mQ now contains Q (n×n), first p columns = Q₁, remaining = Q₂

	; Step 2: Particular solution x_p = Q₁ · R_c⁻ᵀ · WR
	; Solve R_cᵀ · t = WR (lower triangular solve since R_c is upper triangular)
	Local $mTemp = _la_duplicate($mState.Vector_WR)
	_blas_trsv($mRc, $mTemp, "U", "N", "T", $iR)  ; R_cᵀ · t = WR → t in $mTemp
	If @error Then Return SetError($ADJ_ERR_SOLVER, @error, False)

	; x_p = Q₁ · t (Q₁ = first p columns of Q)
	Local $mXp = _blas_createVector($iN)
	_blas_gemv($mQ, $mTemp, $mXp, 1, 0, "N", 1, 1, $iN, $iR)  ; Q₁(n×p) · t(p) → x_p(n)

	; Step 3: Project design matrix and contradiction vector into nullspace of RA
	; Extract Q₂ (n×(n-p)) = columns p+1..n of Q (column-major: offset = n*p elements)
	Local $iNfree = $iN - $iR  ; dimension of free parameter space
	Local $iElemSize = ($mQ.datatype = "DOUBLE") ? 8 : 4
	Local $mQ2 = _blas_createMatrix($iN, $iNfree)
	_lp_lacpy($mQ.ptr + $iN * $iR * $iElemSize, $mQ2, "X", $iN, $iNfree, $iN)

	; Ā = A · Q₂ (m×(n-p))
	Local $mAbar = _blas_createMatrix($iM, $iNfree)
	_blas_gemm($mState.Matrix_A, $mQ2, $mAbar, 1, 0, "N", "N", $iM, $iNfree, $iN)

	; Jacobi-Equilibration of projected Ā (nFree columns)
	; Save pre-equilibration Ā for SVD Qxx, then equilibrate
	If $mState.scaling Then
		If $mState.solver = $ADJ_SOLVER_SVD Then $mState.Abar_orig = _la_duplicate($mAbar)
		__adj_applyEquilibration($mState, $mAbar, $iM, $iNfree)
	EndIf

	; w̄ = w - A · x_p
	Local $mWbar = _la_duplicate($mState.Vector_W)
	_blas_gemv($mState.Matrix_A, $mXp, $mWbar, -1, 1, "N", 1, 1, $iM, $iN)  ; w̄ = w + (-1)·A·x_p

	; Step 3.5: Marquardt D-update and LM augmentation — append √λ·D rows to Ā
	; Damping is applied to the free parameters z in the nullspace of RA,
	; which is equivalent to damping ‖x - x_p‖ (the unconstrained part of x).
	Local $iM_data = $iM  ; save original data row count for r2sum
	If MapExists($mState, "fLMLambda") Then
		Local $fLambda_lse = $mState.fLMLambda

		; compute D_new from Āᵀ·Ā in reduced space (nFree × nFree)
		Local $mC_lse = _blas_createMatrix($iNfree, $iNfree)
		_blas_syrk($mAbar, $mC_lse, 1, 0, "U", "T", $iNfree, $iM_data, $iM_data, $iNfree)
		Local $mDnew_lse = _blas_copy($mC_lse, 0, $mC_lse.rows + 1, 0, 1, $iNfree)
		__adj_applyMarquardtFloorAndSqrt($mDnew_lse, $iNfree)

		; update D = max(D_old, D_new) — D lives in reduced space for LSE
		If MapExists($mState, "LM_D") Then
			__adj_updateMarquardtD($mState.LM_D, $mDnew_lse, $iNfree)
		Else
			$mState.LM_D = $mDnew_lse
		EndIf

		; gradient g = Āᵀ · w̄ in reduced space (before augmentation/DGELSY)
		Local $mGrad_lse = _blas_createVector($iNfree)
		_blas_gemv($mAbar, $mWbar, $mGrad_lse, 1, 0, "T", 1, 1, $iM_data, $iNfree, $iM_data)
		$mState.LM_gradient = $mGrad_lse

		; augment Ā: copy into larger matrix, then fill √λ·D diagonal (Marquardt)
		$iM = $iM_data + $iNfree
		Local $mAbar_aug = _blas_createMatrix($iM, $iNfree)
		_lp_lacpy($mAbar, $mAbar_aug, "X", $iM_data, $iNfree, $iM_data, $iM)
		__adj_fillLMaugmentation($mAbar_aug.ptr, $mState.LM_D, $fLambda_lse, $iM_data, $iM + 1, $iNfree)
		$mAbar = $mAbar_aug

		; augment w̄: extend with zeros (already zero-initialized by _blas_createVector)
		Local $mWbar_aug = _blas_createVector($iM)
		_lp_lacpy($mWbar, $mWbar_aug, "X", $iM_data, 1)
		$mWbar = $mWbar_aug
	EndIf

	; Step 4: Solve reduced unconstrained problem via DGELSY
	; The LDB parameter must be max(m, n-p) for DGELSY
	Local $iLDB_LSE = $iM > $iNfree ? $iM : $iNfree
	If $mWbar.size < $iLDB_LSE Then
		; extend w̄ buffer to accommodate solution if n-p > m
		Local $mWbarExt = _blas_createVector($iLDB_LSE)
		_lp_lacpy($mWbar, $mWbarExt, "X", $iM, 1)
		$mWbar = $mWbarExt
	EndIf

	Local $iLpErr_LSE, $iRank_LSE
	If $mState.solver = $ADJ_SOLVER_SVD Then
		Local $fRCOND_LSE = $mState.solverRCOND
		Local $tSV_LSE = _lp_gelsd($mAbar, $mWbar, "N", Default, Default, Default, Default, Default, $fRCOND_LSE)
		$iLpErr_LSE = @error
		$iRank_LSE = @extended
		If $iLpErr_LSE Then Return SetError($ADJ_ERR_SOLVER, $iLpErr_LSE, False)
		$mState.tSingularValues = $tSV_LSE
	Else
		Local $tJPVT_LSE = _lp_gelsy($mAbar, $mWbar, 1, $iM, $iNfree, $iM, $iLDB_LSE)
		$iLpErr_LSE = @error
		$iRank_LSE = @extended
		If $iLpErr_LSE Then Return SetError($ADJ_ERR_SOLVER, $iLpErr_LSE, False)
		$mState.JPVT = $tJPVT_LSE
	EndIf

	; rank deficiency warning (continue but flag for statistics)
	If $iRank_LSE < $iNfree Then $mState.rankDeficient = True
	$mState.Q2 = $mQ2
	$mState.Abar = $mAbar

	; Step 5: Back-transform to original parameter space
	; z = first (n-p) elements of w̄ (DGELSY solution)
	; Δx = α·x_p + Q₂·z   where α = 1/(1+λ) when LM active, else 1
	; LM damping on z (via augmented Ā) only affects the free parameter space.
	; For models where z ≈ 0 (e.g. identity f(x)=x), the entire correction is x_p.
	; Scaling x_p by α = 1/(1+λ) ensures LM also damps the particular solution.
	Local $fAlpha_xp = MapExists($mState, "fLMLambda") ? (1.0 / (1.0 + $mState.fLMLambda)) : 1.0

	; cache z (reduced step) for Marquardt gain ratio (before back-transformation)
	If MapExists($mState, "fLMLambda") Then
		$mState.LM_step = _blas_copy($mWbar, 0, 1, 0, 1, $iNfree)
	EndIf

	; reverse equilibration on z̃ (reduced solution) BEFORE Q₂ back-transform
	; (LM_step already cached above in equilibrated space for gain ratio)
	If MapExists($mState, "EquilibrationScale") Then
		__adj_reverseEquilibration($mState, $mWbar, $iNfree)
	EndIf

	$mXd = _la_duplicate($mXp)
	_blas_scal($mXd, $fAlpha_xp)  ; x_p → α·x_p
	_blas_gemv($mQ2, $mWbar, $mXd, 1, 1, "N", 1, 1, $iN, $iNfree)  ; Δx = α·x_p + Q₂·z
	$mState.xd = $mXd

	; vᵀ·P·v (residuals in positions (n-p)+1 to m_data of the DGELSY solution vector)
	; Use iM_data (not iM) to exclude LM augmentation rows from the residual norm
	$mState.r2sum = _blas_nrm2($mWbar, $iNfree, 1, $iM_data - $iNfree)^2
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_solveCLS
; Description ...: CLS/WCLS/GCLS: minimum-norm B·v = w via DGELSY/DGELSD
; Syntax.........: __adj_solveCLS($mState)
; Parameters ....: $mState      - [ByRef] Solver state map
; Return values .: Success      - True (residuals and r2sum stored in $mState)
;                  Failure      - SetError($ADJ_ERR_SOLVER) on LAPACK error
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: Back-transforms from whitened space. xd = Null (no parameters).
; Related .......: __adj_dispatchLinearSolver
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __adj_solveCLS(ByRef $mState)
	Local $mW

	; Conditional Least Squares: g(y) = 0
	; ══════════════════════════════════════════════════════════════
	; CLS: underdetermined minimum-norm problem B·v = w, min ‖v‖₂
	; ══════════════════════════════════════════════════════════════
	; For CLS there are no parameters (nParams=0). The condition equations
	; g(y)=0 are linearized to B·v = -w (sign already handled by __adj_computeContradiction).
	; DGELSY finds the minimum-norm solution v with rank detection.

	Local $mB_cls = $mState.Matrix_B  ; m_eq × p_obs (whitened for WCLS)
	Local $iMeq   = $mB_cls.rows       ; number of condition equations
	Local $iPobs  = $mB_cls.cols        ; number of observations

	; save original B before solver overwrites it (needed for statistics in Phase 5)
	$mState.B_orig = _la_duplicate($mB_cls)

	; DGELSY requires RHS vector to have max(m,n) rows
	Local $iLDB_CLS = $iMeq > $iPobs ? $iMeq : $iPobs
	$mW = $mState.Vector_W
	Local $mRHS_CLS
	If $mW.size < $iLDB_CLS Then
		$mRHS_CLS = _blas_createVector($iLDB_CLS)
		_lp_lacpy($mW, $mRHS_CLS, "X", $iMeq, 1)
	Else
		$mRHS_CLS = _la_duplicate($mW)
	EndIf

	; solve B·v = w with rank-revealing solver (consistent with OLS/LSE/GLM solvers)
	Local $mBcopy_CLS = _la_duplicate($mB_cls)
	Local $iLpErr_CLS, $iRank_CLS
	If $mState.solver = $ADJ_SOLVER_SVD Then
		Local $fRCOND_CLS = $mState.solverRCOND
		Local $tSV_CLS = _lp_gelsd($mBcopy_CLS, $mRHS_CLS, "N", Default, Default, Default, Default, Default, $fRCOND_CLS)
		$iLpErr_CLS = @error
		$iRank_CLS = @extended
		If $iLpErr_CLS Then Return SetError($ADJ_ERR_SOLVER, $iLpErr_CLS, False)
		$mState.tSingularValues = $tSV_CLS
	Else
		Local $tJPVT_CLS = _lp_gelsy($mBcopy_CLS, $mRHS_CLS, 1, $iMeq, $iPobs, $iMeq, $iLDB_CLS)
		$iLpErr_CLS = @error
		$iRank_CLS = @extended
		If $iLpErr_CLS Then Return SetError($ADJ_ERR_SOLVER, $iLpErr_CLS, False)
		$mState.JPVT = $tJPVT_CLS
	EndIf

	; rank deficiency warning (continue but flag for statistics)
	If $iRank_CLS < $iMeq Then $mState.rankDeficient = True

	; extract solution v (first iPobs entries of the DGELSY output)
	Local $mV_CLS = _blas_createVector($iPobs)
	_lp_lacpy($mRHS_CLS, $mV_CLS, "X", $iPobs, 1)

	; back-transform for weighted case: v = σ · ṽ (undo whitening B̃ = B·diag(σ))
	If StringRegExp($mState.modelType, "^W") Then
		$mState.r = _la_mul($mState.Vector_ObsStdDev, $mV_CLS, False)
	Else
		$mState.r = $mV_CLS
	EndIf

	; r2sum: overwritten by accumulation in __adj_solveIteration, but set initial value
	$mState.r2sum = _blas_nrm2($mV_CLS)^2
	$mState.xd = Null  ; CLS has no parameters
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_solveGLM
; Description ...: GLM/WGLM/GGLM solver via Cholesky of M = B̃B̃ᵀ + DGELSY/DGELSD
; Syntax.........: __adj_solveGLM($mState)
; Parameters ....: $mState      - [ByRef] Solver state map
; Return values .: Success      - True (xd, residuals and r2sum stored in $mState)
;                  Failure      - SetError($ADJ_ERR_*) on LAPACK/Cholesky error
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: Includes null-space for restrictions, LM augmentation, equilibration.
;                  Residuals via ṽ = B̃ᵀ·M⁻¹·(w-A·Δx).
; Related .......: __adj_dispatchLinearSolver
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __adj_solveGLM(ByRef $mState)
	Local $mXd, $mW

	; General Gauss-Markov Linear Model: f(x,y) = 0
	; ══════════════════════════════════════════════════════════════
	; Cholesky-Transformation + DGELSY
	; ══════════════════════════════════════════════════════════════
	; Transforms GLM (w = AΔx + B̃ṽ, min ‖ṽ‖₂) into standard OLS
	; via Cholesky of M = B̃B̃ᵀ, then solves with DGELSY.
	;
	; 1. M = B̃ · B̃ᵀ                                  ← _blas_syrk
	; 2. Cholesky: M = L · Lᵀ                         ← _lp_potrf
	; 3. Ã = L⁻¹ · A                                  ← _blas_trsm
	; 4. w̃ = L⁻¹ · w                                  ← _blas_trsv
	; 5. LM augmentation: append √λ·I rows to Ã, zeros to w̃
	; 6. DGELSY(Ã, w̃) → Δx                            ← _lp_gelsy
	; 7. vᵀPv from DGELSY residual
	; 8. ṽ = B̃ᵀ · M⁻¹ · (w - A·Δx), then v = σ·ṽ    ← residuals

	Local $mA = $mState.Matrix_A  ; m_eq × n_par (whitened)
	Local $mB = $mState.Matrix_B  ; m_eq × p_obs (whitened)
	Local $iMeq  = $mA.rows        ; number of equations
	Local $iNpar = $mA.cols         ; number of parameters
	Local $iPobs = $mB.cols         ; number of stochastic variables (observations + restrictions)

	; save original matrices before solver overwrites them (needed for Qxx, redundancy matrix, residuals)
	$mState.A_orig = _la_duplicate($mA)
	$mState.B_orig = _la_duplicate($mB)
	Local $mW_saved = _la_duplicate($mState.Vector_W)

	; Steps 1+2: Cholesky factorization of M = B̃ · B̃ᵀ
	; Check for cached Cholesky factor from __adj_initLMDamping (first LM iteration)
	Local $mM, $fANORM_M
	If MapExists($mState, "CholeskyL") Then
		$mM       = $mState.CholeskyL
		$fANORM_M = $mState.CholeskyANORM
		MapRemove($mState, "CholeskyL")
		MapRemove($mState, "CholeskyANORM")
	Else
		; Step 1: M = B̃ · B̃ᵀ (m_eq × m_eq symmetric positive definite)
		$mM = _blas_createMatrix($iMeq, $iMeq)
		_blas_syrk($mB, $mM, 1, 0, "L", "N", $iMeq, $iPobs)

		; compute 1-norm of M before Cholesky overwrites it (needed for condition estimation)
		$fANORM_M = _lp_lansy($mM, "1", "L", $iMeq)

		; Step 2: Cholesky M = L · Lᵀ
		_lp_potrf($mM, "L")
		If @error Then Return SetError($ADJ_ERR_NOT_POS_DEF, @error, False)
	EndIf

	; Step 2b: Condition check — κ(M) > 10⁸ indicates ill-conditioned B̃
	Local $fRCOND = _lp_pocon($mM, $fANORM_M, "L", $iMeq)
	If Not @error And $fRCOND > 0 And (1.0 / $fRCOND) > 1e8 Then Return SetError($ADJ_ERR_SOLVER, Int(1.0 / $fRCOND), False)

	; Step 3: Ã = L⁻¹ · A (forward substitution on copy)
	Local $mAtilde = _la_duplicate($mA)
	_blas_trsm($mM, $mAtilde, 1, "L", "L", "N", "N", $iMeq, $iNpar)

	; Save original Ã before DGELSY/LM-augmentation overwrites it (needed for Phase 5 Qvv)
	$mState.Atilde_orig = _la_duplicate($mAtilde)

	; Step 4: w̃ = L⁻¹ · w (forward substitution on copy)
	$mW = _la_duplicate($mState.Vector_W)
	_blas_trsv($mM, $mW, "L", "N", "N", $iMeq)

	; Override r2sum0 with ‖w̃‖² = wᵀM⁻¹w for correct cost baseline.
	; Must be done EVERY iteration (not just first): after relinearization,
	; the cost at Δx=0 is ‖w̃_current‖², which differs from the previous DGELSY residual.
	; This ensures gain ratio consistency: fActual = ‖w̃‖² - ‖w̃-Ã·Δx‖² vs fPred from gradient.
	$mState.r2sum0 = _blas_nrm2($mW, 0, 1, $iMeq)^2

	; ══════════════════════════════════════════════════════════════
	; Step 4b: Null-space projection for restrictions (analog to LSE)
	; ══════════════════════════════════════════════════════════════
	; min ‖Ã·Δx - w̃‖₂  subject to C·Δx = WR
	; where C = Matrix_RA, WR = Vector_WR (contains -g(x))
	;
	; 1. QR of Cᵀ → Q = [Q₁|Q₂], R_c
	; 2. Particular solution: x_p = Q₁ · R_c⁻ᵀ · WR
	; 3. Project: Ā = Ã·Q₂, w̄ = w̃ - Ã·x_p
	; 4. DGELSY(Ā, w̄) → y
	; 5. Δx = α·x_p + Q₂·y  (α = 1/(1+λ) when LM, else 1)
	Local $iNpar_orig = $iNpar  ; save original nParams for back-transform and residuals
	Local $bHasRestrictions = ($mState.nRestrictions > 0)
	Local $mQ2_glm, $mXp_glm

	If $bHasRestrictions Then
		Local $iR_glm = $mState.nRestrictions
		Local $iNfree_glm = $iNpar - $iR_glm

		; guard: restrictions must not fully determine all parameters
		If $iNfree_glm <= 0 Then Return SetError($ADJ_ERR_RESTR_OVERDETERMINED, $iR_glm, False)

		; Step 4b.1: QR factorization of Cᵀ (nParams × nRestrictions)
		Local $mRAt_glm = _la_transpose($mState.Matrix_RA)  ; Cᵀ: nParams × nRestrictions
		Local $tTau_glm = _lp_geqrf($mRAt_glm, $iNpar, $iR_glm)
		If @error Then Return SetError($ADJ_ERR_SOLVER, @error, False)

		; Extract R_c (upper iR × iR triangular) — LDA of source is iNpar (rows of Cᵀ), not iR
		Local $mRc_glm = _blas_createMatrix($iR_glm, $iR_glm)
		_lp_lacpy($mRAt_glm, $mRc_glm, "U", $iR_glm, $iR_glm, $iNpar)

		; Reconstruct Q (nParams × nParams orthogonal)
		Local $mQ_glm = _blas_createMatrix($iNpar, $iNpar)
		_lp_lacpy($mRAt_glm, $mQ_glm, "X", $iNpar, $iR_glm, $iNpar, $iNpar)
		_lp_orgqr($mQ_glm, $tTau_glm, $iNpar, $iNpar, $iR_glm)
		If @error Then Return SetError($ADJ_ERR_SOLVER, @error, False)

		; Step 4b.2: Particular solution x_p = Q₁ · R_c⁻ᵀ · WR
		Local $mTemp_glm = _la_duplicate($mState.Vector_WR)
		_blas_trsv($mRc_glm, $mTemp_glm, "U", "N", "T", $iR_glm)
		If @error Then Return SetError($ADJ_ERR_SOLVER, @error, False)

		$mXp_glm = _blas_createVector($iNpar)
		_blas_gemv($mQ_glm, $mTemp_glm, $mXp_glm, 1, 0, "N", 1, 1, $iNpar, $iR_glm)

		; Step 4b.3: Extract Q₂ (nParams × nFree) = columns iR+1..nParams of Q
		Local $iElemSize_glm = ($mQ_glm.datatype = "DOUBLE") ? 8 : 4
		$mQ2_glm = _blas_createMatrix($iNpar, $iNfree_glm)
		_lp_lacpy($mQ_glm.ptr + $iNpar * $iR_glm * $iElemSize_glm, $mQ2_glm, "X", $iNpar, $iNfree_glm, $iNpar)

		; Ā = Ã · Q₂ (iMeq × nFree)
		Local $mAbar_glm = _blas_createMatrix($iMeq, $iNfree_glm)
		_blas_gemm($mAtilde, $mQ2_glm, $mAbar_glm, 1, 0, "N", "N", $iMeq, $iNfree_glm, $iNpar)

		; w̄ = w̃ - Ã · x_p
		_blas_gemv($mAtilde, $mXp_glm, $mW, -1, 1, "N", 1, 1, $iMeq, $iNpar)

		; Override r2sum0 with ‖w̄‖² (cost baseline in projected space)
		$mState.r2sum0 = _blas_nrm2($mW, 0, 1, $iMeq)^2

		; Replace Ã with Ā and nParams with nFree for subsequent LM/DGELSY
		$mAtilde = $mAbar_glm
		$iNpar = $iNfree_glm  ; DGELSY now operates on nFree columns

		; Store Q₂ for Qxx back-transform (Phase 5)
		$mState.Q2 = $mQ2_glm

		; save Ā before solver overwrites it (needed for SVD Qxx)
		If $mState.solver = $ADJ_SOLVER_SVD Then
			$mState.Abar_orig = _la_duplicate($mAtilde)
		EndIf
	EndIf

	; Jacobi-Equilibration of the final working matrix (Ã or Ā when restrictions present)
	; At this point: $mAtilde = Ã (no restr.) or Ā = Ã·Q₂ (with restr.)
	;                $iNpar = nParams (no restr.) or nFree (with restr.)
	; Applied AFTER Atilde_orig/Abar_orig saves (pre-equilibration for SVD Qxx)
	; Applied BEFORE D-computation so Marquardt-D̃ captures true curvature
	If $mState.scaling Then
		__adj_applyEquilibration($mState, $mAtilde, $iMeq, $iNpar)
	EndIf

	; Step 5: Marquardt D-update and LM augmentation — append √λ·D rows to Ã (or Ā)
	; Must happen after Cholesky (zero B rows would make M singular)
	Local $iMeq_aug = $iMeq
	If MapExists($mState, "fLMLambda") Then
		Local $fLambda_glm = $mState.fLMLambda

		; compute D_new from Ãᵀ·Ã (or Āᵀ·Ā when restrictions present)
		Local $mC_glm = _blas_createMatrix($iNpar, $iNpar)
		_blas_syrk($mAtilde, $mC_glm, 1, 0, "U", "T", $iNpar, $iMeq, $iMeq, $iNpar)
		Local $mDnew_glm = _blas_copy($mC_glm, 0, $mC_glm.rows + 1, 0, 1, $iNpar)
		__adj_applyMarquardtFloorAndSqrt($mDnew_glm, $iNpar)

		; update D = max(D_old, D_new) (D_old set by initializeLMdampingfactor)
		If MapExists($mState, "LM_D") Then
			__adj_updateMarquardtD($mState.LM_D, $mDnew_glm, $iNpar)
		Else
			$mState.LM_D = $mDnew_glm
		EndIf

		; gradient g = Ãᵀ · w̃ (or Āᵀ · w̄, before augmentation)
		Local $mGrad_glm = _blas_createVector($iNpar)
		_blas_gemv($mAtilde, $mW, $mGrad_glm, 1, 0, "T", 1, 1, $iMeq, $iNpar, $iMeq)
		$mState.LM_gradient = $mGrad_glm

		; augment Ã/Ā: copy into larger matrix, then fill √λ·D diagonal (Marquardt)
		$iMeq_aug = $iMeq + $iNpar
		Local $mAtilde_aug = _blas_createMatrix($iMeq_aug, $iNpar)
		_lp_lacpy($mAtilde, $mAtilde_aug, "X", $iMeq, $iNpar, $iMeq, $iMeq_aug)
		__adj_fillLMaugmentation($mAtilde_aug.ptr, $mState.LM_D, $fLambda_glm, $iMeq, $iMeq_aug + 1, $iNpar)
		$mAtilde = $mAtilde_aug

		; augment w̃/w̄: extend with zeros (already zero-initialized by _blas_createVector)
		Local $mW_aug = _blas_createVector($iMeq_aug)
		_lp_lacpy($mW, $mW_aug, "X", $iMeq, 1)
		$mW = $mW_aug
	EndIf

	; Step 6: DGELSY(Ã, w̃) → Δx  or  DGELSY(Ā, w̄) → y
	Local $iLDB_GLM = ($iMeq_aug > $iNpar) ? $iMeq_aug : $iNpar
	If $mW.size < $iLDB_GLM Then
		Local $mWext = _blas_createVector($iLDB_GLM)
		_lp_lacpy($mW, $mWext, "X", $iMeq_aug, 1)
		$mW = $mWext
	EndIf

	Local $iLpErr_GLM, $iRank_GLM
	If $mState.solver = $ADJ_SOLVER_SVD Then
		Local $fRCOND_GLM = $mState.solverRCOND
		Local $tSV_GLM = _lp_gelsd($mAtilde, $mW, "N", Default, Default, Default, Default, Default, $fRCOND_GLM)
		$iLpErr_GLM = @error
		$iRank_GLM = @extended
		If $iLpErr_GLM Then Return SetError($ADJ_ERR_SOLVER, $iLpErr_GLM, False)
		$mState.tSingularValues = $tSV_GLM
	Else
		Local $tJPVT_GLM = _lp_gelsy($mAtilde, $mW, 1, $iMeq_aug, $iNpar, $iMeq_aug, $iLDB_GLM)
		$iLpErr_GLM = @error
		$iRank_GLM = @extended
		If $iLpErr_GLM Then Return SetError($ADJ_ERR_SOLVER, $iLpErr_GLM, False)
		$mState.JPVT = $tJPVT_GLM
	EndIf

	; rank deficiency warning
	If $iRank_GLM < $iNpar Then $mState.rankDeficient = True
	$mState.CholeskyM = $mM

	If $bHasRestrictions Then
		; Store projected Ā (DGELSY-factored) for Qxx — Q₂ already stored above
		$mState.Abar = $mAtilde

		; Extract y (first nFree elements of DGELSY solution)
		Local $mY_glm = _blas_createVector($iNpar)
		_lp_lacpy($mW, $mY_glm, "X", $iNpar, 1)

		; Cache y (reduced step) for LM gain ratio (before back-transformation)
		If MapExists($mState, "fLMLambda") Then
			$mState.LM_step = _blas_copy($mY_glm, 0, 1, 0, 1, $iNpar)
		EndIf

		; reverse equilibration on y (reduced step) BEFORE Q₂ back-transform
		; (LM_step already cached above in equilibrated space for gain ratio)
		__adj_reverseEquilibration($mState, $mY_glm, $iNpar)

		; Back-transform: Δx = α·x_p + Q₂·y  (α = 1/(1+λ) when LM, else 1)
		Local $fAlpha_xp_glm = MapExists($mState, "fLMLambda") ? (1.0 / (1.0 + $mState.fLMLambda)) : 1.0
		$mXd = _la_duplicate($mXp_glm)
		_blas_scal($mXd, $fAlpha_xp_glm)  ; x_p → α·x_p
		_blas_gemv($mQ2_glm, $mY_glm, $mXd, 1, 1, "N", 1, 1, $iNpar_orig, $iNpar)
		$mState.xd = $mXd
	Else
		$mState.Atilde = $mAtilde
		$mXd = _blas_createVector($iNpar)
		_lp_lacpy($mW, $mXd, "X", $iNpar, 1)
		$mState.xd = $mXd
		; reverse equilibration: dx = S⁻¹ · dx̃
		__adj_reverseEquilibration($mState, $mState.xd, $iNpar)
	EndIf

	; Step 7: DGELSY residual norm
	; For GN: iMeq_aug = iMeq → pure data residual ‖Ã·Δx - w̃‖²
	; For LM: includes penalty rows; actual cost is computed in LM loop after param update
	$mState.r2sum = _blas_nrm2($mW, $iNpar, 1, $iMeq_aug - $iNpar)^2

	; Step 8: Observation residuals
	; Minimum-norm solution of B̃·ṽ = w - A·Δx: ṽ = B̃ᵀ · M⁻¹ · d
	; d = w_saved - A_orig · Δx  (use original nParams dimension)
	_blas_gemv($mState.A_orig, $mXd, $mW_saved, -1, 1, "N", 1, 1, $iMeq, $iNpar_orig)
	; L⁻¹ · d  (forward substitution)
	_blas_trsv($mM, $mW_saved, "L", "N", "N", $iMeq)
	; L⁻ᵀ · (L⁻¹ · d)  = M⁻¹ · d  (back substitution)
	_blas_trsv($mM, $mW_saved, "L", "N", "T", $iMeq)
	; ṽ = B̃ᵀ · (M⁻¹ · d)
	Local $mResiduals = _blas_createVector($iPobs)
	_blas_gemv($mState.B_orig, $mW_saved, $mResiduals, 1, 0, "T", 1, 1, $iMeq, $iPobs)

	; back-transform residuals from whitened space
	If StringRegExp($mState.modelType, "^G(?!LM)") Then
		; GGLM/GCLS: v = L · ṽ where Σₗₗ = L·Lᵀ (TRMV)
		_blas_trmv($mState.CovCholeskyL, $mResiduals, "L", "N", "N", $mState.nObs)
		$mState.r = $mResiduals
	ElseIf StringRegExp($mState.modelType, "^W") Then
		; WGLM/WCLS: v = σ · ṽ (diagonal case, element-wise)
		$mState.r = _la_mul($mState.Vector_ObsStdDev, $mResiduals, False)
	Else
		$mState.r = $mResiduals
	EndIf
EndFunc

#EndRegion ; Linear solvers

#Region Iteration finalization

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_computeResiduals
; Description ...: Computes observation residuals for OLS/LSE via contradiction vector
; Syntax.........: __adj_computeResiduals($mSystem, $mState)
; Parameters ....: $mSystem     - [ByRef] Adjustment system map
;                  $mState      - [ByRef] Solver state map
; Return values .: None
; Author ........: AspirinJunkie
; Remarks .......: GLM/CLS: returns immediately (solver computes residuals).
; ===============================================================================================================================
Func __adj_computeResiduals(ByRef $mSystem, ByRef $mState)
	; For GLM/CLS: observation residuals v are already computed by the solver (v = B̃ᵀ·M⁻¹·d).
	; The contradiction vector f(x̂,ŷ) ≈ 0 at convergence is NOT the observation residual.
	If StringRegExp($mState.modelType, '(CLS|GLM)$') Then Return

	; calculate residuals (= contraction vector with current parameters)
	__adj_computeContradiction($mSystem, $mState)

	; copy to state map
	$mState.r = _blas_duplicate($mState.Vector_W)
EndFunc


; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_initLMDamping
; Description ...: Initializes Levenberg-Marquardt damping
; Syntax.........: __adj_initLMDamping($mState)
; Parameters ....: $mState      - [ByRef] Solver state map
; Return values .: Success      - True (fLMLambda, LM_D stored in $mState)
;                  Failure      - None (returns silently if Cholesky fails for GLM)
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: D = √diag(JᵀPJ), λ₀ = τ·max(diag(JᵀPJ)). OLS: fills augmentation directly.
;                  LSE: defers to reduced space. GLM: pre-computes Cholesky M.
; Related .......: __adj_solveIteration, __adj_fillLMaugmentation
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __adj_initLMDamping(ByRef $mState)
	; Marquardt-Skalierung: λ₀ = τ · max(diag(JᵀPJ)), D = √diag(JᵀPJ)
	; Augmentierung mit √λ·D statt √λ·I (skalierungsinvariant)

	Local $mA, $mC, $iRows, $fLambda, $nP = $mState.nParams

	Switch $mState.modelType
		Case "OLS", "WLS", "GLS"
			; After whitening: A is L⁻¹·A (GLS), diag(1/σ)·A (WLS), or unchanged (OLS)
			; N = Aᵀ·A = JᵀPJ (normal equation matrix in whitened space)
			$mA = $mState.Matrix_A
			$iRows = $mA.rows - $mA.cols  ; nObs (total rows minus LM augmentation rows)
			$mC = _blas_createMatrix($nP, $nP, $mA.datatype)
			_blas_syrk($mA, $mC, 1, 0, "U", "T", $nP, $iRows, $mA.rows, $nP)

			; Marquardt scaling: D = √diag(N) with floor
			Local $mD = _blas_copy($mC, 0, $mC.rows + 1, 0, 1, $nP)
			__adj_applyMarquardtFloorAndSqrt($mD, $nP)
			$mState.LM_D = $mD

			; λ₀ = τ · max(diag(N))
			$fLambda = _blas_amax($mC, 0, $mC.rows + 1, $nP) * $mState.fLMTau
			$mState.fLMLambda = $fLambda

			; fill augmentation diagonal with √λ · Dᵢ (Marquardt)
			Local $iOffset = $mState.nFormulas
			__adj_fillLMaugmentation($mA.ptr, $mD, $fLambda, $iOffset, $mA.rows + 1, $nP)

			; gradient g = Aᵀ · w (whitened, only first nObs rows)
			Local $mGrad = _blas_createVector($nP)
			_blas_gemv($mA, $mState.Vector_W, $mGrad, 1, 0, "T", 1, 1, $iRows, $nP, $mA.rows)
			$mState.LM_gradient = $mGrad

		Case "LSE", "WLSE", "GLSE"
			; N = Aᵀ · A (A is already whitened → gives JᵀPJ)
			; LM augmentation is applied to the reduced system Ā in __adj_dispatchLinearSolver,
			; not directly to A (since the nullspace projection operates on A).
			; D for LSE will be computed in reduced space inside __adj_dispatchLinearSolver.
			$mA = $mState.Matrix_A
			$iRows = $mState.nFormulas
			$mC = _blas_createMatrix($mA.cols, $mA.cols, $mA.datatype)
			_blas_syrk($mA, $mC, 1, 0, "U", "T", $mA.cols, $iRows, $mA.rows, $mA.cols)

			$fLambda = _blas_amax($mC, 0, $mC.rows + 1, $mA.cols) * $mState.fLMTau
			$mState.fLMLambda = $fLambda
			; Note: LM_D, LM_gradient are set in __adj_dispatchLinearSolver for LSE

		Case "GLM", "WGLM", "GGLM"
			; Effective normal equation matrix: N = Ãᵀ·Ã where Ã = M⁻½·A
			; M = B̃B̃ᵀ = L·Lᵀ (Cholesky), B̃ already whitened
			$mA = $mState.Matrix_A
			Local $mB_lm = $mState.Matrix_B
			Local $iMeq_lm = $mA.rows, $iNpar_lm = $mA.cols

			; M = B̃ · B̃ᵀ
			Local $mM_lm = _blas_createMatrix($iMeq_lm, $iMeq_lm)
			_blas_syrk($mB_lm, $mM_lm, 1, 0, "L", "N", $iMeq_lm, $mB_lm.cols)

			; 1-norm of M (needed for condition estimation in solver)
			Local $fANORM_lm = _lp_lansy($mM_lm, "1", "L", $iMeq_lm)

			; Cholesky M = L · Lᵀ
			_lp_potrf($mM_lm, "L")
			If @error Then
				Return ; M not SPD — cannot initialize LM
			EndIf

			; Ã = L⁻¹ · A (on a copy, since A is needed unchanged for the solver)
			Local $mAtilde_lm = _la_duplicate($mA)
			_blas_trsm($mM_lm, $mAtilde_lm, 1, "L", "L", "N", "N", $iMeq_lm, $iNpar_lm)

			; N = Ãᵀ · Ã (normal equation matrix of the transformed system)
			$mC = _blas_createMatrix($iNpar_lm, $iNpar_lm)
			_blas_syrk($mAtilde_lm, $mC, 1, 0, "U", "T", $iNpar_lm, $iMeq_lm, $iMeq_lm, $iNpar_lm)

			; Marquardt scaling: D = √diag(N) with floor
			Local $mD_glm = _blas_copy($mC, 0, $mC.rows + 1, 0, 1, $iNpar_lm)
			__adj_applyMarquardtFloorAndSqrt($mD_glm, $iNpar_lm)
			; When restrictions present, LM_D will be recomputed in reduced space inside __adj_solveGLM
			If $mState.nRestrictions = 0 Then $mState.LM_D = $mD_glm

			; λ₀ = τ · max(diag(N))
			$fLambda = _blas_amax($mC, 0, $mC.rows + 1, $iNpar_lm) * $mState.fLMTau
			$mState.fLMLambda = $fLambda

			; cache Cholesky factor L and ANORM for reuse in solver (first iteration only)
			$mState.CholeskyL     = $mM_lm
			$mState.CholeskyANORM = $fANORM_lm

			; gradient g = Ãᵀ · w̃ — computed in __adj_dispatchLinearSolver (Ã not yet available here after first solve)
			; Note: LM_gradient for GLM is set in __adj_dispatchLinearSolver
	EndSwitch

EndFunc


; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_updateParameters
; Description ...: Applies parameter update x₁ = x₀ + Δx to model params and state vectors
; Syntax.........: __adj_updateParameters($mSystem, $mState)
; Parameters ....: $mSystem     - [ByRef] Adjustment system map
;                  $mState      - [ByRef] Solver state map
; Return values .: None
; Author ........: AspirinJunkie
; ===============================================================================================================================
Func __adj_updateParameters(ByRef $mSystem, ByRef $mState)

	Local $mX0 = $mState.Vector_x0
	Local $mXd = $mState.xd
	Local $mX1 = _blas_duplicate($mXd)
	_blas_axpy($mX0, $mX1)

	; overwrite map of parameter values (for jacobian and contract vector) with new values
	Local $mModel = $mSystem.model
	Local $mIdxParams = $mState.idxParams, $mParams = $mModel.params
	For $sKey In MapKeys($mParams)
		$mParams[$sKey] = __adj_vecGet($mX1, $mIdxParams[$sKey])
	Next
	$mModel.params = $mParams
	$mSystem.model = $mModel ;! MAP WRITE-BACK

	; save x1-array
	$mState.x1 = $mX1
	$mState.Vector_x0 = _blas_duplicate($mX1)
EndFunc


; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_updateObservations
; Description ...: Stub for future observation updates (currently empty)
; Syntax.........: __adj_updateObservations($mState)
; Parameters ....: $mState      - [ByRef] Solver state map
; Return values .: None
; Author ........: AspirinJunkie
; ===============================================================================================================================
Func __adj_updateObservations(ByRef $mState)
	#forceref $mState
EndFunc

#EndRegion ; Iteration finalization

#Region Levenberg-Marquardt

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_updateLMDamping
; Description ...: Computes gain ratio and updates LM damping parameter λ
; Syntax.........: __adj_updateLMDamping($mState)
; Parameters ....: $mState      - [ByRef] Solver state map
; Return values .: Success      - New λ value. @extended = True if step accepted, False if rejected
;                  Failure      - None
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: ρ = (F(x)-F(x+h))/predicted. Nielsen update: λ×⅓ if ρ>0.75, λ×2 if ρ<0.25.
; Related .......: __adj_solveNonlinear
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __adj_updateLMDamping(ByRef $mState)

	; calculate the gain ratio rho with Marquardt scaling (Madsen/Nielsen/Tingleff 2004)
	; pred = hᵀ·g + λ·‖D·h‖² where h = step, g = gradient, D = Marquardt scaling
	Local $fR2sum     = $mState.r2sum
	Local $fR2sum_old = $mState.r2sum0
	Local $fLambda    = $mState.fLMLambda
	Local $fActual    = $fR2sum_old - $fR2sum

	; determine step vector h: reduced step z for LSE, full Δx for other models
	Local $mStep = MapExists($mState, "LM_step") ? $mState.LM_step : $mState.xd
	Local $iStepN = $mStep.size

	; compute predicted decrease: fPred = hᵀ·g + λ·‖D·h‖²
	Local $fPred = 0
	If MapExists($mState, "LM_gradient") And MapExists($mState, "LM_D") Then
		; gradient term: hᵀ · g
		Local $fHtG = _blas_dot($mStep, $mState.LM_gradient, 0, 0, 1, 1, $iStepN)

		; scaled step norm: ‖D·h‖² = Σ Dᵢ² · hᵢ²
		Local $fDh2 = 0
		Local $tD_gr = ($mState.LM_D).struct, $tH_gr = $mStep.struct
		For $j = 1 To $iStepN
			Local $fDiHi = DllStructGetData($tD_gr, 1, $j) * DllStructGetData($tH_gr, 1, $j)
			$fDh2 += $fDiHi * $fDiHi
		Next

		$fPred = $fHtG + $fLambda * $fDh2
	EndIf

	; guard: if predicted decrease is non-positive, fall back to simple Levenberg formula
	If $fPred <= 0 Then
		Local $fh2 = _blas_nrm2($mStep, 0, 1, $iStepN)^2
		$fPred = $fLambda * $fh2
	EndIf

	Local $fRho = ($fPred > 0) ? ($fActual / $fPred) : 0

	; λ update according to Nielsen (1999): ×⅓ if ρ > 0.75, ×2 if ρ < 0.25
	Local $fFactor
	If $fRho > 0.75 Then
		$fFactor = 1 / 3
	ElseIf $fRho < 0.25 Then
		$fFactor = 2.0
	Else
		$fFactor = 1.0
	EndIf

	; accept step if ρ > 0, reject otherwise
	Return SetExtended($fRho > 0 ? True : False, _Max(1e-7, _Min(1e7, $fLambda * $fFactor)))
EndFunc


; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_applyMarquardtFloorAndSqrt
; Description ...: Applies floor and sqrt to Marquardt scaling: D_i = √max(D_i, ε·max(D))
; Syntax.........: __adj_applyMarquardtFloorAndSqrt($mD, $iN)
; Parameters ....: $mD          - BLAS vector with raw diagonal values
;                  $iN          - Number of elements
; Return values .: None (in-place)
; Author ........: AspirinJunkie
; ===============================================================================================================================
Func __adj_applyMarquardtFloorAndSqrt($mD, $iN)
	Local $fMaxD = _blas_amax($mD, 0, 1, $iN)
	Local $fFloor = $fMaxD * $__ADJ_MARQUARDT_REL_FLOOR
	If $fFloor < $__ADJ_MARQUARDT_ABS_FLOOR Then $fFloor = $__ADJ_MARQUARDT_ABS_FLOOR

	Local $tD = $mD.struct
	For $j = 1 To $iN
		Local $fVal = DllStructGetData($tD, 1, $j)
		DllStructSetData($tD, 1, Sqrt($fVal < $fFloor ? $fFloor : $fVal), $j)
	Next
EndFunc


; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_fillLMaugmentation
; Description ...: Fills LM diagonal with √λ·D_i values
; Syntax.........: __adj_fillLMaugmentation($pMatrix, $mD, $fLambda, $iOffset, $iStride, $iN)
; Parameters ....: $pMatrix     - Matrix pointer or struct
;                  $mD          - Marquardt D vector
;                  $fLambda     - Damping parameter
;                  $iOffset     - Start offset
;                  $iStride     - Diagonal stride
;                  $iN          - Number of elements
; Return values .: None
; Author ........: AspirinJunkie
; ===============================================================================================================================
Func __adj_fillLMaugmentation($pMatrix, $mD, $fLambda, $iOffset, $iStride, $iN)
	Local $mScaled = _la_duplicate($mD)
	_blas_scal($mScaled, Sqrt($fLambda))
	_blas_copy($mScaled, 0, 1, $iOffset, $iStride, $iN, $pMatrix, False)
EndFunc

#EndRegion ; Levenberg-Marquardt

#Region Scaling and helper functions

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_applyEquilibration
; Description ...: Jacobi-Equilibration: normalizes column norms
; Syntax.........: __adj_applyEquilibration($mState, $mMatrix, $iNormRows, $iCols)
; Parameters ....: $mState      - [ByRef] Solver state map (stores EquilibrationScale)
;                  $mMatrix     - [ByRef] Matrix to equilibrate (modified in-place)
;                  $iNormRows   - Number of rows used for norm computation
;                  $iCols       - Number of columns to equilibrate
; Return values .: Success      - True (matrix scaled, EquilibrationScale updated)
;                  Failure      - None
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: Max-update strategy: S = max(S_old, S_new). Stores S in $mState.EquilibrationScale.
; Related .......: __adj_reverseEquilibration
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __adj_applyEquilibration(ByRef $mState, ByRef $mMatrix, $iNormRows, $iCols)
	If Not $mState.scaling Then Return

	Local $fEpsMach = 2.2204460492503131e-16  ; machine epsilon for double
	Local $iLDA = $mMatrix.rows

	; compute column 2-norms: S_j = ‖M_*,j‖₂ (only observation rows, not LM augmentation)
	Local $mS = _blas_createVector($iCols)
	Local $tS = $mS.struct
	For $j = 0 To $iCols - 1
		Local $fNorm = _blas_nrm2($mMatrix, $j * $iLDA, 1, $iNormRows)
		; floor to prevent division by zero (null column = unobservable parameter)
		If $fNorm < $fEpsMach Then $fNorm = $fEpsMach
		DllStructSetData($tS, 1, $fNorm, $j + 1)
	Next

	; max-update: S = max(S_old, S_new) — monotonically non-decreasing within VKS iteration
	If MapExists($mState, "EquilibrationScale") Then
		Local $mSOld = $mState.EquilibrationScale
		Local $tSOld = $mSOld.struct
		For $j = 1 To $iCols
			If DllStructGetData($tS, 1, $j) < DllStructGetData($tSOld, 1, $j) Then
				DllStructSetData($tS, 1, DllStructGetData($tSOld, 1, $j), $j)
			EndIf
		Next
	EndIf

	$mState.EquilibrationScale = $mS

	; scale columns in-place: M_*,j /= S_j (only observation rows, NOT LM augmentation rows)
	For $j = 0 To $iCols - 1
		Local $fInvS = 1.0 / DllStructGetData($tS, 1, $j + 1)
		_blas_scal($mMatrix, $fInvS, $j * $iLDA, 1, $iNormRows)
	Next
EndFunc


; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_reverseEquilibration
; Description ...: Back-transforms vector from equilibrated to original space: v_j /= S_j
; Syntax.........: __adj_reverseEquilibration($mState, $mVector, $iN)
; Parameters ....: $mState      - [ByRef] Solver state map (reads EquilibrationScale)
;                  $mVector     - [ByRef] Vector to back-transform (modified in-place)
;                  $iN          - Number of elements
; Return values .: None
; Author ........: AspirinJunkie
; ===============================================================================================================================
Func __adj_reverseEquilibration(ByRef $mState, ByRef $mVector, $iN)
	If Not MapExists($mState, "EquilibrationScale") Then Return

	Local $tS = $mState.EquilibrationScale.struct
	Local $tV = $mVector.struct
	For $j = 1 To $iN
		DllStructSetData($tV, 1, _
			DllStructGetData($tV, 1, $j) / DllStructGetData($tS, 1, $j), $j)
	Next
EndFunc


; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_updateMarquardtD
; Description ...: Element-wise max-update: D = max(D_old, D_new)
; Syntax.........: __adj_updateMarquardtD($mD, $mDnew, $iN)
; Parameters ....: $mD          - Current Marquardt D vector (modified in-place)
;                  $mDnew       - New D vector to compare against
;                  $iN          - Number of elements
; Return values .: None (in-place)
; Author ........: AspirinJunkie
; ===============================================================================================================================
Func __adj_updateMarquardtD($mD, $mDnew, $iN)
	Local $tD = $mD.struct, $tDn = $mDnew.struct
	For $j = 1 To $iN
		Local $fOld = DllStructGetData($tD, 1, $j)
		Local $fNew = DllStructGetData($tDn, 1, $j)
		If $fNew > $fOld Then DllStructSetData($tD, 1, $fNew, $j)
	Next
EndFunc

#EndRegion ; Scaling and helper functions

#Region Numerical differentiation

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_evalDerivative
; Description ...: Returns partial derivative using analytical formula or numerical fallback
; Syntax.........: __adj_evalDerivative($mDerivExec, $sSymbol, $sFunc, $mPar, $mObs, $bObs, $sDeriveMethod)
; Parameters ....: $mDerivExec    - Map of executable derivative strings (or False)
;                  $sSymbol       - Variable name
;                  $sFunc         - Executable function string
;                  $mPar          - Parameter values map
;                  $mObs          - Observation values map
;                  $bObs          - True if deriving w.r.t. observation
;                  $sDeriveMethod - Numerical method
; Return values .: Float — the derivative value
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......:
; Related .......: __adj_derivate1D, __adj_computeJacobians
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __adj_evalDerivative($mDerivExec, $sSymbol, $sFunc, $mPar, $mObs, $bObs, $sDeriveMethod)
	If IsMap($mDerivExec) And MapExists($mDerivExec, $sSymbol) Then Return Execute($mDerivExec[$sSymbol])
	Return __adj_derivate1D($sFunc, $sSymbol, Default, $mPar, $mObs, $bObs, $sDeriveMethod)
EndFunc


; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_derivate1D
; Description ...: Numerical differentiation with multiple methods
; Syntax.........: __adj_derivate1D($sFunc, $sParam [, $fValue [, $mPar [, $mObs [, $bObs [, $sMethod [, $fH [, $iIterationLimit [, $fInitialH_Ridder]]]]]]]])
; Parameters ....: $sFunc              - Executable function string
;                  $sParam             - Variable name to differentiate w.r.t.
;                  $fValue             - [optional] Point (Default from mPar/mObs)
;                  $mPar               - Parameter values map
;                  $mObs               - Observation values map
;                  $bObs               - True for observation
;                  $sMethod            - "Central"/"Forward"/"Backward"/"Ridder"/"Higham" etc.
;                  $fH                 - [optional] Step size
;                  $iIterationLimit    - Max iterations for Ridder/Higham
;                  $fInitialH_Ridder   - Initial step for Ridder method
; Return values .: Float — derivative value
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: Supports Central 1-4th order, Forward, Backward, Ridder (Richardson
;                  extrapolation), Higham (iterative refinement). Scales step size to
;                  parameter magnitude.
; Related .......: __adj_evalDerivative
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __adj_derivate1D($sFunc, $sParam, $fValue = Default, $mPar = Default, $mObs = Default, $bObs = False, $sMethod = "Central", $fH = Default, $iIterationLimit = 10, $fInitialH_Ridder = 1e-5)
	Local $fA, $fB, _
		  $fStep, $fRet

	If IsKeyword($fValue) = 1 Then $fValue = $bObs ? $mObs[$sParam] : $mPar[$sParam]

	; scale step sizes to parameter/observation magnitude
	; Smooth blend h = stepBase · √(x² + stepBase²):
	;   |x| ≫ stepBase  →  h ≈ stepBase · |x|       (relative scaling, avoids catastrophic cancellation for large x)
	;   |x| ≪ stepBase  →  h ≈ stepBase²            (floored, avoids over-stepping for tiny x where the old
	;                                                 max(|x|,1) used h ≈ stepBase, way too coarse for x ≪ 1)
	If IsKeyword($fH) = 1 Then
		; optimal stepBase: eps^(1/3) for central differences, eps^(1/2) for forward/backward
		Local $fStepBase
		Switch $sMethod
			Case "Forward", "Backward"
				$fStepBase = Sqrt($f_LA_DBL_EPS)
			Case Else
				$fStepBase = $f_LA_DBL_STEP
		EndSwitch
		$fH = $fStepBase * Sqrt($fValue * $fValue + $fStepBase * $fStepBase)
	EndIf
	; Ridder/Higham initial step keeps the older max(|x|,1) scaling — they refine adaptively from there
	$fInitialH_Ridder = $fInitialH_Ridder * (Abs($fValue) > 1 ? Abs($fValue) : 1)

	If $bObs Then
		; Function value at which to derive
		If IsKeyword($fValue) = 1 Then $fValue = $mObs[$sParam]

		Switch $sMethod
			Case "Central", "Central1"
				$mObs[$sParam] = $fValue - $fH
				$fA = Execute($sFunc)
				$mObs[$sParam] = $fValue + $fH
				$fB = Execute($sFunc)

				Return ($fB - $fA) / (2 * $fH)

			Case "Central2"
				$mObs[$sParam] = $fValue - 2 * $fH
				$fRet = Execute($sFunc) / 12
				$mObs[$sParam] = $fValue - $fH
				$fRet -= 2/3 * Execute($sFunc)
				$mObs[$sParam] = $fValue + $fH
				$fRet += 2/3 * Execute($sFunc)
				$mObs[$sParam] = $fValue + 2 * $fH
				$fRet -= Execute($sFunc) / 12

				Return $fRet / $fH

			Case "Central3"
				$mObs[$sParam] = $fValue - 3 * $fH
				$fRet = - Execute($sFunc) / 60
				$mObs[$sParam] = $fValue - 2 * $fH
				$fRet += 3 * Execute($sFunc) / 20
				$mObs[$sParam] = $fValue - $fH
				$fRet -= 3 * Execute($sFunc) / 4
				$mObs[$sParam] = $fValue + $fH
				$fRet += 3 * Execute($sFunc) / 4
				$mObs[$sParam] = $fValue + 2 * $fH
				$fRet -= 3 * Execute($sFunc) / 20
				$mObs[$sParam] = $fValue + 3 * $fH
				$fRet += Execute($sFunc) / 60

				Return $fRet / $fH

			Case "Central4"
				$mObs[$sParam] = $fValue - 4 * $fH
				$fRet = Execute($sFunc) / 280
				$mObs[$sParam] = $fValue - 3 * $fH
				$fRet -= 4 * Execute($sFunc) / 105
				$mObs[$sParam] = $fValue - 2 * $fH
				$fRet += Execute($sFunc) / 5
				$mObs[$sParam] = $fValue - $fH
				$fRet -= 4 * Execute($sFunc) / 5
				$mObs[$sParam] = $fValue + $fH
				$fRet += 4 * Execute($sFunc) / 5
				$mObs[$sParam] = $fValue + 2 * $fH
				$fRet -= Execute($sFunc) / 5
				$mObs[$sParam] = $fValue + 3 * $fH
				$fRet += 4 * Execute($sFunc) / 105
				$mObs[$sParam] = $fValue + 4 * $fH
				$fRet -= Execute($sFunc) / 280

				Return $fRet / $fH

			Case "Forward"
				$fA = Execute($sFunc)
				$mObs[$sParam] = $fValue + $fH
				$fB = Execute($sFunc)

				Return ($fB - $fA) / $fH

			Case "Backward"
				$fB = Execute($sFunc)
				$mObs[$sParam] = $fValue - $fH
				$fA = Execute($sFunc)

				Return ($fB - $fA) / $fH

			Case "Ridder"
				Local $fD1, $fD2, $fD3, $fErr, $fErrOld = 1e100

				$fStep = $fInitialH_Ridder

				; initial difference quotient
				$mObs[$sParam] = $fValue + $fStep
				$fA = Execute($sFunc)
				$mObs[$sParam] = $fValue - $fStep
				$fB = Execute($sFunc)
				$fD1 = ($fA - $fB) / (2 * $fStep)

				For $i = 0 To $iIterationLimit
					; half step size
					$fStep /= 2

					; calculate new quotient with halfed step
					$mObs[$sParam] = $fValue + $fStep
					$fA = Execute($sFunc)
					$mObs[$sParam] = $fValue - $fStep
					$fB = Execute($sFunc)
					$fD2 = ($fA - $fB) / (2 * $fStep)

					; extrapolate the new value:
					$fD3 = (4 * $fD2 - $fD1) / 3

					; estimate the error
					$fErr = Abs($fD3 - $fD2) / 3  ; = |D3 - D2| / (2^O - 1)); "O" is typically 2 for the central difference

					; leave if not converge
					If $fErr > $fErrOld Then Return SetError(2, $fD1)
					$fErrOld = $fErr

					$fD1 = $fD3

					; leave if target accuracy is reached
					If $fErr < $fH Then Return $fD1
				Next

			Case "Higham"
				; initial parameters
				Local $fTargetError = $f_LA_DBL_EPS^(6/7) ; 6/7 correct digits for double type
				Local $fError = 1e100
				$fStep = $fH

				; first derivation
				$mObs[$sParam] = $fValue + $fStep
				$fA = Execute($sFunc)
				$mObs[$sParam] = $fValue - $fStep
				$fB = Execute($sFunc)
				Local $fD = ($fA - $fB) / (2 * $fStep)

				; iterative refinement
				Local $fDprev = $fD
				For $i = 0 To $iIterationLimit
					$fStep = $fStep / 2

					$mObs[$sParam] = $fValue + $fStep
					$fA = Execute($sFunc)
					$mObs[$sParam] = $fValue - $fStep
					$fB = Execute($sFunc)
					$fD = ($fA - $fB) / (2 * $fStep)

					$fError = Abs($fD - $fDprev)
					$fDprev = $fD

					If $fError > $fTargetError Then ExitLoop
				Next

				Return $fD

			Case Else
				Return SetError(1, 0, Null)

		EndSwitch

	Else
		; Function value at which to derive
		If IsKeyword($fValue) = 1 Then $fValue = $mPar[$sParam]

		Switch $sMethod
			Case "Central", "Central1"
				$mPar[$sParam] = $fValue - $fH
				$fA = Execute($sFunc)
				$mPar[$sParam] = $fValue + $fH
				$fB = Execute($sFunc)

				Return ($fB - $fA) / (2 * $fH)

			Case "Central2"
				$mPar[$sParam] = $fValue - 2 * $fH
				$fRet = Execute($sFunc) / 12
				$mPar[$sParam] = $fValue - $fH
				$fRet -= 2/3 * Execute($sFunc)
				$mPar[$sParam] = $fValue + $fH
				$fRet += 2/3 * Execute($sFunc)
				$mPar[$sParam] = $fValue + 2 * $fH
				$fRet -= Execute($sFunc) / 12

				Return $fRet / $fH

			Case "Central3"
				$mPar[$sParam] = $fValue - 3 * $fH
				$fRet = - Execute($sFunc) / 60
				$mPar[$sParam] = $fValue - 2 * $fH
				$fRet += 3 * Execute($sFunc) / 20
				$mPar[$sParam] = $fValue - $fH
				$fRet -= 3 * Execute($sFunc) / 4
				$mPar[$sParam] = $fValue + $fH
				$fRet += 3 * Execute($sFunc) / 4
				$mPar[$sParam] = $fValue + 2 * $fH
				$fRet -= 3 * Execute($sFunc) / 20
				$mPar[$sParam] = $fValue + 3 * $fH
				$fRet += Execute($sFunc) / 60

				Return $fRet / $fH

			Case "Central4"
				$mPar[$sParam] = $fValue - 4 * $fH
				$fRet = Execute($sFunc) / 280
				$mPar[$sParam] = $fValue - 3 * $fH
				$fRet -= 4 * Execute($sFunc) / 105
				$mPar[$sParam] = $fValue - 2 * $fH
				$fRet += Execute($sFunc) / 5
				$mPar[$sParam] = $fValue - $fH
				$fRet -= 4 * Execute($sFunc) / 5
				$mPar[$sParam] = $fValue + $fH
				$fRet += 4 * Execute($sFunc) / 5
				$mPar[$sParam] = $fValue + 2 * $fH
				$fRet -= Execute($sFunc) / 5
				$mPar[$sParam] = $fValue + 3 * $fH
				$fRet += 4 * Execute($sFunc) / 105
				$mPar[$sParam] = $fValue + 4 * $fH
				$fRet -= Execute($sFunc) / 280

				Return $fRet / $fH

			Case "Forward"
				$fA = Execute($sFunc)
				$mPar[$sParam] = $fValue + $fH
				$fB = Execute($sFunc)

				Return ($fB - $fA) / $fH

			Case "Backward"
				$fB = Execute($sFunc)
				$mPar[$sParam] = $fValue - $fH
				$fA = Execute($sFunc)

				Return ($fB - $fA) / $fH

			Case "Ridder"
				$fErrOld = 1e100

				$fStep = $fInitialH_Ridder

				; initial difference quotient
				$mPar[$sParam] = $fValue + $fStep
				$fA = Execute($sFunc)
				$mPar[$sParam] = $fValue - $fStep
				$fB = Execute($sFunc)
				$fD1 = ($fA - $fB) / (2 * $fStep)

				For $i = 0 To $iIterationLimit
					; half step size
					$fStep /= 2

					; calculate new quotient with halfed step
					$mPar[$sParam] = $fValue + $fStep
					$fA = Execute($sFunc)
					$mPar[$sParam] = $fValue - $fStep
					$fB = Execute($sFunc)
					$fD2 = ($fA - $fB) / (2 * $fStep)

					; extrapolate the new value:
					$fD3 = (4 * $fD2 - $fD1) / 3

					; estimate the error
					$fErr = Abs($fD3 - $fD2) / 3  ; = |D3 - D2| / (2^O - 1)); "O" is typically 2 for the central difference

					; leave if not converge
					If $fErr > $fErrOld Then Return SetError(2, $fD1)
					$fErrOld = $fErr

					$fD1 = $fD3

					; leave if target accuracy is reached
					If $fErr < $fH Then Return $fD1
				Next

			Case "Higham"
				; initial parameters
				$fTargetError = $f_LA_DBL_EPS^(6/7) ; 6/7 correct digits for double type
				$fError = 1e100
				$fStep = $fH

				; first derivation
				$mPar[$sParam] = $fValue + $fStep
				$fA = Execute($sFunc)
				$mPar[$sParam] = $fValue - $fStep
				$fB = Execute($sFunc)
				$fD = ($fA - $fB) / (2 * $fStep)

				; iterative refinement
				$fDprev = $fD
				For $i = 0 To $iIterationLimit
					$fStep = $fStep / 2

					$mPar[$sParam] = $fValue + $fStep
					$fA = Execute($sFunc)
					$mPar[$sParam] = $fValue - $fStep
					$fB = Execute($sFunc)
					$fD = ($fA - $fB) / (2 * $fStep)

					$fError = Abs($fD - $fDprev)
					$fDprev = $fD

					If $fError > $fTargetError Then ExitLoop
				Next

				Return $fD

			Case Else
				Return SetError(1, 0, Null)

		EndSwitch
	EndIf
EndFunc

#EndRegion ; Numerical differentiation

