#include-once

; #INDEX# =======================================================================================================================
; Title .........: AdjustmentRobust
; AutoIt Version : 3.3.16.1
; Description ...: Robust estimation via IRLS (Iteratively Reweighted Least Squares) for geodetic
;                  adjustment. Supports L1, Huber, Hampel, Biweight, BIBER, and ModifiedM estimators.
; Author(s) .....: AspirinJunkie
; ===============================================================================================================================


#Region Configuration

; #FUNCTION# ====================================================================================================================
; Name...........: _adj_robustDefaults
; Description ...: Returns default tuning parameters for a robust estimator
; Syntax.........: _adj_robustDefaults($sEstimator)
; Parameters ....: $sEstimator  - Estimator name: "L1", "Huber", "Hampel", "Biweight", "BIBER",
;                                 or "ModifiedM"
; Return values .: Success      - Map with estimator-specific tuning parameters:
;                               |"L1"       - no tuning parameters
;                               |"Huber"    - .c = 1.345
;                               |"Hampel"   - .a = 1.7, .b = 3.4, .c = 8.5
;                               |"Biweight" - .c = 4.685
;                               |"BIBER"    - .c = 3.5
;                               |"ModifiedM"- .c = 1.5
;                               +All include .scale = "MAD" and .outlierThreshold = 2.5
;                  Failure      - SetError($ADJ_ERR_INPUT) for unknown estimator
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: Pass the returned Map as $mConfig.robustParams when configuring robust estimation.
;                  Scale methods: "MAD" (default), "s0", "apriori", "fixed".
; Related .......: _adj_solve, __adj_robustIRLS
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func _adj_robustDefaults($sEstimator)
	Local $mParams[]
	Switch $sEstimator
		Case "L1"
			; no tuning parameters
		Case "Huber"
			$mParams.c = 1.345
		Case "Hampel"
			$mParams.a = 1.7
			$mParams.b = 3.4
			$mParams.c = 8.5
		Case "Biweight"
			$mParams.c = 4.685
		Case "BIBER"
			$mParams.c = 3.5
		Case "ModifiedM"
			$mParams.c = 1.5
		Case Else
			Return SetError($ADJ_ERR_INPUT, 0, False)
	EndSwitch
	$mParams.scale = "MAD"
	$mParams.outlierThreshold = 2.5
	Return $mParams
EndFunc

#EndRegion ; Configuration

#Region IRLS iteration loop

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_robustIRLS
; Description ...: Main IRLS loop for robust estimation (Phase 1 of _adj_solve when robust is active)
; Syntax.........: __adj_robustIRLS(ByRef $mSystem, ByRef $mState)
; Parameters ....: $mSystem     - [ByRef] The adjustment system
;                  $mState      - [ByRef] Solver state
; Return values .: Success      - Implicit (results stored in $mState: robustWeights, robustIterations,
;                                 robustScale, robustConverged)
;                  Failure      - SetError($ADJ_ERR_INPUT) for unknown estimator or unsupported model,
;                                 SetError($ADJ_ERR_*) on solver failure
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: Steps: (1) initial non-robust solve, (2) compute redundancy for BIBER/ModifiedM,
;                  (3) save original weights and upgrade model type to weighted variant,
;                  (4) iterate: scale estimation → weight update → convergence check →
;                  whitening update → re-solve, (5) store results in state.
;                  Generalized models with full covariance matrix not yet supported.
; Related .......: __adj_robustScale, __adj_updateRobustWeights, __adj_robustWeight,
;                  __adj_computeRedundancyDiag
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __adj_robustIRLS(ByRef $mSystem, ByRef $mState)
	Local $mConfig = $mSystem.config

	; validate estimator name
	Switch $mConfig.robust
		Case "L1", "Huber", "Hampel", "Biweight", "BIBER", "ModifiedM"
			; ok
		Case Else
			Return SetError($ADJ_ERR_INPUT, 0, False)
	EndSwitch

	; auto-fill default tuning parameters if user did not provide any
	; (otherwise __adj_robustWeight crashes accessing $mParams.c on Null)
	If IsKeyword($mConfig.robustParams) = 2 Or $mConfig.robustParams = Null Then
		$mConfig.robustParams = _adj_robustDefaults($mConfig.robust)
		$mSystem.config = $mConfig
	EndIf

	; guard: generalized models with full covariance not yet supported
	If $mState.hasCovariances Then
		Return SetError($ADJ_ERR_INPUT, 7, False)
	EndIf

	; 1. Initial (non-robust) solution
	__adj_solveNonlinear($mSystem, $mState)
	If @error Then Return SetError(@error, @extended, False)
	__adj_computeResiduals($mSystem, $mState)

	; 2. Compute redundancy numbers once (for BIBER/ModifiedM only)
	Local $aRedundancy = __adj_computeRedundancyDiag($mSystem, $mState)

	; 3. Save original weights and upgrade model type to weighted variant
	;    IRLS re-weights observations → whitening must be applied in subsequent iterations.
	;    For initially-unweighted models (OLS/LSE/CLS/GLM) we must:
	;    a) initialize .weight = 1.0 on every observation (may not exist yet)
	;    b) upgrade the model type so __adj_applyWhitening does not skip
	Local $mModel = $mSystem.model
	Local $mObservations = $mModel.obs
	For $sObsName In MapKeys($mState.idxObs)
		Local $mObs = $mObservations[$sObsName]
		If Not MapExists($mObs, "weight") Then $mObs.weight = 1.0 / ($mObs.stdDev ^ 2)
		$mObs._weightOriginal = $mObs.weight
		$mObservations[$sObsName] = $mObs
	Next
	$mModel.obs = $mObservations
	$mSystem.model = $mModel

	; upgrade unweighted model types to weighted variants for whitening support
	Switch $mState.modelType
		Case "OLS"
			$mState.modelType = "WLS"
		Case "LSE"
			$mState.modelType = "WLSE"
		Case "CLS"
			$mState.modelType = "WCLS"
		Case "GLM"
			$mState.modelType = "WGLM"
	EndSwitch

	; Snapshot the a-priori 1/σᵢ vector BEFORE the IRLS loop starts overwriting
	; Vector_ObsInvStdDev with √weight (which then folds in _robustWeight from prior
	; iterations).  __adj_robustScale must read from this snapshot, otherwise the MAD
	; computation includes the running robust weights — observations get doubly down-
	; weighted, the scale collapses, and ever more points are flagged as outliers.
	$mState.Vector_ObsInvStdDev_apriori = _blas_duplicate($mState.Vector_ObsInvStdDev)

	; 4. IRLS iteration loop
	Local $bConverged = False
	Local $fScale = 0.0
	Local $iIterations = 0

	Local $sScaleMethod = "MAD"
	If $mConfig.robustParams <> Null Then
		If MapExists($mConfig.robustParams, "scale") Then $sScaleMethod = $mConfig.robustParams.scale
	EndIf

	For $iIRLS = 1 To $mConfig.robustMaxIter
		; robust scale estimation
		$fScale = __adj_robustScale($mState, $sScaleMethod, $mConfig.robustParams)

		; guard: MAD = 0 → fallback to s0
		If $fScale < 1e-15 Then
			Local $iDOF = $mState.nObs - $mState.nParams + $mState.nRestrictions
			$fScale = ($iDOF > 0) ? Sqrt($mState.r2sum / $iDOF) : 0
			If $fScale < 1e-15 Then ExitLoop  ; residuals are zero → abort
		EndIf

		; update weights + check convergence
		$bConverged = __adj_updateRobustWeights($mSystem, $mState, $fScale, $aRedundancy)
		$iIterations = $iIRLS
		If $bConverged Then ExitLoop

		; propagate new weights to whitening vectors
		__adj_updateWhiteningVectors($mSystem, $mState)

		; reset iteration state for fresh GN/LM convergence (mirrors VCE reset)
		If MapExists($mState, "nIterations")      Then MapRemove($mState, "nIterations")
		If MapExists($mState, "r_accumulated")     Then MapRemove($mState, "r_accumulated")
		If MapExists($mState, "LM_D")              Then MapRemove($mState, "LM_D")
		If MapExists($mState, "LM_gradient")       Then MapRemove($mState, "LM_gradient")
		If MapExists($mState, "LM_step")           Then MapRemove($mState, "LM_step")
		If MapExists($mState, "EquilibrationScale") Then MapRemove($mState, "EquilibrationScale")

		; full GN/LM solve with updated weights
		__adj_solveNonlinear($mSystem, $mState)
		If @error Then Return SetError(@error, @extended, False)
		__adj_computeResiduals($mSystem, $mState)
	Next

	; 5. Store IRLS results in state
	$mState.robustIterations = $iIterations
	$mState.robustScale = $fScale
	$mState.robustConverged = $bConverged

	; 6. Store robust weights map for results output
	Local $mRobustWeights[]
	Local $mObs2 = ($mSystem.model).obs
	For $sObsName In MapKeys($mState.idxObs)
		$mRobustWeights[$sObsName] = ($mObs2[$sObsName])._robustWeight
	Next
	$mState.robustWeights = $mRobustWeights
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_updateRobustWeights
; Description ...: Updates observation weights based on robust weight function and checks convergence
; Syntax.........: __adj_updateRobustWeights(ByRef $mSystem, ByRef $mState, $fScale, $aRedundancy)
; Parameters ....: $mSystem     - [ByRef] The adjustment system
;                  $mState      - [ByRef] Solver state (provides residuals and index maps)
;                  $fScale      - Robust scale estimate σ̂ (from __adj_robustScale)
;                  $aRedundancy - Array of redundancy values (for BIBER/ModifiedM) or Null
; Return values .: True if converged (max relative weight change < robustConvergence threshold)
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: For leverage-based estimators (BIBER, ModifiedM), the standardized residual
;                  includes √r_i: u_i = v_i / (σ̂ · σ_i · √r_i).
;                  For other estimators: u_i = v_i / (σ̂ · σ_i).
;                  Observations with near-zero redundancy are skipped (weight unchanged).
; Related .......: __adj_robustWeight, __adj_robustIRLS
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __adj_updateRobustWeights(ByRef $mSystem, ByRef $mState, $fScale, $aRedundancy)
	Local $mConfig = $mSystem.config
	Local $sEstimator = $mConfig.robust
	Local $mParams = $mConfig.robustParams
	Local $fConvergenceThreshold = $mConfig.robustConvergence
	Local $bLeverage = ($sEstimator = "BIBER" Or $sEstimator = "ModifiedM")

	Local $mModel = $mSystem.model
	Local $mObservations = $mModel.obs
	Local $mIdxObs = $mState.idxObs

	Local $fMaxChange = 0.0

	For $sObsName In MapKeys($mIdxObs)
		Local $iIdx = $mIdxObs[$sObsName]
		Local $mObs = $mObservations[$sObsName]

		; get residual
		Local $fV = __adj_vecGet($mState.r, $iIdx)

		; get a-priori standard deviation
		Local $fSigma = $mObs.stdDev

		; compute standardized residual u_i
		Local $fU
		If $bLeverage And $aRedundancy <> Null Then
			Local $fR = $aRedundancy[$iIdx]
			If $fR < 1e-10 Then
				; not controllable — skip robust weighting
				$mObs.weight = $mObs._weightOriginal
				$mObservations[$sObsName] = $mObs
				ContinueLoop
			EndIf
			$fU = $fV / ($fScale * $fSigma * Sqrt($fR))
		Else
			$fU = $fV / ($fScale * $fSigma)
		EndIf

		; compute robust weight
		Local $fW = __adj_robustWeight($fU, $sEstimator, $mParams)

		; track convergence (relative weight change — handles L1 unbounded weights gracefully)
		Local $fPrevW = MapExists($mObs, "_robustWeight") ? $mObs._robustWeight : 1.0
		Local $fChange = Abs($fW - $fPrevW) / (1.0 + Abs($fPrevW))
		If $fChange > $fMaxChange Then $fMaxChange = $fChange

		; store robust weight and update observation weight
		$mObs._robustWeight = $fW
		$mObs.weight = $mObs._weightOriginal * $fW
		$mObservations[$sObsName] = $mObs
	Next

	; write back (Map value-type!)
	$mModel.obs = $mObservations
	$mSystem.model = $mModel

	Return ($fMaxChange < $fConvergenceThreshold)
EndFunc

#EndRegion ; IRLS iteration loop

#Region Scale estimation

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_robustScale
; Description ...: Computes robust scale estimate σ̂ using the specified method
; Syntax.........: __adj_robustScale(ByRef $mState, $sMethod[, $mParams = Null])
; Parameters ....: $mState      - [ByRef] Solver state with residuals and observation data
;                  $sMethod     - Scale method: "MAD", "s0", "apriori", or "fixed"
;                  $mParams     - [optional] Robust params map (needed for "fixed": .fixedScale)
; Return values .: Scale estimate σ̂ (Float). Returns 0.0 if MAD is zero.
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: "MAD" = median(|v_i/σ_i|) / 0.6745 (most common).
;                  "s0" = √(vᵀPv/f). "apriori" = 1.0. "fixed" = user-specified value.
;                  Caller must handle MAD=0 fallback (e.g. switch to s0).
; Related .......: __adj_robustMedian, __adj_robustIRLS
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __adj_robustScale(ByRef $mState, $sMethod, $mParams = Null)
	; "apriori": sigma_hat = 1 (use a-priori stochastic model as-is, no scale estimation)
	; Used by BIBER (Wicki) and similar approaches where the stochastic model is trusted.
	If $sMethod = "apriori" Then Return 1.0

	; "fixed": use a user-specified fixed scale value (e.g. Koch uses sigma=0.5)
	If $sMethod = "fixed" Then
		If $mParams <> Null And MapExists($mParams, "fixedScale") Then Return $mParams.fixedScale
		Return 1.0
	EndIf

	If $sMethod = "s0" Then
		Local $iDOF = $mState.nObs - $mState.nParams + $mState.nRestrictions
		If $iDOF <= 0 Then Return 1.0
		Return Sqrt($mState.r2sum / $iDOF)
	EndIf

	; MAD: sigma_hat = median(|v_i / sigma_i|) / 0.6745
	Local $iN = $mState.nObs

	; Use the a-priori 1/σᵢ snapshot taken by __adj_robustIRLS — the live
	; Vector_ObsInvStdDev gets overwritten with √weight (which folds in _robustWeight)
	; on every iteration, so reading it here would doubly down-weight residuals.
	Local $tInvSd = MapExists($mState, "Vector_ObsInvStdDev_apriori") _
		? $mState.Vector_ObsInvStdDev_apriori.struct _
		: $mState.Vector_ObsInvStdDev.struct

	; element-wise |v_i| * (1/sigma_i)
	Local $mAbsStdResid = _blas_duplicate($mState.r)
	Local $tR = $mAbsStdResid.struct
	For $i = 1 To $iN
		DllStructSetData($tR, 1, Abs(DllStructGetData($tR, 1, $i) * DllStructGetData($tInvSd, 1, $i)), $i)
	Next

	Local $fMedian = __adj_robustMedian($mAbsStdResid, $iN)
	If $fMedian < 1e-15 Then Return 0.0  ; caller handles MAD=0 fallback
	Return $fMedian / 0.6745
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_robustMedian
; Description ...: Computes median of a BLAS vector using LAPACK sort
; Syntax.........: __adj_robustMedian($mVec, $iN)
; Parameters ....: $mVec        - BLAS vector map (with .struct, .rows)
;                  $iN          - Number of elements
; Return values .: Median value (Float)
; Author ........: AspirinJunkie
; ===============================================================================================================================
Func __adj_robustMedian($mVec, $iN)
	; copy vector (lasrt sorts in-place)
	Local $mCopy = _blas_duplicate($mVec)
	; sort ascending
	_lp_lasrt($mCopy, "I", 0, $iN)
	; extract median
	Local $tData = $mCopy.struct
	If Mod($iN, 2) = 1 Then
		Return DllStructGetData($tData, 1, Int($iN / 2) + 1)
	Else
		Local $fA = DllStructGetData($tData, 1, $iN / 2)
		Local $fB = DllStructGetData($tData, 1, $iN / 2 + 1)
		Return ($fA + $fB) / 2.0
	EndIf
EndFunc

#EndRegion ; Scale estimation

#Region Weight functions

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_robustWeight
; Description ...: Returns robust weight w(u) for a standardized residual u
; Syntax.........: __adj_robustWeight($fU, $sEstimator, $mParams)
; Parameters ....: $fU          - Standardized residual
;                  $sEstimator  - Estimator name: "L1", "Huber", "BIBER", "ModifiedM", "Hampel",
;                                 or "Biweight"
;                  $mParams     - Tuning parameter map (from _adj_robustDefaults)
; Return values .: Weight w(u) >= 0:
;                  - L1: 1/max(|u|, ε) capped at maxWeight (default 1000)
;                  - Huber/BIBER/ModifiedM: min(1, c/|u|)
;                  - Hampel: three-part redescending (1 → a/|u| → linear taper → 0)
;                  - Biweight: (1-(u/c)²)² for |u| ≤ c, else 0
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: BIBER and ModifiedM use the same Huber weight function — the difference is in
;                  the standardized residual computation (with/without √r_i), handled by the caller.
; Related .......: __adj_updateRobustWeights, _adj_robustDefaults
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __adj_robustWeight($fU, $sEstimator, $mParams)
	Local $fAbsU = Abs($fU)
	Switch $sEstimator
		Case "L1"
			; w = 1 / max(|u|, eps) — singularity guard, capped at maxWeight
			Local $fW_L1 = 1.0 / (($fAbsU > 1e-10) ? $fAbsU : 1e-10)
			Local $fMaxW = MapExists($mParams, "maxWeight") ? $mParams.maxWeight : 1000.0
			Return ($fW_L1 > $fMaxW) ? $fMaxW : $fW_L1

		Case "Huber", "BIBER", "ModifiedM"
			If $fAbsU <= $mParams.c Then Return 1.0
			Return $mParams.c / $fAbsU

		Case "Hampel"
			Local $fA = $mParams.a, $fB = $mParams.b, $fC = $mParams.c
			If $fAbsU <= $fA Then Return 1.0
			If $fAbsU <= $fB Then Return $fA / $fAbsU
			If $fAbsU <= $fC Then Return $fA * ($fC - $fAbsU) / (($fC - $fB) * $fAbsU)
			Return 0.0

		Case "Biweight"
			If $fAbsU > $mParams.c Then Return 0.0
			Local $fRatio = $fU / $mParams.c
			Return (1.0 - $fRatio * $fRatio) * (1.0 - $fRatio * $fRatio)
	EndSwitch
EndFunc

#EndRegion ; Weight functions

#Region Redundancy for leverage-based estimators

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_computeRedundancyDiag
; Description ...: Computes redundancy diagonal r_i from initial solution for leverage-based estimators
; Syntax.........: __adj_computeRedundancyDiag(ByRef $mSystem, ByRef $mState)
; Parameters ....: $mSystem     - [ByRef] The adjustment system
;                  $mState      - [ByRef] Solver state
; Return values .: Success      - Array of redundancy values (for BIBER/ModifiedM)
;                  Not needed   - Null (for non-leverage estimators)
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: Only computed for BIBER and ModifiedM estimators. Calls __adj_computeStatistics
;                  once to get Qxx/Qvv/redundancy, then extracts the diagonal into a simple array.
; Related .......: __adj_robustIRLS, __adj_computeStatistics
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __adj_computeRedundancyDiag(ByRef $mSystem, ByRef $mState)
	Local $mConfig = $mSystem.config
	; only needed for BIBER/ModifiedM
	If $mConfig.robust <> "BIBER" And $mConfig.robust <> "ModifiedM" Then Return Null

	; compute full statistics (Qxx, Qvv, redundancy diagonal)
	__adj_computeStatistics($mSystem, $mState)
	If @error Then Return Null

	; extract redundancy diagonal from results
	Local $mResults = $mSystem.results
	If Not MapExists($mResults, "redundancyDiag") Then Return Null

	; return as a simple AutoIt array for fast access in IRLS loop
	Local $iN = $mState.nObs
	Local $mRedVec = $mResults.redundancyDiag
	Local $aRed[$iN]
	For $i = 0 To $iN - 1
		$aRed[$i] = __adj_vecGet($mRedVec, $i)
	Next
	Return $aRed
EndFunc

#EndRegion ; Redundancy for leverage-based estimators
