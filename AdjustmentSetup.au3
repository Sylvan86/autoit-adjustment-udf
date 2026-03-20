#include-once

; #INDEX# =======================================================================================================================
; Title .........: AdjustmentSetup
; AutoIt Version : 3.3.16.1
; Description ...: Model preparation for geodetic adjustment: validation, formula parsing, index maps,
;                  model classification, weight initialization, and matrix allocation.
; Author(s) .....: AspirinJunkie
; ===============================================================================================================================

#Region Model preparation (orchestrator)

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_prepareModel
; Description ...: Orchestrates model preparation: validation, formula parsing, index building,
;                  classification, weight initialization, and matrix allocation.
; Syntax.........: __adj_prepareModel(ByRef $mSystem)
; Parameters ....: $mSystem     - [ByRef] The adjustment system
; Return values .: Success      - Implicit
;                  Failure      - SetError($ADJ_ERR_INPUT/$ADJ_ERR_DOF/$ADJ_ERR_NOT_POS_DEF)
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: On VCE re-runs (count > 1) only re-initializes weights; all other steps are skipped.
; Related .......: __adj_validateInput, __adj_parseFormulas, __adj_buildIndexMaps,
;                  __adj_classifyModel, __adj_initWeights, __adj_allocateMatrices
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __adj_prepareModel(ByRef $mSystem)
	; increment VCE-outer-loop run counter
	If Not MapExists($mSystem, "_prepareRunCount") Then
		$mSystem["_prepareRunCount"] = 1
	Else
		$mSystem["_prepareRunCount"] += 1
	EndIf

	; weights (needed every run for VKS)
	__adj_initWeights($mSystem)

	; everything below only on first run
	If $mSystem._prepareRunCount > 1 Then Return

	__adj_validateInput($mSystem)
	If @error Then Return SetError(@error, @extended, False)
	__adj_parseFormulas($mSystem)
	__adj_buildIndexMaps($mSystem)
	__adj_classifyModel($mSystem)
	If @error Then Return SetError(@error, @extended, False)
	__adj_allocateMatrices($mSystem)
	If @error Then Return SetError(@error, @extended, False)
EndFunc

#EndRegion ; Model preparation

#Region Validation and weight initialization

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_validateInput
; Description ...: Validates that the system has at least one formula or restriction
; Syntax.........: __adj_validateInput(ByRef $mSystem)
; Parameters ....: $mSystem     - [ByRef] The adjustment system
; Return values .: Failure      - SetError($ADJ_ERR_INPUT) if no formulas and no restrictions
; Author ........: AspirinJunkie
; ===============================================================================================================================
Func __adj_validateInput(ByRef $mSystem)
	Local $mModel = $mSystem.model
	If $mModel.nFormulas = 0 And $mModel.nRestrictions = 0 Then Return SetError($ADJ_ERR_INPUT, 0, False)
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_initWeights
; Description ...: Initializes variance components and observation weights
; Syntax.........: __adj_initWeights(ByRef $mSystem)
; Parameters ....: $mSystem     - [ByRef] The adjustment system
; Return values .: None (modifies $mSystem in-place)
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: Called every VKS run. First run: computes a-priori σ₀ per group as RMS of stdDevs
;                  and sets pᵢ = 1/σᵢ². Subsequent runs: weights already exist (VKS updates them).
; Related .......: __adj_prepareModel, __adj_estimateVCE
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __adj_initWeights(ByRef $mSystem)
	Local $mModel = $mSystem.model
	If $mModel._isWeighted Then
		Local $mObservations = $mModel.obs, $sObsName

		If Not MapExists($mSystem, "VarianceComponents") Then
			; first run: compute a-priori variance components σ₀,g as RMS of group standard deviations
			Local $mVarComps[], $mVarCompsCount[], $sVarComp, $mObs
			For $sObsName in MapKeys($mObservations)
				$mObs = $mObservations[$sObsName]
				$sVarComp = $mObs.varComp

				If Not MapExists($mVarComps, $sVarComp) Then
					$mVarComps[$sVarComp] = $mObs.stdDev^2
					$mVarCompsCount[$sVarComp] = 1
				Else
					$mVarComps[$sVarComp] += $mObs.stdDev^2
					$mVarCompsCount[$sVarComp] += 1
				EndIf
			Next

			For $sVarComp in MapKeys($mVarComps)
				$mVarComps[$sVarComp] = Sqrt($mVarComps[$sVarComp] / $mVarCompsCount[$sVarComp])
			Next
			$mSystem.VarianceComponents = $mVarComps
		EndIf

		; determine weights of the observations pᵢ = 1/σᵢ²
		; weights are only set once here (first run); VKS-loop updates them via σ̂²_k factor
		Local $mVarComps = $mSystem.VarianceComponents, $sVarComp, $mObs
		For $sObsName in MapKeys($mObservations)
			$mObs = $mObservations[$sObsName]
			If Not MapExists($mObs, "weight") Then
				$mObs.weight = 1 / ($mObs.stdDev)^2
				$mObservations[$sObsName] = $mObs
			EndIf
		Next
		$mModel.obs = $mObservations
		$mSystem.model = $mModel ;! MAP WRITE-BACK
	EndIf
EndFunc

#EndRegion ; Validation and weight initialization

#Region Formula processing

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_parseFormulas
; Description ...: Parses formulas and restrictions: replaces fixed parameters with values, builds
;                  executable strings for Execute(), handles observation variable substitution, and
;                  builds derivative executable strings. Sets _isNonLinear flag.
; Syntax.........: __adj_parseFormulas(ByRef $mSystem)
; Parameters ....: $mSystem     - [ByRef] The adjustment system
; Return values .: None (modifies $mSystem.model in-place)
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: Parameter names in formulas are replaced with $mPar['NAME'] map access.
;                  Observation references (#OBS) become $mObs['OBS']. Fixed parameter names are
;                  substituted with their numeric values. Analytical derivative formulas undergo
;                  the same transformation.
; Related .......: __adj_prepareModel, __adj_evalDerivative
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __adj_parseFormulas(ByRef $mSystem)
	Local $mModel = $mSystem.model
	Local $iF, $iP, $aParams, $iC, $sFormula, $fValue, $aObs, $sExecutable, $sObsName
	Local $mDerivRaw, $mDerivExec, $sDKey, $sDeriv, $sFixedKey
	Local $aFormulas     = $mModel.formulas
	Local $aRestrictions = $mModel.restrictions
	Local $mFixed        = $mModel.fixed
	Local $mObservations = $mModel.obs
	Local $bIsNonLinear  = False

	; redim formulas
	If UBound($aFormulas) <> $mModel.nFormulas Then Redim $aFormulas[$mModel.nFormulas]

	; adjust formulas
	Local $mFormula
	For $iF = 0 To UBound($aFormulas) - 1
		$mFormula    = $aFormulas[$iF]
		$aParams     = $mFormula.params
		$sFormula    = $mFormula.formula

		; set system flag non-linear if at least 1 function is nonlinear
		If (NOT $mFormula.isLinear) AND (NOT $bIsNonLinear) Then $bIsNonLinear = True

		; replace fixed parameters with their values
		$iC = 0
		For $iP = 0 To UBound($aParams) - 1
			If MapExists($mFixed, $aParams[$iP]) Then
				; replace fixed param name with value in formula
				$sFormula = StringRegExpReplace($sFormula, "(*UCP)(\$mPar\['\Q" & $aParams[$iP] & "\E'\](*ACCEPT)|(?>\$|\#)\w+(*SKIP)(?!)|\b\Q" & $aParams[$iP] & "\E\b(?!\h*\())", StringFormat("%.20g", $mFixed[$aParams[$iP]]))
			Else
				If $iP <> $iC Then $aParams[$iC] = $aParams[$iP]
				$iC += 1
			EndIf
		Next
		Redim $aParams[$iC]

		; replace parameter variables with their map call
		$sExecutable = UBound($aParams) > 0 ? StringRegExpReplace($sFormula, '(*UCP)(?>"(?>[^"]++|"")*+"|(?>[^'']++|'''')*+'')(*SKIP)(?!)|(\#\w+(*SKIP)(?!)|\b[[:alpha:]]\w*\b(?!\h*\())', "$mPar['$1']") : $sFormula

		; handle the observation value
		$fValue = $mFormula.value

		If $mModel._hasObsFunctions Or $mModel._hasObsRestrictions Then
			If IsString($fValue) Then
				$sExecutable &= " - $mObs['" & $fValue & "']"

				; add to observation array
				$aObs = $mFormula.obs
				Redim $aObs[UBound($aObs) + 1]
				$aObs[UBound($aObs) - 1] = $fValue
				$mFormula.obs = $aObs
			ElseIf $fValue <> 0.0 Then
				$sExecutable &= $fValue < 0 ? StringFormat(" + %.20g", Abs($fValue)) : StringFormat(" - %.20g", $fValue)
			EndIf

			; replace observation variables with their map call
			$sExecutable = StringRegExpReplace($sExecutable, '(*UCP)(\#(\b[[:alpha:]]\w*)\b(?!\h*\())', "$mObs['$2']")

		Else
			If IsString($fValue) Then $fValue = ($mObservations[$fValue]).value
			If $fValue <> 0.0 Then $sExecutable &= $fValue < 0 ? StringFormat(" + %.20g", Abs($fValue)) : StringFormat(" - %.20g", $fValue)
		EndIf

		$mFormula.params     = $aParams
		$mFormula.formula    = $sFormula
		$mFormula.executable = $sExecutable

		; build executable strings for analytical derivatives
		If MapExists($mFormula, "derivatives") Then
			$mDerivRaw = $mFormula.derivatives
			Dim $mDerivExec[]
			For $sDKey In MapKeys($mDerivRaw)
				$sDeriv = $mDerivRaw[$sDKey]
				; replace fixed parameters with their values
				For $sFixedKey In MapKeys($mFixed)
					$sDeriv = StringRegExpReplace($sDeriv, "(*UCP)((?>\$|\#)\w+(*SKIP)(?!)|\b\Q" & $sFixedKey & "\E\b(?!\h*\())", StringFormat("%.20g", $mFixed[$sFixedKey]))
				Next
				; replace parameter variables with their map call
				If UBound($aParams) > 0 Then $sDeriv = StringRegExpReplace($sDeriv, '(*UCP)(?>"(?>[^"]++|"")*+"|(?>[^'']++|'''')*+'')(*SKIP)(?!)|(\#\w+(*SKIP)(?!)|\b[[:alpha:]]\w*\b(?!\h*\())', "$mPar['$1']")
				; replace observation variables with their map call
				$sDeriv = StringRegExpReplace($sDeriv, '(*UCP)(\#(\b[[:alpha:]]\w*)\b(?!\h*\())', "$mObs['$2']")
				$mDerivExec[$sDKey] = $sDeriv
			Next
			$mFormula.derivativesExec = $mDerivExec
		EndIf

		$aFormulas[$iF]      = $mFormula
	Next
	$mModel.formulas = $aFormulas


	; redim restrictions
	If UBound($aRestrictions) <> $mModel.nRestrictions Then Redim $aRestrictions[$mModel.nRestrictions]

	; adjust restrictions
	Local $mRestriction
	For $iF = 0 To UBound($aRestrictions) - 1
		$mRestriction = $aRestrictions[$iF]
		$aParams      = $mRestriction.params
		$sFormula     = $mRestriction.formula

		; set system flag non-linear if at least 1 function is nonlinear
		If (NOT $mRestriction.isLinear) AND (NOT $bIsNonLinear) Then $bIsNonLinear = True

		; replace fixed parameters with their values
		$iC = 0
		For $iP = 0 To UBound($aParams) - 1
			If MapExists($mFixed, $aParams[$iP]) Then
				; replace fixed param name with value in formula
				$sFormula = StringRegExpReplace($sFormula, "(*UCP)(\$mPar\['\Q" & $aParams[$iP] & "\E'\](*ACCEPT)|(?>\$|\#)\w+(*SKIP)(?!)|\b\Q" & $aParams[$iP] & "\E\b(?!\h*\())", StringFormat("%.20g", $mFixed[$aParams[$iP]]))
			Else
				If $iP <> $iC Then $aParams[$iC] = $aParams[$iP]
				$iC += 1
			EndIf
		Next
		Redim $aParams[$iC]

		; replace parameter variables with their map call
		$sExecutable = UBound($aParams) > 0 ? StringRegExpReplace($sFormula, '(*UCP)(?>"(?>[^"]++|"")*+"|(?>[^'']++|'''')*+'')(*SKIP)(?!)|(\#\w+(*SKIP)(?!)|\b[[:alpha:]]\w*\b(?!\h*\())', "$mPar['$1']") : $sFormula

		; replace observation variables with their map call
		$sExecutable = StringRegExpReplace($sExecutable, '(*UCP)(\#(\b[[:alpha:]]\w*)\b(?!\h*\())', "$mObs['$2']")

		; handle the observation value
		$fValue = $mRestriction.value

		; convert f(X,Y) = L --> f(X,Y) - L = 0
		If IsString($fValue) Then $fValue = ($mObservations[$fValue]).value
		If $fValue <> 0.0 Then $sExecutable &= $fValue < 0 ? StringFormat(" + %.20g", Abs($fValue)) : StringFormat(" - %.20g", $fValue)

		$mRestriction.params     = $aParams
		$mRestriction.formula    = $sFormula
		$mRestriction.executable = $sExecutable

		; build executable strings for analytical derivatives
		If MapExists($mRestriction, "derivatives") Then
			$mDerivRaw = $mRestriction.derivatives
			Dim $mDerivExec[]
			For $sDKey In MapKeys($mDerivRaw)
				$sDeriv = $mDerivRaw[$sDKey]
				; replace fixed parameters with their values
				For $sFixedKey In MapKeys($mFixed)
					$sDeriv = StringRegExpReplace($sDeriv, "(*UCP)((?>\$|\#)\w+(*SKIP)(?!)|\b\Q" & $sFixedKey & "\E\b(?!\h*\())", StringFormat("%.20g", $mFixed[$sFixedKey]))
				Next
				; replace parameter variables with their map call
				If UBound($aParams) > 0 Then $sDeriv = StringRegExpReplace($sDeriv, '(*UCP)(?>"(?>[^"]++|"")*+"|(?>[^'']++|'''')*+'')(*SKIP)(?!)|(\#\w+(*SKIP)(?!)|\b[[:alpha:]]\w*\b(?!\h*\())', "$mPar['$1']")
				; replace observation variables with their map call
				$sDeriv = StringRegExpReplace($sDeriv, '(*UCP)(\#(\b[[:alpha:]]\w*)\b(?!\h*\())', "$mObs['$2']")
				$mDerivExec[$sDKey] = $sDeriv
			Next
			$mRestriction.derivativesExec = $mDerivExec
		EndIf

		$aRestrictions[$iF]      = $mRestriction
	Next
	$mModel.restrictions = $aRestrictions
	$mModel._isNonLinear = $bIsNonLinear

	$mSystem.model = $mModel ;! MAP WRITE-BACK
EndFunc

#EndRegion ; Formula processing

#Region Index maps and model classification

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_buildIndexMaps
; Description ...: Creates 0-based index maps for parameters and observations
; Syntax.........: __adj_buildIndexMaps(ByRef $mSystem)
; Parameters ....: $mSystem     - [ByRef] The adjustment system
; Return values .: None (stores idxParams and idxObs in $mSystem.model)
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: For VKS + covariances + multiple groups: sorts observations by group for
;                  block-diagonal Σₗₗ and stores group offsets in $mModel.vceGroupOffsets.
; Related .......: __adj_prepareModel, __adj_allocateMatrices
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __adj_buildIndexMaps(ByRef $mSystem)
	Local $mModel = $mSystem.model

	; create map with indices of the parameters
	Local $iIdx = 0, $mIdxParams = $mModel.params, $sKey
	For $sKey In MapKeys($mIdxParams)
		$mIdxParams[$sKey] = $iIdx
		$iIdx += 1
	Next

	; create map with indices of the observations
	Local $mIdxObs = $mModel.obs

	; if VKS + covariances + multiple groups: sort observations by VKS group
	; → makes Σₗₗ block-diagonal for efficient block operations
	If $mModel._hasCovar And $mModel._multiVarComp Then
		Local $mObservations = $mModel.obs
		Local $mGroupKeys[]  ; group → array of obs names

		; collect obs names per group
		For $sKey In MapKeys($mObservations)
			Local $sGroup = ($mObservations[$sKey]).varComp
			If Not MapExists($mGroupKeys, $sGroup) Then
				Local $aTmp[1] = [$sKey]
				$mGroupKeys[$sGroup] = $aTmp
			Else
				Local $aKeys = $mGroupKeys[$sGroup]
				ReDim $aKeys[UBound($aKeys) + 1]
				$aKeys[UBound($aKeys) - 1] = $sKey
				$mGroupKeys[$sGroup] = $aKeys
			EndIf
		Next

		; assign contiguous indices per group + store offsets
		Local $mGroupOffsets[]
		$iIdx = 0
		For $sGroup In MapKeys($mGroupKeys)
			Local $aKeys2 = $mGroupKeys[$sGroup]
			Local $mOff[]
			$mOff.start = $iIdx
			$mOff.count = UBound($aKeys2)
			$mGroupOffsets[$sGroup] = $mOff
			For $j = 0 To UBound($aKeys2) - 1
				$mIdxObs[$aKeys2[$j]] = $iIdx
				$iIdx += 1
			Next
		Next
		$mModel.vceGroupOffsets = $mGroupOffsets
	Else
		; default: sequential index assignment (original behavior)
		$iIdx = 0
		For $sKey In MapKeys($mIdxObs)
			$mIdxObs[$sKey] = $iIdx
			$iIdx += 1
		Next
	EndIf

	; store in model sub-map (will be copied to state in __adj_allocateMatrices)
	$mModel.idxParams = $mIdxParams
	$mModel.idxObs    = $mIdxObs
	$mSystem.model = $mModel ;! MAP WRITE-BACK
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_classifyModel
; Description ...: Determines adjustment model type and validates degrees of freedom
; Syntax.........: __adj_classifyModel(ByRef $mSystem)
; Parameters ....: $mSystem     - [ByRef] The adjustment system
; Return values .: Success      - Implicit (stores model type in $mModel.AdjustmentModel)
;                  Failure      - SetError($ADJ_ERR_DOF, $iF) if underdetermined (f < 0)
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: Model types: OLS/WLS/GLS (indirect), LSE/WLSE/GLSE (constrained),
;                  CLS/WCLS/GCLS (conditional), GLM/WGLM/GGLM (Gauss-Helmert).
;                  Classification based on: _hasCovar, _isWeighted, _hasObsFunctions,
;                  _hasObsRestrictions, _hasParamRestrictions, nParams.
; Related .......: __adj_prepareModel
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __adj_classifyModel(ByRef $mSystem)
	Local $mModel = $mSystem.model
; set system model type for better further handling
	; "OLS": Ordinary Least Squares (f(x) without restrictions)
	; "LSE": Linear Equality-Constrained Least Squares (f(x) with restrictions g(x))
	; "CLS": Conditional Least Squares (only restrictions g(y) - no x!) — special case of GLM
	; "GLM": General Gauss-Markov Linear Model (f(x,y) with restrictions g(x|y)) — Gauss-Helmert model
; weighted variants (with variances for the observations but not covariances)
	; "WLS": Weighted Least Squares (WLS) — Gauss-Markov model
	; "WLSE": Weighted Linear Equality-Constrained Least Squares (WLSE) — weighted variant of LSE
	; "WCLS": Weighted Conditional Least Squares (WCLS) — weighted variant of CLS
	; "WGLM": Weighted General Gauss-Markov Linear Model (WGLM) — weighted Gauss-Helmert model
; generalized variants (with variances and covariances for the observations)
	; "GLS": Generalized Least Squares (GLS) — generalized variant of OLS
	; "GLSE": Generalized Linear Equality-Constrained Least Squares (GLSE) — generalized variant of LSE
	; "GCLS": Generalized Conditional Least Squares (GCLS) — generalized variant of CLS
	; "GGLM": Generalized Gauss-Markov Linear Model (GGLM) — generalized Gauss-Helmert model
	Local $sAdjModel
	If $mModel._hasCovar Then
		; generalized variants (with covariances between observations)
		If $mModel._hasObsFunctions Or $mModel._hasObsRestrictions Then
			$sAdjModel = $mModel.nParams = 0 ? "GCLS" : "GGLM"
		Else
			$sAdjModel = $mModel._hasParamRestrictions ? "GLSE" : "GLS"
		EndIf
	ElseIf $mModel._isWeighted Then
		If $mModel._hasObsFunctions Or $mModel._hasObsRestrictions Then
			$sAdjModel = $mModel.nParams = 0 ? "WCLS" : "WGLM"
		Else
			$sAdjModel = $mModel._hasParamRestrictions ? "WLSE" : "WLS"
		EndIf
	Else
		If $mModel._hasObsFunctions Or $mModel._hasObsRestrictions Then
			$sAdjModel = $mModel.nParams = 0 ? "CLS" : "GLM"
		Else
			$sAdjModel = $mModel._hasParamRestrictions ? "LSE" : "OLS"
		EndIf
	EndIf
	$mModel.AdjustmentModel = $sAdjModel

	; check: degrees of freedom f >= 0 (underdetermined systems are not solvable)
	Local $iF_check = 0
	If StringRegExp($sAdjModel, '^[OWG]LS$') Then
		$iF_check = $mModel.nObs - $mModel.nParams
	ElseIf StringRegExp($sAdjModel, 'LSE$') Then
		$iF_check = $mModel.nObs - $mModel.nParams + $mModel.nRestrictions
	ElseIf StringRegExp($sAdjModel, 'CLS$') Then
		$iF_check = $mModel.nFormulas
	ElseIf StringRegExp($sAdjModel, 'GLM$') Then
		$iF_check = $mModel.nFormulas + $mModel.nRestrictions - $mModel.nParams
	EndIf

	$mSystem.model = $mModel ;! MAP WRITE-BACK

	If $iF_check < 0 Then Return SetError($ADJ_ERR_DOF, $iF_check, False)
EndFunc

#EndRegion ; Index maps and model classification

#Region Matrix allocation

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_allocateMatrices
; Description ...: Allocates all matrices and vectors based on model type and creates $mSystem.state
; Syntax.........: __adj_allocateMatrices(ByRef $mSystem)
; Parameters ....: $mSystem     - [ByRef] The adjustment system (must have .model and .config set)
; Return values .: Success      - Implicit (creates $mSystem.state sub-map)
;                  Failure      - SetError($ADJ_ERR_NOT_POS_DEF) if Σₗₗ not positive definite
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: Allocates model-specific matrices: A (Jacobian params), B (Jacobian obs),
;                  W (contradiction), P (weights), RA/WR (restrictions), stdDev vectors.
;                  For generalized models: builds Σₗₗ covariance matrix and Cholesky factor.
;                  Creates the .state sub-map that is passed ByRef through the solver hot-path.
; Related .......: __adj_prepareModel, __adj_classifyModel
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __adj_allocateMatrices(ByRef $mSystem)
	Local $mModel        = $mSystem.model
	Local $sAdjModel     = $mModel.AdjustmentModel
	Local $mObservations = $mModel.obs
	Local $mIdxObs       = $mModel.idxObs
	Local $iIdx, $sKey, $mObs

	; create memory for the matrices and vectors
	Local $nARows = $mModel.nFormulas
	; enhance if Levenberg-Marquardt algorithm should be used (currently only appliable for LS/LSE variants not for CLS/GLM)
	Local $mConfig = $mSystem.config
	If $mConfig.algorithm = "LM" AND StringRegExp($sAdjModel, '(OLS|WLS|GLS|LSE)$') Then $nARows += $mModel.nParams

	; ── build .state sub-map for ByRef hot-path ──
	Local $mState[]
	$mState.modelType      = $sAdjModel
	$mState.isNonlinear    = $mModel._isNonLinear
	$mState.hasCovariances = $mModel._hasCovar
	$mState.rankDeficient  = False
	$mState.nObs           = $mModel.nObs
	$mState.nParams        = $mModel.nParams
	$mState.nFormulas      = $mModel.nFormulas
	$mState.nRestrictions  = $mModel.nRestrictions
	$mState.idxObs         = $mIdxObs
	$mState.idxParams      = $mModel.idxParams

	Switch $sAdjModel
		Case "OLS", "WLS", "GLS" ; Ordinary Least Squares: f(x) = y
			$mState.Matrix_A  = _blas_createMatrix($nARows, $mModel.nParams)
			$mState.Vector_W  = _blas_createVector($nARows)

		Case "LSE", "WLSE", "GLSE" ; Linear Equality-Constrained Least Squares: f(x) = y under restrictions g(x) = z
			$mState.Matrix_RA = _blas_createMatrix($mModel.nRestrictions, $mModel.nParams)
			$mState.Matrix_A  = _blas_createMatrix($nARows, $mModel.nParams)
			$mState.Vector_W  = _blas_createVector($nARows)
			$mState.Vector_WR = _blas_createVector($mModel.nRestrictions)

		Case "CLS", "WCLS", "GCLS" ; Conditional Least Squares: g(y) = 0
		    $mState.Matrix_B = _blas_createMatrix($nARows, $mModel.nObs)
			$mState.Matrix_RB = _blas_createMatrix($mModel.nRestrictions, $mModel.nObs)
		    $mState.Vector_W  = _blas_createVector($nARows)
			$mState.Vector_WR = _blas_createVector($mModel.nRestrictions)

		Case "GLM", "WGLM", "GGLM" ; General Gauss-Markov Linear Model: f(x,y) = 0
			; formulas only in A/B/W, restrictions separate in RA/WR (null-space projection)
		    $mState.Matrix_A = _blas_createMatrix($nARows, $mModel.nParams)
		    $mState.Matrix_B = _blas_createMatrix($nARows, $mModel.nObs)
			$mState.Vector_W = _blas_createVector($nARows)
			If $mModel.nRestrictions > 0 Then
				$mState.Matrix_RA = _blas_createMatrix($mModel.nRestrictions, $mModel.nParams)
				$mState.Vector_WR = _blas_createVector($mModel.nRestrictions)
			EndIf
	EndSwitch

	; create vector of a-priori standard deviations for the observations
	Local $mStdDevs    = _blas_createVector(UBound($mIdxObs)), $tStdDevs    = $mStdDevs.struct
	Local $mInvStdDevs = _blas_createVector(UBound($mIdxObs)), $tStdDevsInv = $mInvStdDevs.struct
	For $sKey In MapKeys($mIdxObs)
		$iIdx = $mIdxObs[$sKey] ; get index
		$mObs = $mObservations[$sKey]
		DllStructSetData($tStdDevs, 1, $mObs.stdDev, $iIdx + 1)
		DllStructSetData($tStdDevsInv, 1, 1 / $mObs.stdDev, $iIdx + 1)
	Next
	$mState.Vector_ObsStdDev    = $mStdDevs
	$mState.Vector_ObsInvStdDev = $mInvStdDevs

	; create memory for the weight vector/matrix
	If $mModel._isWeighted Then
		Select
			Case StringRight($sAdjModel, 3) = "GLM"
				Local $mVecP = _blas_createVector($mModel.nObs)
				$mState.Vector_P = $mVecP
			Case StringLeft($sAdjModel, 1) = "W" ; weighted case only
				$mState.Vector_P = _blas_createVector($mModel.nObs)
			Case Else ; general case
				$mState.Matrix_P = _blas_createMatrix($mModel.nObs, $mModel.nObs)
		EndSelect
	EndIf

	; build covariance matrix Σₗₗ and Cholesky factor for generalized models (GLS/GLSE/GGLM/GCLS)
	If $mModel._hasCovar Then
		Local $iNObs = $mModel.nObs
		Local $mSigma = _blas_createMatrix($iNObs, $iNObs)
		Local $tSigma = $mSigma.struct

		; diagonal: σᵢ² from observation stdDevs
		For $sKey In MapKeys($mIdxObs)
			$iIdx = $mIdxObs[$sKey]
			$mObs = $mObservations[$sKey]
			DllStructSetData($tSigma, 1, $mObs.stdDev ^ 2, $iIdx * $iNObs + $iIdx + 1) ; column-major: (i,i) = i*N+i+1
		Next

		; off-diagonal: covariances from _adj_addCovariance() calls
		Local $mCovs = $mModel.covariances
		For $sKey In MapKeys($mCovs)
			Local $aParts = StringSplit($sKey, "|", 2) ; flag 2 = no count element
			Local $sO1 = $aParts[0], $sO2 = $aParts[1]
			If $sO1 = $sO2 Then ContinueLoop ; diagonal already set from stdDev
			Local $iI = $mIdxObs[$sO1], $iJ = $mIdxObs[$sO2]
			; column-major: element (i,j) at position j*N+i+1
			DllStructSetData($tSigma, 1, $mCovs[$sKey], $iJ * $iNObs + $iI + 1)
			DllStructSetData($tSigma, 1, $mCovs[$sKey], $iI * $iNObs + $iJ + 1) ; symmetric
		Next

		$mState.Matrix_Sigma = $mSigma

		; Cholesky factorization: Σₗₗ = L · Lᵀ
		Local $mCholeskyL = _la_duplicate($mSigma)
		_lp_potrf($mCholeskyL, "L")
		If @error Then Return SetError($ADJ_ERR_NOT_POS_DEF, @error, False) ; Σₗₗ not positive definite
		$mState.CovCholeskyL = $mCholeskyL

		; save originals for VKS block-scaling (rebuild from scratch each iteration)
		$mState.Matrix_Sigma_Original = _la_duplicate($mSigma)
		$mState.CovCholeskyL_Original = _la_duplicate($mCholeskyL)
	EndIf

	; copy VKS group offsets from model to state (if present)
	If MapExists($mModel, "vceGroupOffsets") Then $mState.vceGroupOffsets = $mModel.vceGroupOffsets

	; create vector of parameter values x0
	Local $mx0, $tx0, $ml0, $mParams = $mModel.params, $iParam, $mIdxParams = $mModel.idxParams, $mResults
	If $mModel.nParams > 0 Then
		$mx0 = _blas_createVector($mModel.nParams)
		$tx0 = $mx0.struct
		For $sKey In MapKeys($mParams)
			$iParam = $mIdxParams[$sKey] + 1
			DllStructSetData($tx0, 1, $mParams[$sKey], $iParam)
		Next
		$mResults        = $mSystem.results
		$mResults["x0"]  = $mx0
		$mSystem.results = $mResults
	EndIf

	; iteration state (nIterations intentionally NOT initialized — MapExists check used as "first run" guard)
	$mState.r2sum       = 0
	$mState.r2sum0      = 0
	If MapExists($mSystem, "VarianceComponents") Then $mState.VarianceComponents = $mSystem.VarianceComponents
	; x0 from results
	If $mModel.nParams > 0 Then $mState.Vector_x0 = $mSystem.results.x0

	$mSystem.state = $mState
EndFunc

#EndRegion ; Matrix allocation

#Region Formula analysis helpers

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_isLinear
; Description ...: Tests whether a formula string is linear in its variables
; Syntax.........: __adj_isLinear($sFunc)
; Parameters ....: $sFunc       - Formula string to analyze
; Return values .: True if linear, False if nonlinear operators/functions detected
; Author ........: AspirinJunkie
; ===============================================================================================================================
Func __adj_isLinear(Const $sFunc)
	Return Not StringRegExp($sFunc, _
		'(*UCP)(?x)' & @LF & _
		'((?(DEFINE)' & @LF & _
		'   (?<Number> -?(?>0|[1-9]\d*)(?>\.\d+)?(?>[eE][-+]?\d+)? | -? 0x[[:xdigit:]]+ )' & @LF & _
		'   (?<String> "(?> [^"]++ | "" )*+"|''(?> [^'']++ | '''' )*+'')' & @LF & _
		'   (?<CommentBlock> ^\h*\# (?>cs|comments-start)(?sU:(?>(?&CommentBlock)|.)*+)^\h*\#(?>ce|comments-end).*+\R )' & @LF & _
		'   (?<SingleComment> (?> (?P>String) (*SKIP)(*FAIL) | \;.*+$))' & @LF & _
		'   (?<Comment> (?P>SingleComment) | (?P>CommentBlock) )' & @LF & _
		'   (?<WS> \h+\_((?P>Comment)|\h)*+\R(?R) | \h+)   ' & @LF & _
		'   (?<FuncName> [[:alnum:]\_]+)' & @LF & _
		'   (?<FuncCall> \g<FuncName> \g<WS>* \( )' & @LF & _
		'   (?<Param> \b[[:alpha:]]\w*\b)' & @LF & _
		'   (?<NLOps> [*\/^])' & @LF & _
		'   (?<Variable> \$\w+)' & @LF & _
		'   (?<NumOrVar> \g<Number>|\g<Variable>) ' & @LF & _
		'   (?<NumFirst> \b\g<NumOrVar>\g<WS>*\*\g<WS>* (?!\g<FuncCall>))' & @LF & _
		'   (?<NumLast> \#?\g<Param>\g<WS>*[\*\/]\g<WS>(?=\g<NumOrVar>)) ' & @LF & _
		')' & @LF & _
		'(?: \g<NumFirst>' & @LF & _
		'   |\g<NumLast> (?!\g<NumFirst>)  ' & @LF & _
		')(?!\g<WS>\g<NLOps>) ' & @LF & _
		'(*SKIP)(?!)' & @LF & _
		'|\g<NLOps>' & @LF & _
		'|\g<FuncCall>' & @LF & _
		')', 0)
EndFunc

#EndRegion ; Formula analysis helpers
