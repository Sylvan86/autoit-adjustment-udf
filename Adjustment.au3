#include-once

; #CONSTANTS# ===================================================================================================================
Global Const $ADJ_ERR_OK             = 0
Global Const $ADJ_ERR_SOLVER         = 1
Global Const $ADJ_ERR_DOF            = 2
Global Const $ADJ_ERR_RANK           = 3
Global Const $ADJ_ERR_NO_CONVERGENCE = 4
Global Const $ADJ_ERR_NOT_POS_DEF    = 5
Global Const $ADJ_ERR_INPUT          = 6
Global Const $ADJ_ERR_RESTR_OVERDETERMINED = 7

Global Const $ADJ_SOLVER_QR  = "QR"
Global Const $ADJ_SOLVER_SVD = "SVD"

Global Const $__ADJ_VCE_MIN_SIGMA2   = 1e-6
Global Const $__ADJ_LM_TAU_DEFAULT   = 1e-2
Global Const $__ADJ_LM_LAMBDA_MIN    = 1e-7
Global Const $__ADJ_LM_LAMBDA_MAX    = 1e7
Global Const $__ADJ_VCE_MAX_ITER     = 15
Global Const $__ADJ_VCE_CONVERGENCE  = 0.05
Global Const $__ADJ_MAX_ITERATIONS   = 100
Global Const $__ADJ_TOLERANCE            = 1e-10
Global Const $__ADJ_STAGNATION_TOL       = 1e-10
Global Const $__ADJ_MARQUARDT_REL_FLOOR  = 1e-10
Global Const $__ADJ_MARQUARDT_ABS_FLOOR  = 1e-30
; ===============================================================================================================================

#include "LinearAlgebra.au3"
#include "AdjustmentSetup.au3"
#include "AdjustmentSolver.au3"
#include "AdjustmentStats.au3"
#include "AdjustmentRobust.au3"
#include "AdjustmentDiagnostics.au3"
#include "AdjustmentDisplay.au3"

; #INDEX# =======================================================================================================================
; Title .........: Adjustment
; AutoIt Version : 3.3.18.0
; Description ...: Geodetic adjustment (least squares) library supporting OLS/WLS/GLS, LSE, GLM, CLS
;                  and their generalized variants. Includes Gauss-Newton/Levenberg-Marquardt solvers,
;                  variance component estimation, robust estimation (IRLS), and statistical diagnostics.
; Author(s) .....: AspirinJunkie
; Dll ...........: libopenblas.dll, ucrtbase.dll
; ===============================================================================================================================

; #CURRENT# =====================================================================================================================
;_adj_addCovariance
;_adj_addFixedParam
;_adj_addFunction
;_adj_addObs
;_adj_addObsFunction
;_adj_addRestriction
;_adj_createSystem
;_adj_defaultConfig
;_adj_defaultDisplayConfig
;_adj_displayResults
;_adj_getErrorMessage
;_adj_getOutliers
;_adj_getResults
;_adj_removeObs
;_adj_robustDefaults
;_adj_setInitialValue
;_adj_solve
; ===============================================================================================================================

; #INTERNAL_USE_ONLY# ===========================================================================================================
;__adj_addDiagVariance
;__adj_allocateMatrices
;__adj_applyEquilibration
;__adj_applyMarquardtFloorAndSqrt
;__adj_applyWhitening
;__adj_betacf
;__adj_betaRI
;__adj_buildIndexMaps
;__adj_chi2Quantile
;__adj_classifyModel
;__adj_computeCofactors
;__adj_computeContradiction
;__adj_computeCorrelation
;__adj_computeDiagnostics
;__adj_computeDOF
;__adj_computeGlobalTest
;__adj_computeJacobians
;__adj_computeQxx
;__adj_computeQxxFromSVD
;__adj_computeRedundancy
;__adj_computeRedundancyDiag
;__adj_computeResiduals
;__adj_computeStatistics
;__adj_computeVtPV
;__adj_createNewSystem
;__adj_dataSnooping
;__adj_derivate1D
;__adj_dispatchLinearSolver
;__adj_displayColLabel
;__adj_displayEnsureDependencies
;__adj_displayGetValue
;__adj_displayRequiredCompute
;__adj_displayTable
;__adj_ensureComputed
;__adj_erfc
;__adj_estimateVCE
;__adj_evalDerivative
;__adj_fillLMaugmentation
;__adj_fillLowerFromUpper
;__adj_gammaLn
;__adj_gammp
;__adj_gammpApprox
;__adj_gcf
;__adj_gser
;__adj_initLMDamping
;__adj_initWeights
;__adj_inverfc
;__adj_invgammp
;__adj_isLinear
;__adj_normCdf
;__adj_normQuantile
;__adj_popeCdf
;__adj_popeQuantile
;__adj_parseDerivativeInput
;__adj_tCdf
;__adj_parseFormulas
;__adj_prepareModel
;__adj_reverseEquilibration
;__adj_robustIRLS
;__adj_runIRLSPhase
;__adj_robustMedian
;__adj_robustScale
;__adj_robustWeight
;__adj_scaleSymmetricMatrix
;__adj_scaleSymmetricMatrixCholesky
;__adj_solveCLS
;__adj_solveGLM
;__adj_solveIteration
;__adj_solveLSE
;__adj_solveNonlinear
;__adj_solveOLS
;__adj_updateLMDamping
;__adj_updateMarquardtD
;__adj_updateObservations
;__adj_updateParameters
;__adj_updateRobustWeights
;__adj_updateWhiteningVectors
;__adj_validateInput
;__adj_vecGet
; ===============================================================================================================================


; ===============================================================================================================================
; Solver functions - three-tier architecture (→ AdjustmentSolver.au3):
;   Tier 1: _adj_solve()                        - Entry point: prepares system, calls VCE loop
;   Tier 2: __adj_estimateVCE()                  - VCE outer loop (Helmert variance component estimation)
;           __adj_solveNonlinear()                - GN/LM iteration loop
;   Tier 3: __adj_solveIteration()               - Single linear iteration step
;           __adj_computeJacobians()              - Jacobian matrix computation
;           __adj_computeContradiction()          - Contradiction vector computation
;           __adj_applyWhitening()                - Weight transformation (Whitening)
;           __adj_dispatchLinearSolver()           - LAPACK solver dispatch
;           __adj_computeResiduals()               - Residual computation
;           __adj_initLMDamping()                  - LM λ₀ initialization
;           __adj_updateParameters()               - Parameter update x = x₀ + Δx
;           __adj_updateObservations()             - Observation update (stub for GLM)
;           __adj_updateLMDamping()                - LM λ adjustment via gain ratio
; ===============================================================================================================================
#Region System creation

; #FUNCTION# ====================================================================================================================
; Name...........: _adj_createSystem
; Description ...: Creates a new empty adjustment system
; Syntax.........: _adj_createSystem()
; Parameters ....: None
; Return values .: Map — empty adjustment system with .model and .results sub-maps
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: The returned system is populated via _adj_addObs, _adj_addObsFunction, etc.
; Related .......: _adj_addObs, _adj_addObsFunction, _adj_addFunction, _adj_solve
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func _adj_createSystem()
	Return __adj_createNewSystem()
EndFunc

#EndRegion ; System creation

#Region Model definition

; #FUNCTION# ====================================================================================================================
; Name...........: _adj_addObs
; Description ...: Adds an observation to the adjustment system
; Syntax.........: _adj_addObs(ByRef $mSystem[, $sSymbol = Default[, $fValue = 0[, $fStdDev = 1[, $sVarComp = "s0"]]]])
; Parameters ....: $mSystem     - [ByRef] The adjustment system (created if not a Map)
;                  $sSymbol     - [optional] Observation name. Default auto-generates "O1", "O2", ...
;                  $fValue      - [optional] Observed value. Default is 0.
;                  $fStdDev     - [optional] A-priori standard deviation. Default is 1.
;                  $sVarComp    - [optional] Variance component group name. Default is "s0".
; Return values .: None (modifies $mSystem in-place)
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: Symbol names are converted to uppercase. If stdDev ≠ 1, the system is flagged
;                  as weighted. If varComp ≠ "s0", VCE (variance component estimation) is enabled.
; Related .......: _adj_addObsFunction, _adj_addCovariance, _adj_createSystem
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func _adj_addObs(ByRef $mSystem, $sSymbol = Default, Const $fValue = 0, Const $fStdDev = 1, Const $sVarComp = "s0")
		; create empty structure for adjustment system if not exist
	If Not IsMap($mSystem) Then $mSystem = __adj_createNewSystem()

	; validate finite numeric input
	If Not IsNumber($fValue) Or $fValue <> $fValue Or $fValue = 1.0 / 0 Or $fValue = -1.0 / 0 Then Return SetError($ADJ_ERR_INPUT, 10, False)
	If Not IsNumber($fStdDev) Or $fStdDev <= 0 Or $fStdDev <> $fStdDev Or $fStdDev = 1.0 / 0 Then Return SetError($ADJ_ERR_INPUT, 11, False)

	; build observation object
	Local $mObs[]
	$mObs.value   = $fValue
	$mObs.stdDev  = $fStdDev
	$mObs.varComp = $sVarComp

	Local $mModel = $mSystem.model

	; flag that a variance component estimate should be performed
	If Not $mModel._multiVarComp And $sVarComp <> "s0" Then $mModel._multiVarComp = True

	; if (at least) one observation has a std deviation <> 1 set system attribute to "weighted"
	If $fStdDev <> 1 Then $mModel._isWeighted = True

	Local $mObsGlobal = $mModel.obs

	; create symbol if not given
	If IsKeyword($sSymbol) = 1 Then
		$sSymbol = "O" & $mModel._lastOID
		$mModel._lastOID += 1
	Else
		$sSymbol = StringUpper($sSymbol)
	EndIf

	; add observation
	$mObsGlobal[$sSymbol] = $mObs
	$mModel.obs            = $mObsGlobal
	$mModel.nObs          += 1

	$mSystem.model = $mModel ;! MAP WRITE-BACK
EndFunc

; #FUNCTION# ====================================================================================================================
; Name...........: _adj_addObsFunction
; Description ...: Adds an observation with a functional model relating observations and parameters
; Syntax.........: _adj_addObsFunction(ByRef $mSystem[, $sSymbol = Default[, $sFormula = ""[, $fValue = 0[, $fStdDev = 1[, $sVarComp = "s0"[, $fParamInitValue = 0.0[, $vDerivatives = Default]]]]]]])
; Parameters ....: $mSystem         - [ByRef] The adjustment system
;                  $sSymbol         - [optional] Observation name. Default auto-generates.
;                  $sFormula        - [optional] Formula string (e.g. "sqrt((X-x1)^2+(Y-y1)^2)")
;                  $fValue          - [optional] Observed value. Default is 0.
;                  $fStdDev         - [optional] A-priori standard deviation. Default is 1.
;                  $sVarComp        - [optional] Variance component group. Default is "s0".
;                  $fParamInitValue - [optional] Initial value for new parameters. Default is 0.0.
;                  $vDerivatives    - [optional] Analytical derivatives: Map or pipe-delimited string
;                                     (e.g. "X=2*X | Y=2*Y"). Default uses numerical differentiation.
; Return values .: None (modifies $mSystem in-place)
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: Parameters are extracted from the formula via regex. Observation references use
;                  # prefix (e.g. #D1). Creates both observation and formula entries.
; Related .......: _adj_addObs, _adj_addFunction, _adj_addRestriction
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func _adj_addObsFunction(ByRef $mSystem, $sSymbol = Default, $sFormula = "", Const $fValue = 0, Const $fStdDev = 1, Const $sVarComp = "s0", Const $fParamInitValue = 0.0, $vDerivatives = Default)
	; create empty structure for adjustment system if not exist
	If Not IsMap($mSystem) Then $mSystem = __adj_createNewSystem()

	; validate finite numeric input
	If Not IsNumber($fValue) Or $fValue <> $fValue Or $fValue = 1.0 / 0 Or $fValue = -1.0 / 0 Then Return SetError($ADJ_ERR_INPUT, 10, False)
	If Not IsNumber($fStdDev) Or $fStdDev <= 0 Or $fStdDev <> $fStdDev Or $fStdDev = 1.0 / 0 Then Return SetError($ADJ_ERR_INPUT, 11, False)

	; build observation object
	Local $mObs[]
	$mObs.value   = $fValue
	$mObs.stdDev  = $fStdDev
	$mObs.varComp = $sVarComp

	Local $mModel = $mSystem.model

	; flag that a variance component estimate should be performed
	If Not $mModel._multiVarComp And $sVarComp <> "s0" Then $mModel._multiVarComp = True

	; if (at least) one observation has a std deviation <> 1 set system attribute to "weighted"
	If $fStdDev <> 1 Then $mModel._isWeighted = True

	$sFormula = StringUpper($sFormula) ; because Map-Keys are case sensitive

	; extract unknown params and assign their approximation value
	Local $aParams = StringRegExp($sFormula, '(*UCP)(?>"(?>[^"]++|"")*+"|(?>[^'']++|'''')*+'')(*SKIP)(?!)|(\#\w+(*SKIP)(?!)|\b[[:alpha:]]\w*\b(?!\h*\())', 3)
	If UBound($aParams) < 1 Or @error Then
		Dim $aParams[0]
	Else
		$mModel._hasParamFunctions = True
		$aParams             = _ArrayUnique($aParams, 0, 0, 1, 0)
	EndIf

	; add parameters to global list of parameters
	Local $sParam, $mParams = $mModel.params, $bChanged = False
	For $sParam In $aParams
		If Not MapExists($mParams, $sParam) Then
			$mParams[$sParam] = $fParamInitValue
			$bChanged         = True
			$mModel.nParams += 1
		EndIf
	Next
	If $bChanged Then $mModel.params = $mParams

	; extract observations inside the formula
	Local $aObs = StringRegExp($sFormula, '(*UCP)(?x)(?>"(?>[^"]++|"")*+"|(?>[^'']++|'''')*+'')(*SKIP)(?!)|\#\K\b[[:alpha:]]\w*\b(?!\h*\()', 3)
	If UBound($aObs) < 1 Or @error Then
		Dim $aObs[0]
	Else
		$mModel._hasObsFunctions = True
		$aObs                            = _ArrayUnique($aObs, 0, 0, 1, 0)
	EndIf

	Local $mObsGlobal = $mModel.obs

	; create symbol if not given
	If IsKeyword($sSymbol) = 1 Then
		$sSymbol = "O" & $mModel._lastOID
		$mModel._lastOID += 1
	Else
		$sSymbol = StringUpper($sSymbol)
	EndIf

	; add observation
	$mObsGlobal[$sSymbol] = $mObs
	$mModel.obs            = $mObsGlobal

	; create formula object
	Local $mFormula[]
	$mFormula.formula    = $sFormula
	$mFormula.params     = $aParams
	$mFormula.obs        = $aObs
	$mFormula.isLinear   = __adj_isLinear($sFormula)
	$mFormula.value      = IsString($sSymbol) ? StringUpper($sSymbol) : $sSymbol

	; parse analytical derivatives
	Local $mDerivParsed = __adj_parseDerivativeInput($vDerivatives)
	If IsMap($mDerivParsed) Then $mFormula.derivatives = $mDerivParsed

	; add formula
	Local $aFormulas = $mModel.formulas
	$mModel.nFormulas     += 1
	If UBound($aFormulas, 1) = 0 Then Redim $aFormulas[1]
	If UBound($aFormulas, 1) < $mModel.nFormulas Then Redim $aFormulas[2 * UBound($aFormulas, 1)]
	$aFormulas[$mModel.nFormulas - 1] = $mFormula
	$mModel.formulas = $aFormulas

	; increment counters
	$mModel.nObs          += 1

	$mSystem.model = $mModel ;! MAP WRITE-BACK
EndFunc

; #FUNCTION# ====================================================================================================================
; Name...........: _adj_addFunction
; Description ...: Adds a condition equation or restriction without an associated observation
; Syntax.........: _adj_addFunction(ByRef $mSystem, $sFormula[, $fValue = 0[, $fParamInitValue = 0.0[, $vDerivatives = Default[, $bRestriction = False]]]])
; Parameters ....: $mSystem         - [ByRef] The adjustment system
;                  $sFormula        - Formula string (e.g. "(#D1*cos(#T1)-XM)^2 + ...")
;                  $fValue          - [optional] Target value of the equation. Default is 0.
;                  $fParamInitValue - [optional] Initial value for new parameters. Default is 0.0.
;                  $vDerivatives    - [optional] Analytical derivatives (Map or pipe string).
;                  $bRestriction    - [optional] If True, adds as restriction instead of formula.
; Return values .: None (modifies $mSystem in-place)
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: For GLM/CLS models: formulas may contain both parameters and observation references.
;                  Called internally by _adj_addRestriction with $bRestriction = True.
; Related .......: _adj_addObsFunction, _adj_addRestriction
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func _adj_addFunction(ByRef $mSystem, $sFormula, Const $fValue = 0, Const $fParamInitValue = 0.0, $vDerivatives = Default, $bRestriction = False)
	; create empty structure for adjustment system if not exist
	If Not IsMap($mSystem) Then $mSystem = __adj_createNewSystem()

	Local $mModel = $mSystem.model

	$sFormula = StringUpper($sFormula) ; because Map-Keys are case sensitive

	; extract unknown params and assign their approximation value
	Local $aParams = StringRegExp($sFormula, '(*UCP)(?>"(?>[^"]++|"")*+"|(?>[^'']++|'''')*+'')(*SKIP)(?!)|(\#\w+(*SKIP)(?!)|\b[[:alpha:]]\w*\b(?!\h*\())', 3)
	If UBound($aParams) < 1 Or @error Then
		Dim $aParams[0]
	Else
		If $bRestriction Then
			$mModel._hasParamRestrictions = True
		Else
			$mModel._hasParamFunctions = True
		EndIf
		$aParams = _ArrayUnique($aParams, 0, 0, 1, 0)
	EndIf

	; add parameters to global list of parameters
	Local $sParam, $mParams = $mModel.params, $bChanged = False
	For $sParam In $aParams
		If Not MapExists($mParams, $sParam) Then
			$mParams[$sParam] = $fParamInitValue
			$bChanged         = True
			$mModel.nParams += 1
		EndIf
	Next
	If $bChanged Then $mModel.params = $mParams

	; extract observations inside the formula
	Local $aObs = StringRegExp($sFormula, '(*UCP)(?x)(?>"(?>[^"]++|"")*+"|(?>[^'']++|'''')*+'')(*SKIP)(?!)|\#\K\b[[:alpha:]]\w*\b(?!\h*\()', 3)
	If UBound($aObs) < 1 Or @error Then
		Dim $aObs[0]
	Else
		If $bRestriction Then
			$mModel._hasObsRestrictions = True
		Else
			$mModel._hasObsFunctions = True
		EndIf
		$aObs = _ArrayUnique($aObs, 0, 0, 1, 0)
	EndIf

	; create formula object
	Local $mFormula[]
	$mFormula.formula    = $sFormula
	$mFormula.params     = $aParams
	$mFormula.obs        = $aObs
	$mFormula.isLinear   = __adj_isLinear($sFormula)
	$mFormula.value      = IsString($fValue) ? StringUpper($fValue) : $fValue

	; parse analytical derivatives
	Local $mDerivParsed = __adj_parseDerivativeInput($vDerivatives)
	If IsMap($mDerivParsed) Then $mFormula.derivatives = $mDerivParsed

	; add formula/restriction
	If $bRestriction Then
		Local $aRestrictions             = $mModel.restrictions
		If Not IsArray($aRestrictions) Then Dim $aRestrictions[0]
		$mModel.nRestrictions          += 1
		If UBound($aRestrictions, 1)     = 0 Then Redim $aRestrictions[1]
		If UBound($aRestrictions, 1) < $mModel.nRestrictions Then Redim $aRestrictions[2 * UBound($aRestrictions, 1)]
		$aRestrictions[$mModel.nRestrictions - 1] = $mFormula
		$mModel.restrictions = $aRestrictions
	Else
		Local $aFormulas             = $mModel.formulas
		$mModel.nFormulas          += 1
		If UBound($aFormulas, 1)     = 0 Then Redim $aFormulas[1]
		If UBound($aFormulas, 1) < $mModel.nFormulas Then Redim $aFormulas[2 * UBound($aFormulas, 1)]
		$aFormulas[$mModel.nFormulas - 1] = $mFormula
		$mModel.formulas = $aFormulas
	EndIf

	$mSystem.model = $mModel ;! MAP WRITE-BACK

EndFunc

; #FUNCTION# ====================================================================================================================
; Name...........: _adj_addRestriction
; Description ...: Adds a restriction (constraint) to the adjustment system
; Syntax.........: _adj_addRestriction(ByRef $mSystem, $sFormula[, $fValue = 0[, $fParamInitValue = 0.0[, $vDerivatives = Default]]])
; Parameters ....: $mSystem         - [ByRef] The adjustment system
;                  $sFormula        - Restriction formula (e.g. "X + Y - 100")
;                  $fValue          - [optional] Target value. Default is 0.
;                  $fParamInitValue - [optional] Initial value for new parameters. Default is 0.0.
;                  $vDerivatives    - [optional] Analytical derivatives (Map or pipe string).
; Return values .: None (modifies $mSystem in-place)
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: Auto-converts single-parameter restrictions (e.g. "X") to _adj_addFixedParam.
;                  Mixed restrictions (parameters + observations) handled correctly in GLM.
; Related .......: _adj_addFixedParam, _adj_addFunction
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func _adj_addRestriction(ByRef $mSystem, $sFormula, Const $fValue = 0, Const $fParamInitValue = 0.0, $vDerivatives = Default)
	; Auto-convert: if formula is just a single parameter name (optionally negated),
	; treat as _adj_addFixedParam instead (mathematically equivalent, simpler handling)
	Local $aMatch = StringRegExp(StringStripWS($sFormula, 3), '(?i)^(-?)([[:alpha:]]\w*)$', 1)
	If Not @error Then
		; single parameter name detected — convert to fixed parameter
		Local $fFixedValue = $aMatch[0] = "-" ? -$fValue : $fValue
		_adj_addFixedParam($mSystem, $aMatch[1], $fFixedValue)
		Return
	EndIf

	; NOTE: Mixed restrictions (parameters + observations in one formula) are correctly handled in GLM.
	; For separate pure parameter and observation restrictions in the same system, pure parameter
	; restrictions are approximated via the GLM pseudo-observation mechanism (weight 10^12) instead
	; of the exact LSE mechanism — sufficiently accurate in practice.

	_adj_addFunction($mSystem, $sFormula, $fValue, $fParamInitValue, $vDerivatives, True)
EndFunc

; #FUNCTION# ====================================================================================================================
; Name...........: _adj_addFixedParam
; Description ...: Fixes a parameter to a known value (removes from adjustment)
; Syntax.........: _adj_addFixedParam(ByRef $mSystem, $sParam, $fValue)
; Parameters ....: $mSystem     - [ByRef] The adjustment system
;                  $sParam      - Parameter name (case-insensitive, stored uppercase)
;                  $fValue      - Fixed value for the parameter
; Return values .: None (modifies $mSystem in-place)
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: The parameter is removed from the adjustable list and stored in .model.fixed.
;                  Formula substitution happens in __adj_parseFormulas during _adj_solve.
; Related .......: _adj_addRestriction, _adj_setInitialValue
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func _adj_addFixedParam(ByRef $mSystem, $sParam, $fValue)
	; create empty structure for adjustment system if not exist
	If Not IsMap($mSystem) Then $mSystem = __adj_createNewSystem()

	$sParam = StringUpper($sParam) ; because Map-Keys are case sensitive

	Local $mModel = $mSystem.model

	; add to fixed
	Local $mFixed = $mModel.fixed
	$mFixed[$sParam] = $fValue
	$mModel.fixed = $mFixed

	; remove from adjustable parameter list
	Local $mParams = $mModel.params
	If MapRemove($mParams, $sParam) Then $mModel.nParams -= 1
	$mModel.params = $mParams

	$mSystem.model = $mModel ;! MAP WRITE-BACK

	; replace with value in the formulas is the job of __adj_prepareModel()
EndFunc

; #FUNCTION# ====================================================================================================================
; Name...........: _adj_addCovariance
; Description ...: Adds a covariance entry between two observations (or overrides a variance)
; Syntax.........: _adj_addCovariance(ByRef $mSystem, $sObs1, $sObs2, $fCovariance)
; Parameters ....: $mSystem     - [ByRef] The adjustment system (observations must already exist)
;                  $sObs1       - Name of the first observation
;                  $sObs2       - Name of the second observation (= $sObs1 for variance override)
;                  $fCovariance - Covariance value σ₁₂ (for $sObs1 = $sObs2: variance σ²)
; Return values .: Success      - True
;                  Failure      - False, @error = $ADJ_ERR_INPUT:
;                               |@extended = 0: system not initialized
;                               |@extended = 2: $sObs1 not found
;                               |@extended = 3: $sObs2 not found
;                               |@extended = 4: cross-group covariance (incompatible with VKS)
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: For $sObs1 = $sObs2: updates the observation's stdDev to √($fCovariance).
;                  Off-diagonal covariances trigger generalized model (GLS/GLSE/GGLM/GCLS).
; Related .......: _adj_addObs
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func _adj_addCovariance(ByRef $mSystem, $sObs1, $sObs2, Const $fCovariance)
	If Not IsMap($mSystem) Then Return SetError($ADJ_ERR_INPUT, 0, False)

	$sObs1 = StringUpper($sObs1)
	$sObs2 = StringUpper($sObs2)

	Local $mModel = $mSystem.model

	; check that both observations exist
	Local $mObservations = $mModel.obs
	If Not MapExists($mObservations, $sObs1) Then Return SetError($ADJ_ERR_INPUT, 2, False)
	If Not MapExists($mObservations, $sObs2) Then Return SetError($ADJ_ERR_INPUT, 3, False)

	; for self-covariance (variance): update the observation's stdDev
	If $sObs1 = $sObs2 Then
		Local $mObs = $mObservations[$sObs1]
		$mObs.stdDev = Sqrt($fCovariance)
		$mObservations[$sObs1] = $mObs
		$mModel.obs = $mObservations
		If Sqrt($fCovariance) <> 1 Then $mModel._isWeighted = True
	Else
		; off-diagonal covariance → generalized model
		$mModel._hasCovar   = True
		$mModel._isWeighted = True
	EndIf

	; check for cross-group covariances (incompatible with VKS)
	If $sObs1 <> $sObs2 Then
		Local $sVC1 = ($mObservations[$sObs1]).varComp
		Local $sVC2 = ($mObservations[$sObs2]).varComp
		If $sVC1 <> $sVC2 Then
			$mModel._hasCrossGroupCovar = True
			$mSystem.model = $mModel ;! MAP WRITE-BACK
			Return SetError($ADJ_ERR_INPUT, 4, False)
		EndIf
	EndIf

	; store covariance with normalized key (sorted names for uniqueness)
	Local $mCov = $mModel.covariances
	Local $sKey = ($sObs1 <= $sObs2) ? ($sObs1 & "|" & $sObs2) : ($sObs2 & "|" & $sObs1)
	$mCov[$sKey] = $fCovariance
	$mModel.covariances = $mCov

	$mSystem.model = $mModel ;! MAP WRITE-BACK

	Return True
EndFunc

; #FUNCTION# ====================================================================================================================
; Name...........: _adj_setInitialValue
; Description ...: Sets the initial (approximate) value for a parameter
; Syntax.........: _adj_setInitialValue(ByRef $mSystem, $sParam, $fValue)
; Parameters ....: $mSystem     - [ByRef] The adjustment system
;                  $sParam      - Parameter name (case-insensitive)
;                  $fValue      - Initial/approximate value
; Return values .: Success      - True
;                  Failure      - SetError($ADJ_ERR_INPUT) if parameter not found
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: Required for nonlinear models where parameters need starting values for iteration.
; Related .......: _adj_addObsFunction, _adj_addFunction
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func _adj_setInitialValue(ByRef $mSystem, $sParam, $fValue)
	Local $mModel = $mSystem.model
	Local $mParams = $mModel.params
	If Not MapExists($mParams, StringUpper($sParam)) Then Return SetError($ADJ_ERR_INPUT, 0, False)
	$mParams[StringUpper($sParam)] = $fValue
	$mModel.params = $mParams
	$mSystem.model = $mModel ;! MAP WRITE-BACK

	Return True
EndFunc

#EndRegion ; Model definition

#Region Configuration

; #FUNCTION# ====================================================================================================================
; Name...........: _adj_defaultConfig
; Description ...: Creates a default configuration map for _adj_solve
; Syntax.........: _adj_defaultConfig([$sAlgorithm = "GN"[, $bVCE = True[, $iMaxIter = 100[, $fTolerance = 1e-10[, $fAlpha = 0.05[, $sSolver = "QR"[, $fSolverRCOND = Default]]]]]]])
; Parameters ....: $sAlgorithm  - [optional] "GN" (Gauss-Newton) or "LM" (Levenberg-Marquardt)
;                  $bVCE        - [optional] Enable Variance Component Estimation. Default True.
;                  $iMaxIter    - [optional] Maximum iterations. Default 100.
;                  $fTolerance  - [optional] Convergence tolerance ‖Δx‖. Default 1e-10.
;                  $fAlpha      - [optional] Significance level for global test. Default 0.05.
;                  $sSolver     - [optional] "QR" (DGELSY) or "SVD" (DGELSD). Default "QR".
;                  $fSolverRCOND - [optional] Rank threshold for solver. Default 1e-5.
; Return values .: Map with solver configuration including:
;                  - .compute sub-map: .qxx (True), .cofactors, .redundancy, .globalTest,
;                  + .diagnostics (all False by default)
;                  - .diagnostics sub-map: .alpha (0.001), .beta (0.20), .alphaSuspect,
;                  + .testBasis ("pope")
;                  - .robust (""), .robustParams, .robustMaxIter (30), .robustConvergence (1e-3)
;                  - .deriveMethod ("Central"), .scaling (True)
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: Compute flags control what statistics are calculated after solving.
;                  Set .robust to estimator name and .robustParams for IRLS estimation.
;                  Available deriveMethod values: "Central".."Central4", "Forward", "Backward",
;                  "Ridder", "Higham".
;                  .lmTau (default 1e-2) controls the LM initial damping λ₀ = lmTau · max(diag(JᵀPJ)).
;                  Increase (e.g. 1.0) when starting values are far from the optimum;
;                  decrease (e.g. 1e-3) only when the linearisation is already very good.
; Related .......: _adj_solve, _adj_robustDefaults
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func _adj_defaultConfig($sAlgorithm = "LM", $bVCE = False, $iMaxIter = $__ADJ_MAX_ITERATIONS, _
		$fTolerance = $__ADJ_TOLERANCE, $fAlpha = 0.05, _
		$sSolver = "QR", $fSolverRCOND = Default)
	; resolve explicit Default keywords to actual defaults
	If IsKeyword($iMaxIter) = 1 Then $iMaxIter = $__ADJ_MAX_ITERATIONS
	If IsKeyword($fTolerance) = 1 Then $fTolerance = $__ADJ_TOLERANCE
	If IsKeyword($fAlpha) = 1 Then $fAlpha = 0.05

	Local $mCfg[]
	$mCfg.algorithm      = $sAlgorithm
	$mCfg.vce            = $bVCE
	$mCfg.maxIterations  = $iMaxIter
	$mCfg.tolerance      = $fTolerance
	$mCfg.alpha          = $fAlpha
	$mCfg.lmTau          = $__ADJ_LM_TAU_DEFAULT
	$mCfg.lmLambdaMin    = $__ADJ_LM_LAMBDA_MIN
	$mCfg.lmLambdaMax    = $__ADJ_LM_LAMBDA_MAX
	$mCfg.vceMaxIter     = $__ADJ_VCE_MAX_ITER
	$mCfg.vceConvergence = $__ADJ_VCE_CONVERGENCE
	$mCfg.deriveMethod   = "Central"
	$mCfg.solver         = $sSolver
	$mCfg.solverRCOND    = $fSolverRCOND
	$mCfg.scaling        = True                      ; Jacobi-Equilibration (column scaling)
	$mCfg.robust            = ""
	$mCfg.robustParams      = Null
	$mCfg.robustMaxIter     = 30
	$mCfg.robustConvergence = 1e-3
	; inner maxIterations for the GN loop during IRLS re-solves.  Since x* barely
	; moves between IRLS iterations (only weights change), a warm-started GN
	; normally converges in 1–2 inner steps; allowing the full maxIterations (100)
	; lets pathological cases (Biweight with re-descending weights) spin for the
	; full cap before stopping. 10 is a safe upper bound.
	$mCfg.robustInnerMaxIter = 10
	; cascade-init (MM-estimator pattern): run Huber first, then redescending. Default ON
	; for non-convex estimators (Biweight, Hampel, BIBER) since Huber's convex global minimum
	; provides a much better starting point than the non-robust solve alone.
	$mCfg.robustCascadeInit    = True
	$mCfg.robustCascadeMaxIter = 15

	; iterative data snooping (mutually exclusive with robust)
	$mCfg.dataSnooping              = False
	$mCfg.dataSnoopingMethod        = "Pope"   ; "Pope" (a posteriori, default) or "Baarda" (a priori σ₀=1)
	$mCfg.dataSnoopingAlpha         = 0.001
	$mCfg.dataSnoopingMaxIter       = 30
	$mCfg.dataSnoopingMinRedundancy = 0.1      ; Baarda's controllability threshold; below = untestable
	$mCfg.dataSnoopingDownweight    = 1e15     ; new σᵢ for snooped observations

	; compute flags — what to calculate after solving
	Local $mCompute[]
	$mCompute.qxx         = True    ; Qxx + sdx (default: on)
	$mCompute.cofactors   = False   ; Qvv, Qŷ
	$mCompute.redundancy  = False   ; diag(R)
	$mCompute.globalTest  = False   ; χ²-Test
	$mCompute.diagnostics = False   ; Baarda |w|, Pope Tτ, p-value, ∇̂, MDB
	$mCfg.compute = $mCompute

	; diagnostics parameters (only used when compute.diagnostics=True or display requests it)
	Local $mDiag[]
	$mDiag.alpha        = 0.001   ; significance level for Baarda test
	$mDiag.beta         = 0.20    ; Type II error for MDB (power = 1 - beta = 0.80)
	$mDiag.alphaSuspect = Default ; default: 10 * alpha (resolved via IsKeyword)
	$mDiag.testBasis    = "pope"   ; "pope" (a posteriori ŝ₀, normal approx.) or "baarda" (a priori σ₀=1)
	$mCfg.diagnostics = $mDiag

	Return $mCfg
EndFunc

#EndRegion ; Configuration

#Region Solving

; #FUNCTION# ====================================================================================================================
; Name...........: _adj_solve
; Description ...: Main entry point: prepares the model and solves the adjustment
; Syntax.........: _adj_solve(ByRef $mSystem[, $mConfig = Default])
; Parameters ....: $mSystem     - [ByRef] The adjustment system (populated via _adj_add* functions)
;                  $mConfig     - [optional] Solver configuration from _adj_defaultConfig().
;                                 Default uses _adj_defaultConfig().
; Return values .: Success      - Implicit (results stored in $mSystem.results)
;                  Failure      - False, @error set:
;                               |$ADJ_ERR_INPUT  - Invalid solver type or cross-group covariances
;                               |$ADJ_ERR_DOF    - Negative degrees of freedom
;                               |$ADJ_ERR_SOLVER - LAPACK solver failure
;                               |$ADJ_ERR_NO_CONVERGENCE - Iteration did not converge
;                               |$ADJ_ERR_NOT_POS_DEF    - Σₗₗ not positive definite
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: Two-phase solve: (1) robust IRLS if configured, (2) VCE or single-pass.
;                  Results accessible via _adj_getResults() or _adj_displayResults().
; Related .......: _adj_defaultConfig, _adj_getResults, _adj_displayResults
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func _adj_solve(ByRef $mSystem, $mConfig = Default)
	If IsKeyword($mConfig) Then $mConfig = _adj_defaultConfig()
	$mSystem.config = $mConfig

	; reject models with cross-group covariances (incompatible with VKS block-diagonal assumption)
	If ($mSystem.model)._hasCrossGroupCovar Then Return SetError($ADJ_ERR_INPUT, 5, False)

	; validate solver type
	If $mConfig.solver <> $ADJ_SOLVER_QR And $mConfig.solver <> $ADJ_SOLVER_SVD Then
		Return SetError($ADJ_ERR_INPUT, 0, False)
	EndIf

	; data snooping is mutually exclusive with robust IRLS (different outlier philosophies;
	; combining them would mask the snooping test statistics).
	If MapExists($mConfig, "dataSnooping") And $mConfig.dataSnooping And $mConfig.robust <> "" _
			And Not (MapExists($mConfig, "_inSnooping") And $mConfig._inSnooping) Then
		Return SetError($ADJ_ERR_INPUT, 6, False)
	EndIf

	; validate snooping method
	If MapExists($mConfig, "dataSnooping") And $mConfig.dataSnooping _
			And $mConfig.dataSnoopingMethod <> "Pope" And $mConfig.dataSnoopingMethod <> "Baarda" Then
		Return SetError($ADJ_ERR_INPUT, 7, False)
	EndIf

	__adj_prepareModel($mSystem)
	If @error Then Return SetError(@error, @extended, False)

	Local $mState = $mSystem.state       ; ONE-TIME extraction
	$mState.solver = $mConfig.solver
	$mState.solverRCOND = (IsKeyword($mConfig.solverRCOND) = 1) ? 1e-5 : $mConfig.solverRCOND
	$mState.scaling = $mConfig.scaling
	$mState.fLMTau = MapExists($mConfig, "lmTau") ? $mConfig.lmTau : $__ADJ_LM_TAU_DEFAULT

	; ── Phase 1: Robust estimation (IRLS) ──
	If $mConfig.robust <> "" Then
		__adj_robustIRLS($mSystem, $mState)
		Local $__iErrR = @error, $__iExtR = @extended
		If $__iErrR Then
			$mSystem.state = $mState
			Return SetError($__iErrR, $__iExtR, False)
		EndIf
		; reset iteration state for clean VKS/final phase
		If MapExists($mState, "nIterations")      Then MapRemove($mState, "nIterations")
		If MapExists($mState, "r_accumulated")     Then MapRemove($mState, "r_accumulated")
		If MapExists($mState, "LM_D")              Then MapRemove($mState, "LM_D")
		If MapExists($mState, "LM_gradient")       Then MapRemove($mState, "LM_gradient")
		If MapExists($mState, "LM_step")           Then MapRemove($mState, "LM_step")
		If MapExists($mState, "EquilibrationScale") Then MapRemove($mState, "EquilibrationScale")
	EndIf

	; ── Phase 2: VKS or single-pass ──
	__adj_estimateVCE($mSystem, $mState)  ; ByRef -- no copy!
	Local $__iErr = @error, $__iExt = @extended
	If $__iErr Then
		$mSystem.state = $mState         ; Write-Back even on error
		Return SetError($__iErr, $__iExt, False)
	EndIf

	$mSystem.state = $mState              ; ONE-TIME Write-Back

	; sync iteration state back to results for backward compat (_adj_getResults reads $mSystem.results)
	Local $mResults = $mSystem.results
	$mResults.r2sum  = $mState.r2sum
	$mResults.r2sum0 = $mState.r2sum0
	If MapExists($mState, "nIterations") Then $mResults.nIterations = $mState.nIterations
	If MapExists($mState, "r")           Then $mResults.r           = $mState.r
	If MapExists($mState, "Vector_x0")   Then $mResults.x0          = $mState.Vector_x0
	If MapExists($mState, "x1")          Then $mResults.x1          = $mState.x1
	If MapExists($mState, "xd")          Then $mResults.xd          = $mState.xd
	$mResults.rankDeficient = $mState.rankDeficient
	$mResults.solver = $mState.solver
	$mResults.conditionNumber = MapExists($mState, "fConditionNumber") ? $mState.fConditionNumber : Default

	; robust estimation results
	If MapExists($mState, "robustWeights") Then
		$mResults.robustWeights    = $mState.robustWeights
		$mResults.robustIterations = $mState.robustIterations
		$mResults.robustScale      = $mState.robustScale
		$mResults.robustConverged  = $mState.robustConverged
	EndIf

	$mSystem.results = $mResults

	; ── Phase 3: Iterative data snooping (only on outermost call) ──
	If MapExists($mConfig, "dataSnooping") And $mConfig.dataSnooping _
			And Not (MapExists($mConfig, "_inSnooping") And $mConfig._inSnooping) Then
		__adj_dataSnooping($mSystem, $mConfig)
		$mSystem.config = $mConfig   ; restore user config (inner calls overwrote with _inSnooping flag)
	EndIf
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_dataSnooping
; Description ...: Iterative outlier detection via Baarda or Pope test. Each iteration:
;                  computes diagnostics, picks the observation with the largest |statistic|
;                  among those with r_i ≥ minRedundancy, compares against a Sidak-corrected
;                  critical value, soft-downweights the worst (σᵢ → downweight) and re-solves.
;                  Stops when max statistic falls below the critical value, when no testable
;                  observation remains, or when the iteration cap is reached.
; Syntax.........: __adj_dataSnooping(ByRef $mSystem, $mConfig)
; Parameters ....: $mSystem     - [ByRef] The adjustment system (must already be solved)
;                  $mConfig     - User's solver config (snooping fields read here)
; Return values .: None — populates $mSystem.results.snoopingReport (Array of Maps with keys
;                  .key, .iteration, .method, .statistic, .critical, .redundancy,
;                  .stdDevOriginal), .snoopingConverged (Bool), .snoopingIterations (Int).
; Author ........: AspirinJunkie
; Remarks .......: Soft-downweighting only — no observation is physically removed.
;                  Original σᵢ is preserved in obs._stdDevOriginal so the user can restore.
; ===============================================================================================================================
Func __adj_dataSnooping(ByRef $mSystem, $mConfig)
	Local $sMethod    = $mConfig.dataSnoopingMethod
	Local $fAlpha     = $mConfig.dataSnoopingAlpha
	Local $iMaxIter   = $mConfig.dataSnoopingMaxIter
	Local $fMinR      = $mConfig.dataSnoopingMinRedundancy
	Local $fDownSigma = $mConfig.dataSnoopingDownweight

	; inner-call config (suppresses recursive snooping)
	Local $mInner = $mConfig
	$mInner._inSnooping = True

	Local $aReport[0]
	Local $bConverged = False
	Local $iIter = 0

	While $iIter < $iMaxIter
		$iIter += 1

		; ensure diagnostics are present (first iteration may need them; subsequent always)
		__adj_ensureComputed($mSystem, "diagnostics")
		If @error Then ExitLoop

		Local $mResults = $mSystem.results
		Local $iF       = $mResults.f
		If $iF <= 0 Then ExitLoop

		; pick statistic map for current method (key-indexed, populated by __adj_computeDiagnostics)
		Local $mStatMap = ($sMethod = "Pope") ? $mResults.popeT : $mResults.baardaW
		Local $mObsAll  = ($mSystem.model).obs
		Local $mIdxObs  = ($mSystem.state).idxObs   ; key → solver-row index
		Local $mVecRed  = $mResults.redundancyDiag  ; BLAS vector, indexed by idxObs

		; effective n = currently active (not yet snooped) observations
		Local $iNeff = 0
		For $sK In MapKeys($mObsAll)
			If Not (MapExists($mObsAll[$sK], "_snooped") And ($mObsAll[$sK])._snooped) Then $iNeff += 1
		Next
		If $iNeff < 1 Then ExitLoop

		; Sidak-corrected critical value (two-sided test)
		Local $fAlphaEff = 1 - (1 - $fAlpha) ^ (1 / $iNeff)
		Local $fP        = 1 - $fAlphaEff / 2
		Local $fCrit
		If $sMethod = "Pope" Then
			$fCrit = ($iF > 1) ? __adj_popeQuantile($fP, $iF) : __adj_normQuantile($fP)
		Else
			$fCrit = __adj_normQuantile($fP)
		EndIf

		; find candidate with maximal |statistic| among controllable observations (r_i ≥ min)
		Local $sBestKey = "", $fBestStat = -1, $fBestR = 0, $fBestSigma = 0
		For $sK In MapKeys($mStatMap)
			Local $mObsK = $mObsAll[$sK]
			If MapExists($mObsK, "_snooped") And $mObsK._snooped Then ContinueLoop
			Local $vStat = $mStatMap[$sK]
			If IsKeyword($vStat) Then ContinueLoop
			Local $vR = __adj_vecGet($mVecRed, $mIdxObs[$sK])
			If $vR < $fMinR Then ContinueLoop
			Local $fAbs = Abs($vStat)
			If $fAbs > $fBestStat Then
				$fBestStat  = $fAbs
				$sBestKey   = $sK
				$fBestR     = $vR
				$fBestSigma = MapExists($mObsK, "_stdDevOriginal") ? $mObsK._stdDevOriginal : $mObsK.stdDev
			EndIf
		Next

		If $sBestKey = "" Then ExitLoop                ; no testable observation
		If $fBestStat < $fCrit Then
			$bConverged = True
			ExitLoop
		EndIf

		; record snooping decision
		Local $mEntry[]
		$mEntry.key            = $sBestKey
		$mEntry.iteration      = $iIter
		$mEntry.method         = $sMethod
		$mEntry.statistic      = $fBestStat
		$mEntry.critical       = $fCrit
		$mEntry.redundancy     = $fBestR
		$mEntry.stdDevOriginal = $fBestSigma
		ReDim $aReport[UBound($aReport) + 1]
		$aReport[UBound($aReport) - 1] = $mEntry

		; downweight the offending observation (soft removal — preserves model structure)
		Local $mModel = $mSystem.model
		Local $mObs   = $mModel.obs
		Local $mObsB  = $mObs[$sBestKey]
		If Not MapExists($mObsB, "_stdDevOriginal") Then $mObsB._stdDevOriginal = $mObsB.stdDev
		$mObsB.stdDev   = $fDownSigma
		$mObsB._snooped = True
		If MapExists($mObsB, "weight")           Then MapRemove($mObsB, "weight")
		If MapExists($mObsB, "_weightOriginal") Then MapRemove($mObsB, "_weightOriginal")
		$mObs[$sBestKey] = $mObsB
		$mModel.obs      = $mObs
		$mSystem.model   = $mModel

		; reset prepared state to force a full rebuild of weights, Σₗₗ, Cholesky, and matrices
		Local $mEmpty[]
		If MapExists($mSystem, "state") Then MapRemove($mSystem, "state")
		$mSystem.results = $mEmpty
		If MapExists($mSystem, "_prepareRunCount") Then MapRemove($mSystem, "_prepareRunCount")
		If MapExists($mSystem, "VarianceComponents") Then MapRemove($mSystem, "VarianceComponents")

		; re-solve (recursion guarded by $mInner._inSnooping)
		_adj_solve($mSystem, $mInner)
		If @error Then ExitLoop
	WEnd

	; attach report to results of the *last* solve
	Local $mResults = $mSystem.results
	$mResults.snoopingReport     = $aReport
	$mResults.snoopingConverged  = $bConverged
	$mResults.snoopingIterations = $iIter
	$mResults.snoopingMethod     = $sMethod
	$mSystem.results = $mResults
EndFunc

#EndRegion ; Solving

#Region Results and post-processing

; #FUNCTION# ====================================================================================================================
; Name...........: _adj_getResults
; Description ...: Extracts and structures all adjustment results into a user-friendly flat Map
; Syntax.........: _adj_getResults(ByRef $mSystem)
; Parameters ....: $mSystem     - [ByRef] The adjustment system after calling _adj_solve()
; Return values .: Success      - Map with flat results structure:
;                               |.modelType, .algorithm, .s0, .f, .vtpv, .nIterations, .rankDeficient
;                               |.x1, .sdx, .xd (Maps: paramName → value)
;                               |.obsValue, .v, .obsAdj, .sdv, .sdyhat, .r (Maps: obsName → value)
;                               |.Qxx, .Qvv, .Qyhat (matrices or Null)
;                               |.globalTest*, .baardaW, .popeT, .pValue, .blunder, .mdb,
;                               +.testDecision (conditional, only when compute flags set)
;                               |.vceConverged, .vceIterations, .vceGroups (VCE only)
;                               |.robustWeights, .robustIterations, .robustScale (robust only)
;                  Failure      - False, @error = $ADJ_ERR_INPUT
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: All parameter/observation results are returned as named Maps (not index-based).
;                  Conditional keys only present when the corresponding compute flag was enabled.
; Related .......: _adj_solve, _adj_displayResults
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func _adj_getResults(ByRef $mSystem)
	Local $mResults = $mSystem.results
	If Not MapExists($mResults, "s0") Then Return SetError($ADJ_ERR_INPUT, 0, False)

	Local $mRet[]

	; --- Metadata ---
	Local $mState = $mSystem.state
	Local $mModel = $mSystem.model
	$mRet.modelType      = $mState.modelType
	$mRet.algorithm      = $mSystem.config.algorithm
	$mRet.s0             = $mResults.s0
	$mRet.f              = $mResults.f
	$mRet.vtpv           = $mResults.vtpv
	$mRet.nIterations    = $mResults.nIterations
	$mRet.rankDeficient  = MapExists($mResults, "rankDeficient") ? $mResults.rankDeficient : False
	$mRet.solver          = MapExists($mResults, "solver") ? $mResults.solver : "QR"
	$mRet.conditionNumber = MapExists($mResults, "conditionNumber") ? $mResults.conditionNumber : Default

	; --- Parameters (flat named maps) ---
	Local $mX1[], $mSdx[], $mXd[]

	If $mModel.nParams > 0 And MapExists($mResults, "x1") Then
		Local $mIdxParams = $mState.idxParams
		Local $mVecX1  = $mResults.x1
		Local $bHasSdx = MapExists($mResults, "sdx")
		Local $mVecSdx = $bHasSdx ? $mResults.sdx : Null
		Local $bHasXd  = MapExists($mResults, "xd")
		Local $mVecXd  = $bHasXd ? $mResults.xd : Null

		For $sKey In MapKeys($mIdxParams)
			Local $iIdxP = $mIdxParams[$sKey]
			$mX1[$sKey]  = __adj_vecGet($mVecX1, $iIdxP)
			If $bHasSdx Then $mSdx[$sKey] = __adj_vecGet($mVecSdx, $iIdxP)
			If $bHasXd  Then $mXd[$sKey]  = __adj_vecGet($mVecXd, $iIdxP)
		Next
	EndIf
	$mRet.x1  = $mX1
	$mRet.sdx = $mSdx
	$mRet.xd  = $mXd

	; --- Observations (flat named maps) ---
	Local $mObsValue[], $mV[], $mObsAdj[], $mSdv[], $mSdyhat[], $mRedundancy[]
	Local $mIdxObs      = $mState.idxObs
	Local $mObservations = $mModel.obs

	Local $bHasR   = MapExists($mResults, "r")
	Local $bHasSdv = MapExists($mResults, "sdv")
	Local $bHasSdy = MapExists($mResults, "sdy")
	Local $bHasRed = MapExists($mResults, "redundancyDiag")

	Local $mVecR   = $bHasR   ? $mResults.r              : Null
	Local $mVecSdv = $bHasSdv ? $mResults.sdv             : Null
	Local $mVecSdy = $bHasSdy ? $mResults.sdy             : Null
	Local $mVecRed = $bHasRed ? $mResults.redundancyDiag  : Null

	For $sKey In MapKeys($mIdxObs)
		Local $iIdxO    = $mIdxObs[$sKey]
		Local $mObs     = $mObservations[$sKey]

		$mObsValue[$sKey] = $mObs.value
		Local $fResidual  = $bHasR ? __adj_vecGet($mVecR, $iIdxO) : 0
		$mV[$sKey]        = $fResidual
		$mObsAdj[$sKey]   = $mObs.value + $fResidual

		If $bHasSdv Then $mSdv[$sKey]       = __adj_vecGet($mVecSdv, $iIdxO)
		If $bHasSdy Then $mSdyhat[$sKey]    = __adj_vecGet($mVecSdy, $iIdxO)
		If $bHasRed Then $mRedundancy[$sKey] = __adj_vecGet($mVecRed, $iIdxO)
	Next
	$mRet.obsValue = $mObsValue
	$mRet.v        = $mV
	$mRet.obsAdj   = $mObsAdj
	$mRet.sdv      = $mSdv
	$mRet.sdyhat   = $mSdyhat
	$mRet.r        = $mRedundancy

	; --- Sanity check warning (Σrᵢ = f) ---
	If MapExists($mResults, "redundancyTraceMismatch") Then
		$mRet.redundancyTraceMismatch = $mResults.redundancyTraceMismatch
	EndIf

	; --- Raw cofactor matrices ---
	$mRet.Qxx         = MapExists($mResults, "Qxx")         ? $mResults.Qxx         : Null
	$mRet.Qvv         = MapExists($mResults, "Qvv")         ? $mResults.Qvv         : Null
	$mRet.Qyhat       = MapExists($mResults, "Qyhat")       ? $mResults.Qyhat       : Null
	$mRet.correlation = MapExists($mResults, "correlation") ? $mResults.correlation : Null

	; --- Global test (from compute-on-demand, if computed) ---
	If MapExists($mResults, "globalTestPassed") Then
		$mRet.globalTestPassed = $mResults.globalTestPassed
		$mRet.globalTestT      = $mResults.globalTestT
		$mRet.globalTestLower  = $mResults.globalTestLower
		$mRet.globalTestUpper  = $mResults.globalTestUpper
		$mRet.globalTestAlpha  = $mResults.globalTestAlpha
	EndIf

	; --- Diagnostics (from compute-on-demand, if computed) ---
	If MapExists($mResults, "baardaW") Then
		$mRet.baardaW        = $mResults.baardaW
		$mRet.popeT          = $mResults.popeT
		$mRet.pValue         = $mResults.pValue
		$mRet.popePValue     = $mResults.popePValue
		$mRet.blunder        = $mResults.blunder
		$mRet.mdb            = $mResults.mdb
		$mRet.testDecision   = $mResults.testDecision
		$mRet.testBasis      = $mResults.testBasis
		$mRet.baardaWarning  = $mResults.baardaWarning
	EndIf

	; --- VCE results (only present when VCE was active) ---
	If MapExists($mResults, "vceConverged") Then
		$mRet.vceConverged  = $mResults.vceConverged
		$mRet.vceIterations = $mResults.vceIterations
		$mRet.vceGroups     = $mResults.vceGroups
	EndIf

	; --- Robust estimation results ---
	If MapExists($mResults, "robustWeights") Then
		$mRet.robustWeights    = $mResults.robustWeights
		$mRet.robustIterations = $mResults.robustIterations
		$mRet.robustScale      = $mResults.robustScale
		$mRet.robustConverged  = $mResults.robustConverged
	EndIf

	; --- Iterative data snooping ---
	If MapExists($mResults, "snoopingReport") Then
		$mRet.snoopingReport     = $mResults.snoopingReport
		$mRet.snoopingConverged  = $mResults.snoopingConverged
		$mRet.snoopingIterations = $mResults.snoopingIterations
		$mRet.snoopingMethod     = $mResults.snoopingMethod
	EndIf

	Return $mRet
EndFunc

; #FUNCTION# ====================================================================================================================
; Name...........: _adj_getOutliers
; Description ...: Returns observation names identified as outliers
; Syntax.........: _adj_getOutliers(ByRef $mSystem[, $fThreshold = Default[, $sMethod = "robustWeight"]])
; Parameters ....: $mSystem     - [ByRef] The adjustment system after _adj_solve()
;                  $fThreshold  - [optional] Detection threshold (Default: method-dependent):
;                               |"robustWeight": 0.5
;                               |"absU": 2.5
;                               |"baarda"/"pope": $mConfig.diagnostics.alpha
;                  $sMethod     - [optional] Detection method:
;                               |"robustWeight" - outlier if robust weight < threshold
;                               |"absU"         - outlier if |u_i| > threshold
;                               |"baarda"       - outlier if Baarda p-value < alpha
;                               |"pope"         - outlier if Pope p-value < alpha
; Return values .: Success      - Array of outlier observation names (may be empty)
;                  Failure      - SetError($ADJ_ERR_INPUT) if required results not available
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: "robustWeight" and "absU" require prior robust estimation (_adj_solve with
;                  .robust configured). "baarda"/"pope" compute diagnostics if not yet available.
; Related .......: _adj_solve, _adj_removeObs, __adj_computeDiagnostics
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func _adj_getOutliers(ByRef $mSystem, $fThreshold = Default, $sMethod = "robustWeight")
	Local $mResults = $mSystem.results

	; guard: robustWeight/absU require prior robust estimation
	If ($sMethod = "robustWeight" Or $sMethod = "absU") And Not MapExists($mResults, "robustWeights") Then
		Return SetError($ADJ_ERR_INPUT, 0, False)
	EndIf

	Local $aOutliers[0]

	Switch $sMethod
		Case "robustWeight"
			Local $mRobustWeights = $mResults.robustWeights
			; resolve default threshold
			If IsKeyword($fThreshold) Then $fThreshold = 0.5
			; outlier if rob.w < threshold
			For $sKey In MapKeys($mRobustWeights)
				If $mRobustWeights[$sKey] < $fThreshold Then
					ReDim $aOutliers[UBound($aOutliers) + 1]
					$aOutliers[UBound($aOutliers) - 1] = $sKey
				EndIf
			Next

		Case "absU"
			If IsKeyword($fThreshold) Then $fThreshold = 2.5
			; outlier if |u_i| > threshold — requires recomputing |u| from residuals
			Local $mRobustWeightsU = $mResults.robustWeights
			Local $mRes = _adj_getResults($mSystem)
			Local $iErr = @error
			If $iErr Then Return SetError($iErr, @extended, False)
			Local $fScale = MapExists($mRes, "robustScale") ? $mRes.robustScale : 1.0
			If $fScale < 1e-15 Then $fScale = 1.0
			Local $mObs = ($mSystem.model).obs
			For $sKey In MapKeys($mRobustWeightsU)
				Local $fV = $mRes.v[$sKey]
				Local $fSigma = ($mObs[$sKey]).stdDev
				Local $fAbsU = Abs($fV) / ($fScale * $fSigma)
				If $fAbsU > $fThreshold Then
					ReDim $aOutliers[UBound($aOutliers) + 1]
					$aOutliers[UBound($aOutliers) - 1] = $sKey
				EndIf
			Next

		Case "baarda", "pope"
			; outlier if Baarda/Pope p-value < alpha
			__adj_ensureComputed($mSystem, "diagnostics")
			If @error Then Return SetError(@error, @extended, False)
			$mResults = $mSystem.results ; re-read after ensureComputed
			Local $fAlphaThreshold = IsKeyword($fThreshold) ? ($mSystem.config).diagnostics.alpha : $fThreshold
			Local $mPValues = $mResults.pValue
			For $sKey In MapKeys($mPValues)
				If Not IsKeyword($mPValues[$sKey]) And $mPValues[$sKey] < $fAlphaThreshold Then
					ReDim $aOutliers[UBound($aOutliers) + 1]
					$aOutliers[UBound($aOutliers) - 1] = $sKey
				EndIf
			Next

		Case Else
			Return SetError($ADJ_ERR_INPUT, 0, False)
	EndSwitch

	Return $aOutliers
EndFunc

; #FUNCTION# ====================================================================================================================
; Name...........: _adj_removeObs
; Description ...: Removes one or more observations from the adjustment system
; Syntax.........: _adj_removeObs(ByRef $mSystem, $vObs)
; Parameters ....: $mSystem     - [ByRef] The adjustment system
;                  $vObs        - Observation name (String) or Array of observation names
; Return values .: Success      - Number of removed observations (Int)
;                  Failure      - SetError($ADJ_ERR_INPUT) for invalid input type
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: Also removes associated formulas (ObsFunction type), covariances, and clears
;                  prepared state. Call _adj_solve() again after removal to recompute.
; Related .......: _adj_addObs, _adj_getOutliers, _adj_solve
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func _adj_removeObs(ByRef $mSystem, $vObs)
	; normalize to array
	Local $aObs
	If IsString($vObs) Then
		Local $aTmp[1] = [$vObs]
		$aObs = $aTmp
	ElseIf IsArray($vObs) Then
		$aObs = $vObs
	Else
		Return SetError($ADJ_ERR_INPUT, 0, 0)
	EndIf

	Local $mModel = $mSystem.model
	Local $mObservations = $mModel.obs
	Local $iRemoved = 0

	For $k = 0 To UBound($aObs) - 1
		Local $sName = StringUpper($aObs[$k])
		If Not MapExists($mObservations, $sName) Then ContinueLoop

		; 1. Remove from obs map
		MapRemove($mObservations, $sName)
		$mModel.nObs -= 1

		; 2. Remove associated formulas (where .value == obsName → ObsFunction type)
		If $mModel.nFormulas > 0 Then
			Local $aFormulas = $mModel.formulas
			Local $aNewFormulas[$mModel.nFormulas]
			Local $iNewCount = 0
			For $j = 0 To $mModel.nFormulas - 1
				If IsString($aFormulas[$j].value) And StringUpper($aFormulas[$j].value) = $sName Then
					; skip — this formula belongs to the removed observation
					ContinueLoop
				EndIf
				$aNewFormulas[$iNewCount] = $aFormulas[$j]
				$iNewCount += 1
			Next
			Local $iRemovedFormulas = $mModel.nFormulas - $iNewCount
			If $iRemovedFormulas > 0 Then
				ReDim $aNewFormulas[$iNewCount]
				$mModel.formulas = $aNewFormulas
				$mModel.nFormulas = $iNewCount
			EndIf
		EndIf

		; 3. Remove covariances involving this observation
		If MapExists($mModel, "covariances") Then
			Local $mCov = $mModel.covariances
			For $sCovKey In MapKeys($mCov)
				If StringInStr($sCovKey, $sName) Then
					MapRemove($mCov, $sCovKey)
				EndIf
			Next
			$mModel.covariances = $mCov
		EndIf

		$iRemoved += 1
	Next

	; 4. Write back model
	$mModel.obs = $mObservations
	$mSystem.model = $mModel

	; 5. Clear prepared state — next _adj_solve() rebuilds everything
	; Includes _prepareRunCount and VarianceComponents so __adj_prepareModel runs full setup again
	; (otherwise it short-circuits and leaves $mSystem.state empty → next solve crashes).
	; results is reset to an empty map (not removed) because __adj_prepareModel writes into it.
	Local $mEmpty[]
	If MapExists($mSystem, "state") Then MapRemove($mSystem, "state")
	$mSystem.results = $mEmpty
	If MapExists($mSystem, "_prepareRunCount") Then MapRemove($mSystem, "_prepareRunCount")
	If MapExists($mSystem, "VarianceComponents") Then MapRemove($mSystem, "VarianceComponents")

	Return $iRemoved
EndFunc

#Region Error messages

; #FUNCTION# ====================================================================================================================
; Name...........: _adj_getErrorMessage
; Description ...: Returns a human-readable error message for an adjustment error code
; Syntax.........: _adj_getErrorMessage($iErr[, $iExt = 0])
; Parameters ....: $iErr         - Error code ($ADJ_ERR_* constant or @error value)
;                  $iExt         - [optional] Extended error code (@extended value). Default is 0.
; Return values .: String with error description
; Author ........: AspirinJunkie
; Related .......: _adj_solve, _adj_addObs, _adj_addObsFunction
; ===============================================================================================================================
Func _adj_getErrorMessage($iErr, $iExt = 0)
	Switch $iErr
		Case $ADJ_ERR_OK
			Return "No error"
		Case $ADJ_ERR_SOLVER
			Return "LAPACK solver error (LAPACK info = " & $iExt & ")"
		Case $ADJ_ERR_DOF
			Return "Invalid degrees of freedom: f = " & $iExt
		Case $ADJ_ERR_RANK
			Return "Rank deficiency detected"
		Case $ADJ_ERR_NO_CONVERGENCE
			Return "No convergence after " & $iExt & " iterations"
		Case $ADJ_ERR_NOT_POS_DEF
			Return "Matrix not positive definite (Cholesky factorization failed)"
		Case $ADJ_ERR_RESTR_OVERDETERMINED
			Return "Restrictions fully determine all parameters (nRestrictions=" & $iExt & " >= nParams)"
		Case $ADJ_ERR_INPUT
			Switch $iExt
				Case 0
					Return "Invalid input"
				Case 2
					Return "First observation not found in system"
				Case 3
					Return "Second observation not found in system"
				Case 4
					Return "Cross-group covariance not supported"
				Case 10
					Return "Observation value is not a finite number"
				Case 11
					Return "Standard deviation must be a positive finite number"
				Case Else
					Return "Invalid input (extended = " & $iExt & ")"
			EndSwitch
	EndSwitch
	Return "Unknown error code: " & $iErr
EndFunc

#EndRegion ; Error messages

#EndRegion ; Results and post-processing

#Region Internal helpers

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_createNewSystem
; Description ...: Creates the internal system map structure with model/results sub-maps
; Syntax.........: __adj_createNewSystem([$fS0 = 1])
; Parameters ....: $fS0         - [optional] Initial σ₀ (currently unused). Default is 1.
; Return values .: Map with .model sub-map (obs, params, fixed, formulas, restrictions, covariances,
;                  counters, flags) and empty .results sub-map
; Author ........: AspirinJunkie
; ===============================================================================================================================
Func __adj_createNewSystem($fS0 = 1)
	#forceref $fS0
	Local $mRet[], $mEmpty[], $aEmpty[0]

	; model sub-map (user input data)
	Local $mModel[]
	$mModel.obs                        = $mEmpty
	$mModel.params                     = $mEmpty
	$mModel.fixed                      = $mEmpty
	$mModel.formulas                   = $aEmpty
	$mModel.restrictions               = $aEmpty
	$mModel.covariances                = $mEmpty
	$mModel.nObs                       = 0
	$mModel.nFormulas                  = 0
	$mModel.nParams                    = 0
	$mModel.nRestrictions              = 0
	$mModel._lastOID                   = 1
	$mModel._isWeighted                = False
	$mModel._hasCovar                  = False
	$mModel._multiVarComp              = False
	$mModel._hasCrossGroupCovar        = False
	$mModel._hasParamFunctions         = False
	$mModel._hasObsFunctions           = False
	$mModel._hasParamRestrictions      = False
	$mModel._hasObsRestrictions        = False
	$mRet.model = $mModel

	; results sub-map
	$mRet.results                    = $mEmpty

	Return $mRet
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_parseDerivativeInput
; Description ...: Parses derivative input (Map or pipe-delimited string) into a normalized Map
; Syntax.........: __adj_parseDerivativeInput($vInput)
; Parameters ....: $vInput      - Map, pipe-delimited string ("X=2*X | Y=2*Y"), or Default
; Return values .: Success      - Map with upper-cased keys and values
;                  No input     - False (if Default, empty string, or unsupported type)
; Author ........: AspirinJunkie
; ===============================================================================================================================
Func __adj_parseDerivativeInput($vInput)
	If IsKeyword($vInput) = 1 Or $vInput = "" Then Return False

	; already a Map — uppercase keys and values
	If IsMap($vInput) Then
		Local $mResult[], $sKey
		For $sKey In MapKeys($vInput)
			$mResult[StringUpper($sKey)] = StringUpper($vInput[$sKey])
		Next
		Return $mResult
	EndIf

	; parse pipe-delimited string: "X=2*X | Y=2*Y"
	If IsString($vInput) Then
		Local $mResult[], $aParts = StringSplit($vInput, "|", 2)
		Local $sPart, $iEq, $sKey
		For $sPart In $aParts
			$sPart = StringStripWS($sPart, 3)
			If $sPart = "" Then ContinueLoop
			$iEq = StringInStr($sPart, "=")
			If $iEq < 2 Then ContinueLoop
			$sKey = StringStripWS(StringLeft($sPart, $iEq - 1), 3)
			Local $sVal = StringStripWS(StringMid($sPart, $iEq + 1), 3)
			If $sKey <> "" And $sVal <> "" Then $mResult[StringUpper($sKey)] = StringUpper($sVal)
		Next
		Return $mResult
	EndIf

	Return False
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_vecGet
; Description ...: Reads a single element from a BLAS vector map (0-based index)
; Syntax.........: __adj_vecGet($mVec, $iIdx)
; Parameters ....: $mVec        - [Const ByRef] BLAS vector map (with .struct)
;                  $iIdx        - 0-based element index
; Return values .: Element value (Double)
; Author ........: AspirinJunkie
; ===============================================================================================================================
Func __adj_vecGet(Const ByRef $mVec, Const $iIdx)
	Return DllStructGetData($mVec.struct, 1, $iIdx + 1)
EndFunc
#EndRegion ; Internal helpers
