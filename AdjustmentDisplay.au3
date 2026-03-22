#include-once

; #INDEX# =======================================================================================================================
; Title .........: AdjustmentDisplay
; AutoIt Version : 3.3.16.1
; Description ...: Configurable display system for geodetic adjustment results: formatted text output
;                  with selectable columns, compute-on-demand integration, and table rendering.
; Author(s) .....: AspirinJunkie
; ===============================================================================================================================


#Region Display API

; #FUNCTION# ====================================================================================================================
; Name...........: _adj_defaultDisplayConfig
; Description ...: Creates a default display configuration map
; Syntax.........: _adj_defaultDisplayConfig()
; Parameters ....: None
; Return values .: Map with display options:
;                  - .showHeader     [Bool] show header section (model, s0, f, vtpv). Default: True
;                  - .showParams     [Bool] show adjusted parameters table. Default: True
;                  - .showObs        [Bool] show observations table. Default: True
;                  - .showGlobalTest [Bool] show global test results. Default: False
;                  - .showVCE        [Bool] show VCE results. Default: False
;                  - .showRobust     [Bool] show robust estimation info. Default: False
;                  - .paramCols      [String] pipe-separated column names for param table
;                  - .obsCols        [String] pipe-separated column names for obs table
;                  - .precision      [Int] number of significant digits. Default: 6
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: Available columns for obsCols: name, value, v, sdv, sdyhat, r, w, T, p, pPope,
;                  blunder, mdb, decision, robW.
;                  Available columns for paramCols: name, value, sdx, xd.
; Related .......: _adj_displayResults
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func _adj_defaultDisplayConfig()
	Local $mDisplay[]
	$mDisplay.showHeader     = True
	$mDisplay.showParams     = True
	$mDisplay.showObs        = True
	$mDisplay.showGlobalTest = False
	$mDisplay.showVCE        = False
	$mDisplay.showRobust     = False
	$mDisplay.paramCols = "name|value|sdx"
	$mDisplay.obsCols   = "name|value|v|sdv"
	$mDisplay.precision = 6
	Return $mDisplay
EndFunc


; #FUNCTION# ====================================================================================================================
; Name...........: _adj_displayResults
; Description ...: Formats adjustment results as a readable string based on display configuration
; Syntax.........: _adj_displayResults(ByRef $mSystem[, $mDisplayConfig = Default])
; Parameters ....: $mSystem         - [ByRef] The adjustment system after calling _adj_solve()
;                  $mDisplayConfig  - [optional] Display config from _adj_defaultDisplayConfig().
;                                     Default uses _adj_defaultDisplayConfig().
; Return values .: Success          - Formatted string with adjustment summary
;                  Failure          - Empty string, @error set
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: Triggers compute-on-demand for requested columns (e.g. requesting "r" column
;                  computes redundancy if not already computed). Columns that fail to compute
;                  display "---".
; Related .......: _adj_defaultDisplayConfig, _adj_getResults, __adj_ensureComputed
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func _adj_displayResults(ByRef $mSystem, $mDisplayConfig = Default)
	If IsKeyword($mDisplayConfig) Then $mDisplayConfig = _adj_defaultDisplayConfig()

	; trigger compute-on-demand for requested columns
	__adj_displayEnsureDependencies($mSystem, $mDisplayConfig)

	; get flattened results
	Local $mRes = _adj_getResults($mSystem)
	If @error Then Return SetError(@error, @extended, "")

	Local $iPrecision = $mDisplayConfig.precision
	Local $sFmtW = "%15." & $iPrecision & "g"
	Local $sRet = ""

	; --- Header ---
	If $mDisplayConfig.showHeader Then
		$sRet &= "===============================================" & @CRLF
		$sRet &= "  Adjustment Results" & @CRLF
		$sRet &= "===============================================" & @CRLF
		$sRet &= "Model:           " & $mRes.modelType & " (" & $mRes.algorithm & ")" & @CRLF
		If MapExists($mRes, "solver") Then
			Local $sSolverInfo = "Solver:          " & $mRes.solver
			If MapExists($mRes, "conditionNumber") And Not IsKeyword($mRes.conditionNumber) Then
				$sSolverInfo &= StringFormat(" (condition number: %.2e)", $mRes.conditionNumber)
			EndIf
			$sRet &= $sSolverInfo & @CRLF
		EndIf
		$sRet &= "Iterations:      " & $mRes.nIterations & @CRLF
		$sRet &= "DOF:             f = " & $mRes.f & @CRLF
		$sRet &= StringFormat("s0:              %-15." & $iPrecision & "g", $mRes.s0) & @CRLF
		$sRet &= StringFormat("vtPv:            %-15." & $iPrecision & "g", $mRes.vtpv) & @CRLF
		If $mRes.rankDeficient Then $sRet &= "WARNING: System is rank deficient!" & @CRLF
		If $mDisplayConfig.showRobust And MapExists($mRes, "robustIterations") Then
			$sRet &= "  Robust:          " & ($mSystem.config).robust & " (" & $mRes.robustIterations & " IRLS iterations"
			$sRet &= ($mRes.robustConverged ? ", converged" : ", NOT converged") & ")" & @CRLF
			$sRet &= "  Robust scale:    " & StringFormat($sFmtW, $mRes.robustScale) & @CRLF
		EndIf
	EndIf

	; --- Global test ---
	If $mDisplayConfig.showGlobalTest And MapExists($mRes, "globalTestPassed") And Not IsKeyword($mRes.globalTestPassed) Then
		Local $sTestResult = $mRes.globalTestPassed ? "PASSED" : "NOT PASSED"
		$sRet &= @CRLF
		$sRet &= StringFormat("Global test:     %s (alpha=%.4f)", $sTestResult, $mRes.globalTestAlpha) & @CRLF
		$sRet &= StringFormat("  T = vtPv = %-15." & $iPrecision & "g" & ",  chi2 interval: [" & $sFmtW & ", " & $sFmtW & "]", _
			$mRes.globalTestT, $mRes.globalTestLower, $mRes.globalTestUpper) & @CRLF
	EndIf

	; --- Baarda warning ---
	If MapExists($mRes, "baardaWarning") And $mRes.baardaWarning Then
		$sRet &= @CRLF
		$sRet &= StringFormat("WARNING: s0 = %.4f deviates significantly from 1. The Baarda test (a priori sigma0=1)", $mRes.s0) & @CRLF
		$sRet &= "         is unreliable under these conditions. Recommendation: testBasis=""pope""" & @CRLF
	EndIf

	; --- Test basis info ---
	If MapExists($mRes, "testBasis") Then
		$sRet &= "Test basis:      " & $mRes.testBasis & @CRLF
	EndIf

	; --- Parameters ---
	If $mDisplayConfig.showParams And UBound($mRes.x1) > 0 Then
		$sRet &= @CRLF
		$sRet &= "-----------------------------------------------" & @CRLF
		$sRet &= "  Adjusted Parameters" & @CRLF
		$sRet &= "-----------------------------------------------" & @CRLF
		$sRet &= __adj_displayTable($mRes, $mDisplayConfig.paramCols, "param", $iPrecision, $mSystem)
	EndIf

	; --- Observations ---
	If $mDisplayConfig.showObs And UBound($mRes.v) > 0 Then
		$sRet &= @CRLF
		$sRet &= "-----------------------------------------------" & @CRLF
		$sRet &= "  Observations" & @CRLF
		$sRet &= "-----------------------------------------------" & @CRLF
		$sRet &= __adj_displayTable($mRes, $mDisplayConfig.obsCols, "obs", $iPrecision, $mSystem)
	EndIf

	; --- VCE results ---
	If $mDisplayConfig.showVCE And MapExists($mRes, "vceConverged") Then
		Local $sFmtN = "%12." & $iPrecision & "g"
		$sRet &= @CRLF
		$sRet &= "-----------------------------------------------" & @CRLF
		$sRet &= "  Variance Component Estimation" & @CRLF
		$sRet &= "-----------------------------------------------" & @CRLF
		$sRet &= "Converged:       " & ($mRes.vceConverged ? "Yes" : "No") & " (" & $mRes.vceIterations & " iterations)" & @CRLF
		For $sKey In MapKeys($mRes.vceGroups)
			Local $mGrpVCE = ($mRes.vceGroups)[$sKey]
			$sRet &= StringFormat("  %-16s s0=%s  sigma2=%s", $sKey, StringFormat($sFmtN, $mGrpVCE.s0), StringFormat($sFmtN, $mGrpVCE.sigma2)) & @CRLF
		Next
	EndIf

	Return $sRet
EndFunc

#EndRegion ; Display API

#Region Display internals

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_displayEnsureDependencies
; Description ...: Checks requested columns and triggers compute-on-demand for dependencies
; Syntax.........: __adj_displayEnsureDependencies(ByRef $mSystem, $mDisplayConfig)
; Parameters ....: $mSystem         - [ByRef] The adjustment system
;                  $mDisplayConfig  - Display configuration map
; Return values .: None (best-effort; failures result in "---" display)
; Author ........: AspirinJunkie
; ===============================================================================================================================
Func __adj_displayEnsureDependencies(ByRef $mSystem, $mDisplayConfig)
	Local $aObsCols = StringSplit($mDisplayConfig.obsCols, "|", 2)
	Local $mNeeded[]
	For $i = 0 To UBound($aObsCols) - 1
		Local $sReq = __adj_displayRequiredCompute($aObsCols[$i])
		If $sReq <> "" Then $mNeeded[$sReq] = True
	Next
	If $mDisplayConfig.showGlobalTest Then $mNeeded["globalTest"] = True
	; best-effort: if ensureComputed fails, display shows "---" for affected columns
	For $sKey In MapKeys($mNeeded)
		__adj_ensureComputed($mSystem, $sKey)
	Next
EndFunc


; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_displayRequiredCompute
; Description ...: Maps column names to compute dependencies for __adj_ensureComputed
; Syntax.........: __adj_displayRequiredCompute($sCol)
; Parameters ....: $sCol        - Column name (e.g. "sdv", "r", "w")
; Return values .: Dependency string ("cofactors", "redundancy", "diagnostics") or ""
; Author ........: AspirinJunkie
; ===============================================================================================================================
Func __adj_displayRequiredCompute($sCol)
	Switch $sCol
		Case "sdv", "sdyhat"
			Return "cofactors"
		Case "r"
			Return "redundancy"
		Case "w", "T", "p", "pPope", "blunder", "mdb", "decision"
			Return "diagnostics"
		Case Else
			Return ""
	EndSwitch
EndFunc


; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_displayColLabel
; Description ...: Returns the display header label for a column name
; Syntax.........: __adj_displayColLabel($sCol)
; Parameters ....: $sCol        - Column name
; Return values .: Display label string
; Author ........: AspirinJunkie
; ===============================================================================================================================
Func __adj_displayColLabel($sCol)
	Switch $sCol
		Case "name"
			Return "Name"
		Case "value"
			Return "Value"
		Case "sdx"
			Return "sd"
		Case "xd"
			Return "dx"
		Case "v"
			Return "v"
		Case "sdv"
			Return "sd(v)"
		Case "sdyhat"
			Return "sd(y)"
		Case "r"
			Return "r"
		Case "w"
			Return "|w|"
		Case "T"
			Return "Tau"
		Case "p"
			Return "p-value"
		Case "pPope"
			Return "p(Pope)"
		Case "blunder"
			Return "nabla"
		Case "mdb"
			Return "MDB"
		Case "decision"
			Return "Decision"
		Case "robW"
			Return "rob.w"
		Case Else
			Return $sCol
	EndSwitch
EndFunc


; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_displayGetValue
; Description ...: Extracts a value from the results map for a given column and key
; Syntax.........: __adj_displayGetValue($mRes, $sCol, $sKey, $sTableType, ByRef $mSystem)
; Parameters ....: $mRes        - Flattened results from _adj_getResults
;                  $sCol        - Column name (e.g. "value", "sdx", "v", "w", "decision")
;                  $sKey        - Observation or parameter name
;                  $sTableType  - "param" or "obs"
;                  $mSystem     - [ByRef] The adjustment system (for additional data)
; Return values .: Success      - Value (number or string)
;                  Not available - Default keyword
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: Handles all column types: value, sdx, xd, v, sdv, sdyhat, r, w, T, p, pPope,
;                  blunder, mdb, decision, robW.
; Related .......: __adj_displayTable
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __adj_displayGetValue($mRes, $sCol, $sKey, $sTableType, ByRef $mSystem)
	#forceref $mSystem
	Switch $sCol
		Case "value"
			If $sTableType = "param" Then
				If MapExists($mRes.x1, $sKey) Then Return ($mRes.x1)[$sKey]
			Else
				If MapExists($mRes.obsValue, $sKey) Then Return ($mRes.obsValue)[$sKey]
			EndIf

		Case "sdx"
			If MapExists($mRes.sdx, $sKey) Then Return ($mRes.sdx)[$sKey]

		Case "xd"
			If MapExists($mRes.xd, $sKey) Then Return ($mRes.xd)[$sKey]

		Case "v"
			If MapExists($mRes.v, $sKey) Then Return ($mRes.v)[$sKey]

		Case "sdv"
			If MapExists($mRes.sdv, $sKey) Then Return ($mRes.sdv)[$sKey]

		Case "sdyhat"
			If MapExists($mRes.sdyhat, $sKey) Then Return ($mRes.sdyhat)[$sKey]

		Case "r"
			If MapExists($mRes.r, $sKey) Then Return ($mRes.r)[$sKey]

		Case "w"
			If MapExists($mRes, "baardaW") And MapExists($mRes.baardaW, $sKey) Then
				Local $fW = ($mRes.baardaW)[$sKey]
				If Not IsKeyword($fW) Then Return $fW
			EndIf

		Case "T"
			If MapExists($mRes, "popeT") And MapExists($mRes.popeT, $sKey) Then
				Local $fT = ($mRes.popeT)[$sKey]
				If Not IsKeyword($fT) Then Return $fT
			EndIf

		Case "p"
			If MapExists($mRes, "pValue") And MapExists($mRes.pValue, $sKey) Then
				Local $fP = ($mRes.pValue)[$sKey]
				If Not IsKeyword($fP) Then Return $fP
			EndIf

		Case "pPope"
			If MapExists($mRes, "popePValue") And MapExists($mRes.popePValue, $sKey) Then
				Local $fPP = ($mRes.popePValue)[$sKey]
				If Not IsKeyword($fPP) Then Return $fPP
			EndIf

		Case "blunder"
			If MapExists($mRes, "blunder") And MapExists($mRes.blunder, $sKey) Then
				Local $fBl = ($mRes.blunder)[$sKey]
				If Not IsKeyword($fBl) Then Return $fBl
			EndIf

		Case "mdb"
			If MapExists($mRes, "mdb") And MapExists($mRes.mdb, $sKey) Then
				Local $fMDB = ($mRes.mdb)[$sKey]
				If Not IsKeyword($fMDB) Then Return $fMDB
			EndIf

		Case "decision"
			If MapExists($mRes, "testDecision") And MapExists($mRes.testDecision, $sKey) Then
				Return ($mRes.testDecision)[$sKey]
			EndIf

		Case "robW"
			If MapExists($mRes, "robustWeights") And MapExists($mRes.robustWeights, $sKey) Then
				Return ($mRes.robustWeights)[$sKey]
			EndIf
	EndSwitch

	Return Default
EndFunc


; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name...........: __adj_displayTable
; Description ...: Renders a table (params or obs) based on column configuration
; Syntax.........: __adj_displayTable($mRes, $sCols, $sTableType, $iPrecision, ByRef $mSystem)
; Parameters ....: $mRes        - Flattened results from _adj_getResults
;                  $sCols       - Pipe-separated column names (e.g. "name|value|sdx")
;                  $sTableType  - "param" or "obs"
;                  $iPrecision  - Number of significant digits
;                  $mSystem     - [ByRef] The adjustment system
; Return values .: Formatted table string
; Author ........: AspirinJunkie
; Modified.......:
; Remarks .......: Name column: 20-char left-aligned. Data columns: 12-char right-aligned.
;                  Values formatted with StringFormat("%12.<precision>g").
; Related .......: __adj_displayGetValue, __adj_displayColLabel
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __adj_displayTable($mRes, $sCols, $sTableType, $iPrecision, ByRef $mSystem)
	Local $aCols = StringSplit($sCols, "|", 2)
	Local $iNCols = UBound($aCols)
	Local $sRet = ""
	Local $sFmt = "%12." & $iPrecision & "g"

	; --- Header line ---
	Local $sHeader = ""
	For $i = 0 To $iNCols - 1
		Local $sLabel = __adj_displayColLabel($aCols[$i])
		If $aCols[$i] = "name" Then
			$sHeader &= StringFormat("%-20s", $sLabel)
		Else
			$sHeader &= StringFormat(" %12s", $sLabel)
		EndIf
	Next
	$sRet &= $sHeader & @CRLF

	; --- Determine row keys ---
	Local $mKeys
	If $sTableType = "param" Then
		$mKeys = $mRes.x1
	Else
		$mKeys = $mRes.v
	EndIf

	; --- Data rows ---
	For $sKey In MapKeys($mKeys)
		Local $sLine = ""
		For $j = 0 To $iNCols - 1
			If $aCols[$j] = "name" Then
				$sLine &= StringFormat("%-20s", $sKey)
			Else
				Local $vVal = __adj_displayGetValue($mRes, $aCols[$j], $sKey, $sTableType, $mSystem)
				If IsKeyword($vVal) Then
					$sLine &= StringFormat(" %12s", "---")
				ElseIf IsString($vVal) Then
					$sLine &= StringFormat(" %12s", $vVal)
				Else
					$sLine &= " " & StringFormat($sFmt, $vVal)
				EndIf
			EndIf
		Next
		$sRet &= $sLine & @CRLF
	Next

	Return $sRet
EndFunc

#EndRegion ; Display internals
