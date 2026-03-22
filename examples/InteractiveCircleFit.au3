#AutoIt3Wrapper_Au3Check_Parameters=-d -w 1 -w 2 -w 3 -w 4 -w 5 -w 6 -w 7

; ===============================================================================
;   Interactive Circle Fit — Least-Squares Adjustment Demo
; ===============================================================================
;
;   Left-click   = add point
;   Right-click  = remove nearest point (within snap distance)
;   Clear button = remove all points
;
;   With 3+ points a best-fit circle is computed via Gauss-Helmert model
;   (implicit condition equations) and drawn together with residual vectors.
;
;   Model:       WGLM  (each point contributes X and Y as observations)
;   Condition:   (#Xi - XM)^2 + (#Yi - YM)^2 - R^2 = 0
;   Parameters:  XM (center X), YM (center Y), R (radius)
;   Algorithm:   Levenberg-Marquardt
;
; ===============================================================================

#include <GUIConstantsEx.au3>
#include <GDIPlus.au3>
#include <WindowsConstants.au3>
#include "..\Adjustment.au3"

; --- Layout ---
Global Const $g_iCanvasW  = 700
Global Const $g_iCanvasH  = 550
Global Const $g_iPanelW   = 230
Global Const $g_iSnapDist = 15 ; [px] snap distance for right-click removal

; --- Point storage ---
Global $g_aPointsX[0]
Global $g_aPointsY[0]

; --- Adjustment results ---
Global $g_mResults = Null
Global $g_bAdjOK   = False

; --- Covariance data for error ellipse (Cxx = s0² · Qxx) ---
Global $g_fCovXX = 0 ; Var(XM)
Global $g_fCovYY = 0 ; Var(YM)
Global $g_fCovXY = 0 ; Cov(XM, YM)

; --- Robust estimation: per-point weights (0 = outlier, 1 = inlier) ---
Global $g_aRobWeights[0]

; ===== GUI ====================================================================
Global $g_hGUI = GUICreate("Kreisausgleichung", _
	$g_iCanvasW + $g_iPanelW, $g_iCanvasH)

; Separator line
GUICtrlCreateGraphic($g_iCanvasW, 0, 1, $g_iCanvasH)
GUICtrlSetBkColor(-1, 0xA0A0A0)

; Panel: title
GUICtrlCreateLabel("Ergebnisse", $g_iCanvasW + 15, 12, 200, 22)
GUICtrlSetFont(-1, 11, 700)

; Panel: results
Global $g_idLblResults = GUICtrlCreateLabel("", _
	$g_iCanvasW + 15, 42, $g_iPanelW - 30, 280)
GUICtrlSetFont(-1, 10, 400, 0, "Consolas")

; Panel: checkboxes
Global $g_idChkResiduals = GUICtrlCreateCheckbox("Residuen anzeigen", _
	$g_iCanvasW + 15, $g_iCanvasH - 180, $g_iPanelW - 30, 20)
Global $g_idChkEllipse = GUICtrlCreateCheckbox("Fehlerellipse (95%)", _
	$g_iCanvasW + 15, $g_iCanvasH - 155, $g_iPanelW - 30, 20)
Global $g_idChkRobust = GUICtrlCreateCheckbox("Robuste Ausgleichung", _
	$g_iCanvasW + 15, $g_iCanvasH - 130, $g_iPanelW - 30, 20)

; Panel: usage hint
GUICtrlCreateLabel( _
	"Linksklick  = Punkt setzen" & @CRLF & _
	"Rechtsklick = Punkt entfernen", _
	$g_iCanvasW + 15, $g_iCanvasH - 100, $g_iPanelW - 30, 40)
GUICtrlSetColor(-1, 0x808080)

; Panel: clear button
Global $g_idBtnClear = GUICtrlCreateButton("Alle entfernen", _
	$g_iCanvasW + 15, $g_iCanvasH - 50, $g_iPanelW - 30, 30)

; ===== GDI+ ===================================================================
_GDIPlus_Startup()

Global $g_hGraphics = _GDIPlus_GraphicsCreateFromHWND($g_hGUI)
Global $g_hBitmap   = _GDIPlus_BitmapCreateFromGraphics($g_iCanvasW, $g_iCanvasH, $g_hGraphics)
Global $g_hBuffer   = _GDIPlus_ImageGetGraphicsContext($g_hBitmap)
_GDIPlus_GraphicsSetSmoothingMode($g_hBuffer, 4) ; AntiAlias

; Pens & Brushes
Global $g_hPenCircle   = _GDIPlus_PenCreate(0xFF2266BB, 2)     ; blue  — fitted circle
Global $g_hPenResidual = _GDIPlus_PenCreate(0xFFDD6644, 1.5)   ; orange — residual lines
Global $g_hPenCenter   = _GDIPlus_PenCreate(0xFFCC2222, 2)     ; red    — center cross
Global $g_hBrushPoint  = _GDIPlus_BrushCreateSolid(0xFF333333) ; dark   — observed points
Global $g_hBrushAdj    = _GDIPlus_BrushCreateSolid(0xFF44AA44) ; green  — adjusted points
Global $g_hBrushHint   = _GDIPlus_BrushCreateSolid(0xFFAAAAAA) ; gray   — hint text
Global $g_hPenOutlier  = _GDIPlus_PenCreate(0xFFBBBBBB, 1.5)  ; light gray — outlier residuals
Global $g_hBrushOutlier = _GDIPlus_BrushCreateSolid(0xFFBBBBBB) ; light gray — outlier points

_GDIPlus_PenSetDashStyle($g_hPenResidual, 2) ; Dot

Global $g_hPenEllipse = _GDIPlus_PenCreate(0xFFAA22CC, 1.5) ; purple — error ellipse
_GDIPlus_PenSetDashStyle($g_hPenEllipse, 1) ; Dash

; Font for canvas hint text
Global $g_hFontFamily = _GDIPlus_FontFamilyCreate("Segoe UI")
Global $g_hFont       = _GDIPlus_FontCreate($g_hFontFamily, 12)
Global $g_hFormat     = _GDIPlus_StringFormatCreate()
_GDIPlus_StringFormatSetAlign($g_hFormat, 1) ; Center

; Small font + brush for residual labels on canvas
Global $g_hFontSmall    = _GDIPlus_FontCreate($g_hFontFamily, 8)
Global $g_hFormatLeft   = _GDIPlus_StringFormatCreate() ; left-aligned (default)
Global $g_hBrushResLabel = _GDIPlus_BrushCreateSolid(0xFFBB5522) ; dark orange

; ===== Message handlers =======================================================
GUIRegisterMsg($WM_LBUTTONDOWN, "_WM_LBUTTONDOWN")
GUIRegisterMsg($WM_RBUTTONDOWN, "_WM_RBUTTONDOWN")
GUIRegisterMsg($WM_PAINT,       "_WM_PAINT")
GUIRegisterMsg($WM_ERASEBKGND,  "_WM_ERASEBKGND")

GUISetState(@SW_SHOW, $g_hGUI)
_UpdateResultLabel()
_Redraw()

; ===== Main loop ==============================================================
While 1
	Switch GUIGetMsg()
		Case $GUI_EVENT_CLOSE
			ExitLoop
		Case $g_idChkResiduals, $g_idChkEllipse
			_Redraw()
		Case $g_idChkRobust
			_RunAdjustment()
			_Redraw()
		Case $g_idBtnClear
			ReDim $g_aPointsX[0]
			ReDim $g_aPointsY[0]
			ReDim $g_aRobWeights[0]
			$g_bAdjOK = False
			$g_mResults = Null
			_UpdateResultLabel()
			_Redraw()
	EndSwitch
WEnd

_Cleanup()

; ===== Windows message callbacks ==============================================

Func _WM_LBUTTONDOWN($hWnd, $iMsg, $wParam, $lParam)
	#forceref $iMsg, $wParam
	If $hWnd <> $g_hGUI Then Return $GUI_RUNDEFMSG

	Local $iX = BitAND($lParam, 0xFFFF)
	Local $iY = BitAND(BitShift($lParam, 16), 0xFFFF)
	If $iX >= $g_iCanvasW Or $iY >= $g_iCanvasH Then Return $GUI_RUNDEFMSG

	; Add point
	Local $iN = UBound($g_aPointsX)
	ReDim $g_aPointsX[$iN + 1]
	ReDim $g_aPointsY[$iN + 1]
	$g_aPointsX[$iN] = $iX
	$g_aPointsY[$iN] = $iY

	_RunAdjustment()
	_Redraw()
	Return $GUI_RUNDEFMSG
EndFunc

Func _WM_RBUTTONDOWN($hWnd, $iMsg, $wParam, $lParam)
	#forceref $iMsg, $wParam
	If $hWnd <> $g_hGUI Then Return $GUI_RUNDEFMSG

	Local $iX = BitAND($lParam, 0xFFFF)
	Local $iY = BitAND(BitShift($lParam, 16), 0xFFFF)
	If $iX >= $g_iCanvasW Or $iY >= $g_iCanvasH Then Return $GUI_RUNDEFMSG

	Local $iIdx = __FindNearestPoint($iX, $iY)
	If $iIdx < 0 Then Return $GUI_RUNDEFMSG

	; Remove point (swap with last, then shrink)
	Local $iLast = UBound($g_aPointsX) - 1
	If $iIdx < $iLast Then
		$g_aPointsX[$iIdx] = $g_aPointsX[$iLast]
		$g_aPointsY[$iIdx] = $g_aPointsY[$iLast]
	EndIf
	ReDim $g_aPointsX[$iLast]
	ReDim $g_aPointsY[$iLast]

	_RunAdjustment()
	_Redraw()
	Return $GUI_RUNDEFMSG
EndFunc

Func _WM_PAINT($hWnd, $iMsg, $wParam, $lParam)
	#forceref $iMsg, $wParam, $lParam
	If $hWnd = $g_hGUI Then
		_GDIPlus_GraphicsDrawImageRect($g_hGraphics, $g_hBitmap, 0, 0, $g_iCanvasW, $g_iCanvasH)
	EndIf
	Return $GUI_RUNDEFMSG
EndFunc

Func _WM_ERASEBKGND($hWnd, $iMsg, $wParam, $lParam)
	#forceref $iMsg, $lParam
	If $hWnd <> $g_hGUI Then Return $GUI_RUNDEFMSG

	; Only fill the panel area — leave the canvas untouched to avoid flicker
	Local $tRect = DllStructCreate("long;long;long;long")
	DllStructSetData($tRect, 1, $g_iCanvasW) ; left
	DllStructSetData($tRect, 2, 0)           ; top
	DllStructSetData($tRect, 3, $g_iCanvasW + $g_iPanelW) ; right
	DllStructSetData($tRect, 4, $g_iCanvasH)              ; bottom
	Local $aBrush = DllCall("user32.dll", "handle", "GetSysColorBrush", "int", 15) ; COLOR_3DFACE
	If Not @error And IsArray($aBrush) Then
		DllCall("user32.dll", "int", "FillRect", "handle", $wParam, "struct*", $tRect, "handle", $aBrush[0])
	EndIf
	Return 1 ; "background erased" — prevents default erase over canvas
EndFunc

; ===== Helper =================================================================

Func __FindNearestPoint($iX, $iY)
	Local $iN = UBound($g_aPointsX)
	If $iN = 0 Then Return -1

	Local $fMinDist = $g_iSnapDist + 1
	Local $iMinIdx  = -1
	For $i = 0 To $iN - 1
		Local $fDist = Sqrt(($g_aPointsX[$i] - $iX) ^ 2 + ($g_aPointsY[$i] - $iY) ^ 2)
		If $fDist < $fMinDist Then
			$fMinDist = $fDist
			$iMinIdx  = $i
		EndIf
	Next
	Return $iMinIdx
EndFunc

; --- Error ellipse drawing ---

; Draws the 95% confidence ellipse for the center point (XM, YM)
; from the 2x2 covariance submatrix elements.
;
; Math:  eigenvalues of Cxx give principal variances,
;        rotation from eigenvectors via atan2.
;        Semi-axes scaled by sqrt(chi²(0.95, 2)) = sqrt(5.991).
Func __DrawErrorEllipse($fCX, $fCY, $fVarX, $fVarY, $fCovXY)
	Local Const $fChi2_95_2 = 5.991 ; chi²(0.95, df=2)
	Local Const $fPI2 = 6.28318530717959
	Local Const $iSegments = 72

	; Eigenvalues of 2x2 covariance matrix
	Local $fTrace = $fVarX + $fVarY
	Local $fDiff  = $fVarX - $fVarY
	Local $fDisc  = Sqrt($fDiff ^ 2 + 4 * $fCovXY ^ 2)
	Local $fLam1  = ($fTrace + $fDisc) / 2 ; larger
	Local $fLam2  = ($fTrace - $fDisc) / 2 ; smaller
	If $fLam2 < 0 Then $fLam2 = 0 ; numerical guard

	; Semi-axes (95% confidence)
	Local $fA = Sqrt($fChi2_95_2 * $fLam1) ; semi-major
	Local $fB = Sqrt($fChi2_95_2 * $fLam2) ; semi-minor

	; Rotation angle of semi-major axis
	Local $fTheta = __ATan2(2 * $fCovXY, $fDiff) / 2

	; Draw ellipse as polyline (72 segments)
	Local $fCosT = Cos($fTheta), $fSinT = Sin($fTheta)
	Local $fPrevX = $fCX + $fA * $fCosT ; t=0
	Local $fPrevY = $fCY + $fA * $fSinT

	For $j = 1 To $iSegments
		Local $fT = $j * $fPI2 / $iSegments
		Local $fCosJ = Cos($fT), $fSinJ = Sin($fT)
		Local $fCurX = $fCX + $fA * $fCosJ * $fCosT - $fB * $fSinJ * $fSinT
		Local $fCurY = $fCY + $fA * $fCosJ * $fSinT + $fB * $fSinJ * $fCosT
		_GDIPlus_GraphicsDrawLine($g_hBuffer, $fPrevX, $fPrevY, $fCurX, $fCurY, $g_hPenEllipse)
		$fPrevX = $fCurX
		$fPrevY = $fCurY
	Next
EndFunc

Func __ATan2($fY, $fX)
	If $fX > 0 Then Return ATan($fY / $fX)
	If $fX < 0 And $fY >= 0 Then Return ATan($fY / $fX) + ACos(-1)
	If $fX < 0 And $fY < 0 Then Return ATan($fY / $fX) - ACos(-1)
	If $fY > 0 Then Return ACos(-1) / 2
	If $fY < 0 Then Return -(ACos(-1) / 2)
	Return 0
EndFunc

; ===== Circle adjustment ======================================================

Func _RunAdjustment()
	$g_bAdjOK = False
	$g_mResults = Null

	Local $iN = UBound($g_aPointsX)
	If $iN < 3 Then
		_UpdateResultLabel()
		Return
	EndIf

	; --- Build system ---
	Local $mSystem = _adj_createSystem()

	; Observations: Xi, Yi per point (unit weight: sigma = 1 px)
	For $i = 0 To $iN - 1
		_adj_addObs($mSystem, "X" & $i, $g_aPointsX[$i], 1.0)
		_adj_addObs($mSystem, "Y" & $i, $g_aPointsY[$i], 1.0)
	Next

	; Condition equations:  (#Xi - XM)^2 + (#Yi - YM)^2 - R^2 = 0
	For $i = 0 To $iN - 1
		_adj_addFunction($mSystem, _
			"(#X" & $i & " - XM)^2 + (#Y" & $i & " - YM)^2 - R^2")
	Next

	; Initial values: centroid + mean distance
	Local $fCX = 0, $fCY = 0
	For $i = 0 To $iN - 1
		$fCX += $g_aPointsX[$i]
		$fCY += $g_aPointsY[$i]
	Next
	$fCX /= $iN
	$fCY /= $iN

	Local $fR = 0
	For $i = 0 To $iN - 1
		$fR += Sqrt(($g_aPointsX[$i] - $fCX) ^ 2 + ($g_aPointsY[$i] - $fCY) ^ 2)
	Next
	$fR /= $iN
	If $fR < 1 Then $fR = 50 ; guard against degenerate configuration

	_adj_setInitialValue($mSystem, "XM", $fCX)
	_adj_setInitialValue($mSystem, "YM", $fCY)
	_adj_setInitialValue($mSystem, "R",  $fR)

	; Solve
	Local $bRobust = BitAND(GUICtrlRead($g_idChkRobust), 1) ? True : False
	Local $mConfig = _adj_defaultConfig($bRobust ? "GN" : "LM", False)
	If $bRobust Then
		$mConfig.robust           = "Huber"
		$mConfig.robustParams     = _adj_robustDefaults("Huber")
		$mConfig.robustMaxIter    = 10     ; default 30 — fewer IRLS cycles
		$mConfig.robustConvergence = 0.01  ; default 1e-3 — looser convergence
		$mConfig.tolerance        = 1e-6   ; default 1e-10 — looser inner tolerance
		$mConfig.deriveMethod     = "Forward" ; default "Central" — halves Jacobian cost
	EndIf
	_adj_solve($mSystem, $mConfig)
	If @error Then
		_UpdateResultLabel()
		Return
	EndIf

	$g_mResults = _adj_getResults($mSystem)
	$g_bAdjOK = True

	; Extract per-point robust weights (min of X and Y weight per point)
	ReDim $g_aRobWeights[$iN]
	If $mConfig.robust <> "" And IsMap($g_mResults.robustWeights) Then
		Local $mRobW = $g_mResults.robustWeights
		For $i = 0 To $iN - 1
			Local $fWx = $mRobW["X" & $i]
			Local $fWy = $mRobW["Y" & $i]
			$g_aRobWeights[$i] = ($fWx < $fWy) ? $fWx : $fWy
		Next
	Else
		For $i = 0 To $iN - 1
			$g_aRobWeights[$i] = 1.0
		Next
	EndIf

	; Extract 2x2 covariance submatrix for XM, YM (needed for error ellipse)
	$g_fCovXX = 0
	$g_fCovYY = 0
	$g_fCovXY = 0
	If $g_mResults.f > 0 Then
		Local $mQxx  = $g_mResults.Qxx
		Local $mIdx  = $mSystem.state.idxParams
		Local $iIxXM = $mIdx["XM"]
		Local $iIxYM = $mIdx["YM"]
		Local $iRows = $mQxx.rows
		Local $fS0sq = $g_mResults.s0 ^ 2
		; Cxx = s0² · Qxx
		$g_fCovXX = $fS0sq * DllStructGetData($mQxx.struct, 1, $iIxXM + $iIxXM * $iRows + 1)
		$g_fCovYY = $fS0sq * DllStructGetData($mQxx.struct, 1, $iIxYM + $iIxYM * $iRows + 1)
		$g_fCovXY = $fS0sq * DllStructGetData($mQxx.struct, 1, $iIxXM + $iIxYM * $iRows + 1)
	EndIf

	_UpdateResultLabel()
EndFunc

; ===== Result display =========================================================

Func _UpdateResultLabel()
	Local $iN = UBound($g_aPointsX)
	Local $sText = "Punkte: " & $iN & @CRLF

	If $iN < 3 Then
		$sText &= @CRLF & "(mind. 3 Punkte)"
		GUICtrlSetData($g_idLblResults, $sText)
		Return
	EndIf

	If Not $g_bAdjOK Then
		$sText &= @CRLF & "Ausgleichung" & @CRLF & "fehlgeschlagen"
		GUICtrlSetData($g_idLblResults, $sText)
		Return
	EndIf

	Local $mX   = $g_mResults.x1
	Local $mSdx = $g_mResults.sdx
	Local $fXM  = $mX["XM"]
	Local $fYM  = $mX["YM"]
	Local $fR   = $mX["R"]
	Local $fS0  = $g_mResults.s0
	Local $iF   = $g_mResults.f

	;~ $sText &= "Modell: " & $g_mResults.modelType & @CRLF
	$sText &= @CRLF

	If $iF > 0 And IsMap($mSdx) Then
		$sText &= StringFormat("XM =%8.1f %s%.2f", $fXM, Chr(177), $mSdx["XM"]) & @CRLF
		$sText &= StringFormat("YM =%8.1f %s%.2f", $fYM, Chr(177), $mSdx["YM"]) & @CRLF
		$sText &= StringFormat("R  =%8.1f %s%.2f", $fR, Chr(177), $mSdx["R"]) & @CRLF
	Else
		$sText &= StringFormat("XM =%8.1f", $fXM) & @CRLF
		$sText &= StringFormat("YM =%8.1f", $fYM) & @CRLF
		$sText &= StringFormat("R  =%8.1f", $fR) & @CRLF
	EndIf
	$sText &= @CRLF

	If $iF > 0 Then
		$sText &= StringFormat("s0 =%8.1f px", $fS0) & @CRLF
	Else
		$sText &= "s0 = ---  (exakt)" & @CRLF
	EndIf
	$sText &= "f  = " & $iF

	; Outlier count (only when robust is active)
	If UBound($g_aRobWeights) > 0 Then
		Local $iOutliers = 0
		For $i = 0 To UBound($g_aRobWeights) - 1
			If $g_aRobWeights[$i] < 0.5 Then $iOutliers += 1
		Next
		If $iOutliers > 0 Then
			$sText &= @CRLF & @CRLF & "Ausreisser: " & $iOutliers
		EndIf
	EndIf

	GUICtrlSetData($g_idLblResults, $sText)
EndFunc

; ===== Drawing ================================================================

Func _Redraw()
	; Clear canvas
	_GDIPlus_GraphicsClear($g_hBuffer, 0xFFF8F8F8)

	Local $iN = UBound($g_aPointsX)

	; Hint text on empty canvas
	If $iN = 0 Then
		Local $tLayout = _GDIPlus_RectFCreate(0, $g_iCanvasH / 2 - 15, $g_iCanvasW, 30)
		_GDIPlus_GraphicsDrawStringEx($g_hBuffer, _
			"Klicke um Punkte zu setzen", $g_hFont, $tLayout, $g_hFormat, $g_hBrushHint)
	EndIf

	; --- Circle, residuals, center (only when adjustment succeeded) ---
	If $g_bAdjOK Then
		Local $mX  = $g_mResults.x1
		Local $fXM = $mX["XM"]
		Local $fYM = $mX["YM"]
		Local $fR  = $mX["R"]

		; Fitted circle
		_GDIPlus_GraphicsDrawEllipse($g_hBuffer, _
			$fXM - $fR, $fYM - $fR, $fR * 2, $fR * 2, $g_hPenCircle)

		; Residual lines + adjusted positions on the circle
		For $i = 0 To $iN - 1
			; Project observed point onto fitted circle (guaranteed to lie on it)
			Local $fDX = $g_aPointsX[$i] - $fXM
			Local $fDY = $g_aPointsY[$i] - $fYM
			Local $fDist = Sqrt($fDX ^ 2 + $fDY ^ 2)
			Local $fAdjX, $fAdjY
			If $fDist > 0 Then
				$fAdjX = $fXM + $fR * $fDX / $fDist
				$fAdjY = $fYM + $fR * $fDY / $fDist
			Else
				$fAdjX = $fXM + $fR
				$fAdjY = $fYM
			EndIf

			; Outlier detection: weight < 0.5 → gray styling
			Local $bOutlier = ($i < UBound($g_aRobWeights) And $g_aRobWeights[$i] < 0.5)

			; Line: observed point -> projected position on circle
			_GDIPlus_GraphicsDrawLine($g_hBuffer, _
				$g_aPointsX[$i], $g_aPointsY[$i], $fAdjX, $fAdjY, _
				$bOutlier ? $g_hPenOutlier : $g_hPenResidual)

			; Adjusted point (green / gray for outlier)
			_GDIPlus_GraphicsFillEllipse($g_hBuffer, $fAdjX - 3, $fAdjY - 3, 6, 6, _
				$bOutlier ? $g_hBrushOutlier : $g_hBrushAdj)

			; Residual label (if checkbox is checked)
			If BitAND(GUICtrlRead($g_idChkResiduals), 1) Then
				; Combined radial residual: positive = outside circle, negative = inside
				Local $fV = $fDist - $fR
				Local $sLabel = "v=" & StringFormat("%.1f", $fV)
				; Position: midpoint of residual line + perpendicular offset
				Local $fMidX = ($g_aPointsX[$i] + $fAdjX) / 2
				Local $fMidY = ($g_aPointsY[$i] + $fAdjY) / 2
				Local $fOffX = 0, $fOffY = 0
				If $fDist > 0 Then
					$fOffX = -($fDY / $fDist) * 10
					$fOffY = ($fDX / $fDist) * 10
				EndIf
				Local $tLabelRect = _GDIPlus_RectFCreate($fMidX + $fOffX, $fMidY + $fOffY - 7, 70, 14)
				_GDIPlus_GraphicsDrawStringEx($g_hBuffer, $sLabel, _
					$g_hFontSmall, $tLabelRect, $g_hFormatLeft, $g_hBrushResLabel)
			EndIf
		Next

		; Center cross
		_GDIPlus_GraphicsDrawLine($g_hBuffer, $fXM - 8, $fYM, $fXM + 8, $fYM, $g_hPenCenter)
		_GDIPlus_GraphicsDrawLine($g_hBuffer, $fXM, $fYM - 8, $fXM, $fYM + 8, $g_hPenCenter)

		; Error ellipse (95% confidence, only when f > 0 and checkbox checked)
		If BitAND(GUICtrlRead($g_idChkEllipse), 1) And $g_fCovXX > 0 Then
			__DrawErrorEllipse($fXM, $fYM, $g_fCovXX, $g_fCovYY, $g_fCovXY)
		EndIf
	EndIf

	; --- Observed points (always on top) ---
	For $i = 0 To $iN - 1
		Local $bIsOutlier = ($g_bAdjOK And $i < UBound($g_aRobWeights) And $g_aRobWeights[$i] < 0.5)
		_GDIPlus_GraphicsFillEllipse($g_hBuffer, _
			$g_aPointsX[$i] - 5, $g_aPointsY[$i] - 5, 10, 10, _
			$bIsOutlier ? $g_hBrushOutlier : $g_hBrushPoint)
	Next

	; Blit buffer to screen
	_GDIPlus_GraphicsDrawImageRect($g_hGraphics, $g_hBitmap, 0, 0, $g_iCanvasW, $g_iCanvasH)
EndFunc

; ===== Cleanup ================================================================

Func _Cleanup()
	_GDIPlus_FontDispose($g_hFontSmall)
	_GDIPlus_FontDispose($g_hFont)
	_GDIPlus_FontFamilyDispose($g_hFontFamily)
	_GDIPlus_StringFormatDispose($g_hFormatLeft)
	_GDIPlus_StringFormatDispose($g_hFormat)

	_GDIPlus_PenDispose($g_hPenCircle)
	_GDIPlus_PenDispose($g_hPenResidual)
	_GDIPlus_PenDispose($g_hPenCenter)
	_GDIPlus_PenDispose($g_hPenEllipse)
	_GDIPlus_PenDispose($g_hPenOutlier)
	_GDIPlus_BrushDispose($g_hBrushOutlier)
	_GDIPlus_BrushDispose($g_hBrushPoint)
	_GDIPlus_BrushDispose($g_hBrushAdj)
	_GDIPlus_BrushDispose($g_hBrushHint)
	_GDIPlus_BrushDispose($g_hBrushResLabel)

	_GDIPlus_GraphicsDispose($g_hBuffer)
	_GDIPlus_BitmapDispose($g_hBitmap)
	_GDIPlus_GraphicsDispose($g_hGraphics)
	_GDIPlus_Shutdown()
EndFunc
