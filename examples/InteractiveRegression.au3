#AutoIt3Wrapper_Au3Check_Parameters=-d -w 1 -w 2 -w 3 -w 4 -w 5 -w 6 -w 7

; ===============================================================================
;   Interactive Regression — Least-Squares Curve Fitting Tool
; ===============================================================================
;
;   Supports:
;     - Normal regression   (OLS/WLS — only Y has error)
;     - Orthogonal regression (TLS — X and Y have equal error)
;     - Deming regression   (X and Y have different constant error)
;     - York regression     (each measurement has individual errors)
;     - Robust estimation   (Huber, Hampel, Biweight, L1, BIBER)
;
;   Formula syntax:  A * #X + B   (# = observation, bare = parameter)
;   Dependent variable is implicit (default "Y"), configurable.
;
; ===============================================================================

#include <GUIConstantsEx.au3>
#include <GDIPlus.au3>
#include <WindowsConstants.au3>
#include <GuiListView.au3>
#include <ListViewConstants.au3>
#include <EditConstants.au3>
#include <StaticConstants.au3>
#include <ComboConstants.au3>
#include <WinAPISys.au3>
#include "..\Adjustment.au3"
; ===============================================================================
;   Inlined _CellGrid.au3 — editable grid of input/label controls
; ===============================================================================
#include <GUIConstantsEx.au3>
#include <WindowsConstants.au3>
#include <WinAPISys.au3>
#include <WinAPIShellEx.au3>
#include <GuiScrollBars.au3>
#include <ScrollBarConstants.au3>
#include <EditConstants.au3>
#include <StaticConstants.au3>
#include <Clipboard.au3>

; ===============================================================================
;   _CellGrid.au3 — Editierbares Grid aus echten Input-/Label-Controls
; ===============================================================================
;   Public:
;     __cg_create($hParent, $iX, $iY, $iW, $iH, $aSchema, _
;                 $sCommitCallback = "", $sBulkCommitCallback = "") → iGrid
;     __cg_setSchema($iGrid, $aSchema)
;     __cg_setRowCount($iGrid, $iRows)
;     __cg_setCell($iGrid, $iRow, $iCol, $sValue)
;     __cg_getCell($iGrid, $iRow, $iCol) → String
;     __cg_focusCell($iGrid, $iRow, $iCol)
;     __cg_getFocus($iGrid) → Map {row, col} or Null
;     __cg_destroy($iGrid)
; ===============================================================================

; ---- Virtual-Key- und Window-Message-Konstanten (Subclass-Proc) ---------------
Global Const $__cg_VK_TAB    = 0x09
Global Const $__cg_VK_RETURN = 0x0D
Global Const $__cg_VK_ESCAPE = 0x1B
Global Const $__cg_VK_SPACE  = 0x20
Global Const $__cg_VK_LEFT   = 0x25
Global Const $__cg_VK_UP     = 0x26
Global Const $__cg_VK_RIGHT  = 0x27
Global Const $__cg_VK_DOWN   = 0x28
Global Const $__cg_WM_KEYDOWN    = 0x0100
Global Const $__cg_WM_CHAR       = 0x0102
Global Const $__cg_WM_SETFOCUS   = 0x0007
Global Const $__cg_WM_KILLFOCUS  = 0x0008
Global Const $__cg_WM_GETDLGCODE = 0x0087
Global Const $__cg_DLGC_WANTTAB     = 0x0002
Global Const $__cg_DLGC_WANTARROWS  = 0x0001
Global Const $__cg_DLGC_WANTCHARS   = 0x0080

Global $g_aCellGrids[0]     ; Registry: Array of Maps (grid instances)

Global Const $__cg_iSubclassID = 0xCE11   ; willkürliche, aber eindeutige ID
Global $__cg_hSubclassProc = 0             ; DllCallback-Handle
Global $__cg_pSubclassProc = 0             ; Funktionspointer
Global $__cg_mCtrlLookup[]                 ; Map: HWND-String → "iGrid|iRow|iCol"

Func __cg_create($hParent, $iX, $iY, $iW, $iH, $aSchema, _
                 $sCommitCallback = "", $sBulkCommitCallback = "")
    Local $mGrid[]
    $mGrid.hParent             = $hParent
    $mGrid.iX                  = $iX
    $mGrid.iY                  = $iY
    $mGrid.iW                  = $iW
    $mGrid.iH                  = $iH
    $mGrid.iRowH               = 22
    $mGrid.iHeaderH            = 22
    $mGrid.iRowCount           = 0
    $mGrid.iFocusRow           = -1
    $mGrid.iFocusCol           = -1
    $mGrid.iScrollOffset       = 0
    $mGrid.sCommitCallback     = $sCommitCallback
    $mGrid.sBulkCommitCallback = $sBulkCommitCallback
    $mGrid.aSchema             = $aSchema

    Local $aEmpty2D[0][0]
    $mGrid.aControls = $aEmpty2D
    $mGrid.aValues   = $aEmpty2D

    ; Header-GUI (fix, oberste $iHeaderH Pixel)
    $mGrid.hHeaderGui = GUICreate("", $iW, $mGrid.iHeaderH, $iX, $iY, _
                                  $WS_CHILD, 0, $hParent)
    GUISetBkColor(0xE0E0E0, $mGrid.hHeaderGui)
    GUISetState(@SW_SHOW, $mGrid.hHeaderGui)

    ; Scroll-Body-GUI (darunter, mit vertikaler Scrollbar)
    Local $iBodyH = $iH - $mGrid.iHeaderH
    $mGrid.hScrollGui = GUICreate("", $iW, $iBodyH, $iX, $iY + $mGrid.iHeaderH, _
                                  BitOR($WS_CHILD, $WS_VSCROLL, $WS_CLIPCHILDREN, $WS_CLIPSIBLINGS), _
                                  0, $hParent)
    GUISetBkColor(0xFFFFFF, $mGrid.hScrollGui)
    GUISetState(@SW_SHOW, $mGrid.hScrollGui)

    ; Header-Zellen gemäß Schema
    __cg_buildHeader($mGrid)

    GUISwitch($hParent)   ; Caller-GUI-Kontext wiederherstellen

    ; Einmalige globale Registrierung des WM_VSCROLL-Handlers
    Static Local $bVscrollRegistered = False
    If Not $bVscrollRegistered Then
        GUIRegisterMsg(0x0115, "__cg_wmVscroll")   ; WM_VSCROLL = 0x0115
        $bVscrollRegistered = True
    EndIf

    ; Einmalige globale Registrierung des WM_MOUSEWHEEL-Handlers
    Static Local $bMwheelRegistered = False
    If Not $bMwheelRegistered Then
        GUIRegisterMsg(0x020A, "__cg_wmMousewheel")   ; WM_MOUSEWHEEL = 0x020A
        $bMwheelRegistered = True
    EndIf

    ; In Registry ablegen
    ReDim $g_aCellGrids[UBound($g_aCellGrids) + 1]
    Local $iGrid = UBound($g_aCellGrids) - 1
    $g_aCellGrids[$iGrid] = $mGrid

    Return $iGrid
EndFunc

Func __cg_destroy($iGrid)
    If $iGrid < 0 Or $iGrid >= UBound($g_aCellGrids) Then Return
    Local $mGrid = $g_aCellGrids[$iGrid]
    If Not IsMap($mGrid) Then Return

    ; Alle Cell-Subclasses aufheben
    Local $aControls = $mGrid.aControls
    Local $aSchema = $mGrid.aSchema
    For $iRow = 0 To $mGrid.iRowCount - 1
        For $iCol = 0 To UBound($aSchema) - 1
            Local $iId = $aControls[$iRow][$iCol]
            If $iId > 0 Then
                Local $hCtrl = GUICtrlGetHandle($iId)
                If $hCtrl <> 0 Then __cg_unsubclassInput($hCtrl)
            EndIf
        Next
    Next

    If $mGrid.hHeaderGui Then GUIDelete($mGrid.hHeaderGui)
    If $mGrid.hScrollGui Then GUIDelete($mGrid.hScrollGui)
    Local $mEmpty[]
    $g_aCellGrids[$iGrid] = $mEmpty
EndFunc

Func __cg_buildHeader(ByRef $mGrid)
    Local $aSchema = $mGrid.aSchema
    Local $hHeader = $mGrid.hHeaderGui
    GUISwitch($hHeader)

    ; Alte Header-Controls entfernen (falls Rebuild)
    If MapExists($mGrid, "aHeaderIds") Then
        Local $aOld = $mGrid.aHeaderIds
        For $i = 0 To UBound($aOld) - 1
            If $aOld[$i] > 0 Then GUICtrlDelete($aOld[$i])
        Next
    EndIf

    If UBound($aSchema) = 0 Then
        Local $aEmpty[0]
        $mGrid.aHeaderIds = $aEmpty
        GUISwitch($mGrid.hParent)
        Return
    EndIf

    Local $aHeaderIds[UBound($aSchema)]
    Local $iX = 0
    For $i = 0 To UBound($aSchema) - 1
        Local $mCol = $aSchema[$i]
        Local $iW = $mCol.width
        If $mCol.type = "separator" Then
            $aHeaderIds[$i] = GUICtrlCreateLabel("", $iX, 0, $iW, $mGrid.iHeaderH, _
                                                 $SS_SUNKEN)
        Else
            $aHeaderIds[$i] = GUICtrlCreateLabel($mCol.name, $iX + 2, 3, _
                                                 $iW - 4, $mGrid.iHeaderH - 6)
            GUICtrlSetFont($aHeaderIds[$i], 9, 700)
        EndIf
        $iX += $iW
    Next
    $mGrid.aHeaderIds = $aHeaderIds

    GUISwitch($mGrid.hParent)
EndFunc

Func __cg_setRowCount($iGrid, $iNewRows)
    If $iGrid < 0 Or $iGrid >= UBound($g_aCellGrids) Then Return
    Local $mGrid = $g_aCellGrids[$iGrid]
    If Not IsMap($mGrid) Then Return

    Local $iOldRows = $mGrid.iRowCount
    If $iNewRows = $iOldRows Then Return

    Local $aSchema   = $mGrid.aSchema
    Local $iCols     = UBound($aSchema)
    Local $aControls = $mGrid.aControls
    Local $aValues   = $mGrid.aValues
    Local $iId
    Local $hCtrl

    GUISwitch($mGrid.hScrollGui)

    If $iNewRows > $iOldRows Then
        ReDim $aControls[$iNewRows][$iCols]
        ReDim $aValues[$iNewRows][$iCols]

        ; Neue Zeilen anlegen
        For $iRow = $iOldRows To $iNewRows - 1
            Local $iX = 0
            Local $iYRel = $iRow * $mGrid.iRowH - $mGrid.iScrollOffset
            For $iCol = 0 To $iCols - 1
                Local $mColDef = $aSchema[$iCol]
                Local $iW = $mColDef.width
                $iId = __cg_createCell($mColDef.type, $iX, $iYRel, _
                                       $iW, $mGrid.iRowH)
                $aControls[$iRow][$iCol] = $iId
                $aValues[$iRow][$iCol]   = ""
                If $mColDef.type = "input" Then
                    $hCtrl = GUICtrlGetHandle($iId)
                    __cg_subclassInput($hCtrl, $iGrid, $iRow, $iCol)
                EndIf
                $iX += $iW
            Next
        Next
    Else
        ; Zeilen entfernen — vor ReDim, um Zugriff auf alte Rows zu behalten
        For $iRow = $iNewRows To $iOldRows - 1
            For $iCol = 0 To $iCols - 1
                $iId = $aControls[$iRow][$iCol]
                If $iId > 0 Then
                    $hCtrl = GUICtrlGetHandle($iId)
                    If $hCtrl <> 0 Then __cg_unsubclassInput($hCtrl)
                    GUICtrlDelete($iId)
                EndIf
            Next
        Next

        ReDim $aControls[$iNewRows][$iCols]
        ReDim $aValues[$iNewRows][$iCols]
    EndIf

    $mGrid.aControls = $aControls
    $mGrid.aValues   = $aValues
    $mGrid.iRowCount = $iNewRows

    __cg_updateScrollRange($mGrid)

    $g_aCellGrids[$iGrid] = $mGrid

    GUISwitch($mGrid.hParent)
EndFunc

Func __cg_createCell($sType, $iX, $iY, $iW, $iH)
    Local $iId
    Switch $sType
        Case "input"
            ; Kein $WS_BORDER — stattdessen dezenter $WS_EX_STATICEDGE (schmaler 1-px-Rand)
            $iId = GUICtrlCreateInput("", $iX, $iY, $iW, $iH, _
                                      $ES_AUTOHSCROLL, $WS_EX_STATICEDGE)
            GUICtrlSetFont($iId, 9, 400, 0, "Consolas")
            Return $iId
        Case "label"
            $iId = GUICtrlCreateLabel("", $iX + 2, $iY + 3, $iW - 4, $iH - 6)
            GUICtrlSetBkColor($iId, 0xF0F0F0)
            GUICtrlSetFont($iId, 9, 400, 0, "Consolas")
            Return $iId
        Case "separator"
            $iId = GUICtrlCreateLabel("", $iX, $iY, $iW, $iH)
            GUICtrlSetBkColor($iId, 0xC0C0C0)
            Return $iId
    EndSwitch
    Return -1
EndFunc

Func __cg_updateScrollRange(ByRef $mGrid)
    Local $iVirt = $mGrid.iRowCount * $mGrid.iRowH
    Local $iBodyH = $mGrid.iH - $mGrid.iHeaderH
    Local $tSI = DllStructCreate($tagSCROLLINFO)
    DllStructSetData($tSI, "cbSize", DllStructGetSize($tSI))
    DllStructSetData($tSI, "fMask", BitOR($SIF_RANGE, $SIF_PAGE))
    DllStructSetData($tSI, "nMin", 0)
    DllStructSetData($tSI, "nMax", ($iVirt > 0) ? $iVirt - 1 : 0)
    DllStructSetData($tSI, "nPage", $iBodyH)
    _GUIScrollBars_SetScrollInfo($mGrid.hScrollGui, $SB_VERT, $tSI)
EndFunc

Func __cg_setCell($iGrid, $iRow, $iCol, $sValue)
    If $iGrid < 0 Or $iGrid >= UBound($g_aCellGrids) Then Return
    Local $mGrid = $g_aCellGrids[$iGrid]
    If Not IsMap($mGrid) Then Return
    If $iRow < 0 Or $iRow >= $mGrid.iRowCount Then Return
    Local $aSchema = $mGrid.aSchema
    If $iCol < 0 Or $iCol >= UBound($aSchema) Then Return

    Local $mColDef = $aSchema[$iCol]
    If $mColDef.type = "separator" Then Return

    Local $aValues = $mGrid.aValues
    Local $aControls = $mGrid.aControls
    $aValues[$iRow][$iCol] = $sValue
    GUICtrlSetData($aControls[$iRow][$iCol], $sValue)

    $mGrid.aValues = $aValues
    $g_aCellGrids[$iGrid] = $mGrid
EndFunc

Func __cg_getCell($iGrid, $iRow, $iCol)
    If $iGrid < 0 Or $iGrid >= UBound($g_aCellGrids) Then Return ""
    Local $mGrid = $g_aCellGrids[$iGrid]
    If Not IsMap($mGrid) Then Return ""
    If $iRow < 0 Or $iRow >= $mGrid.iRowCount Then Return ""
    Local $aSchema = $mGrid.aSchema
    If $iCol < 0 Or $iCol >= UBound($aSchema) Then Return ""
    Local $mColDef = $aSchema[$iCol]
    If $mColDef.type = "separator" Then Return ""

    ; Input-Zellen: aktueller Control-Inhalt (evtl. ungespeicherte User-Eingabe)
    Local $aControls = $mGrid.aControls
    If $mColDef.type = "input" Then
        Return GUICtrlRead($aControls[$iRow][$iCol])
    EndIf
    ; Label-Zellen: .aValues-Spiegel
    Local $aValues = $mGrid.aValues
    Return $aValues[$iRow][$iCol]
EndFunc

Func __cg_setSchema($iGrid, $aNewSchema)
    If $iGrid < 0 Or $iGrid >= UBound($g_aCellGrids) Then Return
    Local $mGrid = $g_aCellGrids[$iGrid]
    If Not IsMap($mGrid) Then Return

    ; Alle bestehenden Datenzellen entsorgen
    Local $iOldRows = $mGrid.iRowCount
    Local $aControls = $mGrid.aControls
    Local $aOldSchema = $mGrid.aSchema
    Local $iId
    If $iOldRows > 0 Then
        For $iRow = 0 To $iOldRows - 1
            For $iCol = 0 To UBound($aOldSchema) - 1
                $iId = $aControls[$iRow][$iCol]
                If $iId > 0 Then
                    Local $hCtrl = GUICtrlGetHandle($iId)
                    If $hCtrl <> 0 Then __cg_unsubclassInput($hCtrl)
                    GUICtrlDelete($iId)
                EndIf
            Next
        Next
    EndIf

    ; Neues Schema eintragen; Header neu bauen; Zellen wieder auf iOldRows aufbauen
    $mGrid.aSchema = $aNewSchema
    Local $aEmpty2D[0][0]
    $mGrid.aControls = $aEmpty2D
    $mGrid.aValues   = $aEmpty2D
    $mGrid.iRowCount = 0
    $mGrid.iFocusRow = -1
    $mGrid.iFocusCol = -1

    __cg_buildHeader($mGrid)
    $g_aCellGrids[$iGrid] = $mGrid

    ; Zeilen wiederherstellen (mit neuen Spalten; Inhalt initial leer)
    If $iOldRows > 0 Then __cg_setRowCount($iGrid, $iOldRows)
EndFunc

Func __cg_focusCell($iGrid, $iRow, $iCol)
    If $iGrid < 0 Or $iGrid >= UBound($g_aCellGrids) Then Return
    Local $mGrid = $g_aCellGrids[$iGrid]
    If Not IsMap($mGrid) Then Return
    If $iRow < 0 Or $iRow >= $mGrid.iRowCount Then Return
    Local $aSchema = $mGrid.aSchema
    If $iCol < 0 Or $iCol >= UBound($aSchema) Then Return
    Local $mColDef = $aSchema[$iCol]
    If $mColDef.type <> "input" Then Return

    Local $aControls = $mGrid.aControls
    Local $iId = $aControls[$iRow][$iCol]
    If $iId <= 0 Then Return
    Local $hCtrl = GUICtrlGetHandle($iId)
    If $hCtrl = 0 Then Return

    $mGrid.iFocusRow = $iRow
    $mGrid.iFocusCol = $iCol

    ; Auto-Scroll falls Zielzelle außerhalb des Viewports
    Local $iBodyH = $mGrid.iH - $mGrid.iHeaderH
    Local $iTargetY = $iRow * $mGrid.iRowH
    Local $iTopY    = $mGrid.iScrollOffset
    Local $iBotY    = $mGrid.iScrollOffset + $iBodyH - $mGrid.iRowH
    Local $iNewOff  = $mGrid.iScrollOffset
    If $iTargetY < $iTopY Then
        $iNewOff = $iTargetY
    ElseIf $iTargetY > $iBotY Then
        $iNewOff = $iTargetY - ($iBodyH - $mGrid.iRowH)
    EndIf
    If $iNewOff <> $mGrid.iScrollOffset Then
        _GUIScrollBars_ScrollWindow($mGrid.hScrollGui, 0, $mGrid.iScrollOffset - $iNewOff)
        $mGrid.iScrollOffset = $iNewOff
        Local $tSI = DllStructCreate($tagSCROLLINFO)
        DllStructSetData($tSI, "cbSize", DllStructGetSize($tSI))
        DllStructSetData($tSI, "fMask", $SIF_POS)
        DllStructSetData($tSI, "nPos", $iNewOff)
        _GUIScrollBars_SetScrollInfo($mGrid.hScrollGui, $SB_VERT, $tSI)
    EndIf

    $g_aCellGrids[$iGrid] = $mGrid
    _WinAPI_SetFocus($hCtrl)
EndFunc

Func __cg_getFocus($iGrid)
    If $iGrid < 0 Or $iGrid >= UBound($g_aCellGrids) Then Return Null
    Local $mGrid = $g_aCellGrids[$iGrid]
    If Not IsMap($mGrid) Then Return Null
    If $mGrid.iFocusRow < 0 Then Return Null
    Local $mResult[]
    $mResult.row = $mGrid.iFocusRow
    $mResult.col = $mGrid.iFocusCol
    Return $mResult
EndFunc

Func __cg_initSubclass()
    If $__cg_pSubclassProc <> 0 Then Return
    $__cg_hSubclassProc = DllCallbackRegister("__cg_subclassProc", _
        "lresult", "hwnd;uint;wparam;lparam;uint_ptr;dword_ptr")
    $__cg_pSubclassProc = DllCallbackGetPtr($__cg_hSubclassProc)
EndFunc

Func __cg_subclassInput($hCtrl, $iGrid, $iRow, $iCol)
    If $hCtrl = 0 Then Return
    __cg_initSubclass()
    _WinAPI_SetWindowSubclass($hCtrl, $__cg_pSubclassProc, $__cg_iSubclassID, 0)
    $__cg_mCtrlLookup[String($hCtrl)] = $iGrid & "|" & $iRow & "|" & $iCol
EndFunc

Func __cg_unsubclassInput($hCtrl)
    If $hCtrl = 0 Then Return
    _WinAPI_RemoveWindowSubclass($hCtrl, $__cg_pSubclassProc, $__cg_iSubclassID)
    MapRemove($__cg_mCtrlLookup, String($hCtrl))
EndFunc

Func __cg_lookupCtrl($hCtrl)
    Local $sKey = String($hCtrl)
    If Not MapExists($__cg_mCtrlLookup, $sKey) Then Return Null
    Local $aParts = StringSplit($__cg_mCtrlLookup[$sKey], "|", 2)
    Local $mResult[]
    $mResult.iGrid = Number($aParts[0])
    $mResult.iRow  = Number($aParts[1])
    $mResult.iCol  = Number($aParts[2])
    Return $mResult
EndFunc

Func __cg_subclassProc($hCtrl, $iMsg, $wParam, $lParam, $iSubID, $pData)
    #forceref $iSubID, $pData
    Local $mLookup = __cg_lookupCtrl($hCtrl)
    If Not IsMap($mLookup) Then
        Return _WinAPI_DefSubclassProc($hCtrl, $iMsg, $wParam, $lParam)
    EndIf
    Local $iGrid = $mLookup.iGrid
    Local $iRow  = $mLookup.iRow
    Local $iCol  = $mLookup.iCol
    Local $iDir

    ; Dialog-Manager mitteilen, dass wir Tab/Pfeile/Chars selbst verarbeiten
    ; (sonst fängt der Dialog-Keyboard-Code Tab ab und routet zur nächsten Control-Chain)
    If $iMsg = $__cg_WM_GETDLGCODE Then
        Local $lDlg = _WinAPI_DefSubclassProc($hCtrl, $iMsg, $wParam, $lParam)
        Return BitOR($lDlg, $__cg_DLGC_WANTTAB, $__cg_DLGC_WANTARROWS, $__cg_DLGC_WANTCHARS)
    EndIf

    If $iMsg = $__cg_WM_KEYDOWN Then
        Switch $wParam
            Case $__cg_VK_TAB
                $iDir = (_WinAPI_GetAsyncKeyState(0x10) <> 0) ? -1 : 1  ; 0x10 = SHIFT
                If __cg_navHorizontalWrap($iGrid, $iRow, $iCol, $iDir) Then Return 0
            Case $__cg_VK_RETURN
                __cg_navEnter($iGrid, $iRow, $iCol)
                Return 0
            Case $__cg_VK_SPACE
                If __cg_navHorizontalWrap($iGrid, $iRow, $iCol, 1) Then Return 0
            Case $__cg_VK_LEFT
                If __cg_navArrowH($iGrid, $iRow, $iCol, -1) Then Return 0
            Case $__cg_VK_RIGHT
                If __cg_navArrowH($iGrid, $iRow, $iCol, +1) Then Return 0
            Case $__cg_VK_UP
                If __cg_navArrowV($iGrid, $iRow, $iCol, -1) Then Return 0
            Case $__cg_VK_DOWN
                If __cg_navArrowV($iGrid, $iRow, $iCol, +1) Then Return 0
            Case $__cg_VK_ESCAPE
                __cg_navEscape($iGrid, $iRow, $iCol)
                Return 0
            Case 0x56   ; 'V'
                If _WinAPI_GetAsyncKeyState(0x11) <> 0 Then  ; 0x11 = VK_CONTROL
                    __cg_commitCell($iGrid, $iRow, $iCol)
                    __cg_pasteAt($iGrid, $iRow, $iCol)
                    Return 0
                EndIf
        EndSwitch
    EndIf

    If $iMsg = $__cg_WM_CHAR Then
        ; VK_SPACE erzeugt WM_CHAR 0x20 — wir haben WM_KEYDOWN oben schon konsumiert,
        ; aber zur Sicherheit unterdrücken wir auch das Char
        If $wParam = $__cg_VK_SPACE Then Return 0
    EndIf

    If $iMsg = $__cg_WM_SETFOCUS Then
        ; Fokus-State aktualisieren + Text komplett selektieren
        Local $mGrid = $g_aCellGrids[$iGrid]
        If Not IsMap($mGrid) Then Return _WinAPI_DefSubclassProc($hCtrl, $iMsg, $wParam, $lParam)
        $mGrid.iFocusRow = $iRow
        $mGrid.iFocusCol = $iCol
        $g_aCellGrids[$iGrid] = $mGrid
        Local $lRes = _WinAPI_DefSubclassProc($hCtrl, $iMsg, $wParam, $lParam)
        DllCall("user32.dll", "lresult", "SendMessageW", "hwnd", $hCtrl, _
                "uint", $EM_SETSEL, "wparam", 0, "lparam", -1)
        Return $lRes
    EndIf

    If $iMsg = $__cg_WM_KILLFOCUS Then
        __cg_commitCell($iGrid, $iRow, $iCol)
        ; Fokus-State NICHT löschen — Navigation hat den neuen State schon gesetzt
    EndIf

    Return _WinAPI_DefSubclassProc($hCtrl, $iMsg, $wParam, $lParam)
EndFunc

; ---- Input-Spalten-Navigation (Helpers) --------------------------------------

; Gibt Index der nächsten Input-Spalte rechts von $iCol in $aSchema zurück,
; oder -1 wenn keine mehr.
Func __cg_nextInputCol(ByRef $aSchema, $iCol)
    Local $m
    For $i = $iCol + 1 To UBound($aSchema) - 1
        $m = $aSchema[$i]
        If $m.type = "input" Then Return $i
    Next
    Return -1
EndFunc

Func __cg_prevInputCol(ByRef $aSchema, $iCol)
    Local $m
    For $i = $iCol - 1 To 0 Step -1
        $m = $aSchema[$i]
        If $m.type = "input" Then Return $i
    Next
    Return -1
EndFunc

Func __cg_firstInputCol(ByRef $aSchema)
    Local $m
    For $i = 0 To UBound($aSchema) - 1
        $m = $aSchema[$i]
        If $m.type = "input" Then Return $i
    Next
    Return -1
EndFunc

Func __cg_lastInputCol(ByRef $aSchema)
    Local $m
    For $i = UBound($aSchema) - 1 To 0 Step -1
        $m = $aSchema[$i]
        If $m.type = "input" Then Return $i
    Next
    Return -1
EndFunc

; ---- Commit — aktuellen Input-Wert nach .aValues spiegeln + Callback feuern --
Func __cg_commitCell($iGrid, $iRow, $iCol)
    Local $mGrid = $g_aCellGrids[$iGrid]
    If Not IsMap($mGrid) Then Return
    Local $aControls = $mGrid.aControls
    Local $aValues   = $mGrid.aValues
    Local $iId = $aControls[$iRow][$iCol]
    If $iId <= 0 Then Return
    Local $sNew = GUICtrlRead($iId)
    Local $sOld = $aValues[$iRow][$iCol]
    If $sNew = $sOld Then Return   ; kein Commit wenn unverändert
    $aValues[$iRow][$iCol] = $sNew
    $mGrid.aValues = $aValues
    $g_aCellGrids[$iGrid] = $mGrid
    If $mGrid.sCommitCallback <> "" Then
        Call($mGrid.sCommitCallback, $iGrid, $iRow, $iCol, $sNew)
        If @error = 0xDEAD Then _
            ConsoleWrite("! CellGrid: commit callback '" & $mGrid.sCommitCallback & "' nicht gefunden" & @CRLF)
    EndIf
EndFunc

; ---- Navigation: Tab/Shift-Tab/Leertaste mit Wrap ----------------------------
; $iDir: +1 = Tab/Leertaste vorwärts, -1 = Shift+Tab rückwärts
; Returns True wenn konsumiert.
Func __cg_navHorizontalWrap($iGrid, $iRow, $iCol, $iDir)
    Local $mGrid = $g_aCellGrids[$iGrid]
    Local $aSchema = $mGrid.aSchema
    Local $iTargetCol = ($iDir > 0) ? __cg_nextInputCol($aSchema, $iCol) : __cg_prevInputCol($aSchema, $iCol)
    Local $iTargetRow = $iRow

    If $iTargetCol < 0 Then
        ; Zeilen-Wrap
        If $iDir > 0 Then
            $iTargetRow = $iRow + 1
            $iTargetCol = __cg_firstInputCol($aSchema)
            If $iTargetRow >= $mGrid.iRowCount Then
                ; Am Grid-Ende: durchlassen (Windows-Tab-Chain)
                Return False
            EndIf
        Else
            $iTargetRow = $iRow - 1
            $iTargetCol = __cg_lastInputCol($aSchema)
            If $iTargetRow < 0 Then Return False
        EndIf
    EndIf
    __cg_commitCell($iGrid, $iRow, $iCol)
    __cg_focusCell($iGrid, $iTargetRow, $iTargetCol)
    Return True
EndFunc

; ---- Navigation: Enter — gleiche Spalte, nächste Zeile, Auto-Append ----------
; Returns True wenn konsumiert.
Func __cg_navEnter($iGrid, $iRow, $iCol)
    Local $mGrid = $g_aCellGrids[$iGrid]
    Local $iTargetRow = $iRow + 1
    Local $bAppend = ($iTargetRow >= $mGrid.iRowCount)
    __cg_commitCell($iGrid, $iRow, $iCol)
    $mGrid = $g_aCellGrids[$iGrid]                ; refresh nach Callback
    If $bAppend Then __cg_setRowCount($iGrid, $mGrid.iRowCount + 1)
    __cg_focusCell($iGrid, $iTargetRow, $iCol)
    Return True
EndFunc

; ---- Navigation: Pfeiltasten (Links/Rechts) — kein Zeilen-Wrap --------------
; $iDir: +1 = rechts, -1 = links. Returns True wenn konsumiert, False wenn am Rand.
Func __cg_navArrowH($iGrid, $iRow, $iCol, $iDir)
    Local $mGrid = $g_aCellGrids[$iGrid]
    Local $aSchema = $mGrid.aSchema
    Local $iTarget = ($iDir > 0) _
        ? __cg_nextInputCol($aSchema, $iCol) _
        : __cg_prevInputCol($aSchema, $iCol)
    If $iTarget < 0 Then Return False
    __cg_commitCell($iGrid, $iRow, $iCol)
    __cg_focusCell($iGrid, $iRow, $iTarget)
    Return True
EndFunc

; ---- Navigation: Pfeiltasten (Hoch/Runter) — kein Wrap ----------------------
; $iDir: +1 = runter, -1 = hoch. Returns True wenn konsumiert, False wenn am Rand.
Func __cg_navArrowV($iGrid, $iRow, $iCol, $iDir)
    Local $mGrid = $g_aCellGrids[$iGrid]
    Local $iTargetRow = $iRow + $iDir
    If $iTargetRow < 0 Or $iTargetRow >= $mGrid.iRowCount Then Return False
    __cg_commitCell($iGrid, $iRow, $iCol)
    __cg_focusCell($iGrid, $iTargetRow, $iCol)
    Return True
EndFunc

; ---- Navigation: Escape — Fokus raus + Focus-State löschen -------------------
Func __cg_navEscape($iGrid, $iRow, $iCol)
    __cg_commitCell($iGrid, $iRow, $iCol)
    _WinAPI_SetFocus(0)   ; Fokus raus
    Local $mGrid = $g_aCellGrids[$iGrid]
    $mGrid.iFocusRow = -1
    $mGrid.iFocusCol = -1
    $g_aCellGrids[$iGrid] = $mGrid
EndFunc

; ---- Scroll-Handling ---------------------------------------------------------
; Finde Grid-Index anhand der Scroll-Body-GUI-HWND
Func __cg_gridFromScrollHwnd($hWnd)
    For $i = 0 To UBound($g_aCellGrids) - 1
        Local $mGrid = $g_aCellGrids[$i]
        If IsMap($mGrid) And $mGrid.hScrollGui = $hWnd Then Return $i
    Next
    Return -1
EndFunc

Func __cg_wmVscroll($hWnd, $iMsg, $wParam, $lParam)
    #forceref $iMsg, $lParam
    Local $iGrid = __cg_gridFromScrollHwnd($hWnd)
    If $iGrid < 0 Then Return $GUI_RUNDEFMSG
    Local $mGrid = $g_aCellGrids[$iGrid]
    Local $iLo = BitAND($wParam, 0xFFFF)   ; SB_* code

    Local $iOld = $mGrid.iScrollOffset
    Local $iNew = $iOld
    Local $iBodyH = $mGrid.iH - $mGrid.iHeaderH
    Local $iVirt  = $mGrid.iRowCount * $mGrid.iRowH
    Local $iMax   = $iVirt - $iBodyH
    If $iMax < 0 Then $iMax = 0

    Switch $iLo
        Case $SB_LINEUP
            $iNew -= $mGrid.iRowH
        Case $SB_LINEDOWN
            $iNew += $mGrid.iRowH
        Case $SB_PAGEUP
            $iNew -= $iBodyH
        Case $SB_PAGEDOWN
            $iNew += $iBodyH
        Case $SB_THUMBTRACK, $SB_THUMBPOSITION
            $iNew = BitShift($wParam, 16)   ; HIWORD = nTrackPos
        Case $SB_TOP
            $iNew = 0
        Case $SB_BOTTOM
            $iNew = $iMax
    EndSwitch
    If $iNew < 0     Then $iNew = 0
    If $iNew > $iMax Then $iNew = $iMax
    If $iNew = $iOld Then Return 0

    _GUIScrollBars_ScrollWindow($hWnd, 0, $iOld - $iNew)
    $mGrid.iScrollOffset = $iNew
    $g_aCellGrids[$iGrid] = $mGrid

    ; Scrollbar-Thumb synchronisieren
    Local $tSI = DllStructCreate($tagSCROLLINFO)
    DllStructSetData($tSI, "cbSize", DllStructGetSize($tSI))
    DllStructSetData($tSI, "fMask", $SIF_POS)
    DllStructSetData($tSI, "nPos", $iNew)
    _GUIScrollBars_SetScrollInfo($hWnd, $SB_VERT, $tSI)

    Return 0
EndFunc

Func __cg_wmMousewheel($hWnd, $iMsg, $wParam, $lParam)
    #forceref $iMsg, $lParam
    Local $iGrid = __cg_gridFromScrollHwnd($hWnd)
    If $iGrid < 0 Then Return $GUI_RUNDEFMSG
    ; High word of wParam = wheel delta (signed, multiple of 120)
    Local $iDelta = BitShift($wParam, 16)
    If $iDelta >= 0x8000 Then $iDelta = $iDelta - 0x10000
    Local $iLines = Int($iDelta / 120) * 3    ; 3 Zeilen pro Rastung
    If $iLines = 0 Then Return 0

    Local $mGrid = $g_aCellGrids[$iGrid]
    Local $iBodyH = $mGrid.iH - $mGrid.iHeaderH
    Local $iVirt  = $mGrid.iRowCount * $mGrid.iRowH
    Local $iMax   = $iVirt - $iBodyH
    If $iMax < 0 Then $iMax = 0
    Local $iOld = $mGrid.iScrollOffset
    Local $iNew = $iOld - $iLines * $mGrid.iRowH
    If $iNew < 0     Then $iNew = 0
    If $iNew > $iMax Then $iNew = $iMax
    If $iNew = $iOld Then Return 0

    _GUIScrollBars_ScrollWindow($hWnd, 0, $iOld - $iNew)
    $mGrid.iScrollOffset = $iNew
    $g_aCellGrids[$iGrid] = $mGrid

    Local $tSI = DllStructCreate($tagSCROLLINFO)
    DllStructSetData($tSI, "cbSize", DllStructGetSize($tSI))
    DllStructSetData($tSI, "fMask", $SIF_POS)
    DllStructSetData($tSI, "nPos", $iNew)
    _GUIScrollBars_SetScrollInfo($hWnd, $SB_VERT, $tSI)
    Return 0
EndFunc

; ---- Multi-Cell-Paste (Ctrl+V) ----------------------------------------------
; Liest Clipboard-Text, splittet in Zeilen (CRLF/LF) und Felder (TAB),
; schreibt ab ($iStartRow, $iStartCol) — Label-/Separator-Spalten werden übersprungen.
; Fehlende Zeilen werden automatisch per __cg_setRowCount angelegt.
Func __cg_pasteAt($iGrid, $iStartRow, $iStartCol)
    Local $sText = _ClipBoard_GetData()
    If @error Or $sText = "" Then Return

    ; Zeilen splitten
    $sText = StringReplace($sText, @CRLF, @LF)
    $sText = StringReplace($sText, @CR, @LF)
    Local $aLines = StringSplit($sText, @LF, 2)   ; flag 2 = no count element
    ; Trailing leere Zeile entfernen
    While UBound($aLines) > 0 And $aLines[UBound($aLines) - 1] = ""
        ReDim $aLines[UBound($aLines) - 1]
    WEnd
    If UBound($aLines) = 0 Then Return

    Local $mGrid = $g_aCellGrids[$iGrid]
    Local $aSchema = $mGrid.aSchema

    ; Zielzeilenzahl sicherstellen
    Local $iNeeded = $iStartRow + UBound($aLines)
    If $iNeeded > $mGrid.iRowCount Then
        __cg_setRowCount($iGrid, $iNeeded)
        $mGrid = $g_aCellGrids[$iGrid]
    EndIf

    Local $iMaxEndRow = $iStartRow + UBound($aLines) - 1
    Local $bOverflow = False

    For $i = 0 To UBound($aLines) - 1
        Local $aFields = StringSplit($aLines[$i], @TAB, 2)
        Local $iRow = $iStartRow + $i
        Local $iCol = $iStartCol
        For $j = 0 To UBound($aFields) - 1
            ; Skip Label- und Separator-Spalten beim Vorwärtslaufen
            While $iCol < UBound($aSchema)
                Local $mColDef = $aSchema[$iCol]
                If $mColDef.type = "input" Then ExitLoop
                $iCol += 1
            WEnd
            If $iCol >= UBound($aSchema) Then
                $bOverflow = True
                ExitLoop
            EndIf
            __cg_setCell($iGrid, $iRow, $iCol, $aFields[$j])
            $iCol += 1
        Next
    Next

    ; Bulk-Commit-Callback einmalig feuern
    $mGrid = $g_aCellGrids[$iGrid]
    If $mGrid.sBulkCommitCallback <> "" Then
        Call($mGrid.sBulkCommitCallback, $iGrid, $iStartRow, $iMaxEndRow)
        If @error = 0xDEAD Then _
            ConsoleWrite("! CellGrid: bulk-commit callback '" & $mGrid.sBulkCommitCallback & "' nicht gefunden" & @CRLF)
    EndIf

    If $bOverflow Then
        ; Hinweis, dass zu viele Spalten geklappt wurden — Caller entscheidet Darstellung
        ConsoleWrite("! Paste: manche Felder überstiegen Input-Spalten-Bereich" & @CRLF)
    EndIf
EndFunc
; =============== end of inlined _CellGrid.au3 =================================

; --- Layout ---
Global Const $g_iWinW       = 1150
Global Const $g_iWinH       = 750
Global Const $g_iHeaderH    = 70
Global Const $g_iUpperH     = 290
Global Const $g_fDataFrac   = 0.68
Global Const $g_iGraphPadT  = 45
Global Const $g_iGraphPadB  = 35
Global Const $g_iGraphPadL  = 55
Global Const $g_iGraphPadR  = 20

; --- Colors ---
Global Const $g_iClrBg      = 0xFF0A1520
Global Const $g_iClrGrid    = 0xFF1E3040
Global Const $g_iClrAxis    = 0xFF556677
Global Const $g_iClrCurve   = 0xFF4A9EFF
Global Const $g_iClrPoint   = 0xFFFF6B6B
Global Const $g_iClrResid   = 0xFFDD6644
Global Const $g_iClrConf    = 0x304A9EFF
Global Const $g_iClrText    = 0xFF99AABB

; --- Regression types ---
Global Const $REG_NORMAL    = 0
Global Const $REG_ORTHO     = 1
Global Const $REG_DEMING    = 2
Global Const $REG_YORK      = 3

; Max number of σ field pairs pre-allocated in the header (Deming: Obs + DepVar).
; More than this won't fit horizontally before the Robust block at x=750 anyway.
Global Const $g_iMaxSdFields = 6

; --- Known math function names (excluded from parameter detection) ---
Global Const $g_sFuncNames  = "|SIN|COS|TAN|ASIN|ACOS|ATAN|SQRT|EXP|LN|LOG10|ABS|PI|"

; --- State: Formula ---
Global $g_sFormula          = ""
Global $g_sDepVar           = "Y"
Global $g_aObsNames[0]
Global $g_aParamNames[0]

; --- State: Data ---
Global $g_mData[]
Global $g_iRowCount         = 0
Global $g_mBackupData[]

; --- State: Regression ---
Global $g_iRegType          = $REG_NORMAL
Global $g_mSdValues[]       ; Map: varName → σ value (persistent across rebuilds)
Global $g_bRobust           = False
Global $g_sEstimator        = "Huber"

; --- State: Results ---
Global $g_mResults          = Null
Global $g_bAdjOK            = False
Global $g_sErrorMsg         = ""

; --- State: Graph ---
Global $g_fXMin = 0, $g_fXMax = 10
Global $g_fYMin = 0, $g_fYMax = 10
Global $g_bAutoAxis         = True
Global $g_bShowResiduals    = False
Global $g_bShowConfBand     = False

; --- GUI handles ---
Global $g_hGUI
Global $g_idInputFormula, $g_idInputDepVar
Global $g_idRadNormal, $g_idRadOrtho, $g_idRadDeming, $g_idRadYork
Global $g_aidSdLabels[$g_iMaxSdFields]    ; pre-allocated label control IDs
Global $g_aidSdInputs[$g_iMaxSdFields]    ; pre-allocated input control IDs
Global $g_asSdNames[$g_iMaxSdFields]      ; variable name per slot ("" = shared σ for Orthogonal)
Global $g_iSdFieldCount = 0               ; number of currently visible σ slots (0..$g_iMaxSdFields)
Global $g_idSdHintLabel                   ; York-mode hint label
Global $g_idChkRobust, $g_idComboEstimator
Global $g_idBtnCsvLoad, $g_idBtnCsvExport, $g_idBtnRowDel, $g_idBtnClear, $g_idBtnReport
Global $g_idBtnRescaleSigma ; "σ · s₀" — rescale a-priori σ with empirical s₀ (global variance factor)
Global $g_sLastImportedCSV = "" ; remembered for default filename in report export
Global $g_mLastStartValues[] ; paramName -> start value string, captured before each solve
                             ; (grid holds adjusted values after _UpdateResults, so report
                             ; would otherwise show results as "start values")
Global $g_idLblStats, $g_idLblError
Global $g_idChkResid, $g_idChkConf
Global $g_idInputXMin, $g_idInputXMax, $g_idInputYMin, $g_idInputYMax
Global $g_idBtnAutoAxis
Global $g_idInputPredX, $g_idLblPredResult  ; prediction controls

; --- GDI+ handles ---
Global $g_hGraphics, $g_hBitmap, $g_hBuffer
Global $g_iGraphX, $g_iGraphY, $g_iGraphW, $g_iGraphH
Global $g_hPenCurve, $g_hPenAxis, $g_hPenGrid, $g_hPenResid
Global $g_hBrushPoint, $g_hBrushConf, $g_hBrushText
Global $g_aRobWeights[0]       ; per-row robust weights [0..1]
Global $g_hBrushOutlier        ; gray brush for outlier points
Global $g_hBrushResLabel       ; brush for residual value labels
Global $g_hFontAxis, $g_hFontFamily, $g_hFormatCenter, $g_hFormatRight

; --- CellGrid handles ---
Global $g_iDataGrid  = -1
Global $g_iParamGrid = -1

; --- Deferred actions (NEVER modify GUI controls inside WM_ message handlers) ---
Global $g_bPendingAddRow   = False  ; "Enter an letzter Zeile"-Signal, falls benötigt
Global $g_bPendingParse     = False  ; deferred _ParseFormula from WM_COMMAND
Global $g_bPendingSolve     = False  ; deferred _TryAutoSolve from WM_COMMAND
Global $g_bPendingPrediction = False ; deferred _UpdatePrediction from WM_COMMAND

; ===== GUI Setup =============================================================
_SetupGUI()

; ===== Message handlers ======================================================
GUIRegisterMsg($WM_NOTIFY, "_WM_NOTIFY")
GUIRegisterMsg($WM_COMMAND, "_WM_COMMAND")
GUIRegisterMsg($WM_PAINT, "_WM_PAINT")
GUIRegisterMsg($WM_ERASEBKGND, "_WM_ERASEBKGND")

GUISetState(@SW_SHOW, $g_hGUI)
_LoadStartupExample()

; ===== Main loop =============================================================
While 1
	; Check if formula/depvar input is focused and Enter pressed
	If _WinAPI_GetAsyncKeyState(0x0D) Then ; VK_RETURN
		If _WinAPI_GetFocus() = GUICtrlGetHandle($g_idInputFormula) Or _
		   _WinAPI_GetFocus() = GUICtrlGetHandle($g_idInputDepVar) Then
			_ParseFormula()
		EndIf
	EndIf

	; Process deferred actions from WM_ handlers (NEVER modify GUI inside message callbacks)
	If $g_bPendingParse Then
		$g_bPendingParse = False
		_ParseFormula()
	EndIf
	If $g_bPendingSolve Then
		$g_bPendingSolve = False
		_TryAutoSolve()
	EndIf
	If $g_bPendingPrediction Then
		$g_bPendingPrediction = False
		_UpdatePrediction()
	EndIf

	Switch GUIGetMsg()
		Case $GUI_EVENT_CLOSE
			ExitLoop
		Case $g_idRadNormal, $g_idRadOrtho, $g_idRadDeming, $g_idRadYork
			_OnRegTypeChanged()
		Case $g_idChkRobust
			$g_bRobust = (BitAND(GUICtrlRead($g_idChkRobust), $GUI_CHECKED) = $GUI_CHECKED)
			GUICtrlSetState($g_idComboEstimator, $g_bRobust ? $GUI_ENABLE : $GUI_DISABLE)
			_TryAutoSolve()
		Case $g_idComboEstimator
			$g_sEstimator = GUICtrlRead($g_idComboEstimator)
			_TryAutoSolve()
		Case $g_idBtnRowDel
			_DeleteSelectedRow()
		Case $g_idBtnClear
			_ClearAllData()
		Case $g_idBtnCsvLoad
			_ImportCSV()
		Case $g_idBtnCsvExport
			_ExportCSV()
		Case $g_idBtnReport
			_ExportResults()
		Case $g_idBtnRescaleSigma
			_RescaleSigmaFromS0()
		Case $g_idBtnAutoAxis
			_AutoFitAxes()
		Case $g_idChkResid
			$g_bShowResiduals = (BitAND(GUICtrlRead($g_idChkResid), $GUI_CHECKED) = $GUI_CHECKED)
			_DrawGraph()
		Case $g_idChkConf
			$g_bShowConfBand = (BitAND(GUICtrlRead($g_idChkConf), $GUI_CHECKED) = $GUI_CHECKED)
			_DrawGraph()
	EndSwitch
WEnd

_Cleanup()

; ===== GUI Setup function ====================================================

Func _SetupGUI()
	$g_hGUI = GUICreate("Interaktive Regression", $g_iWinW, $g_iWinH, -1, -1, _
		BitOR($WS_MINIMIZEBOX, $WS_CAPTION, $WS_SYSMENU))

	; --- Header bar ---
	Local $iX = 10
	Local $iY = 8

	GUICtrlCreateLabel("Formel:", $iX, $iY + 3, 42, 20)
	$g_idInputFormula = GUICtrlCreateInput("", $iX + 44, $iY, 280, 22)
	GUICtrlSetFont(-1, 10, 400, 0, "Consolas")

	$iX = 345
	GUICtrlCreateLabel("Abh.Var.:", $iX, $iY + 3, 55, 20)
	$g_idInputDepVar = GUICtrlCreateInput("Y", $iX + 57, $iY, 35, 22)
	GUICtrlSetFont(-1, 10, 400, 0, "Consolas")

	$iX = 455
	$g_idRadNormal = GUICtrlCreateRadio("Normal", $iX, $iY, 65, 20)
	GUICtrlSetState(-1, $GUI_CHECKED)
	$g_idRadOrtho  = GUICtrlCreateRadio("Orthogonal", $iX + 68, $iY, 85, 20)
	$g_idRadDeming = GUICtrlCreateRadio("Deming", $iX + 158, $iY, 65, 20)
	$g_idRadYork   = GUICtrlCreateRadio("York", $iX + 228, $iY, 55, 20)

	; Globaler Varianzfaktor: skaliert alle σ-Inputs mit s₀ (erhält das
	; σ-Verhältnis). Parameter und sdx ändern sich dabei nicht; der Chi²-
	; Globaltest passt anschließend per Konstruktion.
	$g_idBtnRescaleSigma = GUICtrlCreateButton(ChrW(963) & " " & ChrW(0x00B7) & " s" & ChrW(0x2080), _
		$iX + 300, $iY - 1, 85, 22)
	GUICtrlSetTip($g_idBtnRescaleSigma, _
		"σ · s₀ — Alle σ-Inputs mit empirischem s₀ skalieren." & @CRLF & _
		"Das Verhältnis der σ-Werte zueinander bleibt erhalten." & @CRLF & _
		"Parameter und sd(x̂) ändern sich nicht." & @CRLF & _
		"Nach Re-Solve liegt s₀ ≈ 1, Chi²-Globaltest passt." & @CRLF & _
		"(Voraussetzung: gültige Ausgleichung liegt vor.)")

	; Pre-allocate σ field controls at on-screen default position, then hide
	; them via direct WinAPI ShowWindow(SW_HIDE). This is Variant B — kept
	; for bisection / debugging; Variant C (child GUI) is the known-good fix.
	$iY = 38
	For $i = 0 To $g_iMaxSdFields - 1
		$g_aidSdLabels[$i] = GUICtrlCreateLabel("", 10, $iY + 3, 40, 20)
		$g_aidSdInputs[$i] = GUICtrlCreateInput("1.0", 10, $iY, 50, 22)
		$g_asSdNames[$i] = ""
		__HideCtrl($g_aidSdLabels[$i])
		__HideCtrl($g_aidSdInputs[$i])
	Next
	$g_idSdHintLabel = GUICtrlCreateLabel(ChrW(963) & " individuell pro Zeile", 10, $iY + 3, 220, 20)
	GUICtrlSetColor($g_idSdHintLabel, 0x888888)
	GUICtrlSetFont($g_idSdHintLabel, 9, 400, 2)
	__HideCtrl($g_idSdHintLabel)

	; Robust controls: second row, right side (σ fields occupy x=10..~950)
	$iX = 950
	$iY = 38
	$g_idChkRobust = GUICtrlCreateCheckbox("Robust:", $iX, $iY, 62, 20)
	$g_idComboEstimator = GUICtrlCreateCombo("", $iX + 65, $iY, 90, 22, $CBS_DROPDOWNLIST)
	GUICtrlSetData(-1, "Huber|Hampel|Biweight|L1|BIBER", "Huber")
	GUICtrlSetState($g_idComboEstimator, $GUI_DISABLE)

	GUICtrlCreateGraphic(0, $g_iHeaderH - 2, $g_iWinW, 1)
	GUICtrlSetBkColor(-1, 0xA0A0A0)

	; --- Upper section: Data table (left) ---
	Local $iDataW = Int($g_iWinW * $g_fDataFrac)
	Local $iUpperY = $g_iHeaderH

	; Buttons above data table
	$g_idBtnCsvLoad   = GUICtrlCreateButton("Laden", 10, $iUpperY + 2, 50, 24)
	$g_idBtnCsvExport = GUICtrlCreateButton("Speichern", 63, $iUpperY + 2, 65, 24)
	$g_idBtnRowDel    = GUICtrlCreateButton("Zeile −", 133, $iUpperY + 2, 55, 24)
	$g_idBtnClear     = GUICtrlCreateButton("Löschen", 193, $iUpperY + 2, 60, 24)
	$g_idBtnReport    = GUICtrlCreateButton("Report…", 258, $iUpperY + 2, 70, 24)

	; Datengrid-Platzhalter — tatsächliches Schema wird von _RebuildDataColumns gesetzt
	Local $aEmptySchema[0]
	$g_iDataGrid = __cg_create($g_hGUI, 10, $iUpperY + 30, $iDataW - 15, $g_iUpperH - 35, _
		$aEmptySchema, "_OnDataCellCommit", "_OnDataBulkCommit")

	; --- Upper section: Parameter panel (right) ---
	Local $iParamX = Int($g_iWinW * $g_fDataFrac) + 5
	Local $iParamW = $g_iWinW - $iParamX - 10

	GUICtrlCreateLabel("Parameter", $iParamX, $iUpperY + 4, 70, 18)
	GUICtrlSetFont(-1, 9, 700)

	Local $aEmptyPSchema[0]
	$g_iParamGrid = __cg_create($g_hGUI, $iParamX, $iUpperY + 24, $iParamW, 120, _
		$aEmptyPSchema, "_OnParamCellCommit", "")

	; Statistics label
	Local $iStatY = $iUpperY + 150
	$g_idLblStats = GUICtrlCreateLabel("", $iParamX, $iStatY, $iParamW, 115)
	GUICtrlSetFont(-1, 9, 400, 0, "Consolas")

	; Error label (red)
	$g_idLblError = GUICtrlCreateLabel("", $iParamX, $iStatY + 118, $iParamW, 30)
	GUICtrlSetFont(-1, 9, 400)
	GUICtrlSetColor(-1, 0xCC2222)

	; --- Lower section: Graph ---
	Local $iGraphTopY = $g_iHeaderH + $g_iUpperH + 5

	$g_idChkResid = GUICtrlCreateCheckbox("Residuen", 10, $iGraphTopY, 75, 20)
	$g_idChkConf  = GUICtrlCreateCheckbox("Konfidenzband", 90, $iGraphTopY, 105, 20)
	GUICtrlSetTip($g_idChkConf, _
		"95 %-Konfidenzband für den Erwartungswert ŷ(x) (zweiseitig, α = 0.05)." & @CRLF & _
		"Bandhalbbreite: ±t · s₀ · √(J·Qxx·Jᵀ), mit J = ∂f/∂θ bei θ = θ̂." & @CRLF & _
		"t-Faktor:  f > 30 → 1.96  |  f ≤ 30 → 2 + 4/f (Approximation, leicht konservativ)." & @CRLF & _
		"Hinweis: Das ist die Unsicherheit des angepassten Modells, NICHT ein" & @CRLF & _
		"Prognoseband für Einzelmessungen (das wäre zusätzlich um σ der Messung breiter).")

	GUICtrlCreateLabel("X:", 220, $iGraphTopY + 2, 15, 18)
	$g_idInputXMin = GUICtrlCreateInput("0", 236, $iGraphTopY, 55, 20)
	GUICtrlCreateLabel("-", 295, $iGraphTopY + 2, 10, 18)
	$g_idInputXMax = GUICtrlCreateInput("10", 308, $iGraphTopY, 55, 20)
	GUICtrlCreateLabel("Y:", 380, $iGraphTopY + 2, 15, 18)
	$g_idInputYMin = GUICtrlCreateInput("0", 396, $iGraphTopY, 55, 20)
	GUICtrlCreateLabel("-", 455, $iGraphTopY + 2, 10, 18)
	$g_idInputYMax = GUICtrlCreateInput("10", 468, $iGraphTopY, 55, 20)
	$g_idBtnAutoAxis = GUICtrlCreateButton("Auto", 535, $iGraphTopY, 42, 20)

	; Prediction row
	Local $iPredY = $iGraphTopY + 24
	GUICtrlCreateLabel("Prädiktion:", 10, $iPredY + 2, 65, 18)
	GUICtrlSetFont(-1, 9, 700)
	$g_idInputPredX = GUICtrlCreateInput("", 78, $iPredY, 80, 20)
	GUICtrlSetFont(-1, 9, 400, 0, "Consolas")
	GUICtrlCreateLabel(ChrW(8594), 162, $iPredY + 1, 15, 18) ; → arrow
	GUICtrlSetFont(-1, 11, 400)
	$g_idLblPredResult = GUICtrlCreateLabel("", 180, $iPredY + 2, 400, 18)
	GUICtrlSetFont(-1, 9, 400, 0, "Consolas")
	GUICtrlSetColor(-1, 0x2266BB)

	; Graph area coordinates
	$g_iGraphX = 0
	$g_iGraphY = $iPredY + 25
	$g_iGraphW = $g_iWinW
	$g_iGraphH = $g_iWinH - $g_iGraphY

	; GDI+ init
	_GDIPlus_Startup()
	$g_hGraphics = _GDIPlus_GraphicsCreateFromHWND($g_hGUI)
	$g_hBitmap   = _GDIPlus_BitmapCreateFromGraphics($g_iGraphW, $g_iGraphH, $g_hGraphics)
	$g_hBuffer   = _GDIPlus_ImageGetGraphicsContext($g_hBitmap)
	_GDIPlus_GraphicsSetSmoothingMode($g_hBuffer, 4) ; AntiAlias

	; Pens & Brushes
	$g_hPenCurve  = _GDIPlus_PenCreate($g_iClrCurve, 2)
	$g_hPenAxis   = _GDIPlus_PenCreate($g_iClrAxis, 1)
	$g_hPenGrid   = _GDIPlus_PenCreate($g_iClrGrid, 1)
	$g_hPenResid  = _GDIPlus_PenCreate($g_iClrResid, 1)
	_GDIPlus_PenSetDashStyle($g_hPenResid, 2) ; Dot
	$g_hBrushPoint = _GDIPlus_BrushCreateSolid($g_iClrPoint)
	$g_hBrushConf  = _GDIPlus_BrushCreateSolid($g_iClrConf)
	$g_hBrushText  = _GDIPlus_BrushCreateSolid($g_iClrText)
	$g_hBrushOutlier = _GDIPlus_BrushCreateSolid(0xFFBBBBBB) ; light gray for outliers
	$g_hBrushResLabel = _GDIPlus_BrushCreateSolid(0xFFDD6644) ; orange, matching residual line color

	; Font for axis labels
	$g_hFontFamily  = _GDIPlus_FontFamilyCreate("Segoe UI")
	$g_hFontAxis    = _GDIPlus_FontCreate($g_hFontFamily, 8)
	$g_hFormatCenter = _GDIPlus_StringFormatCreate()
	_GDIPlus_StringFormatSetAlign($g_hFormatCenter, 1)
	$g_hFormatRight  = _GDIPlus_StringFormatCreate()
	_GDIPlus_StringFormatSetAlign($g_hFormatRight, 2)
EndFunc   ;==>_SetupGUI

; ===== Regression type switching ==============================================

Func _OnRegTypeChanged()
	; Determine selected type
	If BitAND(GUICtrlRead($g_idRadNormal), $GUI_CHECKED) Then
		$g_iRegType = $REG_NORMAL
	ElseIf BitAND(GUICtrlRead($g_idRadOrtho), $GUI_CHECKED) Then
		$g_iRegType = $REG_ORTHO
	ElseIf BitAND(GUICtrlRead($g_idRadDeming), $GUI_CHECKED) Then
		$g_iRegType = $REG_DEMING
	ElseIf BitAND(GUICtrlRead($g_idRadYork), $GUI_CHECKED) Then
		$g_iRegType = $REG_YORK
	EndIf

	_RebuildDataColumns()
	_RebuildSdFields()
	_TryAutoSolve()
EndFunc   ;==>_OnRegTypeChanged

Func _ReadSdValues()
	For $i = 0 To $g_iSdFieldCount - 1
		Local $sKey = ($g_asSdNames[$i] = "") ? "__shared__" : $g_asSdNames[$i]
		Local $fVal = Number(GUICtrlRead($g_aidSdInputs[$i]))
		If $fVal <= 0 Then $fVal = 1.0
		$g_mSdValues[$sKey] = $fVal
	Next
EndFunc   ;==>_ReadSdValues

; Global variance factor: σ_neu = σ_alt · s₀ für alle sichtbaren σ-Inputs.
; Das Verhältnis der σ-Werte bleibt erhalten → Parameter und sdx invariant;
; nur der Chi²-Globaltest rutscht nach s₀ ≈ 1. Bei York liegen die σ_i in
; den sd_-Spalten des Grids, nicht in den Header-Inputs — dort gibt's statt
; Autoskalierung nur einen Hinweis.
Func _RescaleSigmaFromS0()
	If Not $g_bAdjOK Or Not IsMap($g_mResults) Then
		GUICtrlSetData($g_idLblError, "Zuerst eine gültige Ausgleichung durchführen")
		Return
	EndIf
	If Not MapExists($g_mResults, "s0") Or $g_mResults.s0 <= 0 Then
		GUICtrlSetData($g_idLblError, "s₀ nicht verfügbar")
		Return
	EndIf

	Local $fS0 = $g_mResults.s0

	If $g_iRegType = $REG_YORK Then
		GUICtrlSetData($g_idLblError, StringFormat( _
			"Bei York liegen die σ pro Messpunkt in den sd_-Spalten — " & _
			"bitte dort manuell mit s₀=%.6g multiplizieren", $fS0))
		Return
	EndIf

	If $g_iSdFieldCount = 0 Then
		GUICtrlSetData($g_idLblError, "Keine σ-Eingabefelder zum Skalieren vorhanden")
		Return
	EndIf

	For $i = 0 To $g_iSdFieldCount - 1
		Local $fCur = Number(GUICtrlRead($g_aidSdInputs[$i]))
		If $fCur <= 0 Then $fCur = 1.0
		Local $fNew = $fCur * $fS0
		GUICtrlSetData($g_aidSdInputs[$i], StringFormat("%.6g", $fNew))
		Local $sKey = ($g_asSdNames[$i] = "") ? "__shared__" : $g_asSdNames[$i]
		$g_mSdValues[$sKey] = $fNew
	Next

	GUICtrlSetData($g_idLblError, "")
	_TryAutoSolve()
EndFunc   ;==>_RescaleSigmaFromS0

; ===== Model building & solving ==============================================

Func _BuildAndSolve()
	$g_bAdjOK = False
	$g_mResults = Null
	$g_sErrorMsg = ""
	ReDim $g_aRobWeights[0]

	If $g_sFormula = "" Or UBound($g_aObsNames) = 0 Then
		$g_sErrorMsg = "Keine Formel eingegeben"
		GUICtrlSetData($g_idLblError, $g_sErrorMsg)
		Return
	EndIf
	If $g_iRowCount < 1 Then
		$g_sErrorMsg = "Keine Messdaten vorhanden"
		GUICtrlSetData($g_idLblError, $g_sErrorMsg)
		Return
	EndIf

	_ReadSdValues()

	Local $mSystem = _adj_createSystem()

	; Read initial values from parameter table; snapshot them so _ExportResults
	; can still show the true start values after _UpdateResults overwrites the
	; grid with adjusted results.
	Local $mStartSnap[]
	For $i = 0 To UBound($g_aParamNames) - 1
		Local $sStart = __cg_getCell($g_iParamGrid, $i, 1)
		Local $fStart = Number($sStart)
		If $sStart = "" Then
			$fStart = 0.5
			$sStart = "0.5"
		EndIf
		$mStartSnap[$g_aParamNames[$i]] = $sStart
		_adj_setInitialValue($mSystem, $g_aParamNames[$i], $fStart)
	Next
	$g_mLastStartValues = $mStartSnap

	; Build model based on regression type
	Local $iValidRows = 0
	Switch $g_iRegType
		Case $REG_NORMAL
			$iValidRows = _BuildNormalModel($mSystem)
		Case $REG_ORTHO, $REG_DEMING
			$iValidRows = _BuildGLMModel($mSystem)
		Case $REG_YORK
			$iValidRows = _BuildYorkModel($mSystem)
	EndSwitch

	If $iValidRows < UBound($g_aParamNames) Then
		$g_sErrorMsg = "Zu wenige Daten (" & $iValidRows & ") fuer " & UBound($g_aParamNames) & " Parameter"
		GUICtrlSetData($g_idLblError, $g_sErrorMsg)
		_UpdateResults()
		Return
	EndIf

	; Configure solver
	$g_sEstimator = GUICtrlRead($g_idComboEstimator)
	Local $mConfig = _adj_defaultConfig("LM", False)
	$mConfig.maxIterations = 50
	Local $mCompute = $mConfig.compute
	$mCompute.globalTest = True
	$mCompute.cofactors = True
	$mCompute.redundancy = True
	$mConfig.compute = $mCompute
	If $g_bRobust Then
		$mConfig.robust = $g_sEstimator
		$mConfig.robustParams = _adj_robustDefaults($g_sEstimator)
	EndIf

	; Solve
	_adj_solve($mSystem, $mConfig)
	Local $iErr = @error
	Local $iExt = @extended

	If $iErr Then
		$g_sErrorMsg = _adj_getErrorMessage($iErr, $iExt)
		If $iErr = $ADJ_ERR_NO_CONVERGENCE Then
			$g_sErrorMsg &= @CRLF & "-> Startwerte anpassen!"
		EndIf
		GUICtrlSetData($g_idLblError, $g_sErrorMsg)
		_UpdateResults()
		Return
	EndIf

	$g_mResults = _adj_getResults($mSystem)
	$g_bAdjOK = True
	$g_sErrorMsg = ""
	GUICtrlSetData($g_idLblError, "")

	; Extract per-observation robust weights
	ReDim $g_aRobWeights[$g_iRowCount]
	If $g_bRobust And IsMap($g_mResults.robustWeights) Then
		Local $mRobW = $g_mResults.robustWeights
		For $i = 0 To $g_iRowCount - 1
			Local $sObsSuffix = String($i + 1)
			Local $fMinW = 1.0
			; Check dependent variable weight
			Local $sDepKey = $g_sDepVar & $sObsSuffix
			If MapExists($mRobW, $sDepKey) Then
				$fMinW = $mRobW[$sDepKey]
			EndIf
			; For GLM models, also check observation weights (take minimum)
			If $g_iRegType <> $REG_NORMAL Then
				For $j = 0 To UBound($g_aObsNames) - 1
					Local $sObsKey = $g_aObsNames[$j] & $sObsSuffix
					If MapExists($mRobW, $sObsKey) Then
						Local $fW = $mRobW[$sObsKey]
						If $fW < $fMinW Then $fMinW = $fW
					EndIf
				Next
			EndIf
			$g_aRobWeights[$i] = $fMinW
		Next
	Else
		For $i = 0 To $g_iRowCount - 1
			$g_aRobWeights[$i] = 1.0
		Next
	EndIf

	_UpdateResults()

	If $g_bAutoAxis Then
		_AutoFitAxes() ; this calls _DrawGraph() internally
	Else
		_DrawGraph()
	EndIf
EndFunc   ;==>_BuildAndSolve

Func _BuildNormalModel(ByRef $mSystem)
	Local $iValid = 0

	For $i = 0 To $g_iRowCount - 1
		Local $sObsSuffix = String($i + 1)
		Local $bSkip = False

		; Check all observation columns have values
		For $j = 0 To UBound($g_aObsNames) - 1
			Local $aColData = $g_mData[$g_aObsNames[$j]]
			If $aColData[$i] = "" Then
				$bSkip = True
				ExitLoop
			EndIf
		Next
		; Check dependent variable
		Local $aDepData = $g_mData[$g_sDepVar]
		If $aDepData[$i] = "" Then $bSkip = True
		If $bSkip Then ContinueLoop

		; Substitute actual X values directly into formula (true OLS — X is exact)
		Local $sSubstFormula = $g_sFormula
		For $j = 0 To UBound($g_aObsNames) - 1
			Local $sObsName = $g_aObsNames[$j]
			Local $aOD = $g_mData[$sObsName]
			$sSubstFormula = StringRegExpReplace($sSubstFormula, _
				"(?i)#" & $sObsName & "\b", "(" & String(Number($aOD[$i])) & ")")
		Next

		; Add observation function: Y_i = f(x_i, params) — no _adj_addObs for X
		Local $fSdDep = 1.0
		If MapExists($g_mSdValues, $g_sDepVar) Then $fSdDep = $g_mSdValues[$g_sDepVar]
		_adj_addObsFunction($mSystem, $g_sDepVar & $sObsSuffix, _
			$sSubstFormula, Number($aDepData[$i]), $fSdDep)

		$iValid += 1
	Next

	Return $iValid
EndFunc   ;==>_BuildNormalModel

Func _BuildGLMModel(ByRef $mSystem)
	Local $iValid = 0

	; σ-Quelle abhängig vom Regression-Typ: Deming nutzt per-Variable σ
	; (ein Feld pro Obs und DepVar), Orthogonal nutzt *nur* den shared-Key.
	; Legacy per-variable Einträge in $g_mSdValues (z.B. von vorher Deming)
	; dürfen den Ortho-Build nicht beeinflussen, sonst greift der
	; σ·s₀-Rescale-Button nicht (shared-σ würde überschrieben).
	Local $bUsePerVariable = ($g_iRegType = $REG_DEMING)

	Local $fSharedSd = 1.0
	If MapExists($g_mSdValues, "__shared__") Then $fSharedSd = $g_mSdValues["__shared__"]

	For $i = 0 To $g_iRowCount - 1
		Local $sObsSuffix = String($i + 1)
		Local $bSkip = False

		; Check all observation columns
		For $j = 0 To UBound($g_aObsNames) - 1
			Local $aColData = $g_mData[$g_aObsNames[$j]]
			If $aColData[$i] = "" Then
				$bSkip = True
				ExitLoop
			EndIf
		Next
		Local $aDepData = $g_mData[$g_sDepVar]
		If $aDepData[$i] = "" Then $bSkip = True
		If $bSkip Then ContinueLoop

		; Add observations with per-variable (Deming) or shared (Orthogonal) σ
		For $j = 0 To UBound($g_aObsNames) - 1
			Local $sObsName = $g_aObsNames[$j]
			Local $aOD = $g_mData[$sObsName]
			Local $fSd = $fSharedSd
			If $bUsePerVariable And MapExists($g_mSdValues, $sObsName) Then $fSd = $g_mSdValues[$sObsName]
			_adj_addObs($mSystem, $sObsName & $sObsSuffix, Number($aOD[$i]), $fSd, "s0")
		Next

		; Dependent variable
		Local $fSdDep = $fSharedSd
		If $bUsePerVariable And MapExists($g_mSdValues, $g_sDepVar) Then $fSdDep = $g_mSdValues[$g_sDepVar]
		_adj_addObs($mSystem, $g_sDepVar & $sObsSuffix, Number($aDepData[$i]), $fSdDep, "s0")

		; Build condition: #Y_i - (formula_with_#X_i) = 0
		Local $sIndexedFormula = $g_sFormula
		For $j = 0 To UBound($g_aObsNames) - 1
			$sIndexedFormula = StringRegExpReplace($sIndexedFormula, _
				"(?i)#" & $g_aObsNames[$j] & "\b", "#" & $g_aObsNames[$j] & $sObsSuffix)
		Next
		_adj_addFunction($mSystem, "#" & $g_sDepVar & $sObsSuffix & " - (" & $sIndexedFormula & ")")

		$iValid += 1
	Next

	Return $iValid
EndFunc   ;==>_BuildGLMModel

Func _BuildYorkModel(ByRef $mSystem)
	Local $iValid = 0

	For $i = 0 To $g_iRowCount - 1
		Local $sObsSuffix = String($i + 1)
		Local $bSkip = False

		; Check all observation columns
		For $j = 0 To UBound($g_aObsNames) - 1
			Local $aColData = $g_mData[$g_aObsNames[$j]]
			If $aColData[$i] = "" Then
				$bSkip = True
				ExitLoop
			EndIf
		Next
		Local $aDepData = $g_mData[$g_sDepVar]
		If $aDepData[$i] = "" Then $bSkip = True
		If $bSkip Then ContinueLoop

		; Add observations with individual σ from sd_ columns
		For $j = 0 To UBound($g_aObsNames) - 1
			Local $sObsName = $g_aObsNames[$j]
			Local $aOD = $g_mData[$sObsName]
			Local $fSd = 1.0
			Local $sSdKey = "sd_" & $sObsName
			If MapExists($g_mData, $sSdKey) Then
				Local $aSdCol = $g_mData[$sSdKey]
				If $i < UBound($aSdCol) And $aSdCol[$i] <> "" Then $fSd = Number($aSdCol[$i])
				If $fSd <= 0 Then $fSd = 1.0
			EndIf
			_adj_addObs($mSystem, $sObsName & $sObsSuffix, Number($aOD[$i]), $fSd, "s0")
		Next

		; Dependent variable with individual σ
		Local $fSdDep = 1.0
		Local $sSdDepKey = "sd_" & $g_sDepVar
		If MapExists($g_mData, $sSdDepKey) Then
			Local $aSdDepCol = $g_mData[$sSdDepKey]
			If $i < UBound($aSdDepCol) And $aSdDepCol[$i] <> "" Then $fSdDep = Number($aSdDepCol[$i])
			If $fSdDep <= 0 Then $fSdDep = 1.0
		EndIf
		_adj_addObs($mSystem, $g_sDepVar & $sObsSuffix, Number($aDepData[$i]), $fSdDep, "s0")

		; Condition equation
		Local $sIndexedFormula = $g_sFormula
		For $j = 0 To UBound($g_aObsNames) - 1
			$sIndexedFormula = StringRegExpReplace($sIndexedFormula, _
				"(?i)#" & $g_aObsNames[$j] & "\b", "#" & $g_aObsNames[$j] & $sObsSuffix)
		Next
		_adj_addFunction($mSystem, "#" & $g_sDepVar & $sObsSuffix & " - (" & $sIndexedFormula & ")")

		$iValid += 1
	Next

	Return $iValid
EndFunc   ;==>_BuildYorkModel

; ===== Formula parsing & column management ===================================

Func _ParseFormula()
	Local $sFormula = GUICtrlRead($g_idInputFormula)
	Local $sDepVar  = StringStripWS(GUICtrlRead($g_idInputDepVar), 3)
	If $sFormula = "" Then Return
	If $sDepVar = "" Then $sDepVar = "Y"
	$g_sDepVar = StringUpper($sDepVar)

	; Extract observation names (# prefix)
	Local $aObs = StringRegExp($sFormula, "#([A-Za-z]\w*)", 3)
	If @error Then
		Local $aEmpty[0]
		$aObs = $aEmpty
	EndIf

	; Remove duplicates and uppercase
	Local $mSeen[]
	Local $aObsClean[0]
	For $i = 0 To UBound($aObs) - 1
		Local $sName = StringUpper($aObs[$i])
		If Not MapExists($mSeen, $sName) Then
			$mSeen[$sName] = True
			ReDim $aObsClean[UBound($aObsClean) + 1]
			$aObsClean[UBound($aObsClean) - 1] = $sName
		EndIf
	Next

	; Extract parameter names (identifiers not preceded by # and not function names)
	Local $sClean = StringRegExpReplace($sFormula, "#[A-Za-z]\w*", "")
	Local $aIdents = StringRegExp($sClean, "\b([A-Za-z]\w*)\b", 3)
	If @error Then
		Local $aEmpty2[0]
		$aIdents = $aEmpty2
	EndIf

	; Filter: remove function names, keep unique
	Local $mParSeen[]
	Local $aParams[0]
	For $i = 0 To UBound($aIdents) - 1
		Local $sId = StringUpper($aIdents[$i])
		If StringInStr($g_sFuncNames, "|" & $sId & "|") Then ContinueLoop
		If $sId = $g_sDepVar Then ContinueLoop
		If Not MapExists($mParSeen, $sId) Then
			$mParSeen[$sId] = True
			ReDim $aParams[UBound($aParams) + 1]
			$aParams[UBound($aParams) - 1] = $sId
		EndIf
	Next

	; Check: dependent variable must not be an observation name
	For $i = 0 To UBound($aObsClean) - 1
		If $aObsClean[$i] = $g_sDepVar Then
			$g_sErrorMsg = "Abh. Variable '" & $g_sDepVar & "' kollidiert mit Beobachtung #" & $g_sDepVar
			GUICtrlSetData($g_idLblError, $g_sErrorMsg)
			Return
		EndIf
	Next

	$g_sFormula = $sFormula
	$g_aObsNames = $aObsClean
	$g_aParamNames = $aParams
	$g_sErrorMsg = ""
	GUICtrlSetData($g_idLblError, "")

	_RebuildDataColumns()
	_RebuildSdFields()
	_RebuildParamTable()
	_TryAutoSolve()
EndFunc   ;==>_ParseFormula

Func _RebuildDataColumns()
	; Backup current data
	Local $aKeys = MapKeys($g_mData)
	For $i = 0 To UBound($aKeys) - 1
		Local $aCol = $g_mData[$aKeys[$i]]
		If UBound($aCol) > 0 Then
			$g_mBackupData[$aKeys[$i]] = $aCol
		EndIf
	Next

	; Spaltenliste bestimmen (Input-Spalten)
	Local $iInputCols = UBound($g_aObsNames) + 1
	If $g_iRegType = $REG_YORK Then $iInputCols += UBound($g_aObsNames) + 1

	Local $aColNames[$iInputCols]
	Local $idx = 0
	For $i = 0 To UBound($g_aObsNames) - 1
		$aColNames[$idx] = $g_aObsNames[$i]
		$idx += 1
	Next
	$aColNames[$idx] = $g_sDepVar
	$idx += 1
	If $g_iRegType = $REG_YORK Then
		For $i = 0 To UBound($g_aObsNames) - 1
			$aColNames[$idx] = "sd_" & $g_aObsNames[$i]
			$idx += 1
		Next
		$aColNames[$idx] = "sd_" & $g_sDepVar
	EndIf

	; Schema bauen (nur Input-Spalten; Result-Spalten folgen in _UpdateDataResultColumns)
	Local $aSchema[$iInputCols]
	For $i = 0 To $iInputCols - 1
		Local $mCol[]
		$mCol.name     = $aColNames[$i]
		$mCol.type     = "input"
		$mCol.width    = 70
		$mCol.editable = True
		$aSchema[$i] = $mCol
	Next

	__cg_setSchema($g_iDataGrid, $aSchema)

	; Datenmodell rebuilden (wie zuvor)
	Local $mNewData[]
	For $i = 0 To $iInputCols - 1
		If MapExists($g_mBackupData, $aColNames[$i]) Then
			$mNewData[$aColNames[$i]] = $g_mBackupData[$aColNames[$i]]
		Else
			Local $aEmpty[$g_iRowCount]
			For $j = 0 To $g_iRowCount - 1
				$aEmpty[$j] = ""
			Next
			$mNewData[$aColNames[$i]] = $aEmpty
		EndIf
	Next
	$g_mData = $mNewData

	; Zeilenanzahl im Grid setzen und Werte spiegeln
	__cg_setRowCount($g_iDataGrid, $g_iRowCount)
	For $iRow = 0 To $g_iRowCount - 1
		For $i = 0 To $iInputCols - 1
			Local $aColData = $g_mData[$aColNames[$i]]
			Local $sVal = ($iRow < UBound($aColData)) ? $aColData[$iRow] : ""
			__cg_setCell($g_iDataGrid, $iRow, $i, $sVal)
		Next
	Next

	_EnsureEmptyTrailingRow()
EndFunc   ;==>_RebuildDataColumns

Func _RebuildParamTable()
	; Preserve existing user-entered values so that _ParseFormula (triggered by
	; Enter in the formula field or a CSV re-import) does not clobber start
	; values the user already set. Parameters with the same name keep their
	; current value; newly added parameters fall back to the default "0.5".
	Local $mOldValues[]
	If $g_iParamGrid >= 0 And $g_iParamGrid < UBound($g_aCellGrids) Then
		Local $mGridOld = $g_aCellGrids[$g_iParamGrid]
		If IsMap($mGridOld) Then
			Local $iOldRows = $mGridOld.iRowCount
			For $i = 0 To $iOldRows - 1
				Local $sOldName  = __cg_getCell($g_iParamGrid, $i, 0)
				Local $sOldValue = __cg_getCell($g_iParamGrid, $i, 1)
				If $sOldName <> "" And $sOldValue <> "" Then
					$mOldValues[StringUpper($sOldName)] = $sOldValue
				EndIf
			Next
		EndIf
	EndIf

	Local $aSchema[3]
	Local $mN[], $mW[], $mS[]
	$mN.name = "Name"
	$mN.type = "label"
	$mN.width = 55
	$mN.editable = False
	$mW.name = "Wert"
	$mW.type = "input"
	$mW.width = 90
	$mW.editable = True
	$mS.name = ChrW(963)
	$mS.type = "label"
	$mS.width = 70
	$mS.editable = False
	$aSchema[0] = $mN
	$aSchema[1] = $mW
	$aSchema[2] = $mS

	__cg_setSchema($g_iParamGrid, $aSchema)
	__cg_setRowCount($g_iParamGrid, UBound($g_aParamNames))
	For $i = 0 To UBound($g_aParamNames) - 1
		__cg_setCell($g_iParamGrid, $i, 0, $g_aParamNames[$i])
		Local $sStartVal = "0.5"
		If MapExists($mOldValues, $g_aParamNames[$i]) Then $sStartVal = $mOldValues[$g_aParamNames[$i]]
		__cg_setCell($g_iParamGrid, $i, 1, $sStartVal)
		__cg_setCell($g_iParamGrid, $i, 2, "")
	Next
EndFunc   ;==>_RebuildParamTable

Func _OnParamCellCommit($iGrid, $iRow, $iCol, $sValue)
	#forceref $iGrid
	; Nur Spalte 1 = "Wert" ist editierbar; Parameter-Änderung löst KEIN Auto-Solve aus
	; (User setzt evtl. mehrere Startwerte hintereinander, bevor neu gerechnet wird)
	; Komma → Punkt normalisieren (deutsche Tastatur / Excel-Paste)
	If StringInStr($sValue, ",") Then
		__cg_setCell($g_iParamGrid, $iRow, $iCol, StringReplace($sValue, ",", "."))
	EndIf
	; Visuelles Feedback: Kurve + Prädiktion mit neuen Näherungswerten neu zeichnen
	_DrawGraph()
	_UpdatePrediction()
EndFunc   ;==>_OnParamCellCommit

; Liest aktuelle Parameter-Werte aus dem Parametergrid als Map (Name → Zahl).
; Komma wird defensiv durch Punkt ersetzt, leere Zellen liefern Default 0.5.
Func _ReadParamValuesFromGrid()
	Local $mVals[]
	If $g_iParamGrid < 0 Then Return $mVals
	For $i = 0 To UBound($g_aParamNames) - 1
		Local $sVal = __cg_getCell($g_iParamGrid, $i, 1)
		$sVal = StringReplace($sVal, ",", ".")
		Local $fVal = ($sVal = "") ? 0.5 : Number($sVal)
		$mVals[$g_aParamNames[$i]] = $fVal
	Next
	Return $mVals
EndFunc   ;==>_ReadParamValuesFromGrid

Func _UpdateResults()
	If Not $g_bAdjOK Or Not IsMap($g_mResults) Then
		GUICtrlSetData($g_idLblStats, "")
		; Clear previous result values
		For $i = 0 To UBound($g_aParamNames) - 1
			__cg_setCell($g_iParamGrid, $i, 1, "")
			__cg_setCell($g_iParamGrid, $i, 2, "")
		Next
		_UpdateDataResultColumns() ; clear result columns since $g_bAdjOK is False
		Return
	EndIf

	; Update parameter table: Wert + σ columns
	Local $mX   = $g_mResults.x1
	Local $mSdx = $g_mResults.sdx
	For $i = 0 To UBound($g_aParamNames) - 1
		Local $sName = $g_aParamNames[$i]
		If IsMap($mX) And MapExists($mX, $sName) Then
			__cg_setCell($g_iParamGrid, $i, 1, StringFormat("%.6g", $mX[$sName]))
		EndIf
		If IsMap($mSdx) And MapExists($mSdx, $sName) Then
			__cg_setCell($g_iParamGrid, $i, 2, StringFormat("%.4g", $mSdx[$sName]))
		Else
			__cg_setCell($g_iParamGrid, $i, 2, "---")
		EndIf
	Next

	; Statistics with Unicode
	Local $sStats = StringFormat("s" & ChrW(0x2080) & " = %.6g", $g_mResults.s0) & @CRLF ; s₀
	$sStats &= "f  = " & $g_mResults.f & @CRLF
	$sStats &= StringFormat("v" & ChrW(0x1D40) & "Pv = %.6g", $g_mResults.vtpv) ; vᵀPv

	; R² (coefficient of determination)
	Local $fR2 = __ComputeR2()
	If $fR2 > -1e30 Then
		$sStats &= @CRLF & StringFormat("R" & ChrW(0x00B2) & " = %.6f", $fR2) ; R²
	EndIf

	; Globaltest (if available)
	If MapExists($g_mResults, "globalTestPassed") Then
		$sStats &= @CRLF & StringFormat(ChrW(0x03C7) & ChrW(0x00B2) & ": T=%.3g [%.3g..%.3g] %s", _
			$g_mResults.globalTestT, $g_mResults.globalTestLower, $g_mResults.globalTestUpper, _
			$g_mResults.globalTestPassed ? ChrW(0x2713) : ChrW(0x2717)) ; χ² ... ✓/✗
	EndIf

	GUICtrlSetData($g_idLblStats, $sStats)

	_UpdateDataResultColumns()
	_UpdatePrediction()
EndFunc   ;==>_UpdateResults

Func _UpdateDataResultColumns()
	; Aktuelles Schema merken (nur Input-Spalten; Result-Spalten kommen immer neu)
	Local $mGrid = $g_aCellGrids[$g_iDataGrid]
	Local $aCurSchema = $mGrid.aSchema
	Local $iInputCount = 0
	Local $mSchCol
	For $i = 0 To UBound($aCurSchema) - 1
		$mSchCol = $aCurSchema[$i]
		If $mSchCol.type = "input" Then $iInputCount += 1
		If $mSchCol.type = "separator" Then ExitLoop
	Next
	; Wenn bereits Result-Spalten dran hängen, Schema auf reine Input-Spalten kürzen
	Local $aInputSchema[$iInputCount]
	Local $idx = 0
	For $i = 0 To UBound($aCurSchema) - 1
		$mSchCol = $aCurSchema[$i]
		If $mSchCol.type = "input" Then
			$aInputSchema[$idx] = $mSchCol
			$idx += 1
		EndIf
	Next

	; Shared locals für beide (Reset- und Fill-)Pfade
	Local $mCol, $sColName, $aData, $iRow, $iCol

	If Not $g_bAdjOK Or Not IsMap($g_mResults) Then
		__cg_setSchema($g_iDataGrid, $aInputSchema)
		; Daten wieder reinschreiben (Schema-Wechsel hat Zeilen geleert)
		For $iRow = 0 To $g_iRowCount - 1
			For $iCol = 0 To $iInputCount - 1
				$mCol = $aInputSchema[$iCol]
				$sColName = $mCol.name
				If MapExists($g_mData, $sColName) Then
					$aData = $g_mData[$sColName]
					If $iRow < UBound($aData) Then
						__cg_setCell($g_iDataGrid, $iRow, $iCol, $aData[$iRow])
					EndIf
				EndIf
			Next
		Next
		Return
	EndIf

	; Ziel-Observation-Liste
	Local $aResObs[0]
	If $g_iRegType = $REG_NORMAL Then
		ReDim $aResObs[1]
		$aResObs[0] = $g_sDepVar
	Else
		ReDim $aResObs[UBound($g_aObsNames) + 1]
		For $i = 0 To UBound($g_aObsNames) - 1
			$aResObs[$i] = $g_aObsNames[$i]
		Next
		$aResObs[UBound($g_aObsNames)] = $g_sDepVar
	EndIf

	; Gesamtschema: Input-Spalten + Separator + 3 Spalten pro Observation
	Local $iTotal = $iInputCount + 1 + UBound($aResObs) * 3
	Local $aNewSchema[$iTotal]
	For $i = 0 To $iInputCount - 1
		$aNewSchema[$i] = $aInputSchema[$i]
	Next
	Local $mSep[]
	$mSep.name = ""
	$mSep.type = "separator"
	$mSep.width = 5
	$mSep.editable = False
	$aNewSchema[$iInputCount] = $mSep
	Local $iPos = $iInputCount + 1
	Local $mResV, $mResS, $mResR
	For $j = 0 To UBound($aResObs) - 1
		Local $sObs = $aResObs[$j]
		Local $mEmpty1[], $mEmpty2[], $mEmpty3[]
		$mResV = $mEmpty1
		$mResS = $mEmpty2
		$mResR = $mEmpty3
		$mResV.name = "v(" & $sObs & ")"
		$mResV.type = "label"
		$mResV.width = 60
		$mResV.editable = False
		$mResS.name = ChrW(963) & "(" & $sObs & ")"
		$mResS.type = "label"
		$mResS.width = 60
		$mResS.editable = False
		$mResR.name = "r(" & $sObs & ")"
		$mResR.type = "label"
		$mResR.width = 50
		$mResR.editable = False
		$aNewSchema[$iPos] = $mResV
		$aNewSchema[$iPos + 1] = $mResS
		$aNewSchema[$iPos + 2] = $mResR
		$iPos += 3
	Next

	__cg_setSchema($g_iDataGrid, $aNewSchema)

	; Input-Spalten-Daten wieder einfüllen
	For $iRow = 0 To $g_iRowCount - 1
		For $iCol = 0 To $iInputCount - 1
			$mCol = $aInputSchema[$iCol]
			$sColName = $mCol.name
			If MapExists($g_mData, $sColName) Then
				$aData = $g_mData[$sColName]
				If $iRow < UBound($aData) Then
					__cg_setCell($g_iDataGrid, $iRow, $iCol, $aData[$iRow])
				EndIf
			EndIf
		Next
	Next

	; Result-Spalten befüllen
	Local $mV = Null, $mSdv = Null, $mR = Null
	If MapExists($g_mResults, "v")   Then $mV   = $g_mResults.v
	If MapExists($g_mResults, "sdv") Then $mSdv = $g_mResults.sdv
	If MapExists($g_mResults, "r")   Then $mR   = $g_mResults.r

	For $iRow = 0 To $g_iRowCount - 1
		Local $sSuffix = String($iRow + 1)
		$iCol = $iInputCount + 1
		For $j = 0 To UBound($aResObs) - 1
			Local $sKey = $aResObs[$j] & $sSuffix
			If IsMap($mV)   And MapExists($mV,   $sKey) Then __cg_setCell($g_iDataGrid, $iRow, $iCol,     StringFormat("%.4g", $mV[$sKey]))
			If IsMap($mSdv) And MapExists($mSdv, $sKey) Then __cg_setCell($g_iDataGrid, $iRow, $iCol + 1, StringFormat("%.4g", $mSdv[$sKey]))
			If IsMap($mR)   And MapExists($mR,   $sKey) Then __cg_setCell($g_iDataGrid, $iRow, $iCol + 2, StringFormat("%.3f", $mR[$sKey]))
			$iCol += 3
		Next
	Next
EndFunc   ;==>_UpdateDataResultColumns

; ===== Row management ========================================================

Func _AddRow()
	If UBound($g_aObsNames) = 0 Then Return ; no formula parsed yet

	$g_iRowCount += 1

	; Add empty values to each data column
	Local $aKeys = MapKeys($g_mData)
	For $i = 0 To UBound($aKeys) - 1
		Local $aCol = $g_mData[$aKeys[$i]]
		ReDim $aCol[$g_iRowCount]
		$aCol[$g_iRowCount - 1] = ""
		$g_mData[$aKeys[$i]] = $aCol
	Next

	__cg_setRowCount($g_iDataGrid, $g_iRowCount)
EndFunc   ;==>_AddRow

Func _EnsureEmptyTrailingRow()
	If UBound($g_aObsNames) = 0 Then Return ; no formula parsed yet

	; Check if we have any rows and the last row has data
	If $g_iRowCount = 0 Then
		_AddRow()
		Return
	EndIf

	; Check if last row is fully empty
	Local $bLastRowEmpty = True
	Local $aKeys = MapKeys($g_mData)
	For $i = 0 To UBound($aKeys) - 1
		; Skip sd_ columns for emptiness check
		If StringLeft($aKeys[$i], 3) = "sd_" Then ContinueLoop
		Local $aCol = $g_mData[$aKeys[$i]]
		If $g_iRowCount - 1 < UBound($aCol) And $aCol[$g_iRowCount - 1] <> "" Then
			$bLastRowEmpty = False
			ExitLoop
		EndIf
	Next

	If Not $bLastRowEmpty Then _AddRow()
EndFunc   ;==>_EnsureEmptyTrailingRow

Func _DeleteSelectedRow()
	Local $mFocus = __cg_getFocus($g_iDataGrid)
	If Not IsMap($mFocus) Then Return
	Local $iIdx = $mFocus.row
	If $iIdx < 0 Or $iIdx >= $g_iRowCount Then Return

	Local $aKeys = MapKeys($g_mData)
	For $k = 0 To UBound($aKeys) - 1
		Local $aCol = $g_mData[$aKeys[$k]]
		For $i = $iIdx To $g_iRowCount - 2
			$aCol[$i] = $aCol[$i + 1]
		Next
		ReDim $aCol[$g_iRowCount - 1]
		$g_mData[$aKeys[$k]] = $aCol
	Next
	$g_iRowCount -= 1
	__cg_setRowCount($g_iDataGrid, $g_iRowCount)
	_TryAutoSolve()
EndFunc   ;==>_DeleteSelectedRow

Func _ClearAllData()
	__cg_setRowCount($g_iDataGrid, 0)
	$g_iRowCount = 0
	Local $aKeys = MapKeys($g_mData)
	For $i = 0 To UBound($aKeys) - 1
		Local $aEmpty[0]
		$g_mData[$aKeys[$i]] = $aEmpty
	Next
	$g_bAdjOK = False
	$g_mResults = Null
	_EnsureEmptyTrailingRow()
EndFunc   ;==>_ClearAllData

; ===== CellGrid commit callbacks =============================================

Func _OnDataCellCommit($iGrid, $iRow, $iCol, $sValue)
	Local $mGrid = $g_aCellGrids[$iGrid]
	Local $aSchema = $mGrid.aSchema
	Local $mColDef = $aSchema[$iCol]
	Local $sColName = $mColDef.name

	; Nur Input-Spalten — Sicherheitsprüfung
	If $mColDef.type <> "input" Then Return

	; Komma → Punkt normalisieren (deutsche Tastatur / Excel-Paste)
	If StringInStr($sValue, ",") Then
		$sValue = StringReplace($sValue, ",", ".")
		__cg_setCell($iGrid, $iRow, $iCol, $sValue)
	EndIf

	; Zielspalte muss in $g_mData existieren
	If Not MapExists($g_mData, $sColName) Then Return
	Local $aCol = $g_mData[$sColName]
	If UBound($aCol) <= $iRow Then ReDim $aCol[$iRow + 1]
	$aCol[$iRow] = $sValue
	$g_mData[$sColName] = $aCol

	_EnsureEmptyTrailingRow()
	_TryAutoSolve()
EndFunc   ;==>_OnDataCellCommit

Func _OnDataBulkCommit($iGrid, $iStartRow, $iEndRow)
	#forceref $iGrid, $iStartRow, $iEndRow
	; Nach Paste: alle geänderten Zellen liegen in .aValues — nach $g_mData spiegeln.
	; Komma → Punkt für jede Zelle normalisieren (deutsche Excel-Paste).
	Local $mGrid = $g_aCellGrids[$g_iDataGrid]
	Local $aSchema = $mGrid.aSchema
	Local $aValues = $mGrid.aValues
	For $iCol = 0 To UBound($aSchema) - 1
		Local $mColDef = $aSchema[$iCol]
		If $mColDef.type <> "input" Then ContinueLoop
		Local $sColName = $mColDef.name
		If Not MapExists($g_mData, $sColName) Then ContinueLoop
		Local $aCol = $g_mData[$sColName]
		If UBound($aCol) < $mGrid.iRowCount Then ReDim $aCol[$mGrid.iRowCount]
		For $iRow = 0 To $mGrid.iRowCount - 1
			Local $sCell = $aValues[$iRow][$iCol]
			If StringInStr($sCell, ",") Then
				$sCell = StringReplace($sCell, ",", ".")
				__cg_setCell($g_iDataGrid, $iRow, $iCol, $sCell)
			EndIf
			$aCol[$iRow] = $sCell
		Next
		$g_mData[$sColName] = $aCol
	Next
	$g_iRowCount = $mGrid.iRowCount
	_EnsureEmptyTrailingRow()
	_TryAutoSolve()
EndFunc   ;==>_OnDataBulkCommit

; ===== CSV Import ============================================================

Func _ExportCSV()
	If UBound($g_aObsNames) = 0 Or $g_iRowCount = 0 Then
		$g_sErrorMsg = "Keine Daten zum Exportieren"
		GUICtrlSetData($g_idLblError, $g_sErrorMsg)
		Return
	EndIf

	Local $sFile = FileSaveDialog("CSV-Datei speichern", @ScriptDir, _
		"CSV-Dateien (*.csv)|Alle Dateien (*.*)", 16, "", $g_hGUI) ; 16 = prompt overwrite
	If @error Then Return

	; Ensure .csv extension
	If StringRight($sFile, 4) <> ".csv" Then $sFile &= ".csv"

	; Build column list: input data columns only (no sd_ for Normal/Ortho/Deming, with sd_ for York)
	Local $aColNames[0]
	; Observation columns
	For $i = 0 To UBound($g_aObsNames) - 1
		ReDim $aColNames[UBound($aColNames) + 1]
		$aColNames[UBound($aColNames) - 1] = $g_aObsNames[$i]
	Next
	; Dependent variable
	ReDim $aColNames[UBound($aColNames) + 1]
	$aColNames[UBound($aColNames) - 1] = $g_sDepVar
	; sd_ columns for York
	If $g_iRegType = $REG_YORK Then
		For $i = 0 To UBound($g_aObsNames) - 1
			ReDim $aColNames[UBound($aColNames) + 1]
			$aColNames[UBound($aColNames) - 1] = "sd_" & $g_aObsNames[$i]
		Next
		ReDim $aColNames[UBound($aColNames) + 1]
		$aColNames[UBound($aColNames) - 1] = "sd_" & $g_sDepVar
	EndIf

	; Build file content
	Local $sOut = ""

	; Metadata as comments (parseable by _ImportCSV)
	$sOut &= "# formula=" & $g_sFormula & @CRLF
	$sOut &= "# depvar=" & $g_sDepVar & @CRLF
	; Parameter values (initial/current values for solver)
	For $i = 0 To UBound($g_aParamNames) - 1
		Local $sParVal = __cg_getCell($g_iParamGrid, $i, 1)
		$sOut &= "# param_" & $g_aParamNames[$i] & "=" & $sParVal & @CRLF
	Next

	; Header
	For $i = 0 To UBound($aColNames) - 1
		If $i > 0 Then $sOut &= ","
		$sOut &= $aColNames[$i]
	Next
	$sOut &= @CRLF

	; Data rows (skip trailing empty row)
	For $iRow = 0 To $g_iRowCount - 1
		; Check if row is completely empty (trailing row)
		Local $bEmpty = True
		For $j = 0 To UBound($aColNames) - 1
			If MapExists($g_mData, $aColNames[$j]) Then
				Local $aCol = $g_mData[$aColNames[$j]]
				If $iRow < UBound($aCol) And $aCol[$iRow] <> "" Then
					$bEmpty = False
					ExitLoop
				EndIf
			EndIf
		Next
		If $bEmpty Then ContinueLoop

		For $j = 0 To UBound($aColNames) - 1
			If $j > 0 Then $sOut &= ","
			If MapExists($g_mData, $aColNames[$j]) Then
				Local $aColD = $g_mData[$aColNames[$j]]
				If $iRow < UBound($aColD) Then $sOut &= $aColD[$iRow]
			EndIf
		Next
		$sOut &= @CRLF
	Next

	FileDelete($sFile)
	FileWrite($sFile, $sOut)
	$g_sErrorMsg = ""
	GUICtrlSetData($g_idLblError, "Gespeichert: " & $sFile)
EndFunc   ;==>_ExportCSV

; Full result report: runs a second _adj_solve with all compute flags enabled
; (diagnostics, reliability, cofactors, redundancy, globalTest) and writes a
; comprehensive German-language text report to a user-chosen file.
Func _ExportResults()
	If $g_sFormula = "" Or UBound($g_aObsNames) = 0 Then
		$g_sErrorMsg = "Keine Formel — Ergebnisexport nicht moeglich"
		GUICtrlSetData($g_idLblError, $g_sErrorMsg)
		Return
	EndIf
	If _CountValidRows() < UBound($g_aParamNames) Then
		$g_sErrorMsg = "Zu wenige gueltige Datenzeilen fuer Ergebnisexport"
		GUICtrlSetData($g_idLblError, $g_sErrorMsg)
		Return
	EndIf

	; Build default filename: {csvbase}_result_YYYYMMDD_HHMM.txt
	Local $sBase = "regression"
	If $g_sLastImportedCSV <> "" Then
		Local $sFName = $g_sLastImportedCSV
		Local $iSep = StringInStr($sFName, "\", 0, -1)
		If $iSep > 0 Then $sFName = StringMid($sFName, $iSep + 1)
		Local $iDot = StringInStr($sFName, ".", 0, -1)
		If $iDot > 0 Then $sFName = StringLeft($sFName, $iDot - 1)
		If $sFName <> "" Then $sBase = $sFName
	EndIf
	Local $sStamp = StringFormat("%04d%02d%02d_%02d%02d", @YEAR, @MON, @MDAY, @HOUR, @MIN)
	Local $sDefault = $sBase & "_result_" & $sStamp & ".txt"

	Local $sFile = FileSaveDialog("Ergebnisreport speichern", @ScriptDir, _
		"Textdateien (*.txt)|Alle Dateien (*.*)", 16, $sDefault, $g_hGUI)
	If @error Then Return
	If StringRight($sFile, 4) <> ".txt" And StringInStr($sFile, ".") = 0 Then $sFile &= ".txt"

	; Read fresh sd values and build a parallel system with full compute flags
	_ReadSdValues()
	Local $mSystem = _adj_createSystem()
	; Prefer the true start values captured before the original solve; fall back
	; to the current grid if no prior solve happened in this session.
	Local $bHaveSnap = (UBound(MapKeys($g_mLastStartValues)) > 0)
	For $i = 0 To UBound($g_aParamNames) - 1
		Local $sName = $g_aParamNames[$i]
		Local $sStart = ""
		If $bHaveSnap And MapExists($g_mLastStartValues, $sName) Then
			$sStart = $g_mLastStartValues[$sName]
		Else
			$sStart = __cg_getCell($g_iParamGrid, $i, 1)
		EndIf
		Local $fStart = Number($sStart)
		If $sStart = "" Then $fStart = 0.5
		_adj_setInitialValue($mSystem, $sName, $fStart)
	Next

	Switch $g_iRegType
		Case $REG_NORMAL
			_BuildNormalModel($mSystem)
		Case $REG_ORTHO, $REG_DEMING
			_BuildGLMModel($mSystem)
		Case $REG_YORK
			_BuildYorkModel($mSystem)
	EndSwitch

	Local $mConfig = _adj_defaultConfig("LM", False)
	$mConfig.maxIterations = 100
	Local $mCompute = $mConfig.compute
	$mCompute.qxx         = True
	$mCompute.cofactors   = True
	$mCompute.redundancy  = True
	$mCompute.globalTest  = True
	$mCompute.diagnostics = True
	$mCompute.reliability = True
	$mConfig.compute = $mCompute
	If $g_bRobust Then
		$mConfig.robust = GUICtrlRead($g_idComboEstimator)
		$mConfig.robustParams = _adj_robustDefaults($mConfig.robust)
	EndIf

	_adj_solve($mSystem, $mConfig)
	Local $iErr = @error, $iExt = @extended
	If $iErr Then
		$g_sErrorMsg = "Report-Solve fehlgeschlagen: " & _adj_getErrorMessage($iErr, $iExt)
		GUICtrlSetData($g_idLblError, $g_sErrorMsg)
		Return
	EndIf

	Local $sReport = _BuildResultsReport($mSystem, $sFile)

	; Write UTF-8 with BOM
	FileDelete($sFile)
	Local $hFile = FileOpen($sFile, 2 + 8 + 128)
	If $hFile = -1 Then
		$g_sErrorMsg = "Datei konnte nicht geschrieben werden: " & $sFile
		GUICtrlSetData($g_idLblError, $g_sErrorMsg)
		Return
	EndIf
	FileWrite($hFile, $sReport)
	FileClose($hFile)

	$g_sErrorMsg = ""
	GUICtrlSetData($g_idLblError, "Report gespeichert: " & $sFile)
EndFunc   ;==>_ExportResults

Func _BuildResultsReport(ByRef $mSystem, $sFilePath)
	Local $mRes = _adj_getResults($mSystem)
	Local $mX1  = $mRes.x1
	Local $mSdx = $mRes.sdx

	Local $CRLF = @CRLF
	Local $sSep  = "=============================================================" & $CRLF
	Local $sSep2 = "-------------------------------------------------------------" & $CRLF

	Local $s = ""
	$s &= $sSep
	$s &= "  Interaktive Regression — Ergebnisreport" & $CRLF
	$s &= $sSep

	; --- 1. Kopfdaten ---
	$s &= StringFormat("Datum            : %04d-%02d-%02d %02d:%02d:%02d", @YEAR, @MON, @MDAY, @HOUR, @MIN, @SEC) & $CRLF
	$s &= "Zieldatei        : " & $sFilePath & $CRLF
	$s &= "Quell-CSV        : " & ($g_sLastImportedCSV = "" ? "(keine — manuell eingegeben)" : $g_sLastImportedCSV) & $CRLF
	$s &= "Gueltige Zeilen  : " & _CountValidRows() & " (von " & $g_iRowCount & " Grid-Zeilen)" & $CRLF
	$s &= $CRLF

	; --- 2. Modell ---
	$s &= $sSep2
	$s &= "  Modell" & $CRLF
	$s &= $sSep2
	Local $sRegType = "Normal (OLS — nur Y hat Fehler)"
	Switch $g_iRegType
		Case $REG_ORTHO
			$sRegType = "Orthogonal (TLS — X und Y haben gleichen Fehler)"
		Case $REG_DEMING
			$sRegType = "Deming (X und Y haben unterschiedlichen konstanten Fehler)"
		Case $REG_YORK
			$sRegType = "York (pro Messpunkt individuelle Fehler)"
	EndSwitch
	$s &= "Regressionstyp   : " & $sRegType & $CRLF
	$s &= "Formel           : " & $g_sFormula & $CRLF
	$s &= "Abh. Variable    : " & $g_sDepVar & $CRLF
	$s &= "Beobachtungen    : " & _ResJoinArr($g_aObsNames, ", ") & $CRLF
	$s &= "Parameter        : " & _ResJoinArr($g_aParamNames, ", ") & $CRLF
	$s &= "Robuste Schaetz. : " & ($g_bRobust ? "Ja (Estimator: " & GUICtrlRead($g_idComboEstimator) & ")" : "Nein") & $CRLF
	$s &= "Standardabw.     :" & $CRLF
	Local $aSdKeys = MapKeys($g_mSdValues)
	For $k = 0 To UBound($aSdKeys) - 1
		Local $sLabel = ($aSdKeys[$k] = "__shared__") ? "σ (geteilt)" : "σ(" & $aSdKeys[$k] & ")"
		$s &= StringFormat("    %-14s = %s", $sLabel, String($g_mSdValues[$aSdKeys[$k]])) & $CRLF
	Next
	$s &= $CRLF

	; --- 3. Startwerte ---
	$s &= $sSep2
	$s &= "  Startwerte (Naeherungswerte, wie in den Solver gegeben)" & $CRLF
	$s &= $sSep2
	Local $bHaveSnapRep = (UBound(MapKeys($g_mLastStartValues)) > 0)
	For $i = 0 To UBound($g_aParamNames) - 1
		Local $sName = $g_aParamNames[$i]
		Local $sSt = ""
		If $bHaveSnapRep And MapExists($g_mLastStartValues, $sName) Then
			$sSt = $g_mLastStartValues[$sName]
		Else
			$sSt = __cg_getCell($g_iParamGrid, $i, 1) ; fallback for first-ever report
		EndIf
		$s &= StringFormat("  %-8s = %s", $sName, $sSt) & $CRLF
	Next
	$s &= $CRLF

	; --- 4. Rohdaten ---
	$s &= $sSep2
	$s &= "  Rohdaten (Eingabe)" & $CRLF
	$s &= $sSep2
	; Build column list: Obs + DepVar, optionally sd-columns for York
	Local $aCols[0]
	For $i = 0 To UBound($g_aObsNames) - 1
		ReDim $aCols[UBound($aCols) + 1]
		$aCols[UBound($aCols) - 1] = $g_aObsNames[$i]
	Next
	ReDim $aCols[UBound($aCols) + 1]
	$aCols[UBound($aCols) - 1] = $g_sDepVar
	If $g_iRegType = $REG_YORK Then
		For $i = 0 To UBound($g_aObsNames) - 1
			ReDim $aCols[UBound($aCols) + 1]
			$aCols[UBound($aCols) - 1] = "sd_" & $g_aObsNames[$i]
		Next
		ReDim $aCols[UBound($aCols) + 1]
		$aCols[UBound($aCols) - 1] = "sd_" & $g_sDepVar
	EndIf
	; Header line
	$s &= StringFormat("  %-4s", "#")
	For $j = 0 To UBound($aCols) - 1
		$s &= StringFormat("  %14s", $aCols[$j])
	Next
	$s &= $CRLF
	; Data lines
	For $iRow = 0 To $g_iRowCount - 1
		Local $bEmpty = True
		For $j = 0 To UBound($aCols) - 1
			If MapExists($g_mData, $aCols[$j]) Then
				Local $aColD = $g_mData[$aCols[$j]]
				If $iRow < UBound($aColD) And $aColD[$iRow] <> "" Then
					$bEmpty = False
					ExitLoop
				EndIf
			EndIf
		Next
		If $bEmpty Then ContinueLoop
		$s &= StringFormat("  %-4d", $iRow + 1)
		For $j = 0 To UBound($aCols) - 1
			Local $sVal = ""
			If MapExists($g_mData, $aCols[$j]) Then
				Local $aColX = $g_mData[$aCols[$j]]
				If $iRow < UBound($aColX) Then $sVal = $aColX[$iRow]
			EndIf
			$s &= StringFormat("  %14s", $sVal)
		Next
		$s &= $CRLF
	Next
	$s &= $CRLF

	; --- 5. Güteewerte der Ausgleichung ---
	$s &= $sSep2
	$s &= "  Guetewerte der Ausgleichung" & $CRLF
	$s &= $sSep2
	$s &= "Algorithmus      : " & $mRes.algorithm & " (" & $mRes.modelType & ")" & $CRLF
	$s &= "Solver           : " & $mRes.solver
	If IsNumber($mRes.conditionNumber) Then $s &= StringFormat("  (Konditionszahl: %.3e)", $mRes.conditionNumber)
	$s &= $CRLF
	$s &= "Iterationen      : " & $mRes.nIterations & $CRLF
	$s &= "Freiheitsgrade f : " & $mRes.f & $CRLF
	$s &= StringFormat("s0               : %.6g", $mRes.s0) & $CRLF
	$s &= StringFormat("vtPv             : %.6g", $mRes.vtpv) & $CRLF
	$s &= "Rangdefizient    : " & ($mRes.rankDeficient ? "JA (Warnung!)" : "nein") & $CRLF
	If MapExists($mRes, "globalTestPassed") Then
		Local $sPassed = $mRes.globalTestPassed ? "PASSED" : "FAILED"
		$s &= StringFormat("Globaltest       : %s  (alpha=%.4f, T=%.6g, chi²-Intervall=[%.6g, %.6g])", _
			$sPassed, $mRes.globalTestAlpha, $mRes.globalTestT, $mRes.globalTestLower, $mRes.globalTestUpper) & $CRLF
	EndIf
	If $g_bRobust And MapExists($mRes, "robustConverged") Then
		$s &= "Robust           : " & ($mRes.robustConverged ? "konvergiert" : "NICHT konvergiert") _
			& "  (" & $mRes.robustIterations & " IRLS-Iterationen, Skala=" & StringFormat("%.6g", $mRes.robustScale) & ")" & $CRLF
	EndIf
	$s &= $CRLF

	; --- 6. Ausgeglichene Parameter ---
	$s &= $sSep2
	$s &= "  Ausgeglichene Parameter" & $CRLF
	$s &= $sSep2
	$s &= StringFormat("  %-8s  %16s  %16s  %s", "Name", "x̂", "sd(x̂)", "95%-Konfidenzintervall") & $CRLF
	For $i = 0 To UBound($g_aParamNames) - 1
		$sName = $g_aParamNames[$i]
		Local $fX  = MapExists($mX1, $sName)  ? $mX1[$sName]  : 0
		Local $fSd = MapExists($mSdx, $sName) ? $mSdx[$sName] : 0
		Local $fLo = $fX - 1.96 * $fSd
		Local $fHi = $fX + 1.96 * $fSd
		$s &= StringFormat("  %-8s  %16.8g  %16.6g  [%.6g, %.6g]", $sName, $fX, $fSd, $fLo, $fHi) & $CRLF
	Next
	$s &= $CRLF

	; --- 7. Ausgeglichene Beobachtungen ---
	$s &= $sSep2
	$s &= "  Ausgeglichene Beobachtungen" & $CRLF
	$s &= $sSep2
	Local $bHasDiag   = MapExists($mRes, "baardaW")
	Local $bHasRobust = $g_bRobust And MapExists($mRes, "robustWeights")
	; Header
	$s &= StringFormat("  %-4s  %-10s  %12s  %12s  %12s  %10s  %8s", "#", "Name", "l_i", "l_i^", "v", "sd(v)", "r")
	If $bHasDiag Then $s &= StringFormat("  %10s  %10s  %10s", "|w|", "T_Pope", "MDB")
	If $bHasRobust Then $s &= StringFormat("  %8s", "w_rob")
	$s &= $CRLF

	Local $mObsVal  = $mRes.obsValue
	Local $mObsAdj  = $mRes.obsAdj
	Local $mV      = $mRes.v
	Local $mSdv    = $mRes.sdv
	Local $mR      = $mRes.r
	Local $mW      = $bHasDiag   ? $mRes.baardaW      : Null
	Local $mTp     = $bHasDiag   ? $mRes.popeT        : Null
	Local $mMdb    = $bHasDiag   ? $mRes.mdb          : Null
	Local $mRobW   = $bHasRobust ? $mRes.robustWeights : Null

	Local $iCounter = 0
	For $iRow = 0 To $g_iRowCount - 1
		Local $sSuffix = String($iRow + 1)
		; GLM/York: list X-observations first, then dep-var
		If $g_iRegType <> $REG_NORMAL Then
			For $j = 0 To UBound($g_aObsNames) - 1
				Local $sKey = StringUpper($g_aObsNames[$j] & $sSuffix)
				If Not MapExists($mObsVal, $sKey) Then ContinueLoop
				$iCounter += 1
				$s &= _ResFormatObsRow($iCounter, $sKey, $mObsVal, $mObsAdj, $mV, $mSdv, $mR, $mW, $mTp, $mMdb, $mRobW, $bHasDiag, $bHasRobust) & $CRLF
			Next
		EndIf
		Local $sDepKey = StringUpper($g_sDepVar & $sSuffix)
		If Not MapExists($mObsVal, $sDepKey) Then ContinueLoop
		$iCounter += 1
		$s &= _ResFormatObsRow($iCounter, $sDepKey, $mObsVal, $mObsAdj, $mV, $mSdv, $mR, $mW, $mTp, $mMdb, $mRobW, $bHasDiag, $bHasRobust) & $CRLF
	Next
	$s &= $CRLF

	; --- 8. Ausgeglichene Formel ---
	$s &= $sSep2
	$s &= "  Ausgeglichene Formel" & $CRLF
	$s &= $sSep2
	$s &= _ResSubstitutedFormula($g_sFormula, $mX1) & $CRLF
	$s &= $CRLF

	Return $s
EndFunc   ;==>_BuildResultsReport

Func _ResFormatObsRow($iIdx, $sKey, ByRef $mV1, ByRef $mV2, ByRef $mV3, ByRef $mV4, ByRef $mV5, _
		ByRef $mV6, ByRef $mV7, ByRef $mV8, ByRef $mV9, $bHasDiag, $bHasRobust)
	; Builds one row of the observation table. Vectors passed ByRef to avoid copies.
	Local $sRow = StringFormat("  %-4d  %-10s  %12.6g  %12.6g  %12.6g  %10.4g  %8.4g", _
		$iIdx, $sKey, _
		IsMap($mV1) And MapExists($mV1, $sKey) ? $mV1[$sKey] : 0, _
		IsMap($mV2) And MapExists($mV2, $sKey) ? $mV2[$sKey] : 0, _
		IsMap($mV3) And MapExists($mV3, $sKey) ? $mV3[$sKey] : 0, _
		IsMap($mV4) And MapExists($mV4, $sKey) ? $mV4[$sKey] : 0, _
		IsMap($mV5) And MapExists($mV5, $sKey) ? $mV5[$sKey] : 0)
	If $bHasDiag Then
		$sRow &= StringFormat("  %10.4g  %10.4g  %10.4g", _
			IsMap($mV6) And MapExists($mV6, $sKey) ? $mV6[$sKey] : 0, _
			IsMap($mV7) And MapExists($mV7, $sKey) ? $mV7[$sKey] : 0, _
			IsMap($mV8) And MapExists($mV8, $sKey) ? $mV8[$sKey] : 0)
	EndIf
	If $bHasRobust Then
		$sRow &= StringFormat("  %8.4g", _
			IsMap($mV9) And MapExists($mV9, $sKey) ? $mV9[$sKey] : 1.0)
	EndIf
	Return $sRow
EndFunc   ;==>_ResFormatObsRow

Func _ResJoinArr(ByRef $aArr, $sSep)
	If Not IsArray($aArr) Or UBound($aArr) = 0 Then Return ""
	Local $s = $aArr[0]
	For $i = 1 To UBound($aArr) - 1
		$s &= $sSep & $aArr[$i]
	Next
	Return $s
EndFunc   ;==>_ResJoinArr

Func _ResSubstitutedFormula($sFormula, ByRef $mX1)
	Local $sOut = $sFormula
	Local $aKeys = MapKeys($mX1)
	For $i = 0 To UBound($aKeys) - 1
		Local $sName = $aKeys[$i]
		Local $sVal = StringFormat("%.10g", $mX1[$sName])
		; wrap negatives in parentheses so operator precedence stays intact
		If StringLeft($sVal, 1) = "-" Then $sVal = "(" & $sVal & ")"
		$sOut = StringRegExpReplace($sOut, "(?i)\b" & $sName & "\b", $sVal)
	Next
	; Strip # prefix from observation placeholders (keep the variable name)
	$sOut = StringRegExpReplace($sOut, "#([A-Za-z]\w*)", "$1")
	Return $sOut
EndFunc   ;==>_ResSubstitutedFormula

Func _ImportCSV()
	Local $sFile = FileOpenDialog("CSV-Datei laden", @ScriptDir, _
		"CSV/Textdateien (*.csv;*.txt)|Alle Dateien (*.*)", 1, "", $g_hGUI)
	If @error Then Return
	$g_sLastImportedCSV = $sFile

	Local $sContent = FileRead($sFile)
	If $sContent = "" Then
		$g_sErrorMsg = "Datei ist leer"
		GUICtrlSetData($g_idLblError, $g_sErrorMsg)
		Return
	EndIf

	; Split into lines (handle CRLF, LF, CR)
	$sContent = StringReplace($sContent, @CRLF, @LF)
	$sContent = StringReplace($sContent, @CR, @LF)
	Local $aLines = StringSplit(StringStripWS($sContent, 2), @LF, 1)

	; Parse metadata comments (# formula=..., # depvar=..., # param_X=...)
	Local $iFirstDataLine = 1
	Local $mParamInit[]  ; collected parameter initial values
	For $i = 1 To $aLines[0]
		Local $sLine = StringStripWS($aLines[$i], 3)
		If StringLeft($sLine, 1) <> "#" Then ExitLoop
		$iFirstDataLine = $i + 1
		; Parse key=value from comment
		Local $aMatch = StringRegExp($sLine, "^#\s*(\w+)\s*=\s*(.+)$", 1)
		If Not @error Then
			Local $sKey = StringLower($aMatch[0])
			Switch $sKey
				Case "formula"
					GUICtrlSetData($g_idInputFormula, $aMatch[1])
				Case "depvar"
					GUICtrlSetData($g_idInputDepVar, $aMatch[1])
				Case Else
					; param_XXXX=value
					If StringLeft($sKey, 6) = "param_" Then
						$mParamInit[StringUpper(StringMid($aMatch[0], 7))] = StringStripWS($aMatch[1], 3)
					EndIf
			EndSwitch
		EndIf
	Next

	; Parse formula if metadata was found (sets up columns + parameter table)
	_ParseFormula()

	; Apply parameter initial values from file (after _ParseFormula created the table)
	For $i = 0 To UBound($g_aParamNames) - 1
		If MapExists($mParamInit, $g_aParamNames[$i]) Then
			__cg_setCell($g_iParamGrid, $i, 1, $mParamInit[$g_aParamNames[$i]])
		EndIf
	Next

	If UBound($g_aObsNames) = 0 Then
		$g_sErrorMsg = "Keine Formel — weder in Datei noch eingegeben"
		GUICtrlSetData($g_idLblError, $g_sErrorMsg)
		Return
	EndIf

	If $iFirstDataLine > $aLines[0] Or $aLines[0] - $iFirstDataLine < 1 Then
		$g_sErrorMsg = "CSV braucht Header + mind. 1 Datenzeile"
		GUICtrlSetData($g_idLblError, $g_sErrorMsg)
		Return
	EndIf

	; Detect delimiter: count semicolons vs commas in header
	Local $sHeader = $aLines[$iFirstDataLine]
	Local $iSemicolons = StringLen($sHeader) - StringLen(StringReplace($sHeader, ";", ""))
	Local $iCommas     = StringLen($sHeader) - StringLen(StringReplace($sHeader, ",", ""))
	Local $sDelim = ($iSemicolons >= $iCommas) ? ";" : ","
	Local $bCommaDecimal = ($sDelim = ";")

	; Parse header
	Local $aHeaders = StringSplit($sHeader, $sDelim, 1)

	; Map CSV columns to data columns
	Local $aColMap[$aHeaders[0]]
	Local $iMatched = 0
	For $i = 1 To $aHeaders[0]
		Local $sH = StringStripWS(StringUpper($aHeaders[$i]), 3)
		$sH = StringReplace($sH, '"', '')
		$sH = StringReplace($sH, "'", "")
		$aColMap[$i - 1] = ""
		If MapExists($g_mData, $sH) Then
			$aColMap[$i - 1] = $sH
			$iMatched += 1
		EndIf
	Next

	If $iMatched = 0 Then
		Local $sExpected = ""
		Local $aKeys = MapKeys($g_mData)
		For $i = 0 To UBound($aKeys) - 1
			$sExpected &= ($i > 0 ? ", " : "") & $aKeys[$i]
		Next
		$g_sErrorMsg = "Keine Spalten erkannt. Erwartet: " & $sExpected
		GUICtrlSetData($g_idLblError, $g_sErrorMsg)
		Return
	EndIf

	; Clear existing data
	_ClearAllData()

	; Parse data lines
	Local $iSkipped = 0
	For $iLine = $iFirstDataLine + 1 To $aLines[0]
		$sLine = StringStripWS($aLines[$iLine], 3)
		If $sLine = "" Then ContinueLoop

		Local $aFields = StringSplit($sLine, $sDelim, 1)

		; Add row
		$g_iRowCount += 1

		For $j = 1 To $aFields[0]
			If $j - 1 >= UBound($aColMap) Then ExitLoop
			If $aColMap[$j - 1] = "" Then ContinueLoop

			Local $sColName = $aColMap[$j - 1]
			Local $sFieldVal = StringStripWS($aFields[$j], 3)
			$sFieldVal = StringReplace($sFieldVal, '"', '')
			If $bCommaDecimal Then $sFieldVal = StringReplace($sFieldVal, ",", ".")

			; Validate numeric
			If $sFieldVal <> "" And Not __IsNumeric($sFieldVal) Then
				$sFieldVal = ""
				$iSkipped += 1
			EndIf

			If MapExists($g_mData, $sColName) Then
				Local $aCol = $g_mData[$sColName]
				If UBound($aCol) < $g_iRowCount Then ReDim $aCol[$g_iRowCount]
				$aCol[$g_iRowCount - 1] = $sFieldVal
				$g_mData[$sColName] = $aCol
			EndIf
		Next

		; Extend columns that had no data in this CSV row
		Local $aAllKeys = MapKeys($g_mData)
		For $j = 0 To UBound($aAllKeys) - 1
			Local $aCheck = $g_mData[$aAllKeys[$j]]
			If UBound($aCheck) < $g_iRowCount Then
				ReDim $aCheck[$g_iRowCount]
				$aCheck[$g_iRowCount - 1] = ""
				$g_mData[$aAllKeys[$j]] = $aCheck
			EndIf
		Next
	Next

	; Rebuild ListView display
	_RebuildDataColumns()
	_EnsureEmptyTrailingRow()

	If $iSkipped > 0 Then
		$g_sErrorMsg = $iSkipped & " ungueltige Werte uebersprungen"
	Else
		$g_sErrorMsg = ""
	EndIf
	GUICtrlSetData($g_idLblError, $g_sErrorMsg)

	_TryAutoSolve()
EndFunc   ;==>_ImportCSV

Func __IsNumeric($sVal)
	Return StringRegExp($sVal, "^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$")
EndFunc   ;==>__IsNumeric

; ===== Graph: Coordinate transforms =========================================

Func _DataToPixelX($fX)
	If $g_fXMax = $g_fXMin Then Return $g_iGraphPadL
	Return $g_iGraphPadL + ($fX - $g_fXMin) / ($g_fXMax - $g_fXMin) * _
		($g_iGraphW - $g_iGraphPadL - $g_iGraphPadR)
EndFunc   ;==>_DataToPixelX

Func _DataToPixelY($fY)
	If $g_fYMax = $g_fYMin Then Return $g_iGraphH - $g_iGraphPadB
	Return ($g_iGraphH - $g_iGraphPadB) - ($fY - $g_fYMin) / ($g_fYMax - $g_fYMin) * _
		($g_iGraphH - $g_iGraphPadT - $g_iGraphPadB)
EndFunc   ;==>_DataToPixelY

; ===== Graph: Main dispatcher ===============================================

Func _DrawGraph()
	If $g_hBuffer = 0 Then Return
	_GDIPlus_GraphicsClear($g_hBuffer, $g_iClrBg)

	; Read axis range from inputs
	$g_fXMin = Number(GUICtrlRead($g_idInputXMin))
	$g_fXMax = Number(GUICtrlRead($g_idInputXMax))
	$g_fYMin = Number(GUICtrlRead($g_idInputYMin))
	$g_fYMax = Number(GUICtrlRead($g_idInputYMax))
	If $g_fXMax <= $g_fXMin Then $g_fXMax = $g_fXMin + 1
	If $g_fYMax <= $g_fYMin Then $g_fYMax = $g_fYMin + 1

	_DrawGridAndAxes()

	If $g_bAdjOK Then
		If $g_bShowConfBand Then _DrawConfidenceBand()
		_DrawFittedCurve()
		If $g_bShowResiduals Then _DrawResiduals()
	EndIf

	_DrawDataPoints()

	; Blit to screen
	_GDIPlus_GraphicsDrawImageRect($g_hGraphics, $g_hBitmap, _
		$g_iGraphX, $g_iGraphY, $g_iGraphW, $g_iGraphH)
EndFunc   ;==>_DrawGraph

; ===== Graph: Grid and axes =================================================

Func _DrawGridAndAxes()
	Local $iPlotL = $g_iGraphPadL
	Local $iPlotR = $g_iGraphW - $g_iGraphPadR
	Local $iPlotT = $g_iGraphPadT
	Local $iPlotB = $g_iGraphH - $g_iGraphPadB

	Local $fXStep = __NiceInterval($g_fXMax - $g_fXMin)
	Local $fYStep = __NiceInterval($g_fYMax - $g_fYMin)

	; Vertical grid lines + X tick labels
	Local $fTick = Ceiling($g_fXMin / $fXStep) * $fXStep
	While $fTick <= $g_fXMax
		Local $fPx = _DataToPixelX($fTick)
		_GDIPlus_GraphicsDrawLine($g_hBuffer, $fPx, $iPlotT, $fPx, $iPlotB, $g_hPenGrid)
		Local $tRect = _GDIPlus_RectFCreate($fPx - 25, $iPlotB + 3, 50, 16)
		_GDIPlus_GraphicsDrawStringEx($g_hBuffer, StringFormat("%.4g", $fTick), _
			$g_hFontAxis, $tRect, $g_hFormatCenter, $g_hBrushText)
		$fTick += $fXStep
	WEnd

	; Horizontal grid lines + Y tick labels
	$fTick = Ceiling($g_fYMin / $fYStep) * $fYStep
	While $fTick <= $g_fYMax
		Local $fPy = _DataToPixelY($fTick)
		_GDIPlus_GraphicsDrawLine($g_hBuffer, $iPlotL, $fPy, $iPlotR, $fPy, $g_hPenGrid)
		Local $tRect2 = _GDIPlus_RectFCreate(2, $fPy - 8, $iPlotL - 5, 16)
		_GDIPlus_GraphicsDrawStringEx($g_hBuffer, StringFormat("%.4g", $fTick), _
			$g_hFontAxis, $tRect2, $g_hFormatRight, $g_hBrushText)
		$fTick += $fYStep
	WEnd

	; Axes
	_GDIPlus_GraphicsDrawLine($g_hBuffer, $iPlotL, $iPlotB, $iPlotR, $iPlotB, $g_hPenAxis)
	_GDIPlus_GraphicsDrawLine($g_hBuffer, $iPlotL, $iPlotT, $iPlotL, $iPlotB, $g_hPenAxis)

	; Axis labels
	Local $sXLabel = (UBound($g_aObsNames) > 0) ? $g_aObsNames[0] : "X"
	Local $sYLabel = ($g_sDepVar <> "") ? $g_sDepVar : "Y"
	Local $tXL = _GDIPlus_RectFCreate($iPlotL, $iPlotB + 18, $iPlotR - $iPlotL, 16)
	_GDIPlus_GraphicsDrawStringEx($g_hBuffer, $sXLabel, $g_hFontAxis, $tXL, $g_hFormatCenter, $g_hBrushText)
	Local $tYL = _GDIPlus_RectFCreate(2, $iPlotT - 16, $iPlotL - 5, 16)
	_GDIPlus_GraphicsDrawStringEx($g_hBuffer, $sYLabel, $g_hFontAxis, $tYL, $g_hFormatRight, $g_hBrushText)
EndFunc   ;==>_DrawGridAndAxes

Func __NiceInterval($fRange)
	If $fRange <= 0 Then Return 1
	Local $fRough = $fRange / 6
	Local $fPow = 10 ^ Floor(Log($fRough) / Log(10))
	Local $fNorm = $fRough / $fPow
	Local $fNice
	If $fNorm < 1.5 Then
		$fNice = 1
	ElseIf $fNorm < 3 Then
		$fNice = 2
	ElseIf $fNorm < 7 Then
		$fNice = 5
	Else
		$fNice = 10
	EndIf
	Return $fNice * $fPow
EndFunc   ;==>__NiceInterval

; ===== Graph: Data points ===================================================

Func _DrawDataPoints()
	If UBound($g_aObsNames) = 0 Then Return

	Local $sXName = $g_aObsNames[0]
	If Not MapExists($g_mData, $sXName) Or Not MapExists($g_mData, $g_sDepVar) Then Return

	Local $aXData = $g_mData[$sXName]
	Local $aYData = $g_mData[$g_sDepVar]
	Local Const $iR = 4

	For $i = 0 To $g_iRowCount - 1
		If $i >= UBound($aXData) Or $i >= UBound($aYData) Then ExitLoop
		If $aXData[$i] = "" Or $aYData[$i] = "" Then ContinueLoop
		Local $fPx = _DataToPixelX(Number($aXData[$i]))
		Local $fPy = _DataToPixelY(Number($aYData[$i]))

		Local $bOutlier = ($g_bAdjOK And $i < UBound($g_aRobWeights) And $g_aRobWeights[$i] < 0.5)

		_GDIPlus_GraphicsFillEllipse($g_hBuffer, $fPx - $iR, $fPy - $iR, $iR * 2, $iR * 2, _
			$bOutlier ? $g_hBrushOutlier : $g_hBrushPoint)
	Next
EndFunc   ;==>_DrawDataPoints

; ===== Graph: Fitted curve ==================================================

Func _DrawFittedCurve()
	If UBound($g_aObsNames) = 0 Then Return
	If UBound($g_aObsNames) <> 1 Then Return ; only for 2D

	; Kurve wird immer mit den AKTUELLEN Parameter-Werten aus dem Grid gezeichnet —
	; auch nach Parameter-Edit (vor Re-Solve) für visuelles Preview-Feedback
	Local $mX = _ReadParamValuesFromGrid()
	Local $sObsName = $g_aObsNames[0]
	Local Const $iSteps = 200

	Local $sEvalBase = __PrepareFormulaForEval($g_sFormula, $mX)

	Local $fStep = ($g_fXMax - $g_fXMin) / $iSteps
	Local $fPrevPx = 0, $fPrevPy = 0, $bFirst = True

	For $j = 0 To $iSteps
		Local $fXVal = $g_fXMin + $j * $fStep
		Local $sEval = StringReplace($sEvalBase, "#" & $sObsName, String($fXVal))
		Local $fYVal = Execute($sEval)
		If @error Or Not IsNumber($fYVal) Then
			$bFirst = True
			ContinueLoop
		EndIf

		Local $fPx = _DataToPixelX($fXVal)
		Local $fPy = _DataToPixelY($fYVal)

		If Not $bFirst Then
			If $fPrevPy > -500 And $fPrevPy < $g_iGraphH + 500 And _
			   $fPy > -500 And $fPy < $g_iGraphH + 500 Then
				_GDIPlus_GraphicsDrawLine($g_hBuffer, $fPrevPx, $fPrevPy, $fPx, $fPy, $g_hPenCurve)
			EndIf
		EndIf
		$fPrevPx = $fPx
		$fPrevPy = $fPy
		$bFirst = False
	Next
EndFunc   ;==>_DrawFittedCurve

Func __PrepareFormulaForEval($sFormula, $mParamValues)
	Local $sEval = $sFormula

	; Sort param names by length descending (prevent partial matches)
	Local $aNames[UBound($g_aParamNames)]
	For $i = 0 To UBound($g_aParamNames) - 1
		$aNames[$i] = $g_aParamNames[$i]
	Next
	For $i = 0 To UBound($aNames) - 2
		For $j = $i + 1 To UBound($aNames) - 1
			If StringLen($aNames[$j]) > StringLen($aNames[$i]) Then
				Local $sTmp = $aNames[$i]
				$aNames[$i] = $aNames[$j]
				$aNames[$j] = $sTmp
			EndIf
		Next
	Next

	; Replace parameters with values using word boundary regex
	For $i = 0 To UBound($aNames) - 1
		Local $sName = $aNames[$i]
		If IsMap($mParamValues) And MapExists($mParamValues, $sName) Then
			$sEval = StringRegExpReplace($sEval, "(?i)\b" & $sName & "\b", _
				"(" & String($mParamValues[$sName]) & ")")
		EndIf
	Next

	; Translate function names to AutoIt equivalents
	$sEval = StringRegExpReplace($sEval, "(?i)\bLN\(", "Log(")
	$sEval = StringRegExpReplace($sEval, "(?i)\bLOG10\(([^)]+)\)", "(Log($1)/Log(10))")
	$sEval = StringRegExpReplace($sEval, "(?i)\bPI\(\)", String(ACos(-1)))

	Return $sEval
EndFunc   ;==>__PrepareFormulaForEval

; ===== Graph: Auto-fit axes =================================================

Func _AutoFitAxes()
	$g_bAutoAxis = True

	If UBound($g_aObsNames) = 0 Then Return

	Local $sXName = $g_aObsNames[0]
	If Not MapExists($g_mData, $sXName) Or Not MapExists($g_mData, $g_sDepVar) Then Return

	Local $aXData = $g_mData[$sXName]
	Local $aYData = $g_mData[$g_sDepVar]

	Local $fMinX = 1e30, $fMaxX = -1e30
	Local $fMinY = 1e30, $fMaxY = -1e30
	Local $iValid = 0

	For $i = 0 To $g_iRowCount - 1
		If $i >= UBound($aXData) Or $i >= UBound($aYData) Then ExitLoop
		If $aXData[$i] = "" Or $aYData[$i] = "" Then ContinueLoop
		Local $fX = Number($aXData[$i])
		Local $fY = Number($aYData[$i])
		If $fX < $fMinX Then $fMinX = $fX
		If $fX > $fMaxX Then $fMaxX = $fX
		If $fY < $fMinY Then $fMinY = $fY
		If $fY > $fMaxY Then $fMaxY = $fY
		$iValid += 1
	Next

	If $iValid < 1 Then Return

	Local $fPadX = ($fMaxX - $fMinX) * 0.1
	Local $fPadY = ($fMaxY - $fMinY) * 0.1
	If $fPadX = 0 Then $fPadX = 1
	If $fPadY = 0 Then $fPadY = 1

	$g_fXMin = $fMinX - $fPadX
	$g_fXMax = $fMaxX + $fPadX
	$g_fYMin = $fMinY - $fPadY
	$g_fYMax = $fMaxY + $fPadY

	GUICtrlSetData($g_idInputXMin, StringFormat("%.4g", $g_fXMin))
	GUICtrlSetData($g_idInputXMax, StringFormat("%.4g", $g_fXMax))
	GUICtrlSetData($g_idInputYMin, StringFormat("%.4g", $g_fYMin))
	GUICtrlSetData($g_idInputYMax, StringFormat("%.4g", $g_fYMax))

	_DrawGraph()
EndFunc   ;==>_AutoFitAxes

; ===== Graph: Residuals ======================================================

Func _DrawResiduals()
	If UBound($g_aObsNames) <> 1 Then Return ; only for 2D

	Local $sXName = $g_aObsNames[0]
	If Not MapExists($g_mData, $sXName) Or Not MapExists($g_mData, $g_sDepVar) Then Return

	Local $aXData = $g_mData[$sXName]
	Local $aYData = $g_mData[$g_sDepVar]
	Local $mX = $g_mResults.x1

	; For GLM models (Ortho/Deming/York): use adjusted observations for target point
	Local $bUseAdjusted = ($g_iRegType <> $REG_NORMAL)
	Local $mObsAdj = Null
	If $bUseAdjusted And MapExists($g_mResults, "obsAdj") Then
		$mObsAdj = $g_mResults.obsAdj
	Else
		$bUseAdjusted = False
	EndIf

	Local $sEvalBase = __PrepareFormulaForEval($g_sFormula, $mX)

	For $i = 0 To $g_iRowCount - 1
		If $i >= UBound($aXData) Or $i >= UBound($aYData) Then ExitLoop
		If $aXData[$i] = "" Or $aYData[$i] = "" Then ContinueLoop

		Local $fXObs = Number($aXData[$i])
		Local $fYObs = Number($aYData[$i])
		Local $sSuffix = String($i + 1)

		Local $fXTarget, $fYTarget

		If $bUseAdjusted Then
			; Get adjusted observations (the point ON the fitted curve)
			Local $sXKey = $sXName & $sSuffix
			Local $sYKey = $g_sDepVar & $sSuffix
			If IsMap($mObsAdj) And MapExists($mObsAdj, $sXKey) And MapExists($mObsAdj, $sYKey) Then
				$fXTarget = $mObsAdj[$sXKey]
				$fYTarget = $mObsAdj[$sYKey]
			Else
				; Fallback: vertical residual
				$fXTarget = $fXObs
				Local $sEval = StringReplace($sEvalBase, "#" & $sXName, String($fXObs))
				$fYTarget = Execute($sEval)
				If Not IsNumber($fYTarget) Then ContinueLoop
			EndIf
		Else
			; Normal regression: vertical residual (evaluate curve at observed X)
			$fXTarget = $fXObs
			Local $sEval2 = StringReplace($sEvalBase, "#" & $sXName, String($fXObs))
			$fYTarget = Execute($sEval2)
			If Not IsNumber($fYTarget) Then ContinueLoop
		EndIf

		; Draw residual line from observed point to target point (on curve)
		Local $fPxObs = _DataToPixelX($fXObs)
		Local $fPyObs = _DataToPixelY($fYObs)
		Local $fPxTarget = _DataToPixelX($fXTarget)
		Local $fPyTarget = _DataToPixelY($fYTarget)
		_GDIPlus_GraphicsDrawLine($g_hBuffer, $fPxObs, $fPyObs, $fPxTarget, $fPyTarget, $g_hPenResid)

		; Draw residual label at midpoint
		Local $fResidual = Sqrt(($fXObs - $fXTarget) ^ 2 + ($fYObs - $fYTarget) ^ 2)
		; Sign: positive if observed Y > target Y
		If $fYObs < $fYTarget Then $fResidual = -$fResidual
		Local $sLabel = StringFormat("%.2g", $fResidual)
		Local $fMidPx = ($fPxObs + $fPxTarget) / 2
		Local $fMidPy = ($fPyObs + $fPyTarget) / 2
		Local $tLabelRect = _GDIPlus_RectFCreate($fMidPx + 5, $fMidPy - 7, 55, 14)
		_GDIPlus_GraphicsDrawStringEx($g_hBuffer, $sLabel, _
			$g_hFontAxis, $tLabelRect, $g_hFormatCenter, $g_hBrushResLabel)
	Next
EndFunc   ;==>_DrawResiduals

; ===== Graph: Confidence band ================================================

Func _DrawConfidenceBand()
	If Not IsMap($g_mResults) Or UBound($g_aObsNames) <> 1 Then Return
	If $g_mResults.f <= 0 Then Return

	Local $mQxx = $g_mResults.Qxx
	If Not IsMap($mQxx) Then Return

	Local $mX = $g_mResults.x1
	Local $fS0 = $g_mResults.s0
	Local $sObsName = $g_aObsNames[0]
	Local $sEvalBase = __PrepareFormulaForEval($g_sFormula, $mX)

	; t-factor (approximate)
	Local $fT = ($g_mResults.f > 30) ? 1.96 : (2.0 + 4.0 / $g_mResults.f)

	Local $iNp = UBound($g_aParamNames)
	Local Const $iSteps = 100

	; Note: Assumes Qxx parameter ordering matches $g_aParamNames order.
	; This holds when parameters are discovered in left-to-right formula order,
	; which is the case for both our parser and the UDF's internal parser.
	; For exact indexing, $mSystem.state.idxParams would be needed (cf. InteractiveCircleFit.au3).
	Local $iQRows = $mQxx.rows
	Local $aQxx[$iNp][$iNp]
	For $i = 0 To $iNp - 1
		For $j = 0 To $iNp - 1
			$aQxx[$i][$j] = DllStructGetData($mQxx.struct, 1, $i + $j * $iQRows + 1)
		Next
	Next

	; Draw band as vertical line segments
	Local $hPenBand = _GDIPlus_PenCreate($g_iClrConf, 2)
	Local $fDx = ($g_fXMax - $g_fXMin) / $iSteps
	Local Const $fH = 1e-6

	For $k = 0 To $iSteps
		Local $fXVal = $g_fXMin + $k * $fDx

		; Compute fitted Y
		Local $sEval = StringReplace($sEvalBase, "#" & $sObsName, String($fXVal))
		Local $fYFit = Execute($sEval)
		If Not IsNumber($fYFit) Then ContinueLoop

		; Compute Jacobian: df/dp_i by central differences
		Local $aJac[$iNp]
		For $i = 0 To $iNp - 1
			Local $mXPlus = __CloneMap($mX)
			Local $mXMinus = __CloneMap($mX)
			$mXPlus[$g_aParamNames[$i]] = $mX[$g_aParamNames[$i]] + $fH
			$mXMinus[$g_aParamNames[$i]] = $mX[$g_aParamNames[$i]] - $fH

			Local $sEvalP = StringReplace(__PrepareFormulaForEval($g_sFormula, $mXPlus), _
				"#" & $sObsName, String($fXVal))
			Local $sEvalM = StringReplace(__PrepareFormulaForEval($g_sFormula, $mXMinus), _
				"#" & $sObsName, String($fXVal))
			Local $fFp = Execute($sEvalP)
			Local $fFm = Execute($sEvalM)
			$aJac[$i] = (IsNumber($fFp) And IsNumber($fFm)) ? ($fFp - $fFm) / (2 * $fH) : 0
		Next

		; sigma^2 = J * Qxx * J^T
		Local $fVar = 0
		For $i = 0 To $iNp - 1
			For $j = 0 To $iNp - 1
				$fVar += $aJac[$i] * $aQxx[$i][$j] * $aJac[$j]
			Next
		Next
		Local $fSigmaY = $fS0 * Sqrt(Abs($fVar))
		Local $fBand = $fT * $fSigmaY

		; Draw vertical band segment
		Local $fPxX = _DataToPixelX($fXVal)
		Local $fPyUp = _DataToPixelY($fYFit + $fBand)
		Local $fPyLo = _DataToPixelY($fYFit - $fBand)

		If $fPyUp > -500 And $fPyLo < $g_iGraphH + 500 Then
			_GDIPlus_GraphicsDrawLine($g_hBuffer, $fPxX, $fPyUp, $fPxX, $fPyLo, $hPenBand)
		EndIf
	Next

	_GDIPlus_PenDispose($hPenBand)
EndFunc   ;==>_DrawConfidenceBand

Func __CloneMap($mSrc)
	Local $mDst[]
	Local $aKeys = MapKeys($mSrc)
	For $i = 0 To UBound($aKeys) - 1
		$mDst[$aKeys[$i]] = $mSrc[$aKeys[$i]]
	Next
	Return $mDst
EndFunc   ;==>__CloneMap

; ===== Auto-solve ===========================================================

Func _TryAutoSolve()
	; Guard: formula must be set
	If $g_sFormula = "" Or UBound($g_aObsNames) = 0 Then Return
	; Guard: need at least n_params valid data rows
	If _CountValidRows() < UBound($g_aParamNames) Then Return
	_BuildAndSolve()
EndFunc   ;==>_TryAutoSolve

Func _CountValidRows()
	If $g_iRowCount = 0 Then Return 0
	Local $iValid = 0
	For $i = 0 To $g_iRowCount - 1
		Local $bComplete = True
		; Check all observation columns
		For $j = 0 To UBound($g_aObsNames) - 1
			If Not MapExists($g_mData, $g_aObsNames[$j]) Then
				$bComplete = False
				ExitLoop
			EndIf
			Local $aCol = $g_mData[$g_aObsNames[$j]]
			If $i >= UBound($aCol) Or $aCol[$i] = "" Then
				$bComplete = False
				ExitLoop
			EndIf
		Next
		; Check dependent variable
		If $bComplete And MapExists($g_mData, $g_sDepVar) Then
			Local $aDepCol = $g_mData[$g_sDepVar]
			If $i >= UBound($aDepCol) Or $aDepCol[$i] = "" Then $bComplete = False
		Else
			$bComplete = False
		EndIf
		If $bComplete Then $iValid += 1
	Next
	Return $iValid
EndFunc   ;==>_CountValidRows

; ===== Helpers ==============================================================

; Load a nonlinear startup example: exponential growth Y = A * exp(B * X)
Func _LoadStartupExample()
	; Logistic growth: Y = A / (1 + exp(-B*(X-C)))
	; X ∈ [0, 10], Y ∈ [0, 10] — similar scales → orthogonal residuals visibly diagonal
	; Nonlinear, 3 parameters, converges reliably with LM
	GUICtrlSetData($g_idInputFormula, "A / (1 + EXP(-B * (#X - C)))")
	GUICtrlSetData($g_idInputDepVar, "Y")
	_ParseFormula()

	; Noisy logistic data (true: A≈10, B≈0.8, C≈5, with ±0.5..1.0 noise)
	Local $aX[] = [0.0, 1.0, 2.0, 3.0, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0]
	Local $aY[] = [0.8, 0.2, 1.8, 1.2, 4.1, 3.5, 5.8, 5.5, 7.8, 7.9, 9.8, 9.1, 10.5]

	; Populate data table
	For $i = 0 To UBound($aX) - 1
		If $i >= $g_iRowCount Then _AddRow()
		Local $aColX = $g_mData["X"]
		$aColX[$i] = String($aX[$i])
		$g_mData["X"] = $aColX
		Local $aColY = $g_mData["Y"]
		$aColY[$i] = String($aY[$i])
		$g_mData["Y"] = $aColY
	Next

	; Set initial values (logistic converges well from these)
	_RebuildDataColumns()
	__cg_setCell($g_iParamGrid, 0, 1, "10.0") ; A ≈ 10 (saturation)
	__cg_setCell($g_iParamGrid, 1, 1, "1.0")  ; B ≈ 0.8 (steepness)
	__cg_setCell($g_iParamGrid, 2, 1, "5.0")  ; C ≈ 5 (midpoint)

	; Auto-fit and solve
	$g_bAutoAxis = True
	_TryAutoSolve()
EndFunc   ;==>_LoadStartupExample

; Compute R² (coefficient of determination) from current results
Func __ComputeR2()
	If Not IsMap($g_mResults) Or UBound($g_aObsNames) = 0 Then Return -1e31
	If Not MapExists($g_mData, $g_sDepVar) Then Return -1e31

	Local $aYData = $g_mData[$g_sDepVar]

	; Compute mean of Y
	Local $fSum = 0, $iN = 0
	For $i = 0 To $g_iRowCount - 1
		If $i >= UBound($aYData) Or $aYData[$i] = "" Then ContinueLoop
		$fSum += Number($aYData[$i])
		$iN += 1
	Next
	If $iN < 2 Then Return -1e31
	Local $fMean = $fSum / $iN

	; SS_tot = sum((y_i - mean)^2), SS_res = sum(v_i^2)
	; For Normal regression: v = y_obs - y_fit, so SS_res = sum(v^2)
	; For GLM: both X and Y are adjusted, R² is less meaningful but still informative
	Local $fSStot = 0, $fSSres = 0
	Local $mV = Null
	If MapExists($g_mResults, "v") Then $mV = $g_mResults.v

	For $i = 0 To $g_iRowCount - 1
		If $i >= UBound($aYData) Or $aYData[$i] = "" Then ContinueLoop
		Local $fY = Number($aYData[$i])
		$fSStot += ($fY - $fMean) ^ 2

		; Get residual for dependent variable
		Local $sKey = $g_sDepVar & String($i + 1)
		If IsMap($mV) And MapExists($mV, $sKey) Then
			$fSSres += $mV[$sKey] ^ 2
		EndIf
	Next

	If $fSStot = 0 Then Return -1e31
	Return 1 - $fSSres / $fSStot
EndFunc   ;==>__ComputeR2

Func _RebuildSdFields()
	Local Const $iY = 38

	; York: hide pairs, show hint
	If $g_iRegType = $REG_YORK Then
		For $i = 0 To $g_iMaxSdFields - 1
			__HideCtrl($g_aidSdLabels[$i])
			__HideCtrl($g_aidSdInputs[$i])
			$g_asSdNames[$i] = ""
		Next
		GUICtrlSetPos($g_idSdHintLabel, 10, $iY + 3, 220, 20)
		__ShowCtrl($g_idSdHintLabel)
		$g_iSdFieldCount = 0
		Return
	EndIf

	__HideCtrl($g_idSdHintLabel)

	; Determine which variables need σ fields
	Local $aVarNames[0]
	Switch $g_iRegType
		Case $REG_NORMAL
			If $g_sDepVar <> "" Then
				ReDim $aVarNames[1]
				$aVarNames[0] = $g_sDepVar
			EndIf
		Case $REG_ORTHO
			ReDim $aVarNames[1]
			$aVarNames[0] = ""
		Case $REG_DEMING
			ReDim $aVarNames[UBound($g_aObsNames) + 1]
			For $i = 0 To UBound($g_aObsNames) - 1
				$aVarNames[$i] = $g_aObsNames[$i]
			Next
			$aVarNames[UBound($g_aObsNames)] = $g_sDepVar
	EndSwitch

	Local $iCount = UBound($aVarNames)
	If $iCount > $g_iMaxSdFields Then $iCount = $g_iMaxSdFields
	$g_iSdFieldCount = $iCount

	Local $iX = 10
	For $i = 0 To $iCount - 1
		$g_asSdNames[$i] = $aVarNames[$i]

		Local $sLabel
		If $aVarNames[$i] = "" Then
			$sLabel = ChrW(963) & ":"
		Else
			$sLabel = ChrW(963) & "_" & $aVarNames[$i] & ":"
		EndIf

		Local $iLblW = 10 + StringLen($sLabel) * 7
		If $iLblW < 25 Then $iLblW = 25

		GUICtrlSetData($g_aidSdLabels[$i], $sLabel)
		GUICtrlSetPos($g_aidSdLabels[$i], $iX, $iY + 3, $iLblW, 20)
		__ShowCtrl($g_aidSdLabels[$i])
		$iX += $iLblW + 2

		Local $sKey = ($aVarNames[$i] = "") ? "__shared__" : $aVarNames[$i]
		Local $sVal = "1.0"
		If MapExists($g_mSdValues, $sKey) Then $sVal = String($g_mSdValues[$sKey])
		GUICtrlSetData($g_aidSdInputs[$i], $sVal)
		GUICtrlSetPos($g_aidSdInputs[$i], $iX, $iY, 50, 22)
		__ShowCtrl($g_aidSdInputs[$i])
		$iX += 55
	Next

	For $i = $iCount To $g_iMaxSdFields - 1
		__HideCtrl($g_aidSdLabels[$i])
		__HideCtrl($g_aidSdInputs[$i])
		$g_asSdNames[$i] = ""
	Next
EndFunc   ;==>_RebuildSdFields

Func __HideCtrl($iCtrlId)
	Local $hCtrl = GUICtrlGetHandle($iCtrlId)
	If $hCtrl = 0 Then Return
	DllCall("user32.dll", "bool", "ShowWindow", "hwnd", $hCtrl, "int", 0) ; SW_HIDE
EndFunc   ;==>__HideCtrl

Func __ShowCtrl($iCtrlId)
	Local $hCtrl = GUICtrlGetHandle($iCtrlId)
	If $hCtrl = 0 Then Return
	DllCall("user32.dll", "bool", "ShowWindow", "hwnd", $hCtrl, "int", 8) ; SW_SHOWNA
EndFunc   ;==>__ShowCtrl

; Predict Y for a given X using fitted parameters + variance propagation for σ_ŷ
Func _UpdatePrediction()
	Local $sInput = GUICtrlRead($g_idInputPredX)
	If $sInput = "" Then
		GUICtrlSetData($g_idLblPredResult, "")
		Return
	EndIf
	If UBound($g_aObsNames) = 0 Then Return

	; Parameter aus Grid lesen — liefert auch bei Preview-Änderungen aktuelle Werte
	Local $mX = _ReadParamValuesFromGrid()
	Local $sEvalBase = __PrepareFormulaForEval($g_sFormula, $mX)

	; Parse input — for multi-dimensional, expect comma-separated values
	Local $aInputVals = StringSplit(StringStripWS($sInput, 3), ",", 2) ; flag 2 = no count element
	If UBound($aInputVals) <> UBound($g_aObsNames) Then
		Local $sExpected = ""
		For $i = 0 To UBound($g_aObsNames) - 1
			$sExpected &= ($i > 0 ? ", " : "") & $g_aObsNames[$i]
		Next
		GUICtrlSetData($g_idLblPredResult, "Eingabe: " & $sExpected)
		Return
	EndIf

	; Substitute observation values into formula
	Local $sEval = $sEvalBase
	For $i = 0 To UBound($g_aObsNames) - 1
		Local $fVal = Number(StringStripWS($aInputVals[$i], 3))
		$sEval = StringReplace($sEval, "#" & $g_aObsNames[$i], String($fVal))
	Next

	; Evaluate function
	Local $fYPred = Execute($sEval)
	If Not IsNumber($fYPred) Then
		GUICtrlSetData($g_idLblPredResult, "Fehler bei Auswertung")
		Return
	EndIf

	; Compute σ_ŷ via variance propagation: σ² = s₀² · J · Qxx · Jᵀ
	; Nur wenn eine gültige Ausgleichung vorliegt (Preview-Modus liefert keine σ).
	Local $sSigma = ""
	Local $mQxx = Null
	If $g_bAdjOK And IsMap($g_mResults) Then $mQxx = $g_mResults.Qxx
	If IsMap($mQxx) And $g_mResults.f > 0 Then
		Local $iNp = UBound($g_aParamNames)
		Local $fS0 = $g_mResults.s0
		Local $iQRows = $mQxx.rows
		Local Const $fH = 1e-6

		; Numerical Jacobian
		Local $aJac[$iNp]
		For $i = 0 To $iNp - 1
			Local $mXPlus = __CloneMap($mX)
			Local $mXMinus = __CloneMap($mX)
			$mXPlus[$g_aParamNames[$i]] = $mX[$g_aParamNames[$i]] + $fH
			$mXMinus[$g_aParamNames[$i]] = $mX[$g_aParamNames[$i]] - $fH

			Local $sEvalP = __PrepareFormulaForEval($g_sFormula, $mXPlus)
			Local $sEvalM = __PrepareFormulaForEval($g_sFormula, $mXMinus)
			For $j = 0 To UBound($g_aObsNames) - 1
				$sEvalP = StringReplace($sEvalP, "#" & $g_aObsNames[$j], String(Number(StringStripWS($aInputVals[$j], 3))))
				$sEvalM = StringReplace($sEvalM, "#" & $g_aObsNames[$j], String(Number(StringStripWS($aInputVals[$j], 3))))
			Next
			Local $fFp = Execute($sEvalP)
			Local $fFm = Execute($sEvalM)
			$aJac[$i] = (IsNumber($fFp) And IsNumber($fFm)) ? ($fFp - $fFm) / (2 * $fH) : 0
		Next

		; σ² = J · Qxx · Jᵀ (read Qxx directly)
		Local $fVar = 0
		For $i = 0 To $iNp - 1
			For $j = 0 To $iNp - 1
				$fVar += $aJac[$i] * DllStructGetData($mQxx.struct, 1, $i + $j * $iQRows + 1) * $aJac[$j]
			Next
		Next
		Local $fSigmaY = $fS0 * Sqrt(Abs($fVar))
		$sSigma = " " & ChrW(177) & " " & StringFormat("%.6g", $fSigmaY) ; ±
	EndIf

	; Display result
	GUICtrlSetData($g_idLblPredResult, _
		$g_sDepVar & " = " & StringFormat("%.6g", $fYPred) & $sSigma)
EndFunc   ;==>_UpdatePrediction

; ===== Cleanup ==============================================================

Func _Cleanup()
	_GDIPlus_PenDispose($g_hPenCurve)
	_GDIPlus_PenDispose($g_hPenAxis)
	_GDIPlus_PenDispose($g_hPenGrid)
	_GDIPlus_PenDispose($g_hPenResid)
	_GDIPlus_BrushDispose($g_hBrushPoint)
	_GDIPlus_BrushDispose($g_hBrushConf)
	_GDIPlus_BrushDispose($g_hBrushText)
	_GDIPlus_BrushDispose($g_hBrushOutlier)
	_GDIPlus_BrushDispose($g_hBrushResLabel)
	_GDIPlus_FontDispose($g_hFontAxis)
	_GDIPlus_FontFamilyDispose($g_hFontFamily)
	_GDIPlus_StringFormatDispose($g_hFormatCenter)
	_GDIPlus_StringFormatDispose($g_hFormatRight)
	_GDIPlus_GraphicsDispose($g_hBuffer)
	_GDIPlus_BitmapDispose($g_hBitmap)
	_GDIPlus_GraphicsDispose($g_hGraphics)
	_GDIPlus_Shutdown()
	GUIDelete($g_hGUI)
EndFunc   ;==>_Cleanup

Func _WM_NOTIFY($hWnd, $iMsg, $wParam, $lParam)
	#forceref $hWnd, $iMsg, $wParam, $lParam
	Return $GUI_RUNDEFMSG
EndFunc   ;==>_WM_NOTIFY

Func _WM_COMMAND($hWnd, $iMsg, $wParam, $lParam)
	#forceref $hWnd, $iMsg, $lParam
	; IMPORTANT: NEVER call functions that create/delete GUI controls here.
	; Only set flags — the main loop processes them.

	Local $iCtrlId = BitAND($wParam, 0xFFFF)  ; LOWORD
	Local $iCode = BitShift($wParam, 16)       ; HIWORD

	; Formula/depvar inputs lost focus → defer parse
	If ($iCtrlId = $g_idInputFormula Or $iCtrlId = $g_idInputDepVar) And $iCode = $EN_KILLFOCUS Then
		$g_bPendingParse = True
	EndIf

	; σ inputs lost focus → defer re-solve (only check active slots)
	For $i = 0 To $g_iSdFieldCount - 1
		If $iCtrlId = $g_aidSdInputs[$i] And $iCode = $EN_KILLFOCUS Then
			$g_bPendingSolve = True
			ExitLoop
		EndIf
	Next

	; Prediction input changed → defer update
	If $iCtrlId = $g_idInputPredX And ($iCode = $EN_KILLFOCUS Or $iCode = $EN_CHANGE) Then
		$g_bPendingPrediction = True
	EndIf

	Return $GUI_RUNDEFMSG
EndFunc   ;==>_WM_COMMAND

Func _WM_PAINT($hWnd, $iMsg, $wParam, $lParam)
	#forceref $hWnd, $iMsg, $wParam, $lParam
	If $hWnd = $g_hGUI And $g_hBitmap <> 0 Then
		_GDIPlus_GraphicsDrawImageRect($g_hGraphics, $g_hBitmap, _
			$g_iGraphX, $g_iGraphY, $g_iGraphW, $g_iGraphH)
	EndIf
	Return $GUI_RUNDEFMSG
EndFunc   ;==>_WM_PAINT

Func _WM_ERASEBKGND($hWnd, $iMsg, $wParam, $lParam)
	#forceref $iMsg, $lParam
	If $hWnd <> $g_hGUI Then Return $GUI_RUNDEFMSG

	; Paint only the non-graph part of the invalidated clip rect, then claim
	; "fully erased" with Return 1. The graph area is left untouched so that
	; WM_PAINT's GDI+ bitmap overlay shows up without a grey flicker frame.
	;
	; Safety margin of 2 px above $g_iGraphY: the graph bitmap includes a
	; thin top border that visually extends slightly above the exact GraphY
	; line (anti-aliasing / sub-pixel rounding), so we stop the fill a couple
	; of pixels earlier to avoid touching that border.
	;
	; The original buggy version filled a FIXED rect (0,0,WinW,GraphY) and
	; returned 1 regardless of the actual clip. When Windows invalidated a
	; small control rect (e.g. on σ-input hide), FillRect was clipped down
	; to that sliver, yet Return 1 claimed "fully erased" — so neighboring
	; ghost pixels from the deleted control never got cleared. Using the
	; real clip box avoids that.
	Local Const $iGraphTopMargin = 2
	Local $iEraseBottom = $g_iGraphY - $iGraphTopMargin

	Local $tClip = DllStructCreate("long Left;long Top;long Right;long Bottom")
	Local $aClip = DllCall("gdi32.dll", "int", "GetClipBox", "handle", $wParam, "struct*", $tClip)
	If Not IsArray($aClip) Then Return $GUI_RUNDEFMSG

	Local $iClipTop    = DllStructGetData($tClip, "Top")
	Local $iClipBottom = DllStructGetData($tClip, "Bottom")

	If $iClipTop >= $iEraseBottom Then Return 1 ; clip entirely inside graph area
	If $iClipBottom > $iEraseBottom Then DllStructSetData($tClip, "Bottom", $iEraseBottom)

	Local $aBrush = DllCall("user32.dll", "handle", "GetSysColorBrush", "int", 15)
	If Not @error And IsArray($aBrush) Then
		DllCall("user32.dll", "int", "FillRect", "handle", $wParam, "struct*", $tClip, "handle", $aBrush[0])
	EndIf
	Return 1
EndFunc   ;==>_WM_ERASEBKGND

