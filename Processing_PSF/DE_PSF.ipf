#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma modulename=DE_PSF
#include ":\Misc_PanelPrograms\Panel Progs"
#include ":\Misc_PanelPrograms\AsylumNaming"
#include "EquilibriumJumpData"
#include "BigPSD"

Static Function ButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	string saveDF
	variable FOffset, Soffset

	switch( ba.eventCode )
		case 2: // mouse up
			saveDF = GetDataFolder(1)
			wave/T parmWave=root:DE_PSF:MenuStuff:ParmWave
			controlinfo de_PSF_popup0
			SetDataFolder s_value
			controlinfo de_PSF_popup1
			wave StartForceWave=$S_value
			wave StartSepWave=$ReplaceString("Force",S_value,"Sep")
			wave StartZWave=$ReplaceString("Force",S_value,"ZSnsr")
			string ForceWaveString=(S_value)
			NewDataFolder/o $("root:DE_PSF:"+ForceWaveString)
			controlinfo/w=PSF_Panel de_PSF_setvar0
			variable startpnt=v_value
			controlinfo/w=PSF_Panel de_PSF_setvar1
			variable endpnt=v_value

			duplicate/o/R=[startpnt,endpnt] StartForceWave ForceWave
			duplicate/o/R=[startpnt,endpnt] StartSepWave SepWave
			duplicate/o/R=[startpnt,endpnt] StartZWave ZWave

			PSFStepping(ZWave,SepWave,ForceWave,ForceWaveString)
			wave FWidth= $(nameofwave(StartForceWave)+"_Width")
			wave Favg= $(nameofwave(StartForceWave)+"_avg")
			wave SWidth=$(nameofwave(StartForceWave)+"_Width")
			wave Savg=$(nameofwave(StartForceWave)+"_Savg")
			DoWindow ForceWidth
			if(V_flag==0)
			
				display/n=ForceWidth FWidth vs Favg
				ModifyGraph/W=ForceWidth mode=3,marker=19,rgb=(58368,6656,7168)
			else 
				killwindow ForceWidth

				display/n=ForceWidth FWidth vs Favg
				ModifyGraph/W=ForceWidth mode=3,marker=19,rgb=(58368,6656,7168)

			endif
			killwaves ForceWave,SepWave,ZWave
			SetDataFolder saveDF
				
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End
Static Function ButtonProc2(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	string saveDF
	variable FOffset, Soffset

	switch( ba.eventCode )
		case 2: // mouse up
			saveDF = GetDataFolder(1)
			wave/T parmWave=root:DE_PSF:MenuStuff:ParmWave
			controlinfo de_PSF_popup0
			SetDataFolder s_value
			controlinfo de_PSF_popup1
			wave ForceWave=$S_value
			SetVariable de_PSF_setvar0 value= _NUM:pcsr(A,"ForcePlot")
			SetVariable de_PSF_setvar1 value= _NUM:pcsr(B,"ForcePlot")

			SetDataFolder saveDF
				
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End
Static Function SVP(ctrlName,varNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable varNum
	String varStr
	String varName
	strswitch(ctrlName)
		case "de_PSF_setvar0":

			break
	endswitch
End		
Static Function PopMenuProc(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa
	Variable popNum = pa.popNum
	String popStr = pa.popStr
	string saveDF
	strswitch(pa.ctrlName)
	
		case "de_PSF_popup1":
			switch( pa.eventCode )
				case 2: // mouse up
					DE_PSF#WaveSelected()
					break
				case -1: // control being killed
					break
			endswitch
			break 

	endswitch
	 

return 0
End

Function WaveSelected()
	string saveDF

	saveDF = GetDataFolder(1)
	wave/T parmWave=root:DE_PSF:MenuStuff:ParmWave
	controlinfo de_PSF_popup0
	SetDataFolder s_value
	controlinfo de_PSF_popup1
	wave ForceWave=$S_value
	
	DoWindow ForcePlot
	if(V_flag==0)
			
		display/n=ForcePlot ForceWave
		ShowInfo
		Cursor/P/W=ForcePlot A, $S_value,0				
		Cursor/P/W=ForcePlot B, $S_value,(numpnts(ForceWave)-1)

		//ModifyGraph/W=ForceWidth mode=3,marker=19,rgb=(58368,6656,7168)
	else 
		killwindow ForcePlot

		display/n=ForcePlot ForceWave
		ShowInfo			
		Cursor/P/W=ForcePlot A, $S_value,0				
		Cursor/P/W=ForcePlot B, $S_value,(numpnts(ForceWave)-1)
		//ModifyGraph/W=ForceWidth mode=3,marker=19,rgb=(58368,6656,7168)

	endif
	
	
	
	SetDataFolder saveDF

end

Window PSFPanel() : Panel

	PauseUpdate; Silent 1		// building window...
	NewPanel/N=PSF_Panel /W=(0,0,300,450)
	NewDataFolder/o root:DE_PSF
	NewDataFolder/o root:DE_PSF:MenuStuff

	DE_PSF#UpdateParmWave()
	Button de_PSF_button0,pos={75,130},size={150,20},proc=DE_PSF#ButtonProc,title="GO!"
	Button de_PSF_button1,pos={200,80},size={80,20},proc=DE_PSF#ButtonProc2,title="AddCursors!"

	PopupMenu de_PSF_popup0,pos={75,2},size={129,21},Title="Folder"
	PopupMenu de_PSF_popup0,mode=1,popvalue="X",value= #"DE_PanelProgs#ListFolders()"
	PopupMenu de_PSF_popup1,pos={75,40},size={129,21},Title="ForceWave",proc=DE_PSF#PopMenuProc
	PopupMenu de_PSF_popup1,mode=1,popvalue="X",value= #"DE_PSF#ListWaves(\"de_PSF_popup0\")"
	SetVariable de_PSF_setvar0,pos={75,75},size={100,20},proc=DE_PSF#SVP,value= _NUM:0,title="StartPoint"
	SetVariable de_PSF_setvar0 value= _NUM:0
	SetVariable de_PSF_setvar1,pos={75,105},size={100,20},proc=DE_PSF#SVP,value= _NUM:0,title="EndPoint"
	SetVariable de_PSF_setvar1 value= _NUM:0
EndMacro

Static Function UpdateParmWave()
	if(exists("root:DE_PSF:MenuStuff:ParmWave")==1)
		wave/t/z Par=root:DE_PSF:MenuStuff:ParmWave
		wave/z Sel=root:DE_PSF:MenuStuff:SelWave
	Else
		make/t/n=(4,2) root:DE_PSF:MenuStuff:ParmWave
		wave/t/z Par=root:DE_PSF:MenuStuff:ParmWave
		make/n=(4,2) root:DE_PSF:MenuStuff:SelWave
		wave/z Sel=root:DE_PSF:MenuStuff:SelWave
		
		Par[0][0]={"Smoothing","Pulling Rate","Distance to Ignore","Distance to Fit"}
		Par[0][1]={"10e-9","50e-9","3e-9","3e-9"}
		Sel[][0]=0
		Sel[][1]=2
	endif

end

Static Function/S ListWaves(ControlStr)
	string ControlStr
	String saveDF

	saveDF = GetDataFolder(1)
	controlinfo $ControlStr
	SetDataFolder s_value
	String list = WaveList("*Force*", ";", "")
	SetDataFolder saveDF
	return list

end

Static function PSFStepping(ZsnsrWave,SepWave,ForceWave,ForceName)
	wave ZsnsrWave,SepWave,ForceWave
	String ForceName
	String SepName=ReplaceString("Force",ForceName,"Sep")
	String FolderString=("root:DE_PSF:"+ForceName)
	DE_EquilJump#GenerateRegions(ZsnsrWave)
	wave MidPoints
	variable p,FOffset,SepOffset,length
	string Fdup,Sdup
	string FWaveName=nameofwave(ForceWave)+"_cut"
	string SWaveName=nameofwave(SepWave)+"_cut"
	string ZWaveName=nameofwave(ZsnsrWave)+"_cut"
	
	variable totaltime= pnt2x(ForceWave,numpnts(ForceWave)-1)
	variable step= totaltime/numpnts(MidPoints)*.75
	print step
	make/free/n=(numpnts(midpoints)) FWidth,SWidth,Favg,Savg,ZWidth,ZAvg
	variable K0=0
	for(p=0;p<numpnts(midpoints);p+=1)
	
		DE_EquilJump#CutoutRegion(SepWave,MidPoints,p,step)
		wave/z CutS=$SWaveName
		//Smooth/S=2 51, CutS
		Make/n=100/free SHIst;DelayUpdate
		Histogram/C/B=1 CutS,SHIst
		CurveFit/W=2/H="1000"/Q/NTHR=0 gauss  SHIst 
		wave w_coef
		SWidth[p]=w_coef[3]
		Savg[p]=w_coef[2]
		DE_EquilJump#CutoutRegion(ForceWave,MidPoints,p,step)
		wave/z CutF=$FWaveName
		//Smooth/S=2 51, CutF
		BIGPSDFUNC(nameofwave(CutF),4,2,1)
		wave FPSD=$(nameofwave(CutF)+"_psd")
		duplicate/o FPSD $(FolderString+":FPSD_"+num2str(p))
				duplicate/o CutF $(FolderString+":CutF_"+num2str(p))

		Make/n=100/free FHIst;DelayUpdate
		Histogram/C/B=1 CutF,FHIst
		CurveFit/W=2/H="1000"/Q/NTHR=0 gauss  FHIst
		wave w_coef
		duplicate/o FHist $(FolderString+":FHIST_"+num2str(p))
		FWidth[p]=w_coef[3]
		Favg[p]=w_coef[2]
		DE_EquilJump#CutoutRegion(ZsnsrWave,MidPoints,p,step)
		wave/z CutZ=$ZWaveName
		//Smooth/S=2 51, CutF
		Make/n=100/free ZHIst;DelayUpdate
		Histogram/C/B=1 CutZ,ZHIst
		CurveFit/W=2/H="1000"/Q/NTHR=0 gauss  ZHIst 
//		duplicate/o CutZ $(FolderString+":CutZ_"+num2str(p))

		wave w_coef
		ZWidth[p]=w_coef[3]
		Zavg[p]=w_coef[2]
	
	endfor

	duplicate/o FWidth $(ForceName+"_Width")
	duplicate/o Favg $(ForceName+"_avg")
	duplicate/o SWidth $(SepName+"_Width")
	duplicate/o SAvg $(SepName+"_Savg")
	wave w_sigma
	killwaves MidPoints,w_coef,W_Sigma,CutS,CutF,Cutz
end


Menu "Equilibrium"

	"Open PSF Panel", PSFPanel()

end