#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma ModuleName = LSViewer

#include "ScanFiltering"

Window LSViewer() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(50,50,1500,700)
	Button  procbutton,pos={900,100},size={150,25},proc= LSViewer#ButtonProc
	Button procbutton title="Process Data"
	Button  plotbutton,pos={1100,100},size={150,25},proc= LSViewer#ButtonProc
	Button plotbutton title="Plot Data"
	PopupMenu BaseWave,pos={750,100},size={150,25}, proc=LSViewer#PopMenuProc,value=LSViewer#AllFiles()
	SetVariable setvar0, pos={700,10},size={150,25},title="Folded",proc=LSViewer#SetVarProc,value=_NUM:0
	SetVariable setvar1, pos={850,10},size={150,25},title="Unfolded",proc=LSViewer#SetVarProc,value=_NUM:0
	SetVariable setvar2, pos={950,10},size={150,25},title="Edge",proc=LSViewer#SetVarProc,value=_NUM:5
	SetVariable setvar3, pos={1100,10},size={150,25},title="StartFilter",proc=LSViewer#SetVarProc,value=_NUM:150
	SetVariable setvar4, pos={700,35},size={150,25},title="EndFilter",proc=LSViewer#SetVarProc,value=_NUM:2
	SetVariable setvar5, pos={850,35},size={150,25},title="FilterStepSize",proc=LSViewer#SetVarProc,value=_NUM:0
	SetVariable setvar6, pos={950,35},size={150,25},title="PFoldStepSize",proc=LSViewer#SetVarProc,value=_NUM:0
	SetVariable setvar7, pos={1100,35},size={150,25},title="PFoldSmoothing",proc=LSViewer#SetVarProc,value=_NUM:0
	SetVariable setvar8, pos={950,60},size={150,25},title="Dq",proc=LSViewer#SetVarProc,value=_NUM:0

EndMacro

Static Function/S AllFiles()
	return WaveList("*", ";","")
end

Static Function MyWindowHook(s)
	STRUCT WMWinHookStruct &s
	
	Variable hookResult = 0	// 0 if we do not handle event, 1 if we handle it.

	switch(s.eventCode)
		case 7:					// Keyboard event
			UpdateThatShit()
			break
	endswitch

	return hookResult	// If non-zero, we handled event and Igor will ignore it.
End

Static function UpdateThatShit()
	setdatafolder root:FilterScan

	string BN="Boltz_"
	string PN="Pfold_"
	variable point=pcsr(A,"LSViewer#Result")
	wave sm
	wave w1=$BN+num2str(sm[point])
	wave w2=$PN+num2str(sm[point])
	duplicate/o w1 CurrentBolt 
	duplicate/o w2 CurrentPfold 
	variable off=CurrentBolt[(numpnts(CurrentBolt))/2]
	CurrentBolt-=off
		 off=CurrentPfold[(numpnts(CurrentPfold))/2]

	CurrentPfold-=off
	wave BoltW,BoltB,PfoldW,PfoldB,PFOldU
	BoltW=CurrentBolt[str2num(stringfromlist(1,note(CurrentBolt)))]
	SetScale/P x pnt2x(CurrentBolt,str2num(stringfromlist(1,note(CurrentBolt)))),1,"", BoltW
	BoltB=CurrentBolt[str2num(stringfromlist(0,note(CurrentBolt)))]
	SetScale/P x pnt2x(CurrentBolt,str2num(stringfromlist(0,note(CurrentBolt)))),1,"", BoltB

	PfoldW=CurrentPfold[str2num(stringfromlist(0,note(CurrentPfold)))]
	SetScale/P x pnt2x(CurrentPfold,str2num(stringfromlist(0,note(CurrentPfold)))),1,"", PfoldW

	PfoldB=CurrentPfold[str2num(stringfromlist(1,note(CurrentPfold)))]
	SetScale/P x pnt2x(CurrentPfold,str2num(stringfromlist(1,note(CurrentPfold)))),1,"", PfoldB

	PFOldU=CurrentPfold[str2num(stringfromlist(2,note(CurrentPfold)))]
	SetScale/P x pnt2x(CurrentPfold,str2num(stringfromlist(2,note(CurrentPfold)))),1,"", PFOldU

	SetAxis/W=LSViewer#LS/A=2 left;DelayUpdate
	SetAxis/W=LSViewer#LS bottom (pnt2x(CurrentPFold,0)-.5),pnt2x(CurrentPFold,numpnts(CurrentPFold)-1)+.5
	setdatafolder root:
end


Static Function ButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			strswitch( ba.ctrlName)
				case "plotbutton":
					PlotData()
			
					break
				case "procbutton":
					ProcData()
					break
			endswitch
		
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Static Function ProcData()
	variable folded,unfolded, edge,StartFilter,EndFilter,FilterSteps,PFoldStep,PFoldSmoothing,Dq
	controlinfo setvar0
	folded=V_Value
	controlinfo setvar1
	unfolded=V_Value
	controlinfo setvar2
	edge=V_Value
	controlinfo setvar3
	StartFilter=V_Value
	controlinfo setvar4
	EndFilter=V_Value

	controlinfo setvar5
	FilterSteps=V_Value
	controlinfo setvar6
	PFoldStep=V_Value
	controlinfo setvar7
	PFoldSmoothing=V_Value
	controlinfo setvar8
	Dq=V_Value

	controlinfo BaseWave
	print S_Value
	wave BaseWave=$S_Value

	LandscapeFiltering#GeneralizedIterativePFold(BaseWave,folded,unfolded, edge,StartFilter,EndFilter,FilterSteps,PFoldStep,PFoldSmoothing,Dq)
end

Static Function PlotData()
	if(datafolderexists("root:FilterScan")==0)
		print "No Data"
		return -1
	endif
	
	setdatafolder root:filterScan
		
	wave Boltzmann_Height, pfold_Height,Smoothing,Boltzmann_Distance,PFold_Distance,Cbar,Cwell,Ka
	
	if(cmpstr(TraceNameList("LSViewer#Result", ";",1 ),"")==0)
		make/o/n=0 BH,PH,Sm

		display/W=(50,50,700,500)/Host=LSViewer/N=Result
		SetWindow LSViewer, hook(MyHook) = LSViewer#MyWindowHook
		Appendtograph/W=LSViewer#Result BH vs Sm
		Appendtograph/W=LSViewer#Result pH vs Sm

		duplicate/o Boltzmann_Height BH
		duplicate/o pfold_Height pH
		duplicate/o Smoothing sm
		make/o/n=0 BX,PX
			
		AppendToGraph/W=LSViewer#Result/L=L1/B=B1 BX vs sm
		AppendToGraph/W=LSViewer#Result/L=L1/B=B1 PX vs sm
		duplicate/o Boltzmann_Distance BX
		duplicate/o PFold_Distance px
			
		make/o/n=0 CB,CW
			
		AppendToGraph/W=LSViewer#Result/L=L2/B=B2 cb vs sm
		AppendToGraph/W=LSViewer#Result/L=L2/B=B2 cw vs sm
		duplicate/o Cbar CB
		duplicate/o Cwell CW
			
		make/o/n=0 KAD
			
		AppendToGraph/W=LSViewer#Result/L=L3/B=B3 KAD vs sm
		duplicate/o Ka, KAD
			
		ModifyGraph axisEnab(left)={0,0.23}
		ModifyGraph axisEnab(l1)={0.27,.48}
		ModifyGraph axisEnab(l2)={0.52,.73}
		ModifyGraph axisEnab(l3)={0.77,1}

		ModifyGraph fSize=11,axisOnTop=1,standoff=0,font="Arial"
		ModifyGraph freePos(L1)={0,B1}
		ModifyGraph freePos(L2)={0,B2}
		ModifyGraph freePos(L3)={0,B3}
		ModifyGraph freePos(B1)=-95
		ModifyGraph freePos(B2)=-175
		ModifyGraph freePos(B3)=-257
		ModifyGraph noLabel(B1)=2
		ModifyGraph noLabel(B2)=2
		ModifyGraph noLabel(B3)=2
		ModifyGraph mode=4,marker=19,msize=2,rgb(Bx)=(58368,6656,7168);DelayUpdate
		ModifyGraph rgb(BH)=(58368,6656,7168);DelayUpdate
		ModifyGraph rgb(Px)=(14848,32256,47104)
		ModifyGraph rgb(Ph)=(14848,32256,47104)
		ModifyGraph rgb(CB)=(65280,32512,0)
		ModifyGraph rgb(CW)=(19712,44800,18944)
		ModifyGraph rgb(KAD)=(0,0,0)

		ModifyGraph lblPos(L1)=50;DelayUpdate
		ModifyGraph lblPos(L2)=50;DelayUpdate
		ModifyGraph lblPos(L3)=50;DelayUpdate


		Label left "\\f01Barrier (kT)"
		Label l1 "\\f01Distance (nm)"
		Label l2 "\\f01Curvature (pN/nm\S2\M)"
		Label l3 "\\f01kA (s\S-1\M)"
		ModifyGraph tick=2,mirror(left)=1,mirror(L1)=1,mirror(L2)=1,mirror(L3)=1,fSize=11
		Label bottom "\\f01Smoothing (points)"

	else
			duplicate/o Boltzmann_Height BH
		duplicate/o pfold_Height pH
		duplicate/o Smoothing sm
				duplicate/o Boltzmann_Distance BX
		duplicate/o PFold_Distance px
			duplicate/o Cbar CB
		duplicate/o Cwell CW
				duplicate/o Ka, KAD

	endif

	make/o/n=(0,0,0) CurrentPfold, CurrentBolt
	make/o/n=1 BoltW,BoltB,PfoldW,PfoldB,PFOldU
	display/W=(750,150,1300,500)/Host=LSViewer/N=LS
	Appendtograph CurrentPfold
	Appendtograph CurrentBolt
	Appendtograph BoltW,BoltB,PfoldW,PfoldB,PFOldU
	ModifyGraph mode(BoltW)=3,marker(BoltW)=19,mode(BoltB)=3,marker(BoltB)=19;DelayUpdate
	ModifyGraph mode(PfoldW)=3,marker(PfoldW)=19,mode(PfoldB)=3,marker(PfoldB)=19
	ModifyGraph mode(PFOldU)=3,marker(PFOldU)=19
	ModifyGraph rgb(CurrentBolt)=(0,0,0),rgb(BoltW)=(0,0,0),rgb(BoltB)=(0,0,0)
	ModifyGraph rgb(CurrentPfold)=(58368,6656,7168),rgb(PfoldW)=(58368,6656,7168);DelayUpdate
	ModifyGraph rgb(PfoldB)=(58368,6656,7168)
	Cursor/P/W=LSViewer#Result  A Bx (numpnts(Bx)-1)
	setdatafolder root:

	//UpdateThatShit()

	//




end

Static Function PopMenuProc(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa

	switch( pa.eventCode )
		case 2: // mouse up
			Variable popNum = pa.popNum
			String popStr = pa.popStr
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Static Function SetVarProc(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End
