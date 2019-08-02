#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma modulename=DE_RuptureRamp
#include "DE_Filtering"
#include "SimpleWLCPrograms"
#include "C:Users:dedwards:src_prh:IgorUtil:IgorCode:Util:OperatingSystemUtil"
#include "DTE_Dudko"
#include "DE_OverlapRamps"
#include "DE_CorrectRupture"
#include "DE_TwoWLCFit"
#include ":\Misc_PanelPrograms\Panel Progs"
#include ":\Misc_PanelPrograms\AsylumNaming"
//#include ":\Processing_Markov\DE_HMM"
Function CorrectPauses(ForceWave,SepWave)
	wave ForceWave,SepWave
	string originalInds=stringbykey("DE_Ind",note(ForceWave),":","\r")
	string originalpauses=stringbykey("DE_PauseLoc",note(ForceWave),":","\r")
	string PauseState=stringbykey("DE_PauseState",note(ForceWave),":","\r")

	
	variable maxcntInds=itemsinlist(originalInds)
	variable maxcntPauses=itemsinlist(originalInds)
	string newpause=""
	string NewPauseState=""

	variable n,m
	for(n=1;n<maxcntInds;n+=2)
		variable target= str2num(StringFromList(n, originalInds))
		for(m=1;m<maxcntPauses;m+=2)
		variable Pause=str2num( StringFromList(m, originalpauses))
		 	if(Pause==target)
		 		newpause+=(StringFromList(m-1, originalpauses))+";"
		 		newpause+=StringFromList(m, originalpauses)+";"
				NewPauseState+=(StringFromList(m-1, PauseState))+";"
				NewPauseState+=(StringFromList(m, PauseState))+";"

		 	endif
		endfor
	endfor
	string NewNote=ReplaceStringbykey("DE_PauseLoc",note(ForceWave),Newpause,":","\r")
	NewNote=ReplaceStringbykey("DE_PauseState",(NewNote),NewPauseState,":","\r")
	note/K ForceWave NewNote
	note/K SepWave NewNote
end
Static Function PythonFitter( UseWave,Method,Threshold,AMount)//Variables demanded by MarkovFit Code
	Wave UseWave//Input wave
	String Method
	Variable Threshold,AMount
	String Destination = "D:\Data\StepData\Test1.ibw"
	String NewHome = "D:\Data\StepData\Shit.txt"

	Save/O/C UseWave as Destination
	String BasePythonCommand = "cmd.exe /c activate & python D:\Devin\Python\StepAttempt\StepAttempt.py "
	String MethodCommand="-method "+ method +" "
	String InputCom="-inputfile "+ Destination+" "
	String OutputCom="-outputfile "+ NewHome+" "
	String SmoothCommand="-smooth "+ num2str(Amount)+" "
	String ThreshCommand="-threshold "+ num2str(threshold)+" "

//
	String PythonCommand=BasePythonCommand+MethodCommand+InputCom+OutputCom+SmoothCommand+ThreshCommand
	ModOperatingSystemUtil#execute_python(PythonCommand)
	LoadWave/O/N/G/D NewHome
	

end



Static Function PopMenuProc(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa
	switch( pa.eventCode )
		case 2: // mouse up

			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Static Function ButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			strswitch( ba.ctrlName)
				case "de_RupRamp_button0":
					//FitHMMButt()
					break
				case  "de_RupRamp_button1":

					FitPython()
					break
				case  "de_RupRamp_button2":
					MakeStateWave()
					break
				case  "de_RupRamp_button3":
					MakeOffsetForceWave()
					break
				case  "de_RupRamp_button4":
					DetermineWLCParms()
					break
				case  "de_RupRamp_button5":
					CorrectRuptures()
					break
				case  "de_RupRamp_button6":
					MakeSecondStateWave()
					break
				case  "de_RupRamp_button7":
					
					MakeSmoothedWave()
					break
				case  "de_RupRamp_button8":
					MakeSmoothedShWave()
					break
				case  "de_RupRamp_button9":
					AutoProcess()
					break
				case  "de_RupRamp_button10":
					AligntoWLC()
					break
				case  "de_RupRamp_button11":
					MakeSmoothedAlignedWave()
					break
			endswitch
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Static Function MakeSmoothedAlignedWave()
	
	string saveDF

	saveDF = GetDataFolder(1)
	controlinfo/W=RupRampPanel de_RupRamp_popup0
	SetDataFolder s_value
	controlinfo/W=RupRampPanel de_RupRamp_popup14
	wave AlignedForceWave=$S_value
	wave AlignedSepWave=$ReplaceString("Force",S_value,"Sep")
	//controlinfo/W=RupRampPanel de_RupRamp_popup4
	//wave ForceWave=$S_value

	//wave SepWave=$ReplaceString("Force",S_value,"Sep")
	
	controlinfo/w=RupRampPanel de_RupRamp_popup15
	string SmType=S_Value
	controlinfo/w=RupRampPanel de_RupRamp_setvar2
	variable SmAmount=v_value

	duplicate/free AlignedForceWave ForceWaveSm
	duplicate/free AlignedSepWave SepWaveSm

	DE_Filtering#FilterForceSep(AlignedForceWave,AlignedSepWave,ForceWaveSm,SepWaveSm,SmType,SmAmount)
	duplicate/o ForceWaveSm $(nameofwave(AlignedForceWave)+"_Sm")
	duplicate/o SepWaveSm $(nameofwave(AlignedSepWave)+"_Sm")

	SetDataFolder saveDF
end

Static Function AligntoWLC()
	string saveDF
	saveDF = GetDataFolder(1)

	controlinfo/W=RupRampPanel de_RupRamp_popup0
	SetDataFolder s_value
	controlinfo/W=RupRampPanel de_RupRamp_popup4
	wave ForceWave=$S_value
	controlinfo/W=RupRampPanel de_RupRamp_popup4
	wave SepWave=$ReplaceString("Force",S_value,"Sep")
//	controlinfo/W=RupRampPanel de_RupRamp_popup3
//	wave UpPoints=$S_value
//	wave DownPoints=$ReplaceString("PntU",S_value,"PntD")
	controlinfo/W=RupRampPanel de_RupRamp_popup5
	wave StateWave=$S_value
	controlinfo/W=RupRampPanel de_RupRamp_popup6
	wave ShiftedForceWave=$S_value
	
	controlinfo/W=RupRampPanel de_RupRamp_popup9
	wave ForceWaveSM=$S_value
	wave SepWaveSm=$ReplaceString("Force",S_value,"Sep")

	controlinfo/W=RupRampPanel de_RupRamp_popup10
	wave ForceWaveSH_SM=$S_value
	controlinfo/W=RupRampPanel de_RupRamp_popup13
	wave WLCParms=$S_value
	Controlinfo/W=RupRampPanel de_RupRamp_setvar3
//
	variable/C slopes=DE_Dudko#ReturnSeparationSlopes(SepWaveSm,StateWave,500)
	variable pointstoignore=floor(v_value/real(slopes)/dimdelta(ForceWaveSH_SM,0))
	make/o/n=0 FoldedForceOut,FoldedSepOut,UnFoldedForceOut,UnFoldedSepOut
	DE_NewDudko#ReturnFoldandUnfold(ForceWaveSH_SM,SepWaveSm,StateWave,pointstoignore,FoldedForceOut,FoldedSepOut,UnFoldedForceOut,UnFoldedSepOut)

	FoldedForceOut*=-1
	UnFoldedForceOut*=-1

	
	make/o/n=0 ResultsShift,ResultsNoShift
	controlinfo/W=RupRampPanel de_RupRamp_popup17
	variable foldedfit
	strswitch(S_Value)
	
		case "Unfolded":
			AlignSingletoWLC(UnFoldedForceOut,UnFoldedSepOut,WLCParms,0,ResultsNoShift)
			AlignSingletoWLC(UnFoldedForceOut,UnFoldedSepOut,WLCParms,1,ResultsShift)
			foldedfit=0
			break
		
		case "Folded":
			AlignSingletoWLC(FoldedForceOut,FoldedSepOut,WLCParms,0,ResultsNoShift)
			AlignSingletoWLC(FoldedForceOut,FoldedSepOut,WLCParms,1,ResultsShift)
			foldedfit=1
			break
		
		case "Both":
			AlignTwotoWLC(FoldedForceOut,FoldedSepOut,UnFoldedForceOut,UnFoldedSepOut,WLCParms,0,ResultsNoShift,  0)
			AlignTwotoWLC(FoldedForceOut,FoldedSepOut,UnFoldedForceOut,UnFoldedSepOut,WLCParms,1,ResultsShift,0)

			foldedfit=2
			break
	
	endswitch
	
	
	controlinfo/W=RupRampPanel de_RupRamp_check0
	variable allowsepshift=v_value
	
	variable forceshiftused,sepshiftused,forceshiftalt,sepshiftalt
	string shiftString
	if(allowsepshift==1)
		forceshiftused=ResultsShift[0]
		sepshiftused=ResultsShift[1]
		forceshiftalt=ResultsNoShift[0]
		sepshiftalt=ResultsNoShift[1]
		shiftString="YES"
		
	elseif(allowsepshift==0)
		forceshiftalt=ResultsShift[0]
		sepshiftalt=ResultsShift[1]
		forceshiftused=ResultsNoShift[0]
		sepshiftused=ResultsNoShift[1]
		shiftString="No"

	endif
	

	duplicate/o ShiftedForceWave $ReplaceString("Force_Shift",nameofwave(ShiftedForceWave),"Force_Align")
	wave AlignedFWave=$ReplaceString("Force_Shift",nameofwave(ShiftedForceWave),"Force_Align")
	AlignedFWave-=forceshiftused
	print num2str(forceshiftused)+";"+num2str(sepshiftused)
	duplicate/o SepWave $ReplaceString("Sep_Adj",nameofwave(SepWave),"Sep_Align")
	wave AlignedSepWave=$ReplaceString("Sep_Adj",nameofwave(SepWave),"Sep_Align")
	AlignedSepWave-=sepshiftused
	
	String NewNote=ReplaceStringByKey("UsedAlignmentFShift",note(AlignedFWave),num2str(forceshiftused),":","\r")
	NewNote=ReplaceStringByKey("UsedAlignmentSShift",NewNote,num2str(sepshiftused),":","\r")
	
	if(FoldedFit==1)
		NewNote=ReplaceStringByKey("Aligned To",NewNote,"Folded",":","\r")

	elseif(FoldedFit==0)
		NewNote=ReplaceStringByKey("Aligned To",NewNote,"UnFolded",":","\r")
	elseif(FoldedFit==2)
		NewNote=ReplaceStringByKey("Aligned To",NewNote,"Both",":","\r")
	endif
	NumericWaveToStringList(WLCParms)
	NewNote=ReplaceStringByKey("WLCParmsforAlign",NewNote,NumericWaveToStringList(WLCParms),":","\r")

	note/K AlignedFWave, NewNote
	note/K AlignedSepWave, NewNote
	
	print num2str(forceshiftalt)+";"+num2str(sepshiftalt)
	NewNote=ReplaceStringByKey("AltAlignmentFShift",note(AlignedFWave),num2str(forceshiftalt),":","\r")
	NewNote=ReplaceStringByKey("AltAlignmentSShift",NewNote,num2str(sepshiftalt),":","\r")
	note/K AlignedFWave, NewNote
	note/K AlignedSepWave, NewNote

end

Static Function AlignSingletoWLC(ForceIn,SepIn,WLCParms,SepShift,Results)
	Wave ForceIn,SepIn,WLCParms,Results
	variable SepShift
//	
	duplicate/o WLCParms W_coef
	string HoldString="11100"
	if(SepShift==1)
	HoldString="11100"
	else 
	HoldString="11101"
	endif
	
	FuncFit/Q/H=HoldString/NTHR=0 WLC_FIT W_coef  ForceIn /X=SepIn/D

	make/free/n=2 TemporaryResults

	TemporaryResults={w_coef[3]-WLCParms[3],w_coef[4]-WLCParms[4]}

	duplicate/o TemporaryResults Results
end

Static Function AlignTwotoWLC(FoldedForce,FoldedSep,UnfoldedForce,UnfoldedSep,WLCParms,SepShift,Results,FixedShift)
	wave FoldedForce,FoldedSep,UnfoldedForce,UnfoldedSep,WLCParms,Results
	variable SepShift,FixedShift
//	FitForcePair(FoldedForce,FoldedSep,UnFoldedForce,UnFoldedSep,ParmOut,FitOut,[Constrained])
	
//	wave FoldedForce,FoldedSep,UnFoldedForce,UnFoldedSep,ParmOut,FitOut,Constrained
	variable 	timerRefNum = startMSTimer


//	//note that WLCGuess should have the format: Lp,Lc1,Lc2,T,Xoff,Foff,
	make/free/n=0 Fout,Sout
	DE_TwoWLCFit#CombineCurves(FoldedForce,FoldedSep,UnFoldedForce,UnFoldedSep,Fout,SOut)
	duplicate/free WLCParms w_coef
w_coef[4]=FixedShift+WLCParms[4]
	//A bit of code to fix the offset to something specific
	string ConStr="111100"
	if(SepShift==1)
		ConStr="111100"
	elseif(SepShift==0)
		ConStr="111110"
	endif

	FuncFit/N/Q/W=2/H=ConStr/NTHR=0 DE_FitTwo W_coef  Fout /X=SOut 
	make/free/n=2 TemporaryResults

	TemporaryResults={-w_coef[5]+WLCParms[5],w_coef[4]-WLCParms[4]}
	make/free/n=0 FWLCFit
	DE_TwoWLCFit#MakeAMultiFit(SOut,w_coef,FWLCFit)
	make/o/n=(Numpnts(FWLCFit),2) FitTest
	FitTest[][1]=FWLCFit[p]
	FitTest[][0]=SOut[p][0]
	duplicate/o TemporaryResults Results
end

Static Function/S NumericWaveToStringList(w)
	Wave w	// numeric wave (if text, use /T here and %s below)
	String list
	wfprintf list, "%g;", w	// semicolon-separated list
	return list
End

Static Function CheckProc(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End


Static Function AutoProcess()

	string saveDF
	saveDF = GetDataFolder(1)
	controlinfo/W=RupRampPanel de_RupRamp_popup0
	SetDataFolder s_value
	
	controlinfo/W=RupRampPanel de_RupRamp_popup4
	if(cmpstr(S_Value,"")==0)
		print "No Force Wave"
		return 0
	endif
	
		controlinfo/W=RupRampPanel de_RupRamp_popup3
	if(cmpstr(S_Value,"")==0)
		print "No Points Wave"
		return 0
	endif
	MakeSmoothedWave()
	controlinfo/W=RupRampPanel de_RupRamp_popup4
	ControlUpdate/w=RupRampPanel de_RupRamp_popup9
	popupmenu de_RupRamp_popup9 win=RupRampPanel,popmatch=S_value+"_Sm"
	
	MakeStateWave()
	ControlUpdate/w=RupRampPanel de_RupRamp_popup5
	popupmenu de_RupRamp_popup9 win=RupRampPanel,popmatch=replacestring("Force",S_value,"_States")

	MakeOffsetForceWave()
	ControlUpdate/w=RupRampPanel de_RupRamp_popup6
	popupmenu de_RupRamp_popup6 win=RupRampPanel,popmatch=ReplaceString("Force",S_value,"Force_Shift")

	MakeSmoothedShWave()
	ControlUpdate/w=RupRampPanel de_RupRamp_popup10
	popupmenu de_RupRamp_popup10 win=RupRampPanel,popmatch=ReplaceString("Force",S_value,"Force_Shift")+"_Sm"
	//	ControlUpdate/w=RupRampPanel de_RupRamp_popup13
	//popupmenu de_RupRamp_popup13 win=RupRampPanel,popmatch="WLCAlign"

	AligntoWLC()
	ControlUpdate/w=RupRampPanel de_RupRamp_popup14
	popupmenu de_RupRamp_popup14 win=RupRampPanel,popmatch=ReplaceString("Force_Adj",S_value,"Force_Align")
	MakeSmoothedAlignedWave()
	ControlUpdate/w=RupRampPanel de_RupRamp_popup16
	popupmenu de_RupRamp_popup16 win=RupRampPanel,popmatch=ReplaceString("Force_Adj",S_value,"Force_Align")+"_Sm"
	
	
	DetermineWLCParms()
	ControlUpdate/w=RupRampPanel de_RupRamp_popup7
	popupmenu de_RupRamp_popup7 win=RupRampPanel,popmatch=ReplaceString("Force",S_value,"_WLCParms")
	
	CorrectRuptures()
	ControlUpdate/w=RupRampPanel de_RupRamp_popup8
	controlinfo/W=RupRampPanel de_RupRamp_popup3

	popupmenu de_RupRamp_popup8 win=RupRampPanel,popmatch=S_Value+"_Mod"

	MakeSecondStateWave()
	SetDataFolder saveDF


end

Static Function MakeSmoothedWave()
	string saveDF

	saveDF = GetDataFolder(1)
	controlinfo/W=RupRampPanel de_RupRamp_popup0
	SetDataFolder s_value
	controlinfo/W=RupRampPanel de_RupRamp_popup4
	wave ForceWave=$S_value
	wave SepWave=$ReplaceString("Force",S_value,"Sep")
	CorrectPauses(ForceWave,SepWave)
	controlinfo/w=RupRampPanel de_RupRamp_popup11
	string SmType=s_value
	controlinfo/w=RupRampPanel de_RupRamp_setvar0
	variable SmAmount=v_value
	duplicate/free ForceWave ForceWaveSm
	duplicate/free SepWave SepWaveSm
	DE_Filtering#FilterForceSep(ForceWave,SepWave,ForceWaveSm,SepWaveSm,SmType,SmAmount)

	duplicate/o ForceWaveSm $(nameofwave(ForceWave)+"_Sm")
	duplicate/o SepWaveSm $(ReplaceString("Force",nameofwave(ForceWave),"Sep")+"_Sm")
	SetDataFolder saveDF

end
Static Function MakeSmoothedShWave()

	string saveDF

	saveDF = GetDataFolder(1)
	controlinfo/W=RupRampPanel de_RupRamp_popup0
	SetDataFolder s_value
	controlinfo/W=RupRampPanel de_RupRamp_popup6
	wave SHForceWave=$S_value
	controlinfo/W=RupRampPanel de_RupRamp_popup4
	wave ForceWave=$S_value

	wave SepWave=$ReplaceString("Force",S_value,"Sep")
	
	controlinfo/w=RupRampPanel de_RupRamp_popup12
	string SmType=S_Value
	controlinfo/w=RupRampPanel de_RupRamp_setvar1
	variable SmAmount=v_value
	
	duplicate/free ForceWave ForceWaveSm
	duplicate/free SepWave SepWaveSm
	
	DE_Filtering#FilterForceSep(SHForceWave,SepWave,ForceWaveSm,SepWaveSm,SmType,SmAmount)

	duplicate/o ForceWaveSm $(nameofwave(SHForceWave)+"_Sm")
	duplicate/o SepWaveSm $(ReplaceString("Force",nameofwave(SHForceWave),"Sep")+"_Sm")
	SetDataFolder saveDF

end
Static Function MakeSecondStateWave()


	string saveDF

	saveDF = GetDataFolder(1)
	controlinfo/W=RupRampPanel de_RupRamp_popup0
	SetDataFolder s_value
	controlinfo/W=RupRampPanel de_RupRamp_popup4
	wave ForceWave=$S_value
	wave SepWave=$ReplaceString("Force",S_value,"Sep")
	controlinfo/W=RupRampPanel de_RupRamp_popup8
	wave UpPoints=$S_value
	wave DownPoints=$ReplaceString("PntU",S_value,"PntD")
	controlinfo/W=RupRampPanel de_RupRamp_popup6
	wave ShForceWave=$S_value
	controlinfo/W=RupRampPanel de_RupRamp_popup14
	wave AlignedForceWave=$S_value
	wave AlignedsepWave=$ReplaceString("Force",S_value,"Sep")	
//	duplicate/free ShForceWave ForceWaveSm
//	duplicate/free SepWave SepWaveSm

	//ForceWaveSm*=-1;

	//make the state key
	make/o/n=0 $ReplaceString("Force",nameofwave(ForceWave),"_2States")
	wave states=$ReplaceString("Force",nameofwave(ForceWave),"_2States")
	DE_DUDKO#MakeSingleStateKey(AlignedForceWave,UpPoints,DownPoints,States)

	SetDataFolder saveDF

end


Static Function CorrectRuptures()

	string saveDF
	saveDF = GetDataFolder(1)
	controlinfo/W=RupRampPanel de_RupRamp_popup0
	SetDataFolder s_value
	controlinfo/W=RupRampPanel de_RupRamp_popup4
	wave ForceWave=$S_value
	wave SepWave=$ReplaceString("Force",S_value,"Sep")
	controlinfo/W=RupRampPanel de_RupRamp_popup3
	wave UpPoints=$S_value
	wave DownPoints=$ReplaceString("PntU",S_value,"PntD")
	controlinfo/W=RupRampPanel de_RupRamp_popup5
	wave StateWave=$S_value
	controlinfo/W=RupRampPanel de_RupRamp_popup14
	wave AlignedForceWave=$S_value
	wave AlignedsepWave=$ReplaceString("Force",S_value,"Sep")
	controlinfo/W=RupRampPanel de_RupRamp_popup7
	wave WLCParms=$S_value
	duplicate/free AlignedForceWave FWSm,SepWaveSm,FWIN
//	duplicate/free ForceWave SepWaveSm
//	controlinfo/W=RupRampPanel de_RupRamp_popup9
//	wave ForceWaveSM=$S_value
//	wave SepWaveSm=$ReplaceString("Force",S_value,"Sep")
		
	controlinfo/W=RupRampPanel de_RupRamp_popup10
	wave ForceWaveSH_SM=$S_value
	DE_Filtering#FilterForceSep(AlignedForceWave,AlignedsepWave,FWSm,SepWaveSm,"TVD",10e-9)
	FWIN*=-1
	FWSm*=-1

	DE_CorrRup#FixPicks(UpPoints,FWIN,FWSm,WLCParms,SepWaveSm)
//	make/o/n=0 $ReplaceString("Force",nameofwave(ForceWave),"_WLCParms")
//	wave Results=$ReplaceString("Force",nameofwave(ForceWave),"_WLCParms")
//	DE_DUDKO#ContourLengthDetermineCombined(ForceWaveSH_SM,SepWaveSm,StateWave,10000,Results)
	SetDataFolder saveDF


end
Static Function FitPython()
	string saveDF
	variable FOffset, Soffset

	saveDF = GetDataFolder(1)
	wave/T parmWave=root:DE_RupRamp:MenuStuff:PyParmWave

	controlinfo de_RupRamp_popup0
	SetDataFolder s_value
	controlinfo de_RupRamp_popup1
	wave ForceWave=$S_value
	wave SepWave=$ReplaceString("Force",S_value,"Sep")
	//	make/o/n=0 RupPntU,RupPntD
	make/o/n=0 ForceWaveS,SepWaveS
	//DE_Filtering#FilterForceSep(ForceWave,SepWave,ForceWaveS,SepWaveS,"SVG",str2num(parmWave[6][1]))
	DE_Filtering#FilterForceSep(ForceWave,SepWave,ForceWaveS,SepWaveS,"TVD",str2num(parmWave[0][1]))
	ForceWaveS*=-1
	wavestats/q ForceWaveS
	FOffset=v_min
	FOffset-=0e-12
	ForceWaveS-=FOffset
	Soffset= SepWaveS[0]-0e-9
	SepWaveS-=Soffset
	//duplicate/o ForceWaveS LC
	//LC=DE_WLC#ContourTransform(ForceWaveS,SepWaveS,.4e-9,298)		
	//LC*=1e9	
	//setscale/P x dimoffset(ForceWave,0), dimdelta(ForceWave,0), "s", LC//ensuring scaling of input and output wave are the same
	//Resample/DOWN=(str2num(parmWave[7][1]))/N=1/WINF=None LC
	Resample/DOWN=(str2num(parmWave[1][1]))/N=1/WINF=None ForceWaveS
	controlinfo de_RupRamp_popup2
	PythonFitter( ForceWaveS,S_Value,str2num(parmWave[3][1]),str2num(parmWave[2][1]))
	wave no0=wave0
	make/o/n=0 UpP,DownP;DE_RuptureRamp#SortSteps(ForceWaveS,no0,UpP,DownP,20)
	DE_RuptureRamp#RefineSteps(ForceWaveS,UpP,DownP)
	duplicate/o UpP RupPntU
	RupPntU=x2pnt(ForceWave,pnt2x(ForceWaveS,UpP))
	duplicate/o DownP RupPntD
	RupPntD=x2pnt(ForceWave,pnt2x(ForceWaveS,DownP))

	killwaves SepWaveS,ForceWaveS,no0,DownP,UpP
	SetDataFolder saveDF
end

Static Function IdentifyAvgForceDuringPauses(ForceWave,ForceAvg)
	wave ForceWave,ForceAvg
	String PauseLoc=stringbykey("DE_PauseLoc",note(ForceWave),":","\r")
	String Ind=stringbykey("DE_Ind",note(ForceWave),":","\r")
	variable Num=itemsinlist(PauseLoc)
	variable n
	make/free/n=(Num/2,2) FreeInfo
	variable/D startpoint,endpoint,midpoint
	for(n=0;n<Num;n+=2)
		
		startpoint=	str2num(stringfromlist(n,PauseLoc))
				endpoint=str2num(stringfromlist(n+1,PauseLoc))
		if(endpoint>numpnts(ForceWave)-1)
			endpoint=(numpnts(ForceWave)-1)
		
		endif
		midpoint=startpoint+(endpoint-startpoint)/2
		wavestats/q/r=[startpoint,endpoint]  ForceWave
		FreeInfo[n/2][0]=pnt2x(ForceWave,midpoint)
		FreeInfo[n/2][1]=v_avg


	endfor
	duplicate/o FreeInfo ForceAvg 

end

Static Function FindAvgForceAfterPoint(pnt,ForceWave)
	wave ForceWave
	variable pnt
	
	
	String PauseLoc=stringbykey("DE_PauseLoc",note(ForceWave),":","\r")
	String Ind=stringbykey("DE_Ind",note(ForceWave),":","\r")
	variable Num=itemsinlist(PauseLoc)
	variable n
	make/free/n=(Num/2,2) FreeInfo
	variable/D rampstart,pausestart,pauseend,midpoint
	for(n=0;n<Num;n+=2)
		if(n==0)
		rampstart=0
		else
		rampstart=str2num(stringfromlist(n-1,Ind))
		endif
		
		
		pausestart=	str2num(stringfromlist(n,PauseLoc))
		pauseend=str2num(stringfromlist(n+1,PauseLoc))
		if(pauseend>numpnts(ForceWave)-1)
			pauseend=(numpnts(ForceWave)-1)
		endif

		if(pnt>rampstart&&pnt<pauseend)
			midpoint=pausestart+(pauseend-pausestart)/2
			wavestats/q/r=[pausestart,pauseend]  ForceWave
			return v_avg
			break
		endif

	endfor
end

Static Function TurnPntsIntoForcesandTime(Forcewave,PntsIn,ForcesOut,TimeOut)

	wave Forcewave,PntsIn,ForcesOut,TimeOut
	
	duplicate/free PntsIn FreeForces FreeTime
	
	FreeForces=ForceWave[PntsIn]
	FreeTime=pnt2x(Forcewave,PntsIn)
	duplicate/o FreeForces ForcesOut
	duplicate/o FreeTime TimeOut
end


Static Function CutOutPauses(ForceWave,SepWave,DicedForce,DicedSep,CutInfo)
	wave ForceWave,SepWave,DicedForce,DicedSep,CutInfo
	String PauseLoc=stringbykey("DE_PauseLoc",note(ForceWave),":","\r")
	String Ind=stringbykey("DE_Ind",note(ForceWave),":","\r")
	variable Num=itemsinlist(PauseLoc)
	variable n
	duplicate/free ForceWave FreeForce
	duplicate/free SepWave FreeSep
	make/free/n=(Num/2,2) FreeInfo
	variable/D startdelete,enddelete,totalDeletions
	for(n=Num-1;n>0;n-=2)
		if(n==Num-1)
		enddelete=numpnts(ForceWave)-1
		startdelete=	str2num(stringfromlist(n-1,PauseLoc))
		else 
		enddelete=str2num(stringfromlist(n,PauseLoc))
		startdelete=	str2num(stringfromlist(n-1,PauseLoc))


		endif
		FreeInfo[(n-1)/2][0]=startdelete
		FreeInfo[(n-1)/2][1]=(enddelete-startdelete)
		deletepoints startdelete, (enddelete-startdelete), FreeForce,FreeSep


	endfor
	for(n=1;n<num-1;n+=2)
		totalDeletions+=FreeInfo[(n-1)/2][1]
		Ind=ChangeStringItem(Ind,-1*totaldeletions,";",n)
		Ind=ChangeStringItem(Ind,-1*totaldeletions,";",n+1)
	endfor
	string NewNote=replacestringbykey("DE_ind",Note(FreeForce),Ind,":","\r")
	note/K FreeForce, NewNote
		note/K FreeSep, NewNote

	duplicate/o FreeForce DicedForce
	duplicate/o FreeSep DicedSep
	duplicate/o FreeInfo CutInfo 
	
end


Static Function/S ChangeStringItem(ListString,ChangeNumber,separator,location)
	String ListString,separator
	variable/D location,ChangeNumber
	variable/D CurrentValue=str2num(stringfromList(location,ListString))
	variable/D newnumber=CurrentValue+changenumber
	return ReplaceListItem(ListString,num2str(newnumber),separator,location)
end


Static Function/S ReplaceListItem(ListString,NewString,separator,location)
	String ListString,NewString,separator
	variable/D location
	string AdjustedString=ListString
	AdjustedString=removelistItem(location,AdjustedString)
		AdjustedString=addlistitem(Newstring,AdjustedString,separator,location)
		return adjustedstring
end

Static Function DetermineWLCParms()
	string saveDF
	saveDF = GetDataFolder(1)
	controlinfo/W=RupRampPanel de_RupRamp_popup0
	SetDataFolder s_value
	controlinfo/W=RupRampPanel de_RupRamp_popup4
	wave ForceWave=$S_value

	controlinfo/W=RupRampPanel de_RupRamp_popup5
	wave StateWave=$S_value


	controlinfo/W=RupRampPanel de_RupRamp_popup16
	wave ForceWaveAlignSM=$S_value
	wave SepWaveAlignSm=$ReplaceString("Force",S_value,"Sep")

	controlinfo/W=RupRampPanel de_RupRamp_popup10
	wave ForceWaveSH_SM=$S_value
	Controlinfo/W=RupRampPanel de_RupRamp_setvar3
	variable/C slopes=DE_Dudko#ReturnSeparationSlopes(SepWaveAlignSm,StateWave,500)
	variable pointstoignore=floor(v_value/real(slopes)/dimdelta(ForceWaveAlignSM,0))
	ForceWaveAlignSM*=-1
	make/o/n=0 $ReplaceString("Force",nameofwave(ForceWave),"_WLCParms")
	wave Results=$ReplaceString("Force",nameofwave(ForceWave),"_WLCParms")
	controlinfo/W=RupRampPanel de_RupRamp_check2
   DE_DUDKO#ContourLengthDetermineCombined(ForceWaveAlignSM,SepWaveAlignSm,StateWave,pointstoignore,Results,v_value,CopyWavesOUt=1)

	ForceWaveAlignSM*=-1
	make/o/n=0 WLC1,WLC2
	DE_DUDKO#MakeWLCs(SepWaveAlignSm,Results,WLC1,WLC2)
end

Static Function MakeOffsetForceWave()
	//handle offsetting both the smoothed wave and the unsmoothed wave, also corrects our estimates of the rupture forces
	string saveDF
	saveDF = GetDataFolder(1)

	controlinfo/W=RupRampPanel de_RupRamp_popup0
	SetDataFolder s_value
	controlinfo/W=RupRampPanel de_RupRamp_popup4
	wave ForceWave=$S_value
	wave SepWave=$ReplaceString("Force",S_value,"Sep")
	controlinfo/W=RupRampPanel de_RupRamp_popup3
	wave UpPoints=$S_value
	wave DownPoints=$ReplaceString("PntU",S_value,"PntD")
	controlinfo/W=RupRampPanel de_RupRamp_popup9
	wave ForceWaveSm=$S_value
	wave SepWaveSm=$ReplaceString("Force",S_value,"Sep")
	ControlInfo/W=RupRampPanel  de_RupRamp_popup5
	wave StateWave=$S_Value
	//duplicate/free ForceWave ForceWaveSm
	//duplicate/free ForceWave SepWaveSm
	//DE_Filtering#FilterForceSep(ForceWave,SepWave,ForceWaveSm,SepWaveSm,"TVD",20e-9)
	//ForceWaveSm*=-1
	duplicate/o ForceWave $ReplaceString("Force_Adj",nameofwave(ForceWave),"Force_Shift")
	wave FwShift=$ReplaceString("Force_Adj",nameofwave(ForceWave),"Force_Shift")	
	
	duplicate/o SepWave $ReplaceString("Force_Adj",nameofwave(ForceWave),"Sep_Shift")
	String Offsets=DE_OverlapRamps#AddForceOffsetstoForceWave(ForceWaveSm,SepWaveSm,StateWave)
	print offsets
	note/K ForceWave,ReplaceStringByKey("DE_FOff", note(ForceWave), Offsets,":","\r" )
	duplicate/free ForceWave FWHold

	DE_OverlapRamps#OffsetEachStepinForceWave(FWHold,FwShift,StateWave)
	DE_OverlapRamps#AddShiftToStates(FwShift,StateWave)
end			
			
Static Function MakeStateWave()
	string saveDF

	saveDF = GetDataFolder(1)
	wave/T parmWave=root:DE_RupRamp:MenuStuff:PyParmWave

	controlinfo/W=RupRampPanel de_RupRamp_popup0
	SetDataFolder s_value
	controlinfo/W=RupRampPanel de_RupRamp_popup9
	wave ForceWaveSm=$S_value
	controlinfo/W=RupRampPanel de_RupRamp_popup4
	wave ForceWave=$S_value
//	wave SepWave=$ReplaceString("Force",S_value,"Sep")
	controlinfo/W=RupRampPanel de_RupRamp_popup3
	wave UpPoints=$S_value
	wave DownPoints=$ReplaceString("PntU",S_value,"PntD")
//
//	duplicate/free ForceWave ForceWaveSm
//	duplicate/free ForceWave SepWaveSm
//
//	DE_Filtering#FilterForceSep(ForceWave,SepWave,ForceWaveSm,SepWaveSm,"TVD",20e-9)
//	ForceWaveSm*=-1;
//	//make the state key
	make/o/n=0 $ReplaceString("Force",nameofwave(ForceWave),"_States")
	wave states=$ReplaceString("Force",nameofwave(ForceWave),"_States")
	DE_DUDKO#MakeSingleStateKey(ForceWaveSm,UpPoints,DownPoints,States)

end
Static Function RefineSteps(ForceWave,UpPnts,DownPnts)
	wave ForceWave,UpPnts,DownPnts
	variable n

	for (n=0;n<numpnts(UpPnts);n+=1)
//for (n=0;n<1;n+=1)
		duplicate/free/r=[UpPnts[n]-0.01/dimdelta(ForceWave,0),UpPnts[n]+.001/dimdelta(ForceWave,0)] ForceWave New
		FindLevels/Q New wavemax(New)
		wave w_findlevels
		UpPnts[n]=x2pnt(ForceWave,W_FindLevels[numpnts(W_FindLevels)-1])
	endfor	
//
//	UpPos=pnt2x(ForceIn, floor(up-1) )
//	UpForce=ForceIn(UpPos)
//	UpSep=SepIn(UpPos)
//	UpPnt=x2pnt(ForceOrig,UpPos)
//
////	
//	FindLevels/Q/EDGE=2 HMMStates .5
//	duplicate/free w_findlevels Down,DownPos,DownForce,DownSep,DownPnt
//	
	for (n=0;n<numpnts(DownPnts);n+=1)
	//for (n=0;n<1;n+=1)
		duplicate/free/r=[DownPnts[n]-.001/dimdelta(ForceWave,0),DownPnts[n]+.01/dimdelta(ForceWave,0)] ForceWave New
		FindLevels/Q New wavemin(New)
		wave w_findlevels
		DownPnts[n]=x2pnt(ForceWave,W_FindLevels[numpnts(W_FindLevels)-1])
	endfor	
//
//	DownPos=pnt2x(ForceIn,floor(down-1))
//	downForce=ForceIn(downPos)
//	DownSep=SepIn(DownPos)
//	DownPnt=x2pnt(ForceOrig,DownPos)
//
//	duplicate/o UpPnt RupPntU
//	duplicate/o DownPnt RupPntD
//	duplicate/o UpForce RupForcesU
//	duplicate/o UpPos RupTimesU
//	duplicate/o DownForce RupForcesD
//	duplicate/o DownPos RupTimesD
	wave W_Findlevels
	killwaves W_FindLevels

end
Static Function SortSteps(ForceWave,Locations,UpPnts,DownPnts,windowsize)
	wave ForceWave,Locations,UpPnts,DownPnts
	variable windowsize
	variable n,q,a,b
	make/free/n=0 Up,Down
	for(n=0;n<numpnts(Locations);n+=1)
		if (n == 0)
			q = min(windowsize, Locations[n+1]-Locations[n])
		elseif(n == (numpnts(Locations)-1))
			q = min(windowsize, Locations[n]-Locations[n-1])
		else
			q =min( min(windowsize,  Locations[n]-Locations[n-1]),  Locations[n+1]- Locations[n])
		endif
		if(mean(ForceWave,pnt2x(ForceWave,Locations[n]-q ),pnt2x(ForceWave,Locations[n]))>mean(ForceWave,pnt2x(ForceWave,Locations[n]),pnt2x(ForceWave,Locations[n]+q)))
			InsertPoints numpnts(Down),1,Down
			Down[numpnts(Down)-1]= Locations[n]
		else
			InsertPoints numpnts(up),1,up
			Up[numpnts(up)-1]= Locations[n]

		endif
		
	endfor
	duplicate/o Up DownPnts	//This looks confusing becuase I messed up my definitions
	duplicate/o Down UpPnts
end

Static Function ListBoxProc(ctrlName,row,col,event) : ListBoxControl
	String ctrlName
	Variable row
	Variable col
	Variable event	//1=mouse down, 2=up, 3=dbl click, 4=cell select with mouse or keys
					//5=cell select with shift key, 6=begin edit, 7=end

	switch(event)

	endswitch				
	
	return 0
End //ListBoxProc


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
Window RuptureRamp_Panel() : Panel

	PauseUpdate; Silent 1		// building window...
	NewPanel/N=RupRampPanel /W=(0,0,600,700)
	NewDataFolder/o root:DE_RupRamp
	NewDataFolder/o root:DE_RupRamp:MenuStuff

	DE_RuptureRamp#UpdateParmWave()
	//Button de_RupRamp_button0,pos={75,80},size={150,20},proc=DE_RuptureRamp#ButtonProc,title="GO!"
	Button de_RupRamp_button1,pos={75,80},size={150,20},proc=DE_RuptureRamp#ButtonProc,title="Pythong!"
	Button de_RupRamp_button2,pos={165,300},size={150,20},proc=DE_RuptureRamp#ButtonProc,title="Make State Wave"
	Button de_RupRamp_button3,pos={320,330},size={150,20},proc=DE_RuptureRamp#ButtonProc,title="OffSet Waves"
	Button de_RupRamp_button4,pos={320,520},size={150,20},proc=DE_RuptureRamp#ButtonProc,title="MakeWLCFits"
	Button de_RupRamp_button5,pos={320,560},size={150,20},proc=DE_RuptureRamp#ButtonProc,title="Correct Ruptures"
	Button de_RupRamp_button6,pos={320,600},size={150,20},proc=DE_RuptureRamp#ButtonProc,title="Remake StateWave"
	Button de_RupRamp_button7,pos={320,240},size={150,20},proc=DE_RuptureRamp#ButtonProc,title="Smooth"
	Button de_RupRamp_button8,pos={320,360},size={150,20},proc=DE_RuptureRamp#ButtonProc,title="Smooth Shifted"
	Button de_RupRamp_button9,pos={165,620},size={150,20},proc=DE_RuptureRamp#ButtonProc,title="DoThisALL"
	Button de_RupRamp_button10,pos={320,420},size={150,20},proc=DE_RuptureRamp#ButtonProc,title="Align To WLC"
	Button de_RupRamp_button11,pos={320,480},size={150,20},proc=DE_RuptureRamp#ButtonProc,title="Smooth Aligned"

	//Button de_RupRamp_button11,pos={100,420},size={150,20},proc=DE_RuptureRamp#ButtonProc,title="Align To Previous"

	PopupMenu de_RupRamp_popup0,pos={75,2},size={129,21},title="Folder",mode=1
	PopupMenu de_RupRamp_popup0,mode=1,popvalue="X",value= #"DE_PanelProgs#ListFolders()"
	PopupMenu de_RupRamp_popup1,pos={75,40},size={129,21},title="Force Wave"
	PopupMenu de_RupRamp_popup1,mode=1,popvalue="X",value= #"DE_RuptureRamp#ListWaves(\"de_RupRamp_popup0\",\"*Force\")"
	PopupMenu de_RupRamp_popup4,pos={20,240},size={129,21},title="Force Wave",mode=0
	PopupMenu de_RupRamp_popup4,mode=1,popvalue="X",value= #"DE_RuptureRamp#ListWaves(\"de_RupRamp_popup0\",\"*Force*Adj\")"
	PopupMenu de_RupRamp_popup2,pos={75,110},size={129,21},title="Method"
	PopupMenu de_RupRamp_popup2,mode=1,popvalue="X",value= "gauss;ms;"

	//ListBox DE_RupRamp_list0,pos={50,110},size={175,150},proc=DE_RuptureRamp#ListBoxProc,listWave=root:DE_RupRamp:MenuStuff:ParmWave
	//ListBox DE_RupRamp_list0,selWave=root:DE_RupRamp:MenuStuff:SelWave,editStyle= 2,userColumnResize= 1,widths={70,40,70,40}
	ListBox DE_RupRamp_list1,pos={50,140},size={175,75},proc=DE_RuptureRamp#ListBoxProc,listWave=root:DE_RupRamp:MenuStuff:PyParmWave
	ListBox DE_RupRamp_list1,selWave=root:DE_RupRamp:MenuStuff:PySelWave,editStyle= 2,userColumnResize= 1,widths={70,40,70,40}
	DrawLine 11,234,292,234
	PopupMenu de_RupRamp_popup3,pos={310,270},size={129,21},title="RuptureUp"
	PopupMenu de_RupRamp_popup3,mode=1,popvalue="X",value= #"DE_RuptureRamp#ListWaves(\"de_RupRamp_popup0\",\"*RupPntU*Adj\")"
	PopupMenu de_RupRamp_popup5,pos={20,330},size={129,21},title="StateWave"
	PopupMenu de_RupRamp_popup5,mode=1,popvalue="X",value= #"DE_RuptureRamp#ListWaves(\"de_RupRamp_popup0\",\"*State*\")"
	PopupMenu de_RupRamp_popup6,pos={20,360},size={129,21},title="Shifted Wave"
	PopupMenu de_RupRamp_popup6,mode=1,popvalue="X",value= #"DE_RuptureRamp#ListWaves(\"de_RupRamp_popup0\",\"*Force*Shift*\")"
	PopupMenu de_RupRamp_popup7,pos={20,560},size={129,21},title="WLC Parms"
	PopupMenu de_RupRamp_popup7,mode=1,popvalue="X",value= #"DE_RuptureRamp#ListWaves(\"de_RupRamp_popup0\",\"*WLCParms*\")"
	PopupMenu de_RupRamp_popup8,pos={20,600},size={129,21},title="New Up Ruptures"
	PopupMenu de_RupRamp_popup8,mode=1,popvalue="X",value= #"DE_RuptureRamp#ListWaves(\"de_RupRamp_popup0\",\"*RupPntU*\")"
	PopupMenu de_RupRamp_popup9,pos={20,270},size={129,21},title="Smoothed Wave"
	PopupMenu de_RupRamp_popup9,mode=1,popvalue="X",value= #"DE_RuptureRamp#ListWaves(\"de_RupRamp_popup0\",\"*Force*Sm\")"
	PopupMenu de_RupRamp_popup10,pos={20,420},size={129,21},title="Smoothed Shifted Wave"
	PopupMenu de_RupRamp_popup10,mode=1,popvalue="X",value= #"DE_RuptureRamp#ListWaves(\"de_RupRamp_popup0\",\"*Force*Shift*Sm\")"
	PopupMenu de_RupRamp_popup11,pos={480,240},size={129,21},title="Type",proc=DE_RuptureRamp#PopMenuProc
	PopupMenu de_RupRamp_popup11,mode=1,popvalue="X",value= "TVD;SVG"
	PopupMenu de_RupRamp_popup12,pos={480,360},size={129,21},title="Type",proc=DE_RuptureRamp#PopMenuProc
	PopupMenu de_RupRamp_popup12,popvalue="X",value= "TVD;SVG",mode=1
	PopupMenu de_RupRamp_popup15,pos={480,480},size={129,21},title="Type",proc=DE_RuptureRamp#PopMenuProc
	PopupMenu de_RupRamp_popup15,popvalue="X",value= "TVD;SVG",mode=1
	PopupMenu de_RupRamp_popup13,pos={50,450},size={129,21},title="Alignment Parms"
	PopupMenu de_RupRamp_popup13,mode=1,popvalue="X",value= #"DE_RuptureRamp#ListWaves(\"de_RupRamp_popup0\",\"*WLCAlign*\")"
	
	PopupMenu de_RupRamp_popup17,pos={270,450},size={129,21},title="Alignment Type"
	PopupMenu de_RupRamp_popup17,mode=1,popvalue="X",value="Both;Unfolded;Folded"
	
	PopupMenu de_RupRamp_popup14,pos={20,480},size={129,21},title="Aligned Wave"
	PopupMenu de_RupRamp_popup14,mode=1,popvalue="X",value= #"DE_RuptureRamp#ListWaves(\"de_RupRamp_popup0\",\"*Force*Align\")"
	PopupMenu de_RupRamp_popup16,pos={20,520},size={129,21},title="Aligned Smoothed Wave"
	PopupMenu de_RupRamp_popup16,mode=1,popvalue="X",value= #"DE_RuptureRamp#ListWaves(\"de_RupRamp_popup0\",\"*Force*Align*Sm\")"
	//PopupMenu de_RupRamp_popup14,pos={100,450},size={129,21},title="Alignmentto:"
	//PopupMenu de_RupRamp_popup14,mode=1,popvalue="X",value= #"DE_RuptureRamp#ListWaves(\"de_RupRamp_popup0\",\"*Force*Align*\")"
	
	
	SetVariable de_RupRamp_setvar0,pos={547.00,240.00},size={50.00,18.00},proc=DE_RuptureRamp#SetVarProc,value=_num:10e-9
	SetVariable de_RupRamp_setvar0,limits={0,inf,0}

	SetVariable de_RupRamp_setvar1,pos={547.00,360.00},size={50.00,18.00},proc=DE_RuptureRamp#SetVarProc,value=_num:10e-9
	SetVariable de_RupRamp_setvar1,limits={0,inf,0}

	SetVariable de_RupRamp_setvar2,pos={547.00,480.00},size={50.00,18.00},proc=DE_RuptureRamp#SetVarProc,value=_num:10e-9
	SetVariable de_RupRamp_setvar2,limits={0,inf,0}
	
	SetVariable de_RupRamp_setvar3,pos={470.00,520.00},size={99.00,18.00},proc=DE_RuptureRamp#SetVarProc,value=_num:3e-9
	SetVariable de_RupRamp_setvar3,limits={0,inf,0},title="Distance"

	CheckBox de_RupRamp_check0 title="All Sep Shift",pos={450,450},size={150,25},proc=DE_RuptureRamp#CheckProc

	//CheckBox de_RupRamp_check1 title="Fit Folded",pos={550,450},size={150,25},proc=DE_RuptureRamp#CheckProc
	CheckBox de_RupRamp_check2 title="Fit Folded",pos={500,550},size={150,25},proc=DE_RuptureRamp#CheckProc

	ControlUpdate/A/W=RupRampPanel
EndMacro

Static Function UpdateParmWave()
	//	if(exists("root:DE_RupRamp:MenuStuff:ParmWave")==1)
	//		wave/t/z Par=root:DE_RupRamp:MenuStuff:ParmWave
	//		wave/z Sel=root:DE_RupRamp:MenuStuff:SelWave
	//	Else
	//		make/t/n=(8,2) root:DE_RupRamp:MenuStuff:ParmWave
	//		wave/t/z Par=root:DE_RupRamp:MenuStuff:ParmWave
	//		make/n=(8,2) root:DE_RupRamp:MenuStuff:SelWave
	//		wave/z Sel=root:DE_RupRamp:MenuStuff:SelWave
	//		
	//		Par[0][0]={"Number of States","Number of Modes","Drift Bound (nm)","Sd. Deviation (nm)","Transition Bound","Iterations","Smoothing","Decimating"}
	//		Par[0][1]={"2","4",".5",".2",".5","3","50e-9","10"}
	//		Sel[][0]=0
	//		Sel[][1]=2
	//	endif
	
	if(exists("root:DE_RupRamp:MenuStuff:PyParmWave")==1)
		wave/t/z Par=root:DE_RupRamp:MenuStuff:ParmWave
		wave/z Sel=root:DE_RupRamp:MenuStuff:SelWave
	Else
		make/t/n=(4,2) root:DE_RupRamp:MenuStuff:PyParmWave
		wave/t/z Par=root:DE_RupRamp:MenuStuff:PyParmWave
		make/n=(4,2) root:DE_RupRamp:MenuStuff:PySelWave
		wave/z Sel=root:DE_RupRamp:MenuStuff:PySelWave
		
		Par[0][0]={"Smoothing","Decimating","GaussParm","Threshhold"}
		Par[0][1]={"50e-9","1","5","0.1"}
		Sel[][0]=0
		Sel[][1]=2
	endif


end

Static Function/S ListWaves(ControlStr,SearchString)
	string ControlStr,SearchString
	String saveDF

	saveDF = GetDataFolder(1)
	controlinfo $ControlStr
	SetDataFolder s_value
	String list = WaveList(SearchString, ";", "")
	SetDataFolder saveDF
	return list

end

Menu "Ramp"
	//SubMenu "Processing"
	"Open RuptureForce", RuptureRamp_Panel()
	"Open Viewer", MultiRampViewer()


	//end
	
end

//Static Function DriftMarkovFitter( UseWave, stateCount, modeCount, timeStep, driftBound, sigmaBound, transitionBound, iterationCount, [RAM, Threads])//Variables demanded by MarkovFit Code
//	Wave UseWave//Input wave
//	Variable stateCount, modeCount, timeStep, driftBound, sigmaBound, iterationCount, RAM, Threads
//	Variable TransitionBound
//	killwaves /z HidMar0, HidMar1, HidMar2, HidMar3, HidMar4, usable//Getting rid of generated waves to generate new ones
//	RAM = paramIsDefault(RAM) ? 4:RAM
//	Threads = paramIsDefault(Threads) ? 1000:Threads
//	killwaves /z Used
//	duplicate /o UseWave Used
//	if(timeStep==0)
//		timestep = 1.0
//	endif
//	Variable hold
//	if(iterationCount==0)
//		iterationCount = 4
//	endif
//	if(RAM == 0)
//		RAM = 4
//	endif
//	if(modeCount ==0)
//		Variable i
//		for(i=0;i<numpnts(Used); i+=1)
//			Used[i] += -driftBound*i
//		endfor
//	endif
//	String InfoPass = "java -Xmx" + num2str(RAM) +"g -jar C:\MarkovFitter\DriftMarkov2.jar C:\MarkovFitter\UseWave.txt " + num2str(stateCount)+" 0 "//infopass exists to hold the command sent to DOS
//	InfoPass = InfoPass + num2str(modeCount)+" "+num2str(timeStep)+" "+num2str(driftBound)+" "+num2str(sigmaBound)+" "+num2str(transitionBound)+" "+num2str(iterationCount)+" "+num2str(Threads)
//	Save/J/W Used as "C:\MarkovFitter\UseWave.txt"//saving the wave that was given to  proper location
//	print(InfoPass)//gives view of command line in case anything is wrong
//	executescripttext InfoPass//sendng command to command line
//	LoadWave/A=HidMar/J/D/W/K=0 "C:MarkovFitter:DriftMarkovOut.txt"//getting waves from location jar tosses them to(waves have base name HidMar
//	//Display UseWave//displaying wave given
//	variable Temp
//	duplicate/o $"HidMar1" usable//while wave1 is created through this code it cannot regonize it so it must be duplicated
//	Temp =dimoffset(UseWave,0)
//	setscale/P x dimoffset(UseWave,0), dimdelta(UseWave,0), "s", usable//ensuring scaling of input and output wave are the same
//	if(modeCount ==0)
//		for(i=0;i<numpnts(UseWave);i+=1)
//			usable[i] += driftBound*i
//		endfor
//	endif
//	//AppendToGraph usable//putting on same graph
//	//ModifyGraph rgb(usable)=(0,0,65280)//changing color so both waves are visible
//	//display $"HidMar2"//displaying simple jump wave
//	killwaves usable, used
//	executescripttext "java -jar C:\MarkovFitter\GetRidOfUseWave.jar"//Eliminates file created earlier to prevent problems on future runs
//end

//Static Function FindStateChanges(HMMStates,ForceOrig,Forcein,SepIn,RupForcesU,RupTimesU,RupPntU,RupForcesD,RupTimesD,RupPntD)
//	wave HMMStates,ForceOrig,Forcein,SepIn,RupForcesU,RupTimesU,RupForcesD,RupTimesD,RupPntU,RupPntD
//	
//	FindLevels/Q/EDGE=1 HMMStates .5
//	wave w_findlevels
//	duplicate/free w_findlevels Up, UpPos,UpForce,UpSep,UpPnt
//	
//	variable n
//	for (n=0;n<numpnts(Up);n+=1)
//	//for (n=0;n<1;n+=1)
//		duplicate/free/r=(Up[n]-.03,Up[n]+.01) ForceIn New
//		FindLevels/Q New wavemax(New)
//		wave w_findlevels
//		Up[n]=x2pnt(ForceIn,W_FindLevels[numpnts(W_FindLevels)-1])
//	endfor	
//
//	UpPos=pnt2x(ForceIn, floor(up-1) )
//	UpForce=ForceIn(UpPos)
//	UpSep=SepIn(UpPos)
//	UpPnt=x2pnt(ForceOrig,UpPos)
//
////	
//	FindLevels/Q/EDGE=2 HMMStates .5
//	duplicate/free w_findlevels Down,DownPos,DownForce,DownSep,DownPnt
//	
//	for (n=0;n<numpnts(Down);n+=1)
//	//for (n=0;n<1;n+=1)
//		duplicate/free/r=(Down[n]-.03,Down[n]+.01) ForceIn New
//		FindLevels/Q New wavemin(New)
//		wave w_findlevels
//		Down[n]=x2pnt(ForceIn,W_FindLevels[numpnts(W_FindLevels)-1])
//	endfor	
//
//	DownPos=pnt2x(ForceIn,floor(down-1))
//	downForce=ForceIn(downPos)
//	DownSep=SepIn(DownPos)
//	DownPnt=x2pnt(ForceOrig,DownPos)
//
//	duplicate/o UpPnt RupPntU
//	duplicate/o DownPnt RupPntD
//	duplicate/o UpForce RupForcesU
//	duplicate/o UpPos RupTimesU
//	duplicate/o DownForce RupForcesD
//	duplicate/o DownPos RupTimesD
//	wave W_Findlevels
//	killwaves W_FindLevels
//	
//end

//Static Function FindStateChangesSimple(HMMStates,ForceOrig,Forcein,SepIn,RupPntU,RupPntD)
//	wave HMMStates,ForceOrig,Forcein,SepIn,RupPntU,RupPntD
//	
//	FindLevels/Q/EDGE=1 HMMStates .5
//	wave w_findlevels
//	duplicate/free w_findlevels Up, UpPos,UpPnt
////	
//	variable n
//	for (n=0;n<numpnts(Up);n+=1)
//	//for (n=0;n<1;n+=1)
//		duplicate/free/r=(Up[n]-.03,Up[n]+.01) ForceIn New
//			FindLevels/Q New wavemax(New)
//		wave w_findlevels
//		Up[n]=x2pnt(ForceIn,W_FindLevels[numpnts(W_FindLevels)-1])
//	endfor	
//
//	UpPos=pnt2x(ForceIn, floor(up-1) )
//	//UpForce=ForceIn(UpPos)
//	//UpSep=SepIn(UpPos)
//	UpPnt=x2pnt(ForceOrig,UpPos)
//
////	
//	FindLevels/Q/EDGE=2 HMMStates .5
//	duplicate/free w_findlevels Down,DownPos,DownPnt
////	
//	for (n=0;n<numpnts(Down);n+=1)
//	//for (n=0;n<1;n+=1)
//		duplicate/free/r=(Down[n]-.03,Down[n]+.01) ForceIn New
//		FindLevels/Q New wavemin(New)
//		wave w_findlevels
//		Down[n]=x2pnt(ForceIn,W_FindLevels[numpnts(W_FindLevels)-1])
//	endfor	
//
//	DownPos=pnt2x(ForceIn,floor(down-1))
//	//downForce=ForceIn(downPos)
//	//DownSep=SepIn(DownPos)
//	DownPnt=x2pnt(ForceOrig,DownPos)
//
//	duplicate/o UpPnt RupPntU
//	duplicate/o DownPnt RupPntD
//	wave W_Findlevels
//	killwaves W_FindLevels
//	
//end

Static Function FitHMMButt()

				//					saveDF = GetDataFolder(1)
					//					wave/T parmWave=root:DE_RupRamp:MenuStuff:ParmWave
					//
					//					controlinfo de_RupRamp_popup0
					//					SetDataFolder s_value
					//					controlinfo de_RupRamp_popup1
					//					wave ForceWave=$S_value
					//					wave SepWave=$ReplaceString("Force",S_value,"Sep")
					//					//			variable FOffset
					//
					//					make/o/n=0 RupPntU,RupPntD
					//		
					//					make/o/n=0 ForceWaveS,SepWaveS
					//					//DE_Filtering#FilterForceSep(ForceWave,SepWave,ForceWaveS,SepWaveS,"SVG",str2num(parmWave[6][1]))
					//					DE_Filtering#FilterForceSep(ForceWave,SepWave,ForceWaveS,SepWaveS,"TVD",str2num(parmWave[6][1]))
					//					ForceWaveS*=-1
					//					wavestats/q ForceWaveS
					//					FOffset=v_min
					//					FOffset-=0e-12
					//					ForceWaveS-=FOffset
					//					Soffset= SepWaveS[0]-0e-9
					//					SepWaveS-=Soffset
					//					duplicate/o ForceWaveS LC
					//					LC=DE_WLC#ContourTransform(ForceWaveS,SepWaveS,.4e-9,298)		
					//					LC*=1e9	
					//					setscale/P x dimoffset(ForceWave,0), dimdelta(ForceWave,0), "s", LC//ensuring scaling of input and output wave are the same
					//					Resample/DOWN=(str2num(parmWave[7][1]))/N=1/WINF=None LC
					//
					//					DriftMarkovFitter(LC, str2num(parmWave[0][1]), str2num(parmWave[1][1]), dimdelta(LC,0),str2num(parmWave[2][1]),str2num(parmWave[3][1]), str2num(parmWave[4][1]), str2num(parmWave[5][1]))
					//
					//					wave HidMar0,HidMar1,HidMar2,HidMar3,HidMar4
					//					setscale/P x dimoffset(LC,0), dimdelta(LC,0), "s", HidMar2//ensuring scaling of input and output wave are the same
					//
					//					make/o/n=0 RupForcesU,RupTimesU,RupForcesD,RupTimesD,RupPntU,RupPntD
					//					FindStateChangesSimple(HidMar2,ForceWave,ForceWaveS,SepWaveS,RupPntU,RupPntD)
					//					RupForcesU+=FOffset
					//					RupForcesD+=FOffset
					//					ForceWaveS+=FOffset
					//
					//					duplicate/o hidmar0 Data
					//					duplicate/o hidmar1 Fit
					//					Killwaves HidMar0,HidMar1,HidMar2,HidMar3,HidMar4
					//					killwaves SepWaveS,ForceWaveS,LC
					//					//killwaves Data,Fit
					//					SetDataFolder saveDF
					
					
					end