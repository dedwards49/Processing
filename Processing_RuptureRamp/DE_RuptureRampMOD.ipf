#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma modulename=DE_RuptureRamp
#include "DE_Filtering"
#include "SimpleWLCPrograms"
#include "C:Users:Perkins Lab:src_prh:IgorUtil:IgorCode:Util:OperatingSystemUtil"
#include "DTE_Dudko"
#include "DE_OverlapRamps"
#include "DE_CorrectRupture"
#include "DE_TwoWLCFit"
#include "DE_CountRates"

#include ":\Misc_PanelPrograms\Panel Progs"
#include ":\Misc_PanelPrograms\AsylumNaming"
//#include ":\Processing_Markov\DE_HMM"
//	SetDataFolder ReturnPanelString("Folder")
//	wave AlignedForceWave=$ReturnPanelString("AlignedForceWave")
//	wave AlignedSepWave=$ReplaceString("Force",nameofwave(AlignedForceWave),"Sep")
//	string SmType=ReturnPanelString("AlignSmType")
//	variable SmAmount=ReturnPanelVal("AlignSmval")
//	string saveDF
//	saveDF = GetDataFolder(1)
//	
//	
//	controlinfo/W=RupRampPanel de_RupRamp_popup0
//	SetDataFolder s_value
//	controlinfo/W=RupRampPanel de_RupRamp_popup1
//	wave ForceWave=$S_value
//	wave SepWave=$ReplaceString("Force",S_value,"Sep")
//	controlinfo/W=RupRampPanel de_RupRamp_popup2
//	wave UpPoints=$S_value
//	wave DownPoints=$ReplaceString("PntU",S_value,"PntD")
//	controlinfo/W=RupRampPanel de_RupRamp_popup3
//	wave ForceWaveSm=$S_value
//	wave SepWaveSm=$ReplaceString("Force",S_value,"Sep")
//	ControlInfo/W=RupRampPanel  de_RupRamp_popup5
//	wave StateWave=$S_Value
//	
Static Function CutStatebySep(ForceIn,SepIn,ForceOut,SepOut,SepMax,SepMin)

	wave ForceIn,SepIn,ForceOut,SepOut
	variable SepMax,SepMin
	duplicate/free ForceIn FreeForce
	duplicate/free SepIn FreeSep
	variable n
	string NewNote
	if(SepMax!=0&&SepMax!=0)
		for(n=0;n<numpnts(FreeSep);n+=1)
			if(FreeSep[n]<SepMin||FreeSep[n]>SepMax)
				FreeSep[n]=NaN
				FreeForce[n]=NaN
			
			endif
			NewNote=Replacestringbykey("DE_SepMax",note(ForceIn),num2str(SepMax),":","\r")
			NewNote=Replacestringbykey("DE_SepMin",NewNote,num2str(SepMin),":","\r")

		endfor
	
	elseif(SepMax!=0)
		for(n=0;n<numpnts(FreeSep);n+=1)
			if(FreeSep[n]>SepMax)
				FreeSep[n]=NaN
				FreeForce[n]=NaN

			endif

		endfor
		NewNote=Replacestringbykey("DE_SepMax",note(ForceIn),num2str(SepMax),":","\r")

	
	elseif(SepMin!=0)
		for(n=0;n<numpnts(FreeSep);n+=1)
			if(FreeSep[n]<SepMin)
				FreeSep[n]=NaN
				FreeForce[n]=NaN

			endif
		
		endfor
		NewNote=Replacestringbykey("DE_SepMin",note(ForceIn),num2str(SepMin),":","\r")

	endif
	duplicate/o FreeSep SepOut
	duplicate/o FreeForce ForceOut

end

Static Function MakeANicePlot(Type)
	String Type
	string saveDF
	saveDF = GetDataFolder(1)

	//Here we just grab a long list of all the waves we need	
	controlinfo/W=RupRampPanel de_RupRamp_popup0
	SetDataFolder s_value

	StrSwitch (Type)
	
		Case "Align":
			dowindow AlignMent
			if(v_flag==1)
				killwindow Alignment
			endif
			wave AlignUnfolded,Alignfolded,AlignFitShift,AlignFitNoShift,AlignFitOriginal
			Display/N=Alignment AlignUnfolded[][0] vs AlignUnfolded[][1]
			Appendtograph/W=Alignment Alignfolded[][0] vs Alignfolded[][1]
			Appendtograph/W=Alignment AlignFitShift[][0] vs AlignFitShift[][1]
			Appendtograph/W=Alignment AlignFitNoShift[][0] vs AlignFitNoShift[][1]
			Appendtograph/W=Alignment AlignFitOriginal[][0] vs AlignFitOriginal[][1]
			ModifyGraph/W=Alignment rgb(AlignUnFolded)=(14135,32382,47288),rgb(AlignFolded)=(19789,44975,19018)
			ModifyGraph/W=Alignment rgb(AlignFitShift)=(29524,1,58982),rgb(AlignFitNoShift)=(58596,6682,7196),rgb(AlignFitOriginal)=(0,0,0)
		
			ModifyGraph/W=Alignment lsize(AlignFitShift)=2,lsize(AlignFitNoShift)=2,lsize(AlignFitOriginal)=2
		break
	
	default:
	endswitch
	SetDataFolder saveDF

end
Static Function MakeAPlainRupWave()
	

	string saveDF
	variable FOffset, Soffset

	saveDF = GetDataFolder(1)

	controlinfo/W=RupRampPanel de_RupRamp_popup0
	SetDataFolder s_value
	controlinfo/W=RupRampPanel de_RupRamp_popup1
	wave ForceWave=$S_value
	
	String DE_Ind=stringbykey("DE_Ind",note(ForceWave),":","\r")
	string DE_Pause=stringbykey("DE_PauseLoc",note(ForceWave),":","\r")
	variable Num=ItemsInList(DE_Ind)/2
	variable n
	variable tenth=str2num(stringfromlist(0,DE_Pause))/10
	make/o/n=(Num) RupPntU,RupPntD
	RupPntU[0]=tenth
	RupPntD[0]=str2num(stringfromlist(0,DE_Pause))-tenth
	for(n=1;n<num;n+=1)
		RupPntU[n]=str2num(stringfromlist(2*n-1,DE_Pause))+tenth
		RupPntD[n]=str2num(stringfromlist(2*n,DE_Pause))-tenth
	
	
	
	endfor

	SetDataFolder saveDF

end

Function PlotFromFolder(FolderString,ShiftString,GraphString)

	string FolderString,ShiftString,GraphString

	wave ForceWave=$(FolderString+":"+stringfromlist(0,DE_CountRates#ListWaves(FolderString,"*Force_Align_sm")))
	wave SepWave=$(FolderString+":"+stringfromlist(0,DE_CountRates#ListWaves(FolderString,"*Sep_Align")))

	variable/c shift=DE_CountRates#CorrectShift(ForceWave,ShiftString)
	variable 	ForceAdj=real(Shift)
	variable 	SepAdj=imag(Shift)
	if(whichListItem(nameofwave(ForceWave),tracenamelist(graphstring,";",1))>0)
	do
		removefromGraph/W=$Graphstring $nameofwave(ForceWave)
	
	while(whichListItem(nameofwave(ForceWave),tracenamelist(graphstring,";",1))>0)
	endif
	appendtograph/W=$Graphstring ForceWave vs SepWave
	ModifyGraph/W=$GraphString offset($nameofwave(ForceWave))={SepAdj,ForceAdj}

end

Function PlotHelp(GraphString,Type)
	String GraphString,Type
	string Traces= tracenamelist(GraphString,";",1)
	variable tot=itemsinlist(Traces)
	variable n
	String UsedAllowed
	variable ForceAdj,SepAdj,ShiftedFOff,ShiftedSOff,UnShiftedFOff,UnShiftedSOff,UsedFoff,UsedSoff,AltFOff,AltSOff
	for(n=0;n<tot;n+=1)
		wave w1=TraceNameToWaveRef(GraphString, stringfromlist(n,Traces) )
		if(strsearch(nameofwave(w1),"Align",0)==-1)
		
		else
				UsedFoff=str2num(Stringbykey("UsedAlignmentFShift",note(w1),":","\r"))
				UsedSOff=str2num(Stringbykey("UsedAlignmentSShift",note(w1),":","\r"))
				AltFOff=str2num(Stringbykey("AltAlignmentFShift",note(w1),":","\r"))
				AltSOff=str2num(Stringbykey("AltAlignmentSShift",note(w1),":","\r"))
			if(str2num(Stringbykey("UsedAlignmentSShift",note(w1),":","\r"))==0)
				UsedAllowed="No"
				ShiftedFOff=str2num(Stringbykey("AltAlignmentFShift",note(w1),":","\r"))
				ShiftedSOff=str2num(Stringbykey("AltAlignmentSShift",note(w1),":","\r"))
				UnShiftedFOff=str2num(Stringbykey("UsedAlignmentFShift",note(w1),":","\r"))
				UnShiftedSOff=str2num(Stringbykey("UsedAlignmentSShift",note(w1),":","\r"))
				
			else
				UsedAllowed="Yes"
				unShiftedFOff=str2num(Stringbykey("AltAlignmentFShift",note(w1),":","\r"))
				unShiftedSOff=str2num(Stringbykey("AltAlignmentSShift",note(w1),":","\r"))
				ShiftedFOff=str2num(Stringbykey("UsedAlignmentFShift",note(w1),":","\r"))
				ShiftedSOff=str2num(Stringbykey("UsedAlignmentSShift",note(w1),":","\r"))
			endif
			
			Strswitch(Type)
				case "Used":
					ForceAdj=0
					SepAdj=0
					break
				case "Shifted":
					ForceAdj=UsedFoff-ShiftedFOff
					SepAdj=UsedSOff-ShiftedSOff
				break	
				case "UnShifted":
					ForceAdj=UsedFoff-unShiftedFOff
					SepAdj=UsedSOff-UnShiftedSOff
				break	
				case "Alt":
					ForceAdj=UsedFoff-AltFOff
					SepAdj=UsedSOff-AltSOff
				break	
				default:
					ForceAdj=0
					SepAdj=0
			endswitch
				ModifyGraph/W=$GraphString offset($nameofwave(w1))={SepAdj,-ForceAdj}
		endif
		
	endfor
	

end

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
					MakeSmoothedWave()
					break
				case  "de_RupRamp_button1":
					MakeStateWave()
					break
				case  "de_RupRamp_button2":
					MakeOffsetForceWave()
					break
				case  "de_RupRamp_button3":
					MakeSmoothedShWave()
					break
				case  "de_RupRamp_button4":
					DisplayAPlot("Shifted")
					break
				case  "de_RupRamp_button5":
					CrunchAbove()
					break
				case  "de_RupRamp_button6":
					AligntoWLC()
					break
				case  "de_RupRamp_button7":
					MakeSmoothedAlignedWave()
					break
				case  "de_RupRamp_button8":
					DisplayAPlot("Overlap")
					break
				case  "de_RupRamp_button9":
					PullShifts()
					break
				case  "de_RupRamp_button10":
					ApplyShifts()
					break
				
				case  "de_RupRamp_button11":
					MakeSmoothedFinalWave()
					break
				case  "de_RupRamp_button12":
					CrunchMiddle()
					break
				case "de_RupRamp_button13":
					DetermineWLCParms()
					break
				case "de_RupRamp_button14":
					CorrectRuptures()
					break
				case "de_RupRamp_button15":
					MakeSecondStateWave()
					break
				case "de_RupRamp_button16":
					CrunchBottom()
					break
				case "de_RupRamp_button17":
					CrunchAll()
					break
				case "de_RupRamp_button18":
					RunASingle()
					break
		
			endswitch
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Static Function PullShifts()

	DoWindow Overlap
	if(V_flag==0)
		return -1
	endif

	String WaveNameString=ReturnPanelString("SmoothedAlignedForce")
	String TracesString=TraceNameList("Overlap", ";",1)
	if(WhichListItem(WaveNameString,TracesString)==-1)
		return 0
	endif
	String Offsets= (StringByKey("offset(x)",TraceInfo("Overlap", WaveNameString, 0 ),"=",";"	))
	variable FOff,SOff
	sscanf Offsets, "{%f,%f}",Soff,FOff
	SetVariable de_RupRamp_setvar6,value=_num:FOff
	SetVariable de_RupRamp_setvar7,value=_num:Soff

end

Static Function ApplyShifts()
	
	string saveDF
	saveDF = GetDataFolder(1)
	SetDataFolder ReturnPanelString("Folder")
	
	Wave AlignedForceWave=$ReturnPanelString("AlignedForceWave")
	Wave AlignedSepWave=$replacestring("Force",nameofwave(AlignedForceWave),"Sep")
	duplicate/free  AlignedForceWave FreeForce
	duplicate/free  AlignedSepWave FreeSep
	variable Foff=ReturnPanelval("Fshift")
	variable Soff=ReturnPanelval("Sshift")
	FreeForce+=Foff
	FreeSep+=Soff
	duplicate/o FreeForce $replaceString("Align",nameofwave(AlignedForceWave),"Final")
	duplicate/o FreeSep $replaceString("Align",nameofwave(AlignedSepWave),"Final")

	SetDataFolder saveDF

end


Static Function DisplayAPLot(StringCon)

	string StringCon
	string saveDF
	saveDF = GetDataFolder(1)
	
	SetDataFolder ReturnPanelString("Folder")
//	wave AlignedSepWave=$ReplaceString("Force",nameofwave(AlignedForceWave),"Sep")
//	string SmType=ReturnPanelString("AlignSmType")
//	variable SmAmount=ReturnPanelVal("AlignSmval")

	strswitch(StringCon)
	
		case "Shifted":

			controlinfo/W=RupRampPanel de_RupRamp_popup7
			wave SmoothShiftedForce=$S_value
			wave SmoothShiftedSep=$ReplaceString("Force",S_value,"Sep")
			DoWindow $StringCon
			if(V_flag==1)
				killwindow $StringCon
				
			endif
			display/N=$StringCon/W=(550,25,850,200) SmoothShiftedForce vs SmoothShiftedSep

			SetDataFolder saveDF
		
			break 
			
		case "Overlap":
			wave SmoothedAlignedForce=$ReturnPanelString("SmoothedAlignedForce")
			wave SmoothedAlignedSep=$ReplaceString("Force",nameofwave(SmoothedAlignedForce),"Sep")

			DoWindow/F $StringCon
			if(V_flag==1)
				
				if(whichlistitem(nameofwave(SmoothedAlignedForce),tracenamelist(StringCon,";",1))!=-1)
				return 0
				
				endif
				appendtograph/W=$StringCon SmoothedAlignedForce vs SmoothedAlignedSep
				
			else
				if(cmpstr("X",ReturnPanelString("OverlapForce"))==0)
					display/N=$StringCon/W=(550,25,850,200) SmoothedAlignedForce vs SmoothedAlignedSep
				else
					String AlignParm=ReturnPanelString("AlignFolder")
					wave OverlapForce=$(AlignParm+ReturnPanelString("OverlapForce"))
					wave OverlapSep=$ (GetWavesDataFolder(OverlapForce, 0 )+":"+ReplaceString("FSm",nameofwave(OverlapForce),"SSm"))
					display/N=$StringCon/W=(550,25,850,200) OverlapForce vs OverlapSep
					ModifyGraph/W=$StringCon rgb($nameofwave(OVerlapforce))=(0,0,0)
					Appendtograph/W=$StringCon SmoothedAlignedForce vs SmoothedAlignedSep

				endif
				
			endif

			break 
	
	endswitch
	SetDataFolder saveDF

end

Static Function RunASingle()

	string saveDF
	saveDF = GetDataFolder(1)
	SetDataFolder ReturnPanelString("Folder")
	wave UpPoints=$ReturnPanelString("RupPnts")
	wave DownPoints=$ReplaceString("PntU",nameofwave(UpPoints),"PntD")
	String Direc=ReturnPanelString("SingleDirec")
	String Type=ReturnPanelString("SingleType")

	if(cmpstr("Up",Direc)==0)
		wave PntWave=UpPoints
	elseif(cmpstr("Down",Direc)==0)
		wave PntWave=DownPoints

	endif	
	variable index=ReturnPanelVal("SingleNum")
	DE_CorrRup#FixPicksSingle(PntWave,index,101,Type)
	SetDataFolder saveDF

end





Static Function/S ReturnPanelString(NameString)
	String NameString
	
	String PanelName="RupRampPanel"
	String ControlName
	StrSwitch(NameString)
		case "Folder":
			ControlName="de_RupRamp_popup0"
			break
		case "ForceWave":
			ControlName="de_RupRamp_popup1"
			break
		case "RupPnts":
			ControlName="de_RupRamp_popup2"
			break
		case "SmoothedForce":
			ControlName="de_RupRamp_popup3"
			break
		case "ForceSmType":
			ControlName="de_RupRamp_popup4"
			break
		case "StateWave":
			ControlName="de_RupRamp_popup5"
			break
		case "ShiftedForce":
			ControlName="de_RupRamp_popup6"
			break
		case "SmoothedShiftedForce":
			ControlName="de_RupRamp_popup7"
			break
		case "ShiftSmType":
			ControlName="de_RupRamp_popup8"
			break
		case "AlignFolder":
			ControlName="de_RupRamp_popup9"
			break
		case "AlignParms":
			ControlName="de_RupRamp_popup10"
			break
		case "AlignType":
			ControlName="de_RupRamp_popup11"
			break		
		case "AlignedForceWave":
			ControlName="de_RupRamp_popup12"
			break
		case "AlignSmtype":
			ControlName="de_RupRamp_popup13"
			break
		case "SmoothedAlignedForce":
			ControlName="de_RupRamp_popup14"
			break
		case "OverlapForce":
			ControlName="de_RupRamp_popup15"
			break
		case "FinalForceWave":
			ControlName="de_RupRamp_popup16"
			break
		case "SmoothedFinalForce":
			ControlName="de_RupRamp_popup17"
			break
		case "FinalSmType":
			ControlName="de_RupRamp_popup18"
			break
		case "WLCParms":
			ControlName="de_RupRamp_popup19"
			break
		case "NewRupPnts":
			ControlName="de_RupRamp_popup20"
			break
		case "SingleDirec":
			ControlName="de_RupRamp_popup21"
			break
		case "SingleType":
			ControlName="de_RupRamp_popup22"
			break
		default:
			return ""
	endswitch



	Controlinfo/W=$PanelName $ControlName
	return S_Value
end

Static Function ReturnPanelVal(ValString)
	String ValString

	String PanelName="RupRampPanel"
	String ControlName
	StrSwitch(ValString)
		case "ForceSmVal":
			ControlName="de_RupRamp_setvar0"
			break
		case "ShiftSmVal":
			ControlName="de_RupRamp_setvar1"
			break
		case "IgnoreDist":
			ControlName="de_RupRamp_setvar2"
			break	
		case "SepMin":
			ControlName="de_RupRamp_setvar3"
			break
		case "StaticOffset":
			ControlName="de_RupRamp_setvar4"
			break
		case "AlignSmval":
			ControlName="de_RupRamp_setvar5"
			break
		case "Fshift":
			ControlName="de_RupRamp_setvar6"
			break
		case "SShift":
			ControlName="de_RupRamp_setvar7"
			break
		case "FinalSmval":
			ControlName="de_RupRamp_setvar8"
			break
		case "SingleNum":
			ControlName="de_RupRamp_setvar9"
			break	

		default:
			return -1
	endswitch
	
	Controlinfo/W=$PanelName $ControlName
	return v_Value
end

Static Function MakeSmoothedFinalWave()
	
	string saveDF
	saveDF = GetDataFolder(1)

	SetDataFolder ReturnPanelString("Folder")
	wave FinalForceWave=$ReturnPanelString("FinalForceWave")
	wave FinalSepWave=$ReplaceString("Force",nameofwave(FinalForceWave),"Sep")
	string SmType=ReturnPanelString("FinalSmType")
	variable SmAmount=ReturnPanelVal("FinalSmval")
	duplicate/free FinalForceWave ForceWaveSm
	duplicate/free FinalSepWave SepWaveSm
	DE_Filtering#FilterForceSep(FinalForceWave,FinalSepWave,ForceWaveSm,SepWaveSm,SmType,SmAmount)
	duplicate/o ForceWaveSm $(nameofwave(FinalForceWave)+"_Sm")
	duplicate/o SepWaveSm $(nameofwave(FinalSepWave)+"_Sm")

	SetDataFolder saveDF
end

Static Function MakeSmoothedAlignedWave()
	
	string saveDF
	saveDF = GetDataFolder(1)

	SetDataFolder ReturnPanelString("Folder")
	wave AlignedForceWave=$ReturnPanelString("AlignedForceWave")
	wave AlignedSepWave=$ReplaceString("Force",nameofwave(AlignedForceWave),"Sep")
	string SmType=ReturnPanelString("AlignSmType")
	variable SmAmount=ReturnPanelVal("AlignSmval")

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

	//Here we just grab a long list of all the waves we need	
	controlinfo/W=RupRampPanel de_RupRamp_popup0
	SetDataFolder s_value
	controlinfo/W=RupRampPanel de_RupRamp_popup1
	wave ForceWave=$S_value
	wave SepWave=$ReplaceString("Force",S_value,"Sep")
	controlinfo/W=RupRampPanel de_RupRamp_popup6
	wave ShiftedForceWave=$S_value
	controlinfo/W=RupRampPanel de_RupRamp_popup5
	wave StateWave=$S_value

	controlinfo/W=RupRampPanel de_RupRamp_popup3
	wave ForceWaveSM=$S_value
	wave SepWaveSm=$ReplaceString("Force",S_value,"Sep")
	controlinfo/W=RupRampPanel de_RupRamp_popup7
	wave ForceWaveSH_SM=$S_value
	controlinfo/W=RupRampPanel de_RupRamp_popup9
	string WLCWaveFolder=S_value
	controlinfo/W=RupRampPanel de_RupRamp_popup10
	string WLCWaveName=S_value
	wave WLCParms=$(WLCWaveFolder+WLCWaveName)
	controlinfo/W=RupRampPanel de_RupRamp_popup11
	String KindofFit=S_Value
	
	Controlinfo/W=RupRampPanel de_RupRamp_setvar2
	variable distancetoignore=v_value
	controlinfo/W=RupRampPanel de_RupRamp_setvar3
	variable SepMin=v_value 
	controlinfo/W=RupRampPanel de_RupRamp_setvar4
	variable fixedshift=v_value
	controlinfo/W=RupRampPanel de_RupRamp_check0
	variable allowsepshift=v_value

	//Here we make all the waves we need
	make/free/n=0 AlignFoldedForce,AlignFoldedSep,AlignUnFoldedForce,AlignUnFoldedSep //Final aligned waves
	make/free/n=0 ResultsShift,ResultsNoShift
	variable forceshiftused,sepshiftused,forceshiftalt,sepshiftalt,foldedfit
	string shiftString
	
	//This gives us a rough estimate of the pulling speed from the separation waves. 
	//we then use that to calculate how many points to ignore.
	variable/C slopes=DE_Dudko#ReturnSeparationSlopes(SepWaveSm,StateWave,500)
	variable pointstoignore=floor(distancetoignore/real(slopes)/dimdelta(ForceWaveSH_SM,0))


	//Pulls out the unfolded and the folded states from the shifted, smoothed force waves. We ignore about 3 nm of separation after ruptures
	variable timers=stopmstimer(-2)

	DE_NewDudko#ReturnFoldandUnfold(ForceWaveSH_SM,SepWaveSm,StateWave,pointstoignore,AlignFoldedForce,AlignFoldedSep,AlignUnFoldedForce,AlignUnFoldedSep)
	//For processing we make all the forces positive
	AlignFoldedForce*=-1
	AlignUnFoldedForce*=-1


	if(SepMin==0)
	else
	 	CutStatebySep(AlignFoldedForce,AlignFoldedSep,AlignFoldedForce,AlignFoldedSep,0,SepMin)
	endif

	AlignFoldedSep+=fixedshift
	AlignUnFoldedSep+=fixedshift
	//This fits the curves in one of several ways.

	strswitch(KindofFit)
	
		case "Unfolded":
		//Here we fit the UNFOLDED state to a WLC first, then fix everything but LC and fit the FOLDED state. 
		//We actually do this twice first allowing only the force-offset to vary (on the fit of the UNFOLDED state),
		//and then allowing both force and sep offset to vary.
			AlignSingletoWLC(AlignUnFoldedForce,AlignUnFoldedSep,WLCParms,0,ResultsNoShift)
			AlignSingletoWLC(AlignUnFoldedForce,AlignUnFoldedSep,WLCParms,1,ResultsShift)
			foldedfit=0
			break
		
		case "Folded":
			//Here we fit the FOLDED state to a WLC first, then fix everything but LC and fit the UNFOLDED state. 
			//We actually do this twice first allowing only the force-offset to vary (on the fit of the FOLDED state),
			//and then allowing both force and sep offset to vary.
			AlignSingletoWLC(AlignFoldedForce,AlignFoldedSep,WLCParms,0,ResultsNoShift)
			AlignSingletoWLC(AlignFoldedForce,AlignFoldedSep,WLCParms,1,ResultsShift)
			foldedfit=1
			break
		
		case "Both":
			//This is a concurrent fit of BOTH folded and unfolded states, forcing the force and sep offset to be the same for
			//both states. We do this twice, first allwing only the force-offset to vary and then allowing both force and sep offset to vary.
			AlignTwotoWLC(AlignFoldedForce,AlignFoldedSep,AlignUnFoldedForce,AlignUnFoldedSep,WLCParms,0,ResultsNoShift,  0)
			AlignTwotoWLC(AlignFoldedForce,AlignFoldedSep,AlignUnFoldedForce,AlignUnFoldedSep,WLCParms,1,ResultsShift,0)
			
			foldedfit=2
			break

	endswitch
	ResultsShift[1]-=fixedshift
	ResultsNoShift[1]-=fixedshift

	//We recorded the used force and sep shift as well as the alternative shift needed for whichever we didn't use
	//we also note whether we allowed the separation to change.
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
	//Here we report the shifts
	print num2str(forceshiftused)+";"+num2str(sepshiftused)

	//Now we make the aligned Force and Sep waves requested and enact the proper shift.
	duplicate/o ShiftedForceWave $ReplaceString("Force_Shift",nameofwave(ShiftedForceWave),"Force_Align")
	wave AlignedFWave=$ReplaceString("Force_Shift",nameofwave(ShiftedForceWave),"Force_Align")
	FastOP AlignedFWave=ShiftedForceWave-(forceshiftused)
	duplicate/o SepWave $ReplaceString("Sep_Adj",nameofwave(SepWave),"Sep_Align")
	wave AlignedSepWave=$ReplaceString("Sep_Adj",nameofwave(SepWave),"Sep_Align")
	FastOP AlignedSepWave=SepWave-(sepshiftused)
	
	//Begin adding all of this to the wave note starting with the used alignment and whether this allowed 
	//the separation to shift
	String NewNote=ReplaceStringByKey("SepShifted",note(AlignedFWave),shiftString,":","\r")
	NewNote=ReplaceStringByKey("UsedAlignmentFShift",NewNote,num2str(forceshiftused),":","\r")
	NewNote=ReplaceStringByKey("UsedAlignmentSShift",NewNote,num2str(sepshiftused),":","\r")

	//Notes which kind of fitting we did	
	if(FoldedFit==1)
		NewNote=ReplaceStringByKey("Aligned To",NewNote,"Folded",":","\r")
	elseif(FoldedFit==0)
		NewNote=ReplaceStringByKey("Aligned To",NewNote,"UnFolded",":","\r")
	elseif(FoldedFit==2)
		NewNote=ReplaceStringByKey("Aligned To",NewNote,"Both",":","\r")
	endif
	
	
	if(SepMin==0)
	
		NewNote=ReplaceStringByKey("DE_sepMin",NewNote,"0",":","\r")

	else
		NewNote=ReplaceStringByKey("DE_sepMin",NewNote,num2str(SepMin-sepshiftused),":","\r")


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

	//Here we save out the separated FOLDED and UNFOLDED waves for comparison
	make/o/n=(numpnts(AlignFoldedForce),2) AlignFolded
	AlignFolded[][0]= AlignFoldedForce[p]
	AlignFolded[][1]= AlignFoldedSep[p]
	make/o/n=(numpnts(AlignUnFoldedForce),2) AlignUnFolded
	AlignUnFolded[][0]=AlignUnFoldedForce[p]
	AlignUnFolded[][1]=AlignUnFoldedSep[p]
	
	//Making a new Fout and SOut 
	make/free/n=0 Fout,SOut
	DE_TwoWLCFit#CombineCurves(AlignFoldedForce,AlignFoldedSep,AlignUnFoldedForce,AlignUnFoldedSep,Fout,SOut)
	
	//Finally we make the fitted waves both for separation shifted and unshifted 
	make/free/n=6 WLCParmsForFit
	make/free/n=0 FWLCFit
	WLCParmsForFit=WLCParms

	WLCParmsForFit[4]+=ResultsNoShift[1]
	WLCParmsForFit[5]-=ResultsNoShift[0]
	DE_TwoWLCFit#MakeAMultiFit(SOut,WLCParmsForFit,FWLCFit,"Simple")
	make/o/n=(Numpnts(FWLCFit),2) AlignFitNoShift,AlignFitShift,AlignFitOriginal

	AlignFitNoShift[][0]=FWLCFit[p]
	AlignFitNoShift[][1]=SOut[p][0]
	WLCParmsForFit=WLCParms
	WLCParmsForFit[4]+=ResultsShift[1]
	WLCParmsForFit[5]-=ResultsShift[0]

		

	DE_TwoWLCFit#MakeAMultiFit(SOut,WLCParmsForFit,FWLCFit,"Simple")

	AlignFitShift[][0]=FWLCFit[p]
	AlignFitShift[][1]=SOut[p][0]
	
	WLCParmsForFit=WLCParms

	DE_TwoWLCFit#MakeAMultiFit(SOut,WLCParmsForFit,FWLCFit,"Simple")

	AlignFitOriginal[][0]=FWLCFit[p]
	AlignFitOriginal[][1]=SOut[p][0]

	//wave w_sigma,w_coef
	//killwaves W_sigma,W_coef

	SetDataFolder saveDF

end

Static Function CutWavestoMaxandMin(ForceWaveIn,SepWavein,ForceWaveOUt,SepWaveOUt,SepLower,SepMax)

	wave ForceWaveIn,SepWavein,ForceWaveOUt,SepWaveOUt
	variable SepLower,SepMax
	
	duplicate/free ForceWaveIn FreeForce
	duplicate/free SepWavein FreeSep
	variable n
	if(SepLower!=0&&SepMax!=0)
		for(n=0;n<numpnts(FreeSep);n+=1)

			if(FreeSep[n]<SepLower||FreeSep[n]>SepMax)
				FreeSep[n]=NaN
				FreeForce[n]=NaN
			endif
	
		endfor
	
	elseif(SepLower!=0)
		for(n=0;n<numpnts(FreeSep);n+=1)

			if(FreeSep[n]<SepLower)
				FreeSep[n]=NaN
				FreeForce[n]=NaN
			endif
	
		endfor
	
	elseif(SepMax!=0)
		for(n=0;n<numpnts(FreeSep);n+=1)

			if(FreeSep[n]>SepMax)
				FreeSep[n]=NaN
				FreeForce[n]=NaN
			endif
	
		endfor
	
	endif
	duplicate/o FreeForce ForceWaveOUt
	duplicate/o FreeSep SepWaveOUt

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
	
	FuncFit/Q/H=HoldString/NTHR=0 WLC_FIT W_coef  ForceIn /X=SepIn
	make/free/n=2 TemporaryResults

	TemporaryResults={w_coef[3]-WLCParms[3],w_coef[4]-WLCParms[4]}
	wave w_sigma
	killwaves W_coef,W_sigma
	duplicate/o TemporaryResults Results
end

Static Function AlignTwotoWLC(FoldedForce,FoldedSep,UnfoldedForce,UnfoldedSep,WLCParms,SepShift,Results,FixedShift)
	wave FoldedForce,FoldedSep,UnfoldedForce,UnfoldedSep,WLCParms,Results
	variable SepShift,FixedShift

	//note that WLCGuess should have the format: Lp,Lc1,Lc2,T,Xoff,Foff,
	make/free/n=0 Fout,Sout
	
	DE_TwoWLCFit#CombineCurves(FoldedForce,FoldedSep,UnFoldedForce,UnFoldedSep,Fout,SOut)
	duplicate/D/free WLCParms w_coef
	w_coef[4]=FixedShift+WLCParms[4]
	//A bit of code to fix the offset to something specific
	string ConStr="111100"
	if(SepShift==1)
		ConStr="111100"
	elseif(SepShift==0)
		ConStr="111110"
	endif
	if(numpnts(Fout)>.5e6)
		variable down=ceil(numpnts(Fout)/0.5e6)
		Duplicate/free Fout,Fout_samp
		Resample/DOWN=(down)/N=1/WINF=Non Fout_samp
		Duplicate/free Sout,Sout_samp
		Resample/DOWN=(down)/N=1/WINF=Non Sout_samp
		FuncFit/N/Q/W=2/H=ConStr DE_FitTwo W_coef  Fout_samp /X=Sout_samp
	else
		FuncFit/N/Q/W=2/H=ConStr DE_FitTwo W_coef  Fout /X=SOut /D

	endif
	make/free/n=2 TemporaryResults
	print w_Coef
	TemporaryResults={-w_coef[5]+WLCParms[5],w_coef[4]-WLCParms[4]}
	//make/free/n=0 FWLCFit
	//DE_TwoWLCFit#MakeAMultiFit(SOut,w_coef,FWLCFit)
	//make/o/n=(Numpnts(FWLCFit),2) AlignFit
	//AlignFit[][0]=FWLCFit[p]
	//AlignFit[][1]=SOut[p][0]
	//wave w_sigma,w_coef
	//killwaves W_sigma,W_coef
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

Static Function CrunchAbove()
	
	string saveDF
	saveDF = GetDataFolder(1)
	SetDataFolder ReturnPanelString("Folder")

	wave ForceWave=$ReturnPanelString("ForceWave")
	if(cmpstr(nameofwave(ForceWave),"")==0)
		print "No Force Wave"
		return 0
	endif
	wave UpPoints=$ReturnPanelString("RupPnts")
	if(cmpstr(nameofwave(UpPoints),"")==0)
		print "No Pnts Wave"
		return 0
	endif
	
	MakeSmoothedWave()
	popupmenu de_RupRamp_popup3 win=RupRampPanel,popmatch=nameofwave(ForceWave)+"_Sm"

	MakeStateWave()
	popupmenu de_RupRamp_popup5 win=RupRampPanel,popmatch=replacestring("Force",nameofwave(ForceWave),"_States")

	MakeOffsetForceWave()
	ControlUpdate/w=RupRampPanel de_RupRamp_popup6
	popupmenu de_RupRamp_popup6 win=RupRampPanel,popmatch=ReplaceString("Force",nameofwave(ForceWave),"Force_Shift")
	
	
	MakeSmoothedShWave()
	ControlUpdate/w=RupRampPanel de_RupRamp_popup7
	popupmenu de_RupRamp_popup7 win=RupRampPanel,popmatch=ReplaceString("Force",nameofwave(ForceWave),"Force_Shift")+"_Sm"
	
	
	DisplayAPLot("Shifted")

	SetDataFolder saveDF


end

Static Function CrunchMiddle()

	string saveDF
	saveDF = GetDataFolder(1)
	SetDataFolder ReturnPanelString("Folder")
	String AlignFolder=ReturnPanelString("AlignFolder")

	wave AlignParms=$(AlignFolder+ReturnPanelString("AlignParms"))
	if(cmpstr(nameofwave(AlignParms),"")==0)
		print "No Align Parms"
		return 0
	endif
	wave ForceWave=$ReturnPanelString("ForceWave")
	if(cmpstr(nameofwave(ForceWave),"")==0)
		print "No Force Wave"
		return 0
	endif

	AligntoWLC()
	ControlUpdate/w=RupRampPanel de_RupRamp_popup12
	popupmenu de_RupRamp_popup12 win=RupRampPanel,popmatch=ReplaceString("Force_Adj",nameofwave(ForceWave),"Force_Align")
	
	MakeSmoothedAlignedWave()
	ControlUpdate/w=RupRampPanel de_RupRamp_popup14
	popupmenu de_RupRamp_popup14 win=RupRampPanel,popmatch=ReplaceString("Force_Adj",nameofwave(ForceWave),"Force_Align")+"_Sm"
	
	DisplayAPlot("Overlap")

	SetDataFolder saveDF
end

Static Function CrunchBottom()

	string saveDF
	saveDF = GetDataFolder(1)
	SetDataFolder ReturnPanelString("Folder")
	wave ForceWave=$ReturnPanelString("ForceWave")
	if(cmpstr(nameofwave(ForceWave),"")==0)
		print "No Force Wave"
		return 0
	endif
	wave UpPoints=$ReturnPanelString("RupPnts")
	if(cmpstr(nameofwave(UpPoints),"")==0)
		print "No Pnts Wave"
		return 0
	endif
	
	ApplyShifts()
	ControlUpdate/w=RupRampPanel de_RupRamp_popup16
	popupmenu de_RupRamp_popup16 win=RupRampPanel,popmatch=ReplaceString("Force",nameofwave(ForceWave),"Force_Final")

	MakeSmoothedFinalWave()
	ControlUpdate/w=RupRampPanel de_RupRamp_popup17
	popupmenu de_RupRamp_popup17 win=RupRampPanel,popmatch=ReplaceString("Force",nameofwave(ForceWave),"Force_Final_Sm")

	DetermineWLCParms()
	ControlUpdate/w=RupRampPanel de_RupRamp_popup19
	popupmenu de_RupRamp_popup19 win=RupRampPanel,popmatch=ReplaceString("Force",nameofwave(ForceWave),"_WLCParms")
	
	CorrectRuptures()
	ControlUpdate/w=RupRampPanel de_RupRamp_popup20
	controlinfo/W=RupRampPanel de_RupRamp_popup3
	popupmenu de_RupRamp_popup20 win=RupRampPanel,popmatch=nameofwave(UpPoints)+"_Mod"

	MakeSecondStateWave()



	SetDataFolder saveDF
end


Static Function CrunchAll()

	CrunchAbove()
	CrunchMiddle()
	CrunchBottom()
end

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
//////
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
	controlinfo/W=RupRampPanel de_RupRamp_popup1
	wave ForceWave=$S_value
	wave SepWave=$ReplaceString("Force",S_value,"Sep")
	CorrectPauses(ForceWave,SepWave)
	controlinfo/w=RupRampPanel de_RupRamp_popup4
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
	wave SHSepWave=$ReplaceString("Force",nameofwave(SHForceWave),"Sep")

	controlinfo/w=RupRampPanel de_RupRamp_popup8
	string SmType=S_Value
	controlinfo/w=RupRampPanel de_RupRamp_setvar1
	variable SmAmount=v_value
	
	duplicate/free SHForceWave ForceWaveSm
	duplicate/free SHSepWave SepWaveSm
	
	DE_Filtering#FilterForceSep(SHForceWave,SHSepWave,ForceWaveSm,SepWaveSm,SmType,SmAmount)

	duplicate/o ForceWaveSm $(nameofwave(SHForceWave)+"_Sm")
	duplicate/o SepWaveSm $(ReplaceString("Force",nameofwave(SHForceWave),"Sep")+"_Sm")
	SetDataFolder saveDF

end

Static Function MakeSecondStateWave()

	string saveDF
	saveDF = GetDataFolder(1)
	SetDataFolder ReturnPanelString("Folder")

	wave ForceWave=$ReturnPanelString("ForceWave")
	wave SepWave=$ReplaceString("Force",nameofwave(ForceWave),"Sep")
	wave UpPoints=$ReturnPanelString("RupPnts")
	wave DownPoints=$ReplaceString("PntU",nameofwave(UpPoints),"PntD")
	wave StateWave=$ReturnPanelString("StateWave")
	wave FinalForceWave=$ReturnPanelString("FinalForceWave")
	wave WLCParms=$ReturnPanelString("WLCParms")

//	duplicate/free ShForceWave ForceWaveSm
//	duplicate/free SepWave SepWaveSm

	//ForceWaveSm*=-1;

	//make the state key
	make/o/n=0 $ReplaceString("Force",nameofwave(ForceWave),"_2States")
	wave states=$ReplaceString("Force",nameofwave(ForceWave),"_2States")
	DE_DUDKO#MakeSingleStateKey(FinalForceWave,UpPoints,DownPoints,States)

	SetDataFolder saveDF

end


Static Function CorrectRuptures()

	string saveDF
	saveDF = GetDataFolder(1)
	SetDataFolder ReturnPanelString("Folder")

	wave ForceWave=$ReturnPanelString("ForceWave")
	wave SepWave=$ReplaceString("Force",nameofwave(ForceWave),"Sep")
	wave UpPoints=$ReturnPanelString("RupPnts")
	wave DownPoints=$ReplaceString("PntU",nameofwave(UpPoints),"PntD")
	wave StateWave=$ReturnPanelString("StateWave")
	wave FinalForceWave=$ReturnPanelString("FinalForceWave")
	wave FinalSepWave=$ReplaceString("Force",nameofwave(FinalForceWave),"Sep")
	wave WLCParms=$ReturnPanelString("WLCParms")
	//string SmType=ReturnPanelString("FinalSmType")
	//variable SmAmount=ReturnPanelVal("FinalSmval")
	




	
//	controlinfo/W=RupRampPanel de_RupRamp_popup14
//	wave AlignedForceWave=$S_value
//	wave AlignedsepWave=$ReplaceString("Force",S_value,"Sep")
	
	duplicate/free FinalForceWave FWSm,SepWaveSm,FWIN
//	duplicate/free ForceWave SepWaveSm
//	controlinfo/W=RupRampPanel de_RupRamp_popup9
//	wave ForceWaveSM=$S_value
//	wave SepWaveSm=$ReplaceString("Force",S_value,"Sep")
		
//	controlinfo/W=RupRampPanel de_RupRamp_popup10
//	wave ForceWaveSH_SM=$S_value
	DE_Filtering#FilterForceSep(FinalForceWave,FinalSepWave,FWSm,SepWaveSm,"TVD",10e-9)
	FWIN*=-1
	FWSm*=-1

	DE_CorrRup#FixPicks(UpPoints,FWIN,FWSm,WLCParms,SepWaveSm)
//	make/o/n=0 $ReplaceString("Force",nameofwave(ForceWave),"_WLCParms")
//	wave Results=$ReplaceString("Force",nameofwave(ForceWave),"_WLCParms")
//	DE_DUDKO#ContourLengthDetermineCombined(ForceWaveSH_SM,SepWaveSm,StateWave,10000,Results)
	SetDataFolder saveDF


end

Static Function PlotWLCParmFit()
	wave WLC_Folded,WLC_Unfolded,WLCFitfolded,WLCFitUnfolded
	if(waveexists(WLC_Folded)==0)
		return 0
	
	endif
	DoWindow WLCParmFit
	if(V_flag==1)
		Killwindow WLCParmFit
	endif
		
	Display/N=WLCParmFit WLC_Folded[][0] vs WLC_Folded[][1]
	Appendtograph/W=WLCParmFit WLC_unFolded[][0] vs WLC_unFolded[][1]
	Appendtograph/W=WLCParmFit WLCFitfolded
	Appendtograph/W=WLCParmFit WLCFitUnfolded
	ModifyGraph/W=WLCParmFit lsize(WLCFitFolded)=1.5,lsize(WLCFitUnfolded)=1.5
	ModifyGraph/W=WLCParmFit rgb(WLC_Folded)=(19789,44975,19018),rgb(WLC_UnFolded)=(14906,32382,47288),rgb(WLCFitFolded)=(0,0,0),rgb(WLCFitUnfolded)=(0,0,0)


end

Static Function PlotAlignFit()
	wave AlignFolded,AlignUnFolded,AlignFit
	if(waveexists(AlignFolded)==0)
		return 0
	
	endif
	DoWindow Aligned
	if(V_flag==1)
		Killwindow Aligned
	endif
		
	Display/N=Aligned AlignFolded[][0] vs AlignFolded[][1]
	Appendtograph/W=Aligned AlignUnFolded[][0] vs AlignUnFolded[][1]
	Appendtograph/W=Aligned AlignFit[][0] vs AlignFit[][1]
	ModifyGraph/W=Aligned lsize(AlignFit)=1.5
	ModifyGraph/W=Aligned rgb(AlignFolded)=(19789,44975,19018),rgb(AlignFolded)=(14906,32382,47288),rgb(AlignFit)=(0,0,0)


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
	SetDataFolder ReturnPanelString("Folder")
	wave StateWave=$ReturnPanelString("StateWave")

	wave FinalForceWave=$ReturnPanelString("FinalForceWave")
	wave FinalSepWave=$ReplaceString("Force",nameofwave(FinalForceWave),"Sep")
	wave ForceWave=$ReturnPanelString("ForceWave")
	controlinfo/W=RupRampPanel de_RupRamp_popup16
	wave ForceWaveFinalSM=$ReturnPanelString("SmoothedFinalForce")
	wave SepWaveFinalSm=$ReplaceString("Force",nameofwave(ForceWaveFinalSM),"Sep")

	controlinfo/W=RupRampPanel de_RupRamp_popup10
	wave ForceWaveSH_SM=$S_value
	
	
	variable IgnoreDist=ReturnPanelVal("IgnoreDist")
	//Controlinfo/W=RupRampPanel de_RupRamp_setvar3
	variable/C slopes=DE_Dudko#ReturnSeparationSlopes(SepWaveFinalSm,StateWave,500)
	variable pointstoignore=floor(IgnoreDist/real(slopes)/dimdelta(ForceWaveFinalSM,0))
	ForceWaveFinalSM*=-1
	make/o/n=0 $ReplaceString("Force",nameofwave(ForceWave),"_WLCParms")
	wave Results=$ReplaceString("Force",nameofwave(ForceWave),"_WLCParms")
	controlinfo/W=RupRampPanel de_RupRamp_check1
   DE_DUDKO#ContourLengthDetermineCombined(ForceWaveFinalSM,SepWaveFinalSm,StateWave,pointstoignore,Results,v_value,CopyWavesOUt=1)

	ForceWaveFinalSM*=-1
	make/o/n=0 WLCFitFolded,WLCFitUnfolded
	DE_DUDKO#MakeWLCs(SepWaveFinalSm,Results,WLCFitFolded,WLCFitUnfolded)
end

Static Function MakeOffsetForceWave()
	//handle offsetting both the smoothed wave and the unsmoothed wave, also corrects our estimates of the rupture forces
	
	
	string saveDF
	saveDF = GetDataFolder(1)
	
	
	controlinfo/W=RupRampPanel de_RupRamp_popup0
	SetDataFolder s_value
	controlinfo/W=RupRampPanel de_RupRamp_popup1
	wave ForceWave=$S_value
	wave SepWave=$ReplaceString("Force",S_value,"Sep")
	controlinfo/W=RupRampPanel de_RupRamp_popup2
	wave UpPoints=$S_value
	wave DownPoints=$ReplaceString("PntU",S_value,"PntD")
	controlinfo/W=RupRampPanel de_RupRamp_popup3
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

	controlinfo/W=RupRampPanel de_RupRamp_popup0
	SetDataFolder s_value
	controlinfo/W=RupRampPanel de_RupRamp_popup3
	wave ForceWaveSm=$S_value
	controlinfo/W=RupRampPanel de_RupRamp_popup1
	wave ForceWave=$S_value
//	wave SepWave=$ReplaceString("Force",S_value,"Sep")
	controlinfo/W=RupRampPanel de_RupRamp_popup2
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
	NewPanel/N=RupRampPanel /W=(0,0,750,775)
	NewDataFolder/o root:DE_RupRamp
	NewDataFolder/o root:DE_RupRamp:MenuStuff
	Button de_RupRamp_button0,pos={10,30},size={150,20},proc=DE_RuptureRamp#ButtonProc,title="Smooth"
	Button de_RupRamp_button1,pos={10,90},size={150,20},proc=DE_RuptureRamp#ButtonProc,title="Make State Wave"
	Button de_RupRamp_button2,pos={10,120},size={150,20},proc=DE_RuptureRamp#ButtonProc,title="OffSet Waves"
	Button de_RupRamp_button3,pos={10,150},size={150,20},proc=DE_RuptureRamp#ButtonProc,title="Smooth Shifted"
	Button de_RupRamp_button4,pos={550,150},size={150,20},proc=DE_RuptureRamp#ButtonProc,title="Plot Smoothed"
	Button de_RupRamp_button5,pos={250,200},size={150,20},proc=DE_RuptureRamp#ButtonProc,title="Crunch Above"

	SetDrawEnv linethick= 2.00;DrawLine 0,230,750,230;SetDrawEnv linethick= 2.00;	DrawLine 0,235,750,235 

	Button de_RupRamp_button6,pos={10,285},size={150,20},proc=DE_RuptureRamp#ButtonProc,title="Align To WLC"
	Button de_RupRamp_button7,pos={10,365},size={150,20},proc=DE_RuptureRamp#ButtonProc,title="Smooth Aligned"
	Button de_RupRamp_button8,pos={600,365},size={150,20},proc=DE_RuptureRamp#ButtonProc,title="Plot Aligned"
	
	SetDrawEnv linethick= 2.00;DrawLine 0,445,750,445;SetDrawEnv linethick= 2.00;	DrawLine 0,450,750,450 
	
	Button de_RupRamp_button9,pos={10,460},size={150,20},proc=DE_RuptureRamp#ButtonProc,title="Pull Shifts"
	Button de_RupRamp_button10,pos={10,500},size={150,20},proc=DE_RuptureRamp#ButtonProc,title="Apply Additional"
	Button de_RupRamp_button11,pos={10,530},size={150,20},proc=DE_RuptureRamp#ButtonProc,title="Smooth Final"
	Button de_RupRamp_button12,pos={250,400},size={150,20},proc=DE_RuptureRamp#ButtonProc,title="Crunch Middle"
	Button de_RupRamp_button13,pos={10,590},size={150,20},proc=DE_RuptureRamp#ButtonProc,title="MakeWLCFits"
	Button de_RupRamp_button14,pos={10,640},size={150,20},proc=DE_RuptureRamp#ButtonProc,title="Correct Ruptures"
	Button de_RupRamp_button15,pos={10,670},size={150,20},proc=DE_RuptureRamp#ButtonProc,title="Remake StateWave"
	Button de_RupRamp_button16,pos={250,670},size={150,20},proc=DE_RuptureRamp#ButtonProc,title="Crunch Bottom"
	Button de_RupRamp_button17,pos={500,620},size={150,50},proc=DE_RuptureRamp#ButtonProc,title="Crunch All"

	SetDrawEnv linethick= 5.00;DrawLine 0,700,750,700;SetDrawEnv linethick= 5.00;	DrawLine 0,710,750,710 

	Button de_RupRamp_button18,pos={450,740},size={150,20},proc=DE_RuptureRamp#ButtonProc,title="Single"

	PopupMenu de_RupRamp_popup0,pos={10,2},size={129,21},title="Folder",mode=1
	PopupMenu de_RupRamp_popup0,mode=1,popvalue="X",value= #"DE_PanelProgs#ListFolders()"
	PopupMenu de_RupRamp_popup1,pos={200,2},size={129,21},title="Force Wave",mode=0
	PopupMenu de_RupRamp_popup1,mode=1,popvalue="X",value= #"DE_RuptureRamp#ListWaves(\"de_RupRamp_popup0\",\"*Force*Adj\")"
	PopupMenu de_RupRamp_popup2,pos={450,2},size={129,21},title="RuptureUp"
	PopupMenu de_RupRamp_popup2,mode=1,popvalue="X",value= #"DE_RuptureRamp#ListWaves(\"de_RupRamp_popup0\",\"*RupPntU*Adj\")"
	PopupMenu de_RupRamp_popup3,pos={200,30},size={129,21},title="Smoothed Wave"
	PopupMenu de_RupRamp_popup3,mode=1,popvalue="X",value= #"DE_RuptureRamp#ListWaves(\"de_RupRamp_popup0\",\"*Force*Sm\")"
	PopupMenu de_RupRamp_popup4,pos={10,55},size={129,21},title="Type",proc=DE_RuptureRamp#PopMenuProc
	PopupMenu de_RupRamp_popup4,mode=1,value= "TVD;SVG"
	PopupMenu de_RupRamp_popup5,pos={200,90},size={129,21},title="StateWave"
	PopupMenu de_RupRamp_popup5,mode=1,popvalue="X",value= #"DE_RuptureRamp#ListWaves(\"de_RupRamp_popup0\",\"*State*\")"
	PopupMenu de_RupRamp_popup6,pos={200,120},size={129,21},title="Shifted Wave"
	PopupMenu de_RupRamp_popup6,mode=1,popvalue="X",value= #"DE_RuptureRamp#ListWaves(\"de_RupRamp_popup0\",\"*Force*Shift\")"
	PopupMenu de_RupRamp_popup7,pos={200,150},size={129,21},title="Smoothed Shifted Wave"
	PopupMenu de_RupRamp_popup7,mode=1,popvalue="X",value= #"DE_RuptureRamp#ListWaves(\"de_RupRamp_popup0\",\"*Force*Shift*Sm\")"
	PopupMenu de_RupRamp_popup8,pos={10,175},size={129,21},title="Type",proc=DE_RuptureRamp#PopMenuProc
	PopupMenu de_RupRamp_popup8,value= "TVD;SVG",mode=1
	PopupMenu de_RupRamp_popup9,pos={10,250},size={129,21},title="Alignment Folder"
	PopupMenu de_RupRamp_popup9,mode=1,popvalue="X",value= #"DE_PanelProgs#ListFolders()"
	PopupMenu de_RupRamp_popup10,pos={300,250},size={129,21},title="Alignment Parms"
	PopupMenu de_RupRamp_popup10,mode=1,popvalue="X",value= #"DE_RuptureRamp#ListWaves(\"de_RupRamp_popup9\",\"*WLCAlign*\")"
	PopupMenu de_RupRamp_popup11,pos={10,315},size={129,21},title="Alignment Type"
	PopupMenu de_RupRamp_popup11,mode=1,value="Both;Unfolded;Folded"
	PopupMenu de_RupRamp_popup12,pos={200,285},size={129,21},title="Aligned Wave"
	PopupMenu de_RupRamp_popup12,mode=1,popvalue="X",value= #"DE_RuptureRamp#ListWaves(\"de_RupRamp_popup0\",\"*Force*Align\")"
	PopupMenu de_RupRamp_popup13,pos={10,390},size={129,21},title="Type",proc=DE_RuptureRamp#PopMenuProc
	PopupMenu de_RupRamp_popup13,value= "TVD;SVG",mode=1
	PopupMenu de_RupRamp_popup14,pos={200,365},size={129,21},title="Aligned Smoothed Wave"
	PopupMenu de_RupRamp_popup14,mode=1,popvalue="X",value= #"DE_RuptureRamp#ListWaves(\"de_RupRamp_popup0\",\"*Force*Align*Sm\")"
	PopupMenu de_RupRamp_popup15,pos={450,390},size={129,21},title="Overlapped Curve"
	PopupMenu de_RupRamp_popup15,mode=1,popvalue="X",value= #"DE_RuptureRamp#ListWaves(\"de_RupRamp_popup9\",\"*FSM*\")"
	PopupMenu de_RupRamp_popup16,pos={200,500},size={129,21},title="Final Wave"
	PopupMenu de_RupRamp_popup16,mode=1,popvalue="X",value= #"DE_RuptureRamp#ListWaves(\"de_RupRamp_popup0\",\"*Force*Final\")"
	PopupMenu de_RupRamp_popup17,pos={200,530},size={129,21},title="Final Smoothed Wave"
	PopupMenu de_RupRamp_popup17,mode=1,popvalue="X",value= #"DE_RuptureRamp#ListWaves(\"de_RupRamp_popup0\",\"*Force*Final*Sm\")"
	PopupMenu de_RupRamp_popup18,pos={10,555},size={129,21},title="Type",proc=DE_RuptureRamp#PopMenuProc
	PopupMenu de_RupRamp_popup18,value= "TVD;SVG",mode=1
	PopupMenu de_RupRamp_popup19,pos={200,590},size={129,21},title="WLC Parms"
	PopupMenu de_RupRamp_popup19,mode=1,popvalue="X",value= #"DE_RuptureRamp#ListWaves(\"de_RupRamp_popup0\",\"*WLCParms*\")"
	PopupMenu de_RupRamp_popup20,pos={200,640},size={129,21},title="New Up Ruptures"
	PopupMenu de_RupRamp_popup20,mode=1,popvalue="X",value= #"DE_RuptureRamp#ListWaves(\"de_RupRamp_popup0\",\"*RupPntU*\")"


	PopupMenu de_RupRamp_popup21,pos={10,740},size={150,20},title="Up or Down"
	PopupMenu de_RupRamp_popup21,mode=1,popvalue="X",value="Up;Down"
	PopupMenu de_RupRamp_popup22,pos={200,740},size={150,20},title="WLC or Line"
	PopupMenu de_RupRamp_popup22,mode=1,popvalue="root:",value="WLC;Line"
	
	SetVariable de_RupRamp_setvar0,pos={105,55},size={50.00,18.00},proc=DE_RuptureRamp#SetVarProc,value=_num:10e-9
	SetVariable de_RupRamp_setvar0,limits={0,inf,0}
	SetVariable de_RupRamp_setvar1,pos={105,175},size={50.00,18.00},proc=DE_RuptureRamp#SetVarProc,value=_num:10e-9
	SetVariable de_RupRamp_setvar1,limits={0,inf,0}
	SetVariable de_RupRamp_setvar2,pos={175,315},size={99.00,18.00},proc=DE_RuptureRamp#SetVarProc,value=_num:3e-9
	SetVariable de_RupRamp_setvar2,limits={0,inf,0},title="Ignore\rDistance"
	SetVariable de_RupRamp_setvar3,pos={395,315},size={150,18.00},proc=DE_RuptureRamp#SetVarProc,value=_num:0
	SetVariable de_RupRamp_setvar3,limits={0,inf,0},title="SepMin"
	SetVariable de_RupRamp_setvar4,pos={560,315},size={125,18.00},proc=DE_RuptureRamp#SetVarProc,value=_num:0
	SetVariable de_RupRamp_setvar4,limits={0,inf,0},title="Static\rOffset"
	SetVariable de_RupRamp_setvar5,pos={105,390},size={50.00,18.00},proc=DE_RuptureRamp#SetVarProc,value=_num:10e-9
	Setvariable de_RupRamp_setvar5,limits={0,inf,0}
	SetVariable de_RupRamp_setvar6,pos={200,460},size={150.00,18.00},proc=DE_RuptureRamp#SetVarProc,value=_num:0
	SetVariable de_RupRamp_setvar6,limits={0,inf,0},title="Custom\rForce Shift"
	SetVariable de_RupRamp_setvar7,pos={400,460},size={150.00,18.00},proc=DE_RuptureRamp#SetVarProc,value=_num:0
	SetVariable de_RupRamp_setvar7,limits={0,inf,0},title="Custom\rSepShift"
	SetVariable de_RupRamp_setvar8,pos={105,555},size={50.00,18.00},proc=DE_RuptureRamp#SetVarProc,value=_num:10e-9
	Setvariable de_RupRamp_setvar8,limits={0,inf,0}
	SetVariable de_RupRamp_setvar9,pos={350,740},size={50.00,18.00},proc=DE_RuptureRamp#SetVarProc,value=_num:0
	SetVariable de_RupRamp_setvar9,limits={0,inf,0}

	CheckBox de_RupRamp_check0 title="All Sep Shift",pos={295,315},size={150,25},proc=DE_RuptureRamp#CheckProc
//	//CheckBox de_RupRamp_check1 title="Fit Folded",pos={550,450},size={150,25},proc=DE_RuptureRamp#CheckProc
	CheckBox de_RupRamp_check1 title="Fit Folded",pos={10,610},size={150,25},proc=DE_RuptureRamp#CheckProc

EndMacro



Static Function/S ListWaves(ControlStr,SearchString)
	string ControlStr,SearchString
	String saveDF

	saveDF = GetDataFolder(1)
	controlinfo $ControlStr
	string Result=s_value
	SetDataFolder Result
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

Static Function CheckPointWave()
	
	controlinfo/W=RupRampPanel de_RupRamp_popup4
	if(cmpstr(S_Value,"")==0)
		print "No Force Wave"
		return 0
	else
		wave ForceWave=$S_Value
	endif
	
	controlinfo/W=RupRampPanel de_RupRamp_popup3
	if(cmpstr(S_Value,"")==0)
		print "No Points Wave"
		return 0
	else
		wave UpPoints=$S_value
		wave DownPoints=$ReplaceString("PntU",S_value,"PntD")
	endif

	variable n
	for(n=0;n<numpnts(UpPoints);n+=1)
	endfor
	for(n=0;n<numpnts(DownPoints);n+=1)
	endfor

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

