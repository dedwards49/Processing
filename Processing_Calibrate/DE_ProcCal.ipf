#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma modulename= DE_ProcCal
 
#include "C:\Devin\Documents\Software\Library\processing\Misc_PanelPrograms\Panel Progs"
#include "C:\Devin\Documents\Software\Library\processing\Misc_PanelPrograms\AsylumNaming"
Constant cKb = 1.3806504e-23;		//Boltzman constant in J / K 2006 CODATA value
Constant ThermalT = 298;		//Estimated temperature for Thermal

Static Function/S ListFolders(NameString)
	String NameString
	string list=DE_PanelProgs#PrintAllFolders_String(NameString)
	return list
End

Static Function TracesinFolder(FolderString)
	String FolderString
	DFREF saveDFR = GetDataFolderDFR()	// Get reference to current data folder

	SetDataFolder FolderString	// Create a new data folder
	variable DefsinFolder=itemsinlist(Wavelist("DeflV_*",";",""))	// Do some work in it
	wave w1=Results
	variable dimennum=dimsize(w1,1)
	if(DefsinFolder!=dimennum)
		return -1
	endif
	SetDataFolder saveDFR					// Restore current data folder

	return DefsinFolder
End
Static Function/S  StringTracesinFolder()
	Controlinfo/W=ProcCal DE_ProcCal_popup0
	
		String FolderString=S_Value
		variable top=TracesinFolder(FolderString)
		variable n
		if(top==0)
		return ""
		endif
		String Result=""
		for(n=0;n<top;n+=1)
			Result+=num2str(n)+";"
		
		endfor
		return Result
End

Static Function InitiateActiveResults()
	ControlInfo/W=ProcCal DE_ProcCal_popup1
	variable tracenumber=v_value-1
	ControlInfo/W=ProcCal DE_ProcCal_popup0
	String FolderString=s_value
	DFREF saveDFR = GetDataFolderDFR()	// Get reference to current data folder
	SetDataFolder FolderString	// Create a new data folder
	wave AllInfo=Results
	wave ActiveResults
	if(WaveExists(ActiveResults)==0)
		duplicate/o AllInfo ActiveResults
	endif
	
	SetDataFolder saveDFR	// Create a new data folder
end

Static Function ResetActiveResults()
	ControlInfo/W=ProcCal DE_ProcCal_popup1
	variable tracenumber=v_value-1
	ControlInfo/W=ProcCal DE_ProcCal_popup0
	String FolderString=s_value
	DFREF saveDFR = GetDataFolderDFR()	// Get reference to current data folder
	SetDataFolder FolderString	// Create a new data folder
	wave AllInfo=Results
	wave ActiveResults
	
	duplicate/o AllInfo ActiveResults
		
	SetDataFolder saveDFR	// Create a new data folder
end


Static Function InitiateOrigThermalFit()
	ControlInfo/W=ProcCal DE_ProcCal_popup1
	variable tracenumber=v_value-1
	ControlInfo/W=ProcCal DE_ProcCal_popup0
	String FolderString=s_value
	DFREF saveDFR = GetDataFolderDFR()	// Get reference to current data folder
	SetDataFolder FolderString	// Create a new data folder
	
			
	wave DeflV=$("DeflV_"+num2str(tracenumber))
	wave ZSnsr=$("ZSnsr_"+num2str(tracenumber))
	wave ThermalTop=$("Thermal_"+num2str(tracenumber)+"_top")
	wave ThermalBottom=$("Thermal_"+num2str(tracenumber)+"_bottom")
	wave AllInfo=Results
	
	make/o/n=5 OrigThermal
	SetDimLabel 0,0,ADC,OrigThermal
	SetDimLabel 0,1,Omega,OrigThermal
	SetDimLabel 0,2,Q,OrigThermal
	SetDimLabel 0,3,WhiteNoise,OrigThermal
	SetDimLabel 0,4,FitWidth,OrigThermal
	

	variable k=AllInfo[%TopSpringConstant][tracenumber]
	variable omega=AllInfo[%TopFrequency][tracenumber]
	variable Q=AllInfo[%TopQValue][tracenumber]
	ThermalParmsToWave(k,Q,omega,OrigThermal)
	variable FitWidthVar=omega*1.95
	OrigThermal[%FitWidth]=FitWidthVar

	variable WN=FindWN(ThermalTop,OrigThermal)
	OrigThermal[%WhiteNoise]=WN
	
	Duplicate/o OrigThermal ActiveThermal
	SetDataFolder saveDFR	// Create a new data folder
end

Static Function ThermalParmsToWave(k,Q,omega,WaveOut)
	variable k,Q,omega
	wave WaveOut
	//w[0] = Amplitude at DC
	//w[1] = Q
	//w[2] = omega0, in Hz
	//x is the variable omega
	
	
	//ref:
	//Short cantilevers for atomic force microscopy.  Walters DA, Cleveland
	//JP, Thomson NH, Hansma PK, Wendman MA, Gurley G, Elings V, Rev. Sci.
	//Instrum. 67 (10) 1996, p. 3583, endnote 39
	k*=1e9
	variable energy = cKb * ThermalT
	variable ADC=sqrt((2*energy/(Pi*omega*abs(Q)*k)))
	
	WaveOut[%ADC]=ADC
	WaveOut[%Q]=Q
	WaveOut[%omega]=omega
	
	//return(2*energy/(Pi*w[2]*abs(w[1])*w[0]^2))
end

Static Function MakeAthermal(ThermalParms,ThermalReturn)

	wave ThermalParms,ThermalReturn

	ControlInfo/W=ProcCal DE_ProcCal_popup1
	variable tracenumber=v_value-1
	ControlInfo/W=ProcCal DE_ProcCal_popup0
	String FolderString=s_value
	DFREF saveDFR = GetDataFolderDFR()	// Get reference to current data folder
	SetDataFolder FolderString	// Create a new data folder
	wave ThermalTop=$("Thermal_"+num2str(tracenumber)+"_top")
	wave ThermalBottom=$("Thermal_"+num2str(tracenumber)+"_bottom")

	variable startpnt=x2pnt(ThermalTop,ThermalParms[%omega]-ThermalParms[%FitWidth]/2)
	variable endpnt=x2pnt(ThermalTop,ThermalParms[%omega]+ThermalParms[%FitWidth]/2)	
	duplicate/Free/r=[startpnt,endpnt] ThermalTop Fit,XWave
	XWave=pnt2x(Fit,p)
	make/free/n=4 SHOParms
	SHOParms={ThermalParms[%Adc],ThermalParms[%Q],ThermalParms[%Omega],ThermalParms[%WhiteNoise]}
	Fit=SHOAmpWhite(SHOParms,Fit,XWave)
	duplicate/o Fit ThermalReturn
	SetDataFolder saveDFR					// Restore current data folder


	//
end

Static Function InitiateCurrentThermal()
	ControlInfo/W=ProcCal DE_ProcCal_popup1
	variable tracenumber=v_value-1
	ControlInfo/W=ProcCal DE_ProcCal_popup0
	String FolderString=s_value
	DFREF saveDFR = GetDataFolderDFR()	// Get reference to current data folder
	SetDataFolder FolderString	// Create a new data folder
	wave AllInfo=Results
	make/o/n=5 CurrentThermal
	SetDimLabel 0,0,ADC,CurrentThermal
	SetDimLabel 0,1,Omega,CurrentThermal
	SetDimLabel 0,2,Q,CurrentThermal
	SetDimLabel 0,3,WhiteNoise,CurrentThermal
	SetDimLabel 0,4,FitWidth,CurrentThermal

	SetDataFolder saveDFR	// Create a new data folder
end

Static Function FindWN(ThermalWave,WaveIn)
	wave ThermalWave,WaveIn
	
	make/free/n=0 FitWave
	Fitwave={WaveIn[%ADC],WaveIn[%Q],WaveIn[%omega],WaveIn[%ADC]/25}

	variable startpnt=x2pnt(ThermalWave,Fitwave[2]/2)
	variable endpnt=x2pnt(ThermalWave,3*Fitwave[2]/2)
	FuncFit/Q/W=2/NTHR=0/H="1110" SHOAmpWhite Fitwave  ThermalWave[startpnt,endpnt]
	wave w_coef
	return Fitwave[3]
	
end

Function SHOAmpWhite(Parm,Output,XWave)	: FitFunc		//this is an all at once fit function
	wave Parm, Output, XWave							//the XWave does nothing at the moment
	
	//Amplitude response for a SHO. "Vibration and Waves" p. 89. This is SHOAmp with white noise added.
	//Parm[0] = Amplitude at DC
	//Parm[1] = Q
	//Parm[2] = omega0, in either radians or Hz
	//Parm[3] = white noise. Units are m/rtHz
	//x is the variable omega
	
	variable DC = Parm[0]							//using local variables is faster than using wave points
	variable QTerm = 1/Parm[1]^2
	variable Omega0 = Parm[2]
	variable WNTerm = Parm[3]^2
	
	Output = Sqrt((DC*Omega0/XWave)^2/ ((Omega0/XWave - XWave/Omega0)^2 + QTerm) + WNTerm)	//calculate the wave
	
	return 0

end //SHOAmpWhite
Static function Thermalk(w)
	wave w
	//returns the thermal spring constant calculated from one of the SHO functions
	//This assumes you have fit a linear Power Spectral Density of the cantilever motion (units of m/rtHz)
	//using a function such as SHOAmp or SHOAmpPinkWhite
	//w[0] = Amplitude at DC
	//w[1] = Q
	//w[2] = omega0, in Hz
	//x is the variable omega
	
	
	//ref:
	//Short cantilevers for atomic force microscopy.  Walters DA, Cleveland
	//JP, Thomson NH, Hansma PK, Wendman MA, Gurley G, Elings V, Rev. Sci.
	//Instrum. 67 (10) 1996, p. 3583, endnote 39
	
	variable energy = cKb * ThermalT
	return(2*energy/(Pi*w[2]*abs(w[1])*w[0]^2))

end //Thermalk
	
Static function InverseThermalk(parm,springConstant,invOLS)
	wave parm
	variable springConstant, invOLS
	
	//calculating the new InvOLS if the spring constant is assumed to be the same
	variable dcMeters = parm[0]
	variable qVar = parm[1]
	variable freq = parm[2]
	variable dcVolts = dcMeters/invOLS
	variable energy = cKb * ThermalT
	
	return sqrt((2*energy)/(springConstant*pi*freq*qVar))/dcVolts

end //InversThermalk


Static Function SelectedANumber()
	ControlInfo/W=ProcCal DE_ProcCal_popup1
	variable tracenumber=v_value-1
	ControlInfo/W=ProcCal DE_ProcCal_popup0
	String FolderString=s_value
	DFREF saveDFR = GetDataFolderDFR()	// Get reference to current data folder
	SetDataFolder FolderString	// Create a new data folder
	
	wave DeflV=$("DeflV_"+num2str(tracenumber))
	wave ZSnsr=$("ZSnsr_"+num2str(tracenumber))
	wave ThermalTop=$("Thermal_"+num2str(tracenumber)+"_top")
	wave ThermalBottom=$("Thermal_"+num2str(tracenumber)+"_bottom")
	wave Select=SelectWave

	checkbox DE_ProcCal_Check0,value=Select[tracenumber]
	duplicate/o DeflV 	 root:DE_ProcCal:YDisp
	duplicate/o ZSnsr 	 root:DE_ProcCal:XDisp
	duplicate/o ThermalTop 	 root:DE_ProcCal:Thermal
	wavestats/Q DeflV
	duplicate/o DeflV root:DE_ProcCal:ZColor
	variable SurfaceTriggerTime=V_maxloc
	wave ZColor=root:DE_ProcCal:ZColor
	ZColor[0,x2pnt(DeflV,SurfaceTriggerTime)]=3
	ZColor[x2pnt(DeflV,SurfaceTriggerTime),]=-7
	
	UpdateDisplayVars()
	
	InitiateOrigThermalFit()
	wave ThermalParms=OrigThermal
	wave OrigThermalFit=root:DE_ProcCal:OThermalFit
	MakeAthermal(ThermalParms,OrigThermalFit)
	SetDataFolder saveDFR					// Restore current data folder
	
end

Static Function UpdateDisplayVars()

	ControlInfo/W=ProcCal DE_ProcCal_popup1
	variable tracenumber=v_value-1
	ControlInfo/W=ProcCal DE_ProcCal_popup0
	String FolderString=s_value
	DFREF saveDFR = GetDataFolderDFR()	// Get reference to current data folder
	SetDataFolder FolderString	// Create a new data folder
	
	wave DeflV=$("DeflV_"+num2str(tracenumber))
	wave ZSnsr=$("ZSnsr_"+num2str(tracenumber))
	wave ThermalTop=$("Thermal_"+num2str(tracenumber)+"_top")
	wave ThermalBottom=$("Thermal_"+num2str(tracenumber)+"_bottom")
	wave Select=SelectWave


	NVAR OInv=root:DE_ProcCal:OrigInvols
	NVAR OTopSpr=root:DE_ProcCal:OrigTopSpring
	NVAR OBotSpr=root:DE_ProcCal:OrigBotSpring

	wave Results=$("Results")
	OInv=Results[%Invols][tracenumber]
	OTopSpr=Results[%TopSpringConstant][tracenumber]*1e12
	OBotSpr=Results[%BotSpringConstant][tracenumber]*1e12
	wave AResult=ActiveResults
	if(WaveExists(AResult)==0)
	
	else
		NVAR CInv=root:DE_ProcCal:CurrentInvols
		NVAR CTopSpr=root:DE_ProcCal:CurrentTopSpring
		NVAR CBotSpr=root:DE_ProcCal:CurrentBotSpring
		CInv=AResult[%Invols][tracenumber]
		CTopSpr=AResult[%TopSpringConstant][tracenumber]*1e12
		CBotSpr=AResult[%BotSpringConstant][tracenumber]*1e12
		
		variable Ostart,Oend,CStart,Cend
		OStart=Results[%CSrPnt1][tracenumber]
		OEnd=Results[%CSrPnt2][tracenumber]
		CStart=AResult[%CSrPnt1][tracenumber]
		CEnd=AResult[%CSrPnt2][tracenumber]
		
		
	endif
	
	//Check to see if I should make a local Invols

end


Static Function MakeSelectionWave()
	ControlInfo/W=ProcCal DE_ProcCal_popup0
	String FolderString=s_value
	variable num=TracesinFolder(FolderString)
	DFREF saveDFR = GetDataFolderDFR()	// Get reference to current data folder
	SetDataFolder FolderString	// Create a new data folder
	wave SelectWave
	if(WaveExists(SelectWave)==0)
		make/o/n=(num) SelectWave
		SelectWave=1
	endif
	SetDataFolder saveDFR					// Restore current data folder
end

Static Function CheckonExisting()
	ControlInfo/W=ProcCal DE_ProcCal_popup1
	variable tracenumber=v_value-1
	ControlInfo/W=ProcCal DE_ProcCal_popup0
	String FolderString=s_value
	DFREF saveDFR = GetDataFolderDFR()	// Get reference to current data folder
	SetDataFolder FolderString	// Create a new data folder
	wave Results
	wave ActiveResults
	
	variable Ostart,Oend,CStart,Cend
	OStart=Results[%CSrPnt1][tracenumber]
	OEnd=Results[%CSrPnt2][tracenumber]
	CStart=ActiveResults[%CSrPnt1][tracenumber]
	CEnd=ActiveResults[%CSrPnt2][tracenumber]

	if(OStart==CStart&&OEnd==CEnd)
			
		make/o/n=0 root:DE_ProcCal:CInvolsFit,root:DE_ProcCal:CInvolsRegionY,root:DE_ProcCal:CInvolsRegionX
	else
		MakeInvolsFit("Existing")
	endif

	SetDataFolder saveDFR					// Restore current data folder

end

Static Function PopupMenuAction(PU_Struct) : PopupMenuControl
	STRUCT WMPopupAction &PU_Struct
	Switch(PU_Struct.eventCode)
		case 2:

			StrSwitch(PU_Struct.ctrlname)
				case "DE_ProcCal_popup0":
					MakeSelectionWave()
					InitiateActiveResults()
					CurrentReport()
					OriginalReport()
				break
				case "DE_ProcCal_popup1":
					SelectedANumber()
					MakeInvolsFit("Original")
					CheckonExisting()

					break
				
	
			endswitch
			break
	
	
	endswitch	
End

Static Function MakeInvolsFit(TypeString)
	String TypeString
	
	ControlInfo/W=ProcCal DE_ProcCal_popup1	
	variable tracenumber=v_value-1
	ControlInfo/W=ProcCal DE_ProcCal_popup0
	String FolderString=s_value
	DFREF saveDFR = GetDataFolderDFR()	// Get reference to current data folder
	SetDataFolder FolderString	// Create a new data folder
	wave DeflV=$("DeflV_"+num2str(tracenumber))
	wave ZSnsr=$("ZSnsr_"+num2str(tracenumber))
	wave ThermalTop=$("Thermal_"+num2str(tracenumber)+"_top")
	wave ThermalBottom=$("Thermal_"+num2str(tracenumber)+"_top")
	wave AllInfo=Results
	wave ActiveResults

	variable startpnt,endpnt,newslope
	StrSwitch(TypeString)
	
		case "Original":
			startpnt=AllInfo[%CsrPnt1][tracenumber]
			endpnt=AllInfo[%CsrPnt2][tracenumber]
			CurveFit/Q/W=2/NTHR=0 line  DeflV[startpnt,endpnt]/X=ZSnsr[startpnt,endpnt] /D
			wave AutoFitWave=$("fit_"+nameofwave(DeflV))
			duplicate/o AutoFitWave root:DE_ProcCal:OInvolsFit
			duplicate/o/R=[startpnt,endpnt] DeflV root:DE_ProcCal:OInvolsRegionY
			duplicate/o/r=[startpnt,endpnt] ZSnsr root:DE_ProcCal:OInvolsRegionX
		break
		
		case "Existing":
			startpnt=ActiveResults[%CsrPnt1][tracenumber]
			endpnt=ActiveResults[%CsrPnt2][tracenumber]
			CurveFit/Q/W=2/NTHR=0 line  DeflV[startpnt,endpnt]/X=ZSnsr[startpnt,endpnt] /D
			wave AutoFitWave=$("fit_"+nameofwave(DeflV))
			duplicate/o AutoFitWave root:DE_ProcCal:CInvolsFit
			duplicate/o/R=[startpnt,endpnt] DeflV root:DE_ProcCal:CInvolsRegionY
			duplicate/o/r=[startpnt,endpnt] ZSnsr root:DE_ProcCal:CInvolsRegionX
		break
		
		case "New":
			string AString=CSRINFO(A,"ProcCal#Invols")
			string BString=CSRINFO(B,"ProcCal#Invols")
			if(cmpstr(AString,"")==0||cmpstr(BString,"")==0)
			print "Incomplete Cursors"
			else
			startpnt=pcsr(A,"ProcCal#Invols")
			endpnt=pcsr(B,"ProcCal#Invols")
			ActiveResults[%CsrPnt1][tracenumber]=startpnt
			ActiveResults[%CsrPnt2][tracenumber]=endpnt
			CurveFit/Q/W=2/NTHR=0 line  DeflV[startpnt,endpnt]/X=ZSnsr[startpnt,endpnt] /D
			wave AutoFitWave=$("fit_"+nameofwave(DeflV))
			newslope=InvolsFromFit(AutoFitWave)/1e-9
			duplicate/o AutoFitWave root:DE_ProcCal:CInvolsFit
			duplicate/o/R=[startpnt,endpnt] DeflV root:DE_ProcCal:CInvolsRegionY
			duplicate/o/r=[startpnt,endpnt] ZSnsr root:DE_ProcCal:CInvolsRegionX
			endif
		break
		
	
	endswitch
	Killwaves AutoFitWave
	wave FitWave=root:DE_ProcCal:InvolsFit
	
	SetDataFolder saveDFR// Restore current data folder
	return newslope
end

Static Function UpdateInvolsActive(ActiveResults,OrigResults,tracenumber,NewInvols)

	wave ActiveResults,OrigResults
	variable tracenumber,NewInvols
	variable oldinvols,conversion

	if(numtype(ActiveResults[%TopSpringConstant][tracenumber])!=0||numtype(ActiveResults[%TopSpringConstant][tracenumber])!=0)
		oldInvols=OrigResults[%Invols][tracenumber]
		ActiveResults[%Invols][tracenumber]=NewInvols
		conversion=(OldInvols/NewInvols)^2
		
		ActiveResults[%TopSpringConstant][tracenumber]=conversion*OrigResults[%TopSpringConstant][tracenumber]
		ActiveResults[%BotSpringConstant][tracenumber]=conversion*OrigResults[%BotSpringConstant][tracenumber]
	else
		oldInvols=ActiveResults[%Invols][tracenumber]
		ActiveResults[%Invols][tracenumber]=NewInvols
		conversion=(OldInvols/NewInvols)^2
		ActiveResults[%TopSpringConstant][tracenumber]*=conversion
		ActiveResults[%BotSpringConstant][tracenumber]*=conversion
	
	endif


end

Static Function ApplyNewInvols(newinvols)
	
	variable newinvols
	ControlInfo/W=ProcCal DE_ProcCal_popup1
	variable tracenumber=v_value-1
	ControlInfo/W=ProcCal DE_ProcCal_popup0
	String FolderString=s_value
	DFREF saveDFR = GetDataFolderDFR()	// Get reference to current data folder
	SetDataFolder FolderString	// Create a new data folder
	wave ActiveResults=ActiveResults
	wave OrigResults=Results

	UpdateInvolsActive(ActiveResults,OrigResults,tracenumber,newinvols)

	SetDataFolder saveDFR					// Restore current data folder


end

Static Function InvolsFromFit(FitWave)
	wave FitWave
	String WaveNote=note(FitWave)
	String APart=stringfromlist(11,WaveNote,"\r")
	String BPart=stringfromlist(12,WaveNote,"\r")
	variable test1,test2
	sscanf BPart,  "  	b	=%f ± %f",test1,test2
	variable ZLVDTSens=1.6517e-06
	variable ResultingSlope=ZLVDTSens/test1
	return ResultingSlope
end

Static Function InvolsFromORes()

	ControlInfo/W=ProcCal DE_ProcCal_popup1
	variable tracenumber=v_value-1
	ControlInfo/W=ProcCal DE_ProcCal_popup0
	String FolderString=s_value
	DFREF saveDFR = GetDataFolderDFR()	// Get reference to current data folder
	SetDataFolder FolderString	// Create a new data folder
	wave AllInfo=Results
	variable invOLS=AllInfo[%Invols][tracenumber]*1e-9
	SetDataFolder saveDFR					// Restore current data folder
	return Invols

end

Window ProcCal_Panel() : Panel

	PauseUpdate; Silent 1		// building window...
	NewPanel/N=ProcCal /W=(0,0,1200,	500)
	NewDataFolder/o root:DE_ProcCal
	NewDataFolder/o root:DE_ProcCal:MenuStuff
	PopupMenu DE_ProcCal_popup0,pos={75,2},size={129,21},title="Folder",mode=1,proc=DE_ProcCal#PopupMenuAction
	PopupMenu DE_ProcCal_popup0,mode=1,popvalue="X",value= #"DE_ProcCal#ListFolders(\"Results\")"
	PopupMenu DE_ProcCal_popup1,pos={250,2},size={129,21},title="TouchOff",mode=1,proc=DE_ProcCal#PopupMenuAction
	PopupMenu DE_ProcCal_popup1,mode=1,popvalue="X",value= #"DE_ProcCal#StringTracesinFolder()"
	Checkbox DE_ProcCal_Check0,pos={400,2},size={129,21},title="Use",mode=0,proc=DE_ProcCAl#CheckBoxProc
	variable/G root:DE_ProcCal:OrigInvols=0,root:DE_ProcCal:OrigTopSpring=0,root:DE_ProcCal:CurrentInvols=0,root:DE_ProcCal:CurrentTopSpring=0
	variable/G root:DE_ProcCal:OrigBotSpring=0,root:DE_ProcCal:CurrentBotSpring=0

	ValDisplay DE_ProcCal_VDisp0,pos={75,425},size={129,21},title="Orig Invols",value=root:DE_ProcCal:OrigInvols
	ValDisplay DE_ProcCal_VDisp1,pos={225,425},size={129,21},title="Orig Top Spring",value=root:DE_ProcCal:OrigTopSpring
	ValDisplay DE_ProcCal_VDisp2,pos={375,425},size={129,21},title="Orig Bot Spring",value=root:DE_ProcCal:OrigBotSpring
	ValDisplay DE_ProcCal_VDisp3,pos={75,450},size={129,21},title="Orig Invols",value=root:DE_ProcCal:CurrentInvols
	ValDisplay DE_ProcCal_VDisp4,pos={225,450},size={129,21},title="Orig Top Spring",value=root:DE_ProcCal:CurrentTopSpring
	ValDisplay DE_ProcCal_VDisp5,pos={375,450},size={129,21},title="Orig Bot Spring",value=root:DE_ProcCal:CurrentBotSpring

	Button DE_ProcCal_button0,pos={75,475},size={129,21},title="Take Invols",proc=DE_ProcCAl#ButtonProc
	//ValDisplay DE_ProcCal_VDisp1,pos={75,500},size={129,21},title="Orig Invols"
	//ValDisplay DE_ProcCal_VDisp2,pos={75,500},size={129,21},title="Orig Invols"

	make/o/n=(7,3)/T root:DE_ProcCal:OSummary
	root:DE_ProcCal:OSummary[][0]={"0","0","0","0","0","0","0"}
	root:DE_ProcCal:OSummary[][1]={"0","0","0","0","0","0","0"}
	root:DE_ProcCal:OSummary[][2]={"","nm/V","pN/nm","pN/nm","pN/nm","Hz",""}
	duplicate/o root:DE_ProcCal:OSummary root:DE_ProcCal:CSummary

	ListBox DE_ProcCal_LB0,pos={867,72},size={150,110},title="Use",mode=0,proc=DE_ProcCAl#ListBoxProc
 	ListBox DE_ProcCal_LB0,listWave= root:DE_ProcCal:OSummary
	ListBox DE_ProcCal_LB1,pos={867,200},size={150,110},title="Use",mode=0,proc=DE_ProcCAl#ListBoxProc
 	ListBox DE_ProcCal_LB1,listWave= root:DE_ProcCal:CSummary
	
	make/o/n=0 root:DE_ProcCal:Thermal, root:DE_ProcCal:YDisp,root:DE_ProcCal:XDisp,root:DE_ProcCal:OThermalFit,root:DE_ProcCal:LineFit, root:DE_ProcCal:ZColor
	make/o/n=0 root:DE_ProcCal:OInvolsRegionY,root:DE_ProcCal:OInvolsRegionX,root:DE_ProcCal:OInvolsFit,root:DE_ProcCal:CInvolsFit
	make/o/n=0 root:DE_ProcCal:CInvolsRegionY,root:DE_ProcCal:CInvolsRegionX,root:DE_ProcCal:CThermalFit
	
	
	display/host=ProcCal/N=InvOLS/W=(10,50,400,400) root:DE_ProcCal:YDisp vs root:DE_ProcCal:XDisp
	ModifyGraph/W=ProcCal#InvOLS zColor($nameofwave(root:DE_ProcCal:YDisp))={ root:DE_ProcCal:ZColor,-10,10,Spectrum,0}
	Appendtograph/W=ProcCal#Invols root:DE_ProcCal:OInvolsRegionY vs root:DE_ProcCal:OInvolsRegionX
	ModifyGraph lsize(OInvolsRegionY)=4,rgb(OInvolsRegionY)=(0,0,0)
	
	Appendtograph/W=ProcCal#Invols root:DE_ProcCal:CInvolsRegionY vs root:DE_ProcCal:CInvolsRegionX
	ModifyGraph lsize(CInvolsRegionY)=4,rgb(CInvolsRegionY)=(29440,0,58880)
	
	Appendtograph/W=ProcCal#Invols root:DE_ProcCal:OInvolsFit
	ModifyGraph lsize(OInvolsFit)=2,rgb(OInvolsFit)=(16384,65280,41216)
	
	Appendtograph/W=ProcCal#Invols root:DE_ProcCal:CInvolsFit
	ModifyGraph lsize(CInvolsFit)=2,rgb(CInvolsFit)=(65280,65280,32768)
	
	display/host=ProcCal/N=Thermal/W=(450,50,850,400) root:DE_ProcCal:Thermal
	ModifyGraph/W=ProcCal#Thermal log=1
	Appendtograph/W=ProcCal#Thermal root:DE_ProcCal:OThermalFit
	Appendtograph/W=ProcCal#Thermal root:DE_ProcCal:CThermalFit
	ModifyGraph/W=ProcCal#Thermal hideTrace(CThermalFit)=1
	ModifyGraph/W=ProcCal#Thermal lsize(OThermalFit)=2,rgb(OThermalFit)=(0,65280,0)
	SetAxis/W=ProcCal#Thermal bottom 100,*
end	

Static Function ListBoxProc(LB_Struct) : ListboxControl
	STRUCT WMListboxAction &LB_Struct

End
Static Function ButtonProc(B_Struct) : ButtonControl
	STRUCT WMButtonAction &B_Struct
	Switch(B_Struct.eventcode)
		case 2:
			StrSwitch(B_Struct.ctrlname)
				case "DE_ProcCal_button0":
					variable newinvols=MakeInvolsFit("New")
					
					ApplyNewInvols(newinvols)
					UpdateDisplayVars()
					CurrentReport()
				break
			
			
			endswitch
		break
	
	endswitch
End

Static Function PopulateAReport(ResultsWave,SelectWave,TextWave)
	wave ResultsWave,SelectWave
	wave/T TextWave
	
	TextWave[0][0]=num2str(sum(SelectWave))
	
	make/free/n=0 AveragesOut1,AveragesOut2
	AverageStats(ResultsWave,SelectWave,AveragesOut1,"Invols")
	wavestats/Q AveragesOut1
	TextWave[1][0]=num2str(v_avg)
	TextWave[1][1]=num2str(v_sdev)
	
	AverageStats(ResultsWave,SelectWave,AveragesOut1,"TopSpringConstant")
	AverageStats(ResultsWave,SelectWave,AveragesOut2,"BotSpringConstant")
	AveragesOut1*=1e12
	AveragesOut2*=1e12
	wavestats/Q AveragesOut1
	TextWave[3][0]=num2str(v_avg)
	TextWave[3][1]=num2str(v_sdev)
	wavestats/Q AveragesOut2
	TextWave[4][0]=num2str(v_avg)
	TextWave[4][1]=num2str(v_sdev)
	Concatenate/NP/o {AveragesOUt1,averagesout2}, Combined
	wavestats/q Combined 
	TextWave[2][0]=num2str(v_avg)
	TextWave[2][1]=num2str(v_sdev)

	AverageStats(ResultsWave,SelectWave,AveragesOut1,"TopFrequency")
	AverageStats(ResultsWave,SelectWave,AveragesOut2,"BotFrequency")
	Concatenate/NP/o {AveragesOUt1,averagesout2}, Combined
	wavestats/q Combined 
	TextWave[5][0]=num2str(v_avg)
	TextWave[5][1]=num2str(v_sdev)
	
	AverageStats(ResultsWave,SelectWave,AveragesOut1,"TopQValue")
	AverageStats(ResultsWave,SelectWave,AveragesOut2,"BotQValue")
	Concatenate/NP/o {AveragesOUt1,averagesout2}, Combined
	wavestats/q Combined 
	TextWave[6][0]=num2str(v_avg)
	TextWave[6][1]=num2str(v_sdev)
	
	killwaves combined
end



Static Function OriginalReport()
	ControlInfo/W=ProcCal DE_ProcCal_popup0
	String FolderString=s_value
	wave OResults=$(FolderString+"Results")
	make/free/n=(dimsize(OResults,1)) SelectWave
	SelectWave=1
	wave/T TextWave=$"root:DE_ProcCAl:OSummary"
	PopulateAReport(OResults,SelectWave,TextWave)

end

Static Function CurrentReport()

	ControlInfo/W=ProcCal DE_ProcCal_popup0
	String FolderString=s_value
	wave OResults=$(FolderString+"ActiveResults")

	wave SelectWave=$(FolderString+"SelectWave")
	wave/T TextWave=$"root:DE_ProcCAl:CSummary"
	PopulateAReport(OResults,SelectWave,TextWave)

end

Static Function AverageStats(ResultsWave,SelectWave,AveragesOut,ParmString)
	wave ResultsWave,SelectWave,AveragesOut
	String ParmString
	variable tot=dimsize(ResultsWave,1)
	variable n
	make/free/n=0 Parm
	for(n=0;n<tot;n+=1)
		if(Selectwave[n]==1)
			insertpoints 0,1, Parm
			Parm[0]=resultsWave[%$ParmString][n]	
		endif
	
	
	endfor
	duplicate/o Parm AveragesOut 
end


Static Function CheckBoxProc(CB_Struct) : CheckBoxControl
	STRUCT WMCheckboxAction &CB_Struct
	switch(CB_Struct.eventCode)
		case 2:
			ControlInfo/W=ProcCal DE_ProcCal_popup1
			variable tracenumber=v_value-1
			ControlInfo/W=ProcCal DE_ProcCal_popup0
			String FolderString=s_value
			DFREF saveDFR = GetDataFolderDFR()	// Get reference to current data folder
			SetDataFolder FolderString	// Create a new data folder
			wave SelectWave
			SelectWave[tracenumber]=CB_Struct.checked
		
			SetDataFolder saveDFR					// Restore current data folder
			CurrentReport()

			break
		
	
	endswitch

End

Menu "Calibrate"
	"Load Calibrate Processing", ProcCal_Panel()
	End
//Basic Code
//ControlInfo/W=ProcCal DE_ProcCal_popup1
//	variable tracenumber=v_value-1
//	ControlInfo/W=ProcCal DE_ProcCal_popup0
//	String FolderString=s_value
//	DFREF saveDFR = GetDataFolderDFR()	// Get reference to current data folder
//	SetDataFolder FolderString	// Create a new data folder
//
//
//	SetDataFolder saveDFR					// Restore current data folder
