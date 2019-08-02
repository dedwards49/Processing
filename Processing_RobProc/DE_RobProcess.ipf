#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma modulename=DE_RobProc
#include "DE_Filtering"
#include ":\Misc_PanelPrograms\Panel Progs"
#include ":\Misc_PanelPrograms\AsylumNaming"
Static Function GenTracesinnm(DefV,ZSnsrV,Def,ZSnsr,Force,Sep)
	wave DefV,ZSnsrV,Def,ZSnsr,Force,Sep
	variable InVols=str2num(stringbykey("Invols",note(DefV),"=","\r"))
	if(InVols>1e-6)
		InVols*=1e-9
	
	endif
	variable k=str2num(stringbykey("K",note(DefV),"=","\r"))
	if(k>1)
		k*=1e-3
	
	endif

	variable LVDTSens=str2num(stringbykey("ZLVDTSens",note(DefV),"=","\r"))
	duplicate/free DefV TestDefl,TestZSnsr,TestForce,TestSep
	FastOp TestDefl=(InVols)*TestDefl
	FastOp TestZSnsr=(LVDTSens)*ZSnsrV
	FastOp TestForce=(k)*TestDefl
	FastOp TestSep=TestDefl-TestZSnsr

	duplicate/o TestDefl Def
	duplicate/o TestZSnsr ZSnsr
	duplicate/o TestForce Force
	duplicate/o TestSep Sep

end

Static Function InitialRampPrep(DefV,ZSnsrV,Def,ZSnsr,Force,Sep)
	wave DefV,ZSnsrV,Def,ZSnsr,Force,Sep
	variable InVols=str2num(stringbykey("Invols",note(DefV),"=","\r"))
	if(InVols>1e-6)
		InVols*=1e-9
	
	endif
	variable k=str2num(stringbykey("K",note(DefV),"=","\r"))
	if(k>1)
		k*=1e-3
	
	endif

	variable LVDTSens=str2num(stringbykey("ZLVDTSens",note(DefV),"=","\r"))
	duplicate/free DefV TestDefl,TestZSnsr,TestForce,TestSep
	FastOp TestDefl=(InVols)*TestDefl
	FastOp TestZSnsr=(LVDTSens)*ZSnsrV
	FastOp TestForce=(k)*TestDefl
	FastOp TestSep=TestDefl-TestZSnsr

	duplicate/o TestDefl Def
	duplicate/o TestZSnsr ZSnsr
	duplicate/o TestForce Force
	duplicate/o TestSep Sep

end

Static Function FindBreaksandPauses(ZSensorSetWave,Settings,OutWave)

	wave ZSensorSetWave,Settings,OutWave
	variable numberofscans=floor(Settings[%TotalTime]/(Settings[%Distance]/Settings[%RetractVElocity]+Settings[%Distance]/Settings[%ApproachVElocity]+Settings[%DwellTime]))
	variable OutTime=Settings[%Distance]/Settings[%RetractVelocity]
	variable InTime=Settings[%Distance]/Settings[%ApproachVelocity]
	variable Pause=Settings[%DwellTime]
	make/free/n=(numberofscans,3) Times
	variable n
	for(n=0;n<numberofscans;n+=1)
		if(n==0)
		Times[n][0]=OutTime
		else
		Times[n][0]=Times[n-1][2]+OutTime
		endif
		Times[n][1]=Times[n][0]+InTime
		Times[n][2]=Times[n][1]+Pause
	endfor
	//Times[numberofscans-1][2]=pnt2x(ZSensorSetWave,numpnts(ZSensorSetWave)-1)
	duplicate/o Times OutWave

//	duplicate/free ZSensorSetWave ZSetDif
//	Differentiate ZSensorSetWave/D=ZSetDif
//	make/free/n=(numberofscans,3) Times
//	variable startslope=ZSetDif[0]
//	variable startsearch
//	variable n
//	for(n=0;n<numberofscans;n+=1)
//		FindValue/S=(startsearch)/T=.1e-2/V=(-1*startslope)  ZSetDif
//		Times[n][0]=pnt2x(ZSensorSetWave,v_value)
//		startsearch=V_value
//		FindValue/S=(startsearch)/T=.1e-2/V=(0)  ZSetDif
//		Times[n][1]=pnt2x(ZSensorSetWave,v_value)
//		startsearch=V_value
//		FindValue/S=(startsearch)/T=.1e-2/V=(startslope)  ZSetDif
//		Times[n][2]=pnt2x(ZSensorSetWave,v_value)
//		startsearch=v_value
//	endfor
//	if(v_value!=-1)
//		return -1
//	endif
//
//	
//	duplicate/o Times OutWave
	return 0
end

Static Function AddSectionNotes(DeflVWave,WaveofLocations)
	wave DeflVWave,WaveofLocations
	variable n=0
	String Locs=""
	String Dirs=""
	String PauseLoc=""
	String PauseState=""

	for(n=0;n<(dimsize(WaveofLocations,0)-1);n+=1)
		
		Locs+=num2istr(x2pnt(DeflVWave,WaveofLocations[n][0]))+";"
		Dirs+="Ret;"
		Locs+=num2istr(x2pnt(DeflVWave,WaveofLocations[n][2]))+";"
		Dirs+="Ext;"
		PauseLoc+=num2istr(x2pnt(DeflVWave,WaveofLocations[n][1]))+";"
		PauseState+="Start;"
		PauseLoc+=num2istr(x2pnt(DeflVWave,WaveofLocations[n][2]))+";"
		PauseState+="End;"
	endfor
	
			Locs+=num2istr(x2pnt(DeflVWave,WaveofLocations[(dimsize(WaveofLocations,0)-1)][0]))+";"
		Dirs+="Ret;"
		Locs+=num2istr(numpnts(DeflVWave)-1)+";"
		Dirs+="Ext;"
		PauseLoc+=num2istr(x2pnt(DeflVWave,WaveofLocations[(dimsize(WaveofLocations,0)-1)][1]))+";"
		PauseState+="Start;"
		PauseLoc+=num2istr(numpnts(DeflVWave)-1)+";"
		PauseState+="End;"
	String NoteString=note(DeflVWave)

	NoteString=ReplaceStringbyKey("DE_Dir",NoteString,Dirs,":","\r")
	NoteString=ReplaceStringbyKey("DE_Ind",NoteString,Locs,":","\r")
	NoteString=ReplaceStringbyKey("DE_PauseLoc",NoteString,PauseLoc,":","\r")
	NoteString=ReplaceStringbyKey("DE_PauseState",NoteString,PauseState,":","\r")
	note/k DeflVWave,NoteString
	
end

Static Function SplitInitial(IDefl,IZsnsr,IDefRetOut)
	wave IDefl,IZsnsr,IDefRetOut
	make/free/n=0 ZSm
	Differentiate IZsnsr/D=ZSm
	make/free/n=0 DiffSm
	DE_Filtering#TVD1D_denoise(ZSm,2e-6,DiffSm)
	FindLevel/Q/edge=2  DiffSm,0
	variable Touch=v_levelx
	FindLevel/Q/B=50/edge=2/R=(v_levelx,)  DiffSm,-5e-9
	variable Leave=v_levelx
	String OutZRet=ReplaceString("Defl",nameofwave(IDefRetOut),"ZSnsr")
	String OutDTowd=ReplaceString("Ret",nameofwave(IDefRetOut),"Towd")
	String OutZTowd=ReplaceString("Defl",OutDTowd,"ZSnsr")
	String OutDExt=ReplaceString("Ret",nameofwave(IDefRetOut),"Ext")
	String OutZExt=ReplaceString("Defl",OutDExt,"ZSnsr")

	duplicate/free/r=(0,Touch) IDefl IdeflExt
	duplicate/free/r=(Touch,Leave) IDefl IdeflTowd
	duplicate/free/r=(Leave,) IDefl IdeflRet
	duplicate/free/r=(0,Touch) IZsnsr IZExt
	duplicate/free/r=(Touch,Leave) IZsnsr IZTowd
	duplicate/free/r=(Leave,) IZsnsr IZRet
	duplicate/o IdeflRet IDefRetOut
	duplicate/o IZRet $OutZRet
	
	duplicate/o IdeflTowd $OutDTowd
	duplicate/o IZTowd $OutZTowd
	
	duplicate/o IdeflExt $OutDExt
	duplicate/o IZExt $OutZExt

end

Static Function ButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	Strswitch(ba.ctrlname)
			string saveDF
			
		case "de_RobProc_button0":
			switch( ba.eventCode )
				case 2: // mouse up
					saveDF = GetDataFolder(1)
					controlinfo/W=RobProcPanel de_RobProc_popup0
					SetDataFolder s_value
					controlinfo/W=RobProcPanel de_RobProc_popup1
					wave InitDefV=$S_value
					controlinfo/W=RobProcPanel de_RobProc_popup2
					wave InitZV=$S_value
					controlinfo/W=RobProcPanel de_RobProc_setvar0
					string BaseName=S_Value
					controlinfo/W=RobProcPanel de_RobProc_setvar1
					variable NewNumber=V_Value
					controlinfo/W=RobProcPanel de_RobProc_popup3
					wave DFV=$S_value	
					controlinfo/W=RobProcPanel de_RobProc_popup4
					wave ZSV=$S_value
					wave ZSet=$replacestring("Sensor",S_Value,"SetPoint")
					controlinfo/W=RobProcPanel de_RobProc_popup5
					wave Settings=$S_value
					String NewNote=note(DFV)
					NewNote=ReplaceStringByKey("BaseSuffix", NewNote, num2str(NewNumber),":","\r")
					note/K InitDefV NewNote
					note/K ZSV NewNote
					NewNote=note(DFV)
					NewNote=ReplaceStringByKey("BaseSuffix", NewNote, num2str(NewNumber),":","\r")
					note/K DFV NewNote
					note/K InitZV NewNote
					make/o/n=0 $DE_Naming#StringCreate(BaseName,NewNumber,"Defl","Ret")
					wave InitDeflRet=$DE_Naming#StringCreate(BaseName,NewNumber,"Defl","Ret")
					//make/o/n=0 $DE_Naming#StringCreate(BaseName,NewNumber,"ZSnsr","Ret")
					//wave InitZSnsrRet=$DE_Naming#StringCreate(BaseName,NewNumber,"ZSnsr","Ret")

					
					
					make/o/n=0 $DE_Naming#StringCreate(BaseName,NewNumber,"Defl","All")
					wave Deflection=$DE_Naming#StringCreate(BaseName,NewNumber,"Defl","All")
					make/o/n=0 $DE_Naming#StringCreate(BaseName,NewNumber,"ZSnsr","All")
					wave ZWave=$DE_Naming#StringCreate(BaseName,NewNumber,"ZSnsr","All")
					make/o/n=0 $DE_Naming#StringCreate(BaseName,NewNumber,"Force","All")
					wave Force=$DE_Naming#StringCreate(BaseName,NewNumber,"Force","All")
					make/o/n=0 $DE_Naming#StringCreate(BaseName,NewNumber,"Sep","All")
					wave Sep=$DE_Naming#StringCreate(BaseName,NewNumber,"Sep","All")
					
					make/free/n=0 Garbage1,Garbage2,IDefl,IZSen
					GenTracesinnm(InitDefV,InitZV,IDefl,IZSen,Garbage1,Garbage2)
					SplitInitial(IDefl,IZSen,InitDeflRet)
					
					
					make/free/n=0 OutWave
					FindBreaksandPauses(ZSet,Settings,OutWave)
					AddSectionNotes(DFV,OutWave)
					GenTracesinnm(DFV,ZSV,Deflection,ZWave,Force,Sep)

//					make/o/n=0 $(Add+replacestring("DEfV",nameofwave(DFV),"Defl"))
//					wave Deflection=$(Add+replacestring("DEfV",nameofwave(DFV),"Defl"))
//					make/o/n=0 $(Add+replacestring("ZSensor",nameofwave(ZSV),"Extension"))
//					wave ZWave=$(Add+replacestring("ZSensor",nameofwave(ZSV),"Extension"))
//					
//					make/o/n=0 $(Add+replacestring("DEfV",nameofwave(DFV),"Force"))
//					wave Force=$(Add+replacestring("DEfV",nameofwave(DFV),"Force"))
//					make/o/n=0 $(Add+replacestring("ZSensor",nameofwave(ZSV),"Sep"))
//					wave Sep=$(Add+replacestring("ZSensor",nameofwave(ZSV),"Sep"))
//					
//					
//						make/o/n=0 FFil,SFil
//					DE_Filtering# FilterForceSep(Force,Sep,FFil,SFil,"TVD",5e-8)
					SetDataFolder saveDF

					break
				case -1: // control being killed
					break
			endswitch
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

	endswitch
End	
Window RobProcPanel() : Panel

	PauseUpdate; Silent 1		// building window...
	NewPanel/N=RobProcPanel /W=(0,0,400,300)
	NewDataFolder/o root:DE_RobProc
	NewDataFolder/o root:DE_RobProc:MenuStuff

//	DE_RobProc#UpdateParmWave()
	Button de_RobProc_button0,pos={75,130},size={150,20},proc=DE_RobProc#ButtonProc,title="GO!"
//	Button de_PSF_button1,pos={200,80},size={80,20},proc=DE_PSF#ButtonProc2,title="AddCursors!"

	PopupMenu de_RobProc_popup0,pos={75,2},size={129,21},Title="Folder"
	PopupMenu de_RobProc_popup0,mode=1,popvalue="X",value= #"DE_PanelProgs#ListFolders()"
	PopupMenu de_RobProc_popup1,pos={75,40},size={129,21},Title="Initial DefV"
	PopupMenu de_RobProc_popup1,mode=1,popvalue="X",value= #"DE_RobProc#ListWaves(\"de_RobProc_popup0\",\"*InitRamp_D*\")"
	PopupMenu de_RobProc_popup2,pos={75,70},size={129,21},Title="Initial ZVolts"
	PopupMenu de_RobProc_popup2,mode=1,popvalue="X",value= #"DE_RobProc#ListWaves(\"de_RobProc_popup0\",\"*InitRamp_Z*\")"

	
	PopupMenu de_RobProc_popup3,pos={220,40},size={129,21},Title="DefVWave"
	PopupMenu de_RobProc_popup3,mode=1,popvalue="X",value= #"DE_RobProc#ListWaves(\"de_RobProc_popup0\",\"*DefV*\")"
	PopupMenu de_RobProc_popup4,pos={220,70},size={129,21},Title="ZsnsVWave"
	PopupMenu de_RobProc_popup4,mode=1,popvalue="X",value= #"DE_RobProc#ListWaves(\"de_RobProc_popup0\",\"*Sensor*\")"
	PopupMenu de_RobProc_popup5,pos={220,95},size={129,21},Title="Settings"
	PopupMenu de_RobProc_popup5,mode=1,popvalue="X",value= #"DE_RobProc#ListWaves(\"de_RobProc_popup0\",\"*Settings*\")"
	SetVariable de_RobProc_setvar0,pos={75,175},size={150,16},proc=DE_robProc#SVP,value= _STR:"Image",title="Base Name"
	SetVariable de_RobProc_setvar1,pos={75,200},size={150,16},proc=DE_robProc#SVP,value= _NUM:0,title="New Number"
	

//	SetVariable de_PSF_setvar0,pos={75,75},size={100,20},proc=DE_PSF#SVP,value= _NUM:0,title="StartPoint"
//	SetVariable de_PSF_setvar0 value= _NUM:0
//	SetVariable de_PSF_setvar1,pos={75,105},size={100,20},proc=DE_PSF#SVP,value= _NUM:0,title="EndPoint"
//	SetVariable de_PSF_setvar1 value= _NUM:0
EndMacro

//Static Function UpdateParmWave()
//	if(exists("root:DE_RobProc:MenuStuff:ParmWave")==1)
//		wave/t/z Par=root:DE_RobProc:MenuStuff:ParmWave
//		wave/z Sel=root:DE_RobProc:MenuStuff:SelWave
//	Else
//		make/t/n=(4,2) root:DE_RobProc:MenuStuff:ParmWave
//		wave/t/z Par=root:DE_RobProc:MenuStuff:ParmWave
//		make/n=(4,2) root:DE_RobProc:MenuStuff:SelWave
//		wave/z Sel=root:DE_RobProc:MenuStuff:SelWave
//		
//		Par[0][0]={"Smoothing","Pulling Rate","Distance to Ignore","Distance to Fit"}
//		Par[0][1]={"10e-9","50e-9","3e-9","3e-9"}
//		Sel[][0]=0
//		Sel[][1]=2
//	endif
//
//end

Static Function/S ListWaves(ControlStr,SearchStr)
	string ControlStr,SearchStr
	String saveDF

	saveDF = GetDataFolder(1)
	controlinfo $ControlStr
	SetDataFolder s_value
	String list = WaveList(SearchStr, ";", "")
	SetDataFolder saveDF
	return list

end

Menu "Ramp"
	"Open RobProcess", RobProcPanel()

end