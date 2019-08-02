#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma modulename=DE_CTFCMenuProc
#include "DE_Filtering"
#include "DE_TimeFuncs"

#include ":\Misc_PanelPrograms\Panel Progs"
#include ":\Misc_PanelPrograms\AsylumNaming"


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Menu "Ramp"
	"Open CTFCMenuProcPanel", CTFCMenuProcPanel()

end


Window CTFCMenuProcPanel() : Panel

	PauseUpdate; Silent 1		// building window...
	NewPanel/N=CTFCProcPanel /W=(0,0,400,300)
	NewDataFolder/o root:DE_CTFCProc
	NewDataFolder/o root:DE_CTFCProc:MenuStuff

	Button de_CTFCProc_button0,pos={75,130},size={150,20},proc=DE_CTFCMenuProc#ButtonProc,title="GO!"
	Button de_CTFCProc_button1,pos={250,130},size={150,20},proc=DE_CTFCMenuProc#ButtonProc,title="Print!"


	PopupMenu de_CTFCProc_popup0,pos={75,2},size={129,21},Title="Folder"
	PopupMenu de_CTFCProc_popup0,mode=1,popvalue="X",value= #"DE_PanelProgs#ListFolders()"
	PopupMenu de_CTFCProc_popup1,pos={75,40},size={129,21},Title="DefV",proc=DE_CTFCMenuProc#PopMenuProc
	PopupMenu de_CTFCProc_popup1,mode=1,popvalue="X",value= #"DE_CTFCMenuProc#ListWaves(\"de_CTFCProc_popup0\",\"*\")"
	PopupMenu de_CTFCProc_popup2,pos={75,70},size={129,21},Title="ZVolts"
	PopupMenu de_CTFCProc_popup2,mode=1,popvalue="X",value= #"DE_CTFCMenuProc#ListWaves(\"de_CTFCProc_popup0\",\"*\")"

	PopupMenu de_CTFCProc_popup3,pos={75,100},size={129,21},Title="Type"
	PopupMenu de_CTFCProc_popup3,mode=1,popvalue="X",value= "Single;Multi;Initial;EquilRamp;Equil;500 kHz;Touch"

	SetVariable de_CTFCProc_setvar0,pos={75,175},size={150,16},value= _STR:"Image",title="Base Name"

	SetVariable de_CTFCProc_setvar1,pos={75,200},size={150,16},value= _NUM:0,title="New Number"
	
	SetVariable de_CTFCProc_setvar2,pos={150,40},size={150,16},value= _STR:"Image",title="Search Str"

EndMacro


Static Function PopMenuProc(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa

	switch( pa.eventCode )
		case 2: // mouse up
			Variable popNum = pa.popNum
			String popStr = pa.popStr
			strswitch(pa.ctrlName)
			
				case "de_CTFCProc_popup1":
					break			
			endswitch
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End



Static Function ButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	Strswitch(ba.ctrlname)
			string saveDF
			
		case "de_CTFCProc_button0":
			switch( ba.eventCode )
				case 2: // mouse up
					//				
					saveDF = GetDataFolder(1)
					controlinfo/W=CTFCProcPanel de_CTFCProc_popup0
					SetDataFolder s_value
					controlinfo/W=CTFCProcPanel de_CTFCProc_popup1
					wave DefV=$S_value
					controlinfo/W=CTFCProcPanel de_CTFCProc_popup2
					wave ZV=$S_value
					controlinfo/W=CTFCProcPanel de_CTFCProc_setvar0
					string BaseName=S_Value
					controlinfo/W=CTFCProcPanel de_CTFCProc_setvar1
					variable NewNumber=V_Value
					controlinfo/W=CTFCProcPanel de_CTFCProc_popup3
					String Type=S_value
					ProcessVoltWaves(Type,DefV,ZV,BaseName,NewNumber)
					

					SetDataFolder saveDF

					break
				case -1: // control being killed
					break

			endswitch
			break
		case "de_CTFCProc_button1":
			switch( ba.eventCode )
				case 2: // mouse up
					saveDF = GetDataFolder(1)
					controlinfo/W=CTFCProcPanel de_CTFCProc_popup0
					SetDataFolder s_value
					controlinfo/W=CTFCProcPanel de_CTFCProc_popup1
					wave DefV=$S_value
					TellRelevantWaves(DefV)
					SetDataFolder saveDF
					break
			endswitch
			
			break
	endswitch
	return 0
End

Static Function/S ListWaves(ControlStr,SearchStr)
	string ControlStr,SearchStr
	String saveDF
	controlinfo/w=CTFCProcPanel de_CTFCProc_setvar2
	SearchStr=S_Value
	saveDF = GetDataFolder(1)
	controlinfo $ControlStr
	SetDataFolder s_value
	String list = WaveList(SearchStr, ";", "")
	SetDataFolder saveDF
	return list

end

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Static Function GenTracesinnm(DefV,ZSnsrV,Def,ZSnsr,Force,Sep,[Rate])
	wave DefV,ZSnsrV,Def,ZSnsr,Force,Sep
	variable Rate
	variable InVols=str2num(stringbykey("Invols",note(DefV),":","\r"))
	if(InVols>1e-6)
		InVols*=1e-9
	
	endif
	variable k=str2num(stringbykey("Spring Constant",note(DefV),":","\r"))
	if(k>1)
		k*=1e-3
	
	endif

	variable LVDTSens=str2num(stringbykey("ZLVDTSens",note(DefV),":","\r"))
	duplicate/free DefV TestDefl,TestZSnsr,TestForce,TestSep
	duplicate/free ZSnsrV ,TestZSnsr,TestSep
	WaveTransform zapNaNs, TestDefl
	WaveTransform zapNaNs, TestZSnsr
	WaveTransform zapNaNs, TestForce
	WaveTransform zapNaNs, TestSep

	FastOp TestDefl=(InVols)*TestDefl
	FastOp TestZSnsr=(LVDTSens)*TestZSnsr
	if(ParamisDefault(Rate))
	
		DE_Filtering#MakeZPositionFinal(TestZSnsr,wavemax(TestZSnsr)-wavemin(TestZSnsr))
	else
		DE_Filtering#MakeZPositionFinal(TestZSnsr,wavemax(TestZSnsr)-wavemin(TestZSnsr),rate=rate)

	endif
	FastOp TestForce=(k)*TestDefl
	FastOp TestSep=TestDefl-TestZSnsr
	duplicate/o TestDefl Def
	duplicate/o TestZSnsr ZSnsr
	duplicate/o TestForce Force
	duplicate/o TestSep Sep

end


Static Function ProcessVoltWaves(Type,DefV,ZVolt,BaseString,WaveNumber)
	String Type,BaseString
	wave DefV,ZVolt
	variable WaveNumber
	String NewNote=note(DefV)
	NewNote=ReplaceStringByKey("BaseSuffix", NewNote, num2str(WaveNumber),":","\r")
	note/K DefV NewNote
	note/K ZVolt NewNote
	make/o/n=0 Garbage1,Garbage2,IDefl,IZSen
	string tracename

	variable cutpoint
	strSwitch(Type)
	
		case "Multi":
			GenTracesinnm(DefV,ZVolt,IDefl,IZSen,Garbage1,Garbage2)//,rate=1e2)
			make/free/n=0 Mbreaks
			FindMultiBreaks(IDefl,IZSen,Mbreaks)
			AddSectionNotes(IDefl,MBreaks)
			duplicate/o IDefl $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"Defl","All")
			duplicate/o IZSen $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"ZSnsr","All")
			break
	
		case "EquilRamp":
			GenTracesinnm(DefV,ZVolt,IDefl,IZSen,Garbage1,Garbage2)
			duplicate/o IDefl $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"Defl","Ret")
			duplicate/o IZSen $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"ZSnsr","Ret")
			break
	
		case "Equil":
			GenTracesinnm(DefV,ZVolt,IDefl,IZSen,Garbage1,Garbage2)
			duplicate/o IDefl $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"Defl","Equil")
			duplicate/o IZSen $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"ZSnsr","Equil")
			break
	
		case "Initial":
			GenTracesinnm(DefV,ZVolt,IDefl,IZSen,Garbage1,Garbage2)
			make/free/n=0 IBreaks
			FindInitialBreaksandPauses(IDefl,IZSen,IBreaks)
			cutpoint=x2pnt(IDefl,IBreaks[0])
			duplicate/o/R=[0,cutpoint-1] IDefl $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"Defl","Ext")
			duplicate/o/R=[0,cutpoint-1] IZSen $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"ZSnsr","Ext")
			duplicate/o/R=[cutpoint,] IDefl $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"Defl","Ret")
			duplicate/o/R=[cutpoint,] IZSen $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"ZSnsr","Ret")
			break
			
		case "Touch":
			GenTracesinnm(DefV,ZVolt,IDefl,IZSen,Garbage1,Garbage2)
			make/free/n=0 IBreaks
			FindTouchBreaksandPauses(IDefl,IZSen,IBreaks)
			cutpoint=x2pnt(IDefl,IBreaks[0])
			duplicate/o/R=[0,cutpoint-1] IDefl $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"Defl","Ext")
			duplicate/o/R=[0,cutpoint-1] IZSen $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"ZSnsr","Ext")
			duplicate/o/R=[cutpoint,] IDefl $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"Defl","Ret")
			duplicate/o/R=[cutpoint,] IZSen $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"ZSnsr","Ret")
			break
	
		case "Single":
			GenTracesinnm(DefV,ZVolt,IDefl,IZSen,Garbage1,Garbage2)
			make/free/n=0 SBreaks
			FindMSingleBreaksandPauses(IDefl,IZSen,SBreaks)
			if(numpnts(Sbreaks)==0)
			print nameofwave(DefV)
			endif
			cutpoint=x2pnt(IDefl,SBreaks[0])
			if( cutpoint<0)
				cutpoint=0
			endif
			duplicate/o/R=[0,cutpoint-1] IDefl $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"Defl","Ret")
			duplicate/o/R=[0,cutpoint-1] IZSen $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"ZSnsr","Ret")
			duplicate/o/R=[cutpoint,] IDefl $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"Defl","Ext")
			duplicate/o/R=[cutpoint,] IZSen $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"ZSnsr","Ext")
			break
		
		case "500kHz":
			Prompt traceName,"Wave",popup,wavelist("*Defl",";","")
			DoPrompt "Pick a wave",traceName
			WAVE w = $traceName
							
			NewNote=note(w)
			NewNote=ReplaceStringByKey("BaseSuffix", NewNote, num2str(WaveNumber),":","\r")
			NewNote=ReplaceStringByKey("Bandwidth", NewNote, num2str(500),":","\r")

			note/K DefV NewNote
			note/K ZVolt NewNote
		
			GenTracesinnm(DefV,ZVolt,IDefl,IZSen,Garbage1,Garbage2,rate=50e3)
			make/o/n=0 MBreaks,Mbreaks2
			FindMultiBreaks(IDefl,IZSen,Mbreaks2,smoothing=501)
			AddSectionNotes(IDefl,Mbreaks2)
			duplicate/o IDefl $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"Defl","Fast")
			duplicate/o IZSen $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"ZSnsr","Fast")
			break		
		default:
			break
	
	endswitch
	killwaves Garbage1,Garbage2,IDefl,IZSen

end

Static Function ProcessCentering(CXRX,CXRZ,CYRY,CYRZ,BaseString,WaveNumber)
	String BaseString
	wave CXRX,CXRZ,CYRY,CYRZ
	variable WaveNumber
	String NewNote=note(CXRX)
	NewNote=ReplaceStringByKey("BaseSuffix", NewNote, num2str(WaveNumber),":","\r")
	note/K CXRX NewNote
	note/K CXRZ NewNote
	note/K CYRY NewNote
	note/K CYRZ NewNote
	variable XLVDTSENS=str2num(stringbykey("XLVDTSens",note(CXRX),":","\r"))
	variable YLVDTSENS=str2num(stringbykey("YLVDTSens",note(CXRX),":","\r"))
	variable ZLVDTSENS=str2num(stringbykey("ZLVDTSens",note(CXRX),":","\r"))

	duplicate/free CXRX CXRXnm
	CXRXnm*=XLVDTSENS
	duplicate/free CYRY CYRYnm
	CYRYnm*=YLVDTSENS
	duplicate/free CXRZ CXRZnm
	CXRZnm*=ZLVDTSENS
	duplicate/free CYRZ CYRZnm
	CYRZnm*=YLVDTSENS


	duplicate/o CXRXnm $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"CenterX","X")
	duplicate/o CXRZnm $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"CenterX","Z")
	duplicate/o CYRYnm $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"CenterY","Y")
	duplicate/o CYRZnm $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"CenterY","Z")

end

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



Static Function FindMultiBreaks(DEFV,ZSensor,OutWave,[smoothing])
	variable smoothing
	wave DEFV,ZSensor,OutWave
	if(paramisdefault(smoothing))
		smoothing=501
	endif
	variable AppVel=str2num(stringbykey("ApproachVelocity",note(DefV),":","\r"))
	variable RetVel=str2num(stringbykey("RetractVelocity",note(DefV),":","\r"))
	variable/D AppPause=str2num(stringbykey("ApproachPause",note(DefV),":","\r"))
	variable/D RetPause=str2num(stringbykey("RetractPause",note(DefV),":","\r"))
	variable Repeats=str2num(stringbykey("Repeats",note(DefV),":","\r"))

	variable/D Distance=wavemax(ZSensor)-WaveMin(ZSensor)
	variable/D RetPoints=Distance/(RetVel*1e-9)/dimdelta(DEFV,0)
	variable/D AppPoints=Distance/(AppVel*1e-9)/dimdelta(DEFV,0)
	variable CycleTime=RetPoints+AppPoints+AppPause/dimdelta(DEFV,0)+RetPause/dimdelta(DEFV,0)
	variable RepeatUnit=numpnts(DefV)/Repeats
	variable endslope,startslope,partramp
	variable startsearch=0
	variable n
	variable tol=1e-9
	//DUPLICATE/FREE ZSensor zsM
	//Smooth/B (smoothing), zsM

	make/free/n=(Repeats,3) Times
	Times=0
	for(n=0;n<Repeats;n+=1)
		if(n==0)
		
			duplicate/free/R=[0,1.1*RepeatUnit] ZSensor zcUT
		else
			duplicate/free/R=[Times[n-1][2],Times[n-1][2]+1.1*RepeatUnit] ZSensor zcUT
		endif
		//	Interpolate2/S=.1E-10/T=3/N=(NUMPNTS(zcUT))/F=1/Y=ZSetSM zcUT;DelayUpdate
		duplicate/free zcUT ZSetDif
		Differentiate zcUT/D=ZSetDif
		

		//startslope=.98*wavemin(ZSetDif,pnt2x(ZSetdif,0),pnt2x(ZSetdif,RepeatUnit))
		startslope=-1*RetVel*1e-9
		//endslope=wavemax(ZSetDif,pnt2x(ZSetdif,0),pnt2x(ZSetdif,RepeatUnit))
		endslope=1*AppVel*1e-9
		duplicate/free/r=[0,100] ZSetDif ZStart

		if(wavemin(ZStart)<startslope)
			partramp=0

		else
			FindLevel/Q/P/edge=1/R=[0,CycleTime] ZSetDif,startslope
			partramp=v_levelx
		endif
		

		if(n==0)
			Times[n][0]=0

		else
			FindLevel/Q/P/edge=2/R=[partramp,0] ZSetDif,0//startslope/100
			Times[n][0]=x2pnt(ZSensor,pnt2x(zcUT,v_levelx))
		endif
		//startsearch=partramp+5
		FindLevel/Q/P/edge=1/R=[partramp] ZSetDif,0
		Times[n][1]=x2pnt(ZSensor,pnt2x(zcUT,v_levelx))
		startsearch=v_levelx+50
		FindLevel/Q/P/edge=1/R=[startsearch] ZSetDif,endslope
		FindLevel/Q/P/edge=2/R=[v_levelx] ZSetDif,0
		//
		Times[n][2]=x2pnt(ZSensor,pnt2x(zcUT,v_levelx))
			
	endfor
	make/free/n=(dimsize(Times,0)-1) Gap
	Gap=Times[p+1][1]-Times[p][1]
	wavestats/q Gap
	duplicate/free Times TimesEven
	TimesEven[0][0]=Times[0][1]
	TimesEven[0][1]=Times[0][2]
	TimesEven[0][2]=Times[1][0]-1
	TimesEven=TimesEven[0][q]+p*v_avg

	duplicate/o TimesEven OutWave
	killwaves zcUT
end


Static Function TellRelevantWaves(WaveIn)
	wave WaveIn
	string InitialRampNumber=stringbykey("Initial Ramp Number",note(waveIn),":","\r")
	Print "Initial Ramp: "+"InitRamp"+InitialRampNumber
	string InitialCenterNumber=stringbykey("Last Centered ",note(waveIn),":","\r")
	Print "Center Wave: "+"CenterWave"+InitialCenterNumber
	print "Sramp Before: "+FindAWave(WaveIn,"Sramp_D*","Before")
	print "Sramp After: "+FindAWave(WaveIn,"Sramp_D*","After")
end

Static Function FindMSingleBreaksandPauses(DEFV,ZSensor,OutWave)

	wave DEFV,ZSensor,OutWave

	duplicate/free ZSensor ZSetSm
	duplicate/free ZSetSm	ZSetDif
	Differentiate ZSetSm/D=ZSetDif

	make/free/n=(1) Times

	FindLevel/R=(5e-3)/Q/edge=1 ZSetDif 0
	if(V_flag==1)
	print nameofwave(DefV)
	endif
	Times[0]=v_levelx
	duplicate/o Times OutWave
	wave W_FindLevels
	killwaves W_FindLevels
	return 0
end


Static Function FindTouchBreaksandPauses(DEFV,ZSensor,OutWave)

	wave DEFV,ZSensor,OutWave

	duplicate/free ZSensor ZSetSm
	duplicate/free ZSetSm	ZSetDif
	Differentiate ZSetSm/D=ZSetDif

	make/free/n=(1) Times

	FindLevel/R=(5e-3)/Q/edge=2 ZSetDif 0
	if(V_flag==1)
	print nameofwave(DefV)
	endif
	Times[0]=v_levelx
	duplicate/o Times OutWave
	wave W_FindLevels
	killwaves W_FindLevels
	return 0
end



Static Function FindInitialBreaksandPauses(DEFV,ZSensor,OutWave)

	wave DEFV,ZSensor,OutWave

	duplicate/free ZSensor ZSetSm
	//	Resample/DOWN=5/N=1/WINF=None ZSetSm;DelayUpdate
	Smooth/B 11, ZSetSm
	duplicate/free ZSetSm	ZSetDif
	Differentiate ZSetSm/D=ZSetDif
	//DE_Filtering#TVD1D_denoise(ZSetDif,1e-9,ZSetDif)
	//	variable startslope=wavemin(ZSetDif)

	make/free/n=(1) Times
	 
	variable startsearch=.01*numpnts(ZSensor)
	variable n
	variable tol=20e-9
	FindValue/S=(startsearch)/T=(tol)/V=(0)  ZSetDif
		
	Times[0]=pnt2x(ZSetSm,v_value)
	duplicate/o Times OutWave
	wave W_FindLevels
	killwaves W_FindLevels
	return 0
end


Static Function AddSectionNotes(DeflVWave,WaveofLocations)
	wave DeflVWave,WaveofLocations
	variable n=0
	String Locs=""
	String Dirs=""
	String PauseLoc=""
	String PauseState=""


	for(n=0;n<(dimsize(WaveofLocations,0)-0);n+=1)
		
		Locs+=num2istr(WaveofLocations[n][0])+";"
		Dirs+="Ret;"
		Locs+=num2istr(WaveofLocations[n][2])+";"
		Dirs+="Ext;"
		PauseLoc+=num2istr(WaveofLocations[n][1])+";"
		PauseState+="Start;"
		PauseLoc+=num2istr(WaveofLocations[n][2])+";"
		PauseState+="End;"
	endfor


	String NoteString=note(DeflVWave)

	NoteString=ReplaceStringbyKey("DE_Dir",NoteString,Dirs,":","\r")
	NoteString=ReplaceStringbyKey("DE_Ind",NoteString,Locs,":","\r")
	NoteString=ReplaceStringbyKey("DE_PauseLoc",NoteString,PauseLoc,":","\r")
	NoteString=ReplaceStringbyKey("DE_PauseState",NoteString,PauseState,":","\r")
	note/k DeflVWave,NoteString
	
end






Function/S FindAWave(InputWave,WaveTypeToFind,Type)
	wave InputWave
	String WaveTypeToFind,Type
	String Result
	if(cmpstr(Type,"Before")!=0&&cmpstr(Type,"After")!=0)
		return "-1"
	endif
	DFREF saveDF = GetDataFolderDFR()
	SetDataFolder GetWavesDataFolder(InputWave, 1)
	variable TotalTime=DE_TimeFuncs#ReturnDateFromWave(InputWave)+DE_TimeFuncs#ReturnTimeFromWave(InputWave)
	String ListofWaves= WaveList(WaveTypeToFind,";","")
	variable n
	make/D/free/n=(itemsinlist(ListofWaves)) Times
	make/Free/T/n=(itemsinlist(ListofWaves)) WaveNameWave

	for(n=0;n<itemsinlist(ListofWaves);n+=1)
		Times[n]=DE_TimeFuncs#ReturnDateFromWave($stringfromlist(n,ListofWaves))+DE_TimeFuncs#ReturnTimeFromWave($stringfromlist(n,ListofWaves))
		WaveNameWave[n]=stringfromlist(n,ListofWaves)
	endfor
	SORT Times WaveNameWave
	SORT Times Times
	if(TotalTime>Times[numpnts(Times)-1])	
		if(cmpstr(Type,"Before")==0)
			Result= WaveNameWave[numpnts(WaveNameWave)-1]
		elseif(cmpstr(Type,"After")==0)
			Result= "-1"
		endif
	
	elseif(TotalTime<Times[0])
		if(cmpstr(Type,"Before")==0)
			Result= "-1"
		elseif(cmpstr(Type,"After")==0)
			Result= WaveNameWave[0]
		endif
	
	else
		FindLevel/Q Times TotalTime
		if(cmpstr(Type,"Before")==0)
			Result= WaveNameWave[floor(v_levelx)]
		elseif(cmpstr(Type,"After")==0)
			Result= WaveNameWave[ceil(v_levelx)]
		endif
	
	endif
	SetDataFolder savedf
	return Result
end

Static Function/S FindAnAssociatedForceWave(WaveIn,Type)
	wave WaveIn
	string type
	String basicWave
	DFREF savedf=GetDataFolderDFR()
	variable n=0
	string FolderName
	
	StrSwitch(Type)
	
		case "Initial":
			basicWave="InitRamp_D"+stringbykey("Initial Ramp Number",note(waveIn),":","\r")
			do
				if(strlen( wavelist(basicwave,";","")) != 0)
					Foldername=""
					break
				endif
				FolderName= GetIndexedObjNameDFR(root:, 4,n )
			
				if (strlen(FolderName) == 0)
					return ""
				endif
				DFREF Folder=$FolderName
				SetDataFolder Folder
				if(strlen( wavelist(basicwave,";","")) != 0)
					SetDataFolder savedf

					break
				endif
				SetDataFolder savedf
				n+=1
			while(1>0)
			wave FinalWave=$(":"+FolderName+":"+basicwave)
			return	GetDataFolder(1)+FindTheRawwave(FinalWave)

			break
		
		case "RampBefore":
			basicWave=FindAWave(WaveIn,"*Force_Ret","Before")
						return	GetDataFolder(1)+basicWave

			break
				
		case "RampAfter":
			basicWave=FindAWave(WaveIn,"*Force_Ret","After")
						return	GetDataFolder(1)+basicWave

			break
				
		default:
			break

	endswitch
end

Static Function DisplayAssociatedWave(InputWave,Type)
wave inputwave
string type
string ForceName=FindAnAssociatedForceWave(inputwave,Type)
string SepName=replaceString("Force",ForceName,"Sep")
wave ForceWave=$ForceName
wave SepWave=$SepName
display ForceWave vs SepWave


end

Static Function/S FindTheRawwave(InputWave)
	wave InputWave

	variable/D TargetTime=DE_TimeFuncs#ReturnDateFromWave(InputWave)+DE_TimeFuncs#ReturnTimeFromWave(InputWave)
	variable/D FulLTime
	String AllWaves=wavelist("*",";","")
	variable n=0
	for(n=0;n<itemsinlist(AllWaves);n+=1)
		wave w1=$stringfromlist(n,AllWaves)
		FulLTime=DE_TimeFuncs#ReturnDateFromWave(w1)+DE_TimeFuncs#ReturnTimeFromWave(w1)
		if(FulLTime==TargetTime)
			return nameofwave(w1)
			break
		
		endif
	
	endfor


end



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Static Function FindMultiBreaksandPausesNoNote(DEFV,ZSensor,OutWave)
//
//	wave DEFV,ZSensor,OutWave
//	duplicate/free ZSensor ZSetSm
//	//Resample/DOWN=51/N=1/WINF=None ZSetSm;DelayUpdate
//	//Smooth/B 41, ZSetSm
//	duplicate/o ZSetSm	ZSetDif
//	Differentiate ZSetSm/D=ZSetDif
//	//DE_Filtering#TVD1D_denoise(ZSetDif,50e-9,ZSetDif)
//	FindLevels/Q ZSetSm mean(ZSetSm)
//	variable numberofscans=v_levelsfound/2
//	if(round(numberofscans)!=numberofscans)
//	print "FUCK"
//	endif
//	variable spacing=floor(numpnts(ZSetSm)/numberofscans)
//	make/free/n=(numberofscans,3) Times
//	variable startslope=wavemin(ZSetDif,pnt2x(ZSetdif,0),pnt2x(ZSetdif,spacing))
//	variable endslope=wavemax(ZSetDif,pnt2x(ZSetdif,0),pnt2x(ZSetdif,spacing))
//
//	variable startsearch=0
//	variable n
//	variable tol=1e-9
//	FindLevel/Q/P/edge=1/R=[startsearch,(n+1)*spacing] ZSetDif,0
//	Times[0][0]=v_levelx
//	startsearch=v_levelx+50
//	FindLevel/Q/P/edge=1/R=[startsearch,(n+1)*spacing] ZSetDif,(endslope)
//
//	startsearch=v_levelx+50
//	FindLevel/Q/P/edge=2/R=[startsearch,(n+1)*spacing] ZSetDif,0
//	Times[0][1]=v_levelx
//	//startsearch=v_levelx
//	//FindLevel/Q/P/edge=1/R=[(n+1)*numpnts(ZSensor)/numberofscans,startsearch] ZSetDif,(0)
//	Times[0][2]=spacing
//		
//	Times=Times[0][q]+p*spacing
//		
//		 
//	if(Times[dimsize(Times,0)-1][2]==numpnts(ZSensor))
//		//Times[dimsize(Times,0)-1][2]-=1
//	endif
//	duplicate/o Times OutWave
//	wave W_FindLevels
//	killwaves W_FindLevels
//	return 0
//end


//Static Function AddSectionNotesfor500kHz(DeflIn,ZSnsrIn,wavetosteal)
//	wave DeflIn,ZSnsrIn,wavetosteal
//	String UseOther
//	Prompt UseOther,"Wanna",popup,"Yes;No;"
//	DoPrompt "Wanna steal the ramp points?",UseOther
//		
//			variable sampling=numpnts(DeflIn)/numpnts(wavetosteal)
//			string Dirs= stringbykey("DE_Dir",note(wavetosteal),":","\r")
//			string Ind= stringbykey("DE_Ind",note(wavetosteal),":","\r")
//			string PauseLoc= stringbykey("DE_PauseLoc",note(wavetosteal),":","\r")
//			string PauseState= stringbykey("DE_PauseState",note(wavetosteal),":","\r")
//			Ind=ReplaceStringListBySampling(Ind,sampling)
//			PauseLoc=ReplaceStringListBySampling(PauseLoc,sampling)
//
//			
//			string NoteString=note(DeflIn)
//			NoteString=ReplaceStringbyKey("DE_Dir",NoteString,Dirs,":","\r")
//			NoteString=ReplaceStringbyKey("DE_Ind",NoteString,Ind,":","\r")
//			NoteString=ReplaceStringbyKey("DE_PauseLoc",NoteString,PauseLoc,":","\r")
//			NoteString=ReplaceStringbyKey("DE_PauseState",NoteString,PauseState,":","\r")
//			note/k DeflIn,NoteString
//end


//Static Function/S ReplaceStringListBySampling(inString,sampling)
//	String inString
//	variable sampling
//	String Out=""
//	variable tot=itemsinlist(inString)
//	variable n=0
//	for(n=0;n<tot;n+=1)
//		 Out+=num2str(sampling*str2num(stringfromlist(n,inString)))+";"
//	endfor
//
//return Out
//end

//Static Function ProcessMulti(DefV,ZVolt,BaseString,WaveNumber)
//	wave DefV,ZVolt
//	string BaseString
//	variable WaveNumber
//	String NewNote=note(DefV)
//	NewNote=ReplaceStringByKey("BaseSuffix", NewNote, num2str(WaveNumber),":","\r")
//	note/K DefV NewNote
//	note/K ZVolt NewNote
//	make/o/n=0 Garbage1,Garbage2,IDefl,IZSen
//	GenTracesinnm(DefV,ZVolt,IDefl,IZSen,Garbage1,Garbage2,rate=10e3)
//	make/free/n=0 MBreaks,Mbreaks2
//	//FindMultiBreaksandPausesNoNote(IDefl,IZSen,MBreaks)
//	FindMultiBreaksandPausesFromNote(IDefl,IZSen,Mbreaks2)
//	AddSectionNotes(IDefl,Mbreaks2)
//	duplicate/o IDefl $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"Defl","All")
//	duplicate/o IZSen $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"ZSnsr","All")
//	killwaves Garbage1,Garbage2,IDefl,IZSen
//
//end

//Static Function Process500kHz(DefV,ZVolt,BaseString,WaveNumber)
//	wave DefV,ZVolt
//	string BaseString
//	variable WaveNumber
//	String NewNote=note(DefV)
//	NewNote=ReplaceStringByKey("BaseSuffix", NewNote, num2str(WaveNumber),":","\r")
//	note/K DefV NewNote
//	note/K ZVolt NewNote
//	make/o/n=0 Garbage1,Garbage2,IDefl,IZSen
//	
//	string tracename
//	Prompt traceName,"Wave",popup,wavelist("*Defl",";","")
//	DoPrompt "Pick a wave",traceName
//	WAVE w = $traceName
//							
//	NewNote=note(w)
//	NewNote=ReplaceStringByKey("BaseSuffix", NewNote, num2str(WaveNumber),":","\r")
//	NewNote=ReplaceStringByKey("Bandwidth", NewNote, num2str(500),":","\r")
//
//	note/K DefV NewNote
//	note/K ZVolt NewNote
//		
//	GenTracesinnm(DefV,ZVolt,IDefl,IZSen,Garbage1,Garbage2,rate=50e3)
//	make/o/n=0 MBreaks,Mbreaks2
//	FindMultiBreaksandPausesFromNote(IDefl,IZSen,Mbreaks2,smoothing=501)
//	AddSectionNotes(IDefl,Mbreaks2)
//	duplicate/o IDefl $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"Defl","Fast")
//	duplicate/o IZSen $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"ZSnsr","Fast")
//	killwaves Garbage1,Garbage2,IDefl,IZSen
//
//end
//
//Static Function ProcessEquilRamp(DefV,ZVolt,BaseString,WaveNumber)
//	wave DefV,ZVolt
//	string BaseString
//	variable WaveNumber
//	String NewNote=note(DefV)
//	NewNote=ReplaceStringByKey("BaseSuffix", NewNote, num2str(WaveNumber),":","\r")
//	note/K DefV NewNote
//	note/K ZVolt NewNote
//	make/o/n=0 Garbage1,Garbage2,IDefl,IZSen
//	GenTracesinnm(DefV,ZVolt,IDefl,IZSen,Garbage1,Garbage2)
//	duplicate/o IDefl $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"Defl","Ret")
//	duplicate/o IZSen $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"ZSnsr","Ret")
//	killwaves Garbage1,Garbage2,IDefl,IZSen
//
//end

//Static Function ProcessEquil(DefV,ZVolt,BaseString,WaveNumber)
//	wave DefV,ZVolt
//	string BaseString
//	variable WaveNumber
//	String NewNote=note(DefV)
//	NewNote=ReplaceStringByKey("BaseSuffix", NewNote, num2str(WaveNumber),":","\r")
//	note/K DefV NewNote
//	note/K ZVolt NewNote
//	make/o/n=0 Garbage1,Garbage2,IDefl,IZSen
//	GenTracesinnm(DefV,ZVolt,IDefl,IZSen,Garbage1,Garbage2)
//	duplicate/o IDefl $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"Defl","Equil")
//	duplicate/o IZSen $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"ZSnsr","Equil")
//	killwaves Garbage1,Garbage2,IDefl,IZSen
//
//end

//Static Function ProcessInitial(DefV,ZVolt,BaseString,WaveNumber)
//	wave DefV,ZVolt
//	string BaseString
//	variable WaveNumber					
//	variable cutpoint
//	String NewNote=note(DefV)
//	NewNote=ReplaceStringByKey("BaseSuffix", NewNote, num2str(WaveNumber),":","\r")
//	note/K DefV NewNote
//	note/K ZVolt NewNote
//	make/o/n=0 Garbage1,Garbage2,IDefl,IZSen
//	
//	
//	GenTracesinnm(DefV,ZVolt,IDefl,IZSen,Garbage1,Garbage2)
//	make/free/n=0 IBreaks
//	FindInitialBreaksandPauses(IDefl,IZSen,IBreaks)
//	cutpoint=x2pnt(IDefl,IBreaks[0])
//	duplicate/o/R=[0,cutpoint-1] IDefl $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"Defl","Ext")
//	duplicate/o/R=[0,cutpoint-1] IZSen $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"ZSnsr","Ext")
//	duplicate/o/R=[cutpoint,] IDefl $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"Defl","Ret")
//	duplicate/o/R=[cutpoint,] IZSen $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"ZSnsr","Ret")
//	killwaves Garbage1,Garbage2,IDefl,IZSen
//
//end

//Static Function ProcessSingle(DefV,ZVolt,BaseString,WaveNumber)
//	wave DefV,ZVolt
//	string BaseString
//	variable WaveNumber					
//	variable cutpoint
//	String NewNote=note(DefV)
//	NewNote=ReplaceStringByKey("BaseSuffix", NewNote, num2str(WaveNumber),":","\r")
//	note/K DefV NewNote
//	note/K ZVolt NewNote
//	make/o/n=0 Garbage1,Garbage2,IDefl,IZSen
//	GenTracesinnm(DefV,ZVolt,IDefl,IZSen,Garbage1,Garbage2)
//	make/free/n=0 SBreaks
//	FindMSingleBreaksandPauses(IDefl,IZSen,SBreaks)
//	cutpoint=x2pnt(IDefl,SBreaks[0])
//	if( cutpoint<0)
//		cutpoint=0
//	endif
//	duplicate/o/R=[0,cutpoint-1] IDefl $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"Defl","Ret")
//	duplicate/o/R=[0,cutpoint-1] IZSen $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"ZSnsr","Ret")
//	duplicate/o/R=[cutpoint,] IDefl $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"Defl","Ext")
//	duplicate/o/R=[cutpoint,] IZSen $"root:"+DE_Naming#StringCreate(BaseString,WaveNumber,"ZSnsr","Ext")
//	killwaves Garbage1,Garbage2,IDefl,IZSen
//
//	
//end
//Static Function FindMultiBreaksandPausesNoNote(DEFV,ZSensor,OutWave)
//
//	wave DEFV,ZSensor,OutWave
//	duplicate/free ZSensor ZSetSm
//	//Resample/DOWN=51/N=1/WINF=None ZSetSm;DelayUpdate
//	//Smooth/B 41, ZSetSm
//	duplicate/o ZSetSm	ZSetDif
//	Differentiate ZSetSm/D=ZSetDif
//	//DE_Filtering#TVD1D_denoise(ZSetDif,50e-9,ZSetDif)
//	FindLevels/Q ZSetSm mean(ZSetSm)
//	variable numberofscans=v_levelsfound/2
//	if(round(numberofscans)!=numberofscans)
//	print "FUCK"
//	endif
//	variable spacing=floor(numpnts(ZSetSm)/numberofscans)
//	make/free/n=(numberofscans,3) Times
//	variable startslope=wavemin(ZSetDif,pnt2x(ZSetdif,0),pnt2x(ZSetdif,spacing))
//	variable endslope=wavemax(ZSetDif,pnt2x(ZSetdif,0),pnt2x(ZSetdif,spacing))
//
//	variable startsearch=0
//	variable n
//	variable tol=1e-9
//	FindLevel/Q/P/edge=1/R=[startsearch,(n+1)*spacing] ZSetDif,0
//	Times[0][0]=v_levelx
//	startsearch=v_levelx+50
//	FindLevel/Q/P/edge=1/R=[startsearch,(n+1)*spacing] ZSetDif,(endslope)
//
//	startsearch=v_levelx+50
//	FindLevel/Q/P/edge=2/R=[startsearch,(n+1)*spacing] ZSetDif,0
//	Times[0][1]=v_levelx
//	//startsearch=v_levelx
//	//FindLevel/Q/P/edge=1/R=[(n+1)*numpnts(ZSensor)/numberofscans,startsearch] ZSetDif,(0)
//	Times[0][2]=spacing
//		
//	Times=Times[0][q]+p*spacing
//		
//		 
//	if(Times[dimsize(Times,0)-1][2]==numpnts(ZSensor))
//		//Times[dimsize(Times,0)-1][2]-=1
//	endif
//	duplicate/o Times OutWave
//	wave W_FindLevels
//	killwaves W_FindLevels
//	return 0
//end
