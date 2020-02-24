#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma modulename=DE_RupPntFinder

Menu "Ramp"
	"Open PntFinder", RupturePntFinder()

	
end

Static Function/S RemovePauses(ForceIn,SepIn,ForceOut,Sepout)
	wave ForceIn,SepIn,ForceOut,Sepout
	string NoteString=note(ForceIn)
	String PauseLocString=stringbykey("DE_PauseLoc",NoteString,":","\r")
	variable last=itemsinlist(PauseLocString)
	variable n,startdelete,enddelete,todelete
	string DeletedString,DeletedStringNew
	DeletedString=""
	duplicate/free ForceIn FreeForce
	duplicate/free SepIn FreeSep

	for(n=last-2;n>=0;n-=2)
	
		startdelete=str2num(stringfromlist(n,PauseLocString))
		enddelete=str2num(stringfromlist(n+1,PauseLocString))
		if(n==(last-2))
			todelete=enddelete-startdelete
			DeletedString=num2str(todelete)


		endif
		deletepoints startdelete,todelete, FreeForce,FreeSep
					DeletedStringNew=num2str(startdelete-todelete*n/2)+";"+DeletedString
		DeletedString=DeletedStringNew
	endfor
	duplicate/o FreeForce ForceOut
	duplicate/o FreeSep Sepout
	return DeletedString
end

Static Function PythonFitter( UseWave,Method,Threshold,AMount)//Variables demanded by MarkovFit Code
	Wave UseWave//Input wave
	String Method
	Variable Threshold,AMount
	String Destination = "C:\Data\StepData\Test1.ibw"
	String NewHome = "C:\Data\StepData\Shit.txt"

	Save/O/C UseWave as Destination
	String BasePythonCommand = "cmd.exe /C activate & python C:\Devin\Python\StepAttempt\StepAttempt.py "
	String MethodCommand="-method "+ method +" "
	String InputCom="-inputfile "+ Destination+" "
	String OutputCom="-outputfile "+ NewHome+" "
	String SmoothCommand="-smooth "+ num2str(Amount)+" "
	String ThreshCommand="-threshold "+ num2str(threshold)+" "

	String PythonCommand=BasePythonCommand+MethodCommand+InputCom+OutputCom+SmoothCommand+ThreshCommand
	print PythonCommand
	ModOperatingSystemUtil#execute_python(PythonCommand)
	LoadWave/O/N/G/D NewHome

end

Static Function FitPython()
	string saveDF
	variable FOffset, Soffset

	saveDF = GetDataFolder(1)
	wave/T parmWave=root:DE_RupturePntFinder:MenuStuff:PyParmWave

	controlinfo/W=RupPntFinder de_RupPnt_popup0
	SetDataFolder s_value
	controlinfo/W=RupPntFinder de_RupPnt_popup1
	wave ForceWave=$S_value
	wave SepWave=$ReplaceString("Force",S_value,"Sep")
	controlinfo/W=RupPntFinder de_RupPnt_check3
	variable nodelay=v_value
	//	make/o/n=0 RupPntU,RupPntD
	make/o/n=0 ForceWaveS,SepWaveS
	//DE_Filtering#FilterForceSep(ForceWave,SepWave,ForceWaveS,SepWaveS,"SVG",str2num(parmWave[6][1]))
	make/free/n=0 ForceOut,Sepout
	if(nodelay==1)
	string deletionString=RemovePauses(ForceWave,SepWave,ForceOut,Sepout)
	DE_Filtering#FilterForceSep(ForceOut,Sepout,ForceWaveS,SepWaveS,"TVD",str2num(parmWave[0][1]))
	else
	DE_Filtering#FilterForceSep(ForceWave,SepWave,ForceWaveS,SepWaveS,"TVD",str2num(parmWave[0][1]))
	endif
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
	controlinfo/W=RupPntFinder de_RupPnt_popup2
	PythonFitter( ForceWaveS,S_Value,str2num(parmWave[3][1]),str2num(parmWave[2][1]))
	wave no0=wave0
	make/o/n=0 UpP,DownP;DE_RuptureRamp#SortSteps(ForceWaveS,no0,UpP,DownP,20)
	DE_RupPntFinder#RefineSteps(ForceWaveS,UpP,DownP)
	duplicate/o UpP RupPntU
	RupPntU=x2pnt(ForceWave,pnt2x(ForceWaveS,UpP))
	duplicate/o DownP RupPntD
	RupPntD=x2pnt(ForceWave,pnt2x(ForceWaveS,DownP))
	string PauseIndices=stringbykey("DE_PauseLoc",note(ForceWave),":","\r")//
	if(nodelay==1)
	RupPntU=CorrectforRemovedPause(RupPntU,PauseIndices,deletionString)
	RupPntD=CorrectforRemovedPause(RupPntD,PauseIndices,deletionString)
	else
	endif

	killwaves SepWaveS,ForceWaveS,no0,DownP,UpP
	SetDataFolder saveDF
end

static Function CorrectforRemovedPause(PointsIn,PauseIndices,deletionString)
	variable PointsIn
	String PauseIndices,deletionString
	variable tot=itemsinlist(deletionString)
	variable deletesize=str2num(stringfromlist(tot-1,deletionString))
	variable n,CurrentPnt,TotAdd
	for(n=0;n<tot;n+=1)
		CurrentPnt=str2num(stringfromlist(n,deletionString))
		if(PointsIn<CurrentPnt)
			return PointsIn+n*deletesize
		else
		endif
	endfor
	return -1
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

Static Function ButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			strswitch( ba.ctrlName)
				case "de_RupPnt_button0":
					FitPython()
					break

			endswitch
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

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

Window RupturePntFinder() 
	PauseUpdate; Silent 1		// building window...
	NewPanel/N=RupPntFinder /W=(0,0,300,175)
	NewDataFolder/o root:DE_RupturePntFinder
	NewDataFolder/o root:DE_RupturePntFinder:MenuStuff

	DE_RupPntFinder#UpdateParmWave()
	Button de_RupPnt_button0,pos={10,33},size={150,20},proc=DE_RupPntFinder#ButtonProc,title="Pythong!"
	PopupMenu de_RupPnt_popup0,pos={10,2},size={129,21},title="Folder",mode=1
	PopupMenu de_RupPnt_popup0,mode=1,popvalue="X",value= #"DE_PanelProgs#ListFolders()"
	PopupMenu de_RupPnt_popup1,pos={157,2},size={129,21},title="Force Wave"
	PopupMenu de_RupPnt_popup1,mode=1,popvalue="X",value= #"DE_RupPntFinder#ListWaves(\"de_RupPnt_popup0\",\"*Force\")"
	PopupMenu de_RupPnt_popup2,pos={208,65},size={129,21},title="Method"
	PopupMenu de_RupPnt_popup2,mode=1,popvalue="X",value= "gauss;ms;"

	ListBox DE_RupPnt_list1,pos={7,70},size={175,75},proc=DE_RupPntFinder#ListBoxProc,listWave=root:DE_RupturePntFinder:MenuStuff:PyParmWave
	ListBox DE_RupPnt_list1,selWave=root:DE_RupturePntFinder:MenuStuff:PySelWave,editStyle= 2,userColumnResize= 1,widths={70,40,70,40}

	CheckBox de_RupPnt_check3 title="No Delays",pos={202,36},size={150,20},proc=DE_RuptureRamp#CheckProc

//	ControlUpdate/A/W=RupRampPanel
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
	
	if(exists("root:DE_RupturePntFinder:MenuStuff:PyParmWave")==1)
		wave/t/z Par=root:DE_RupturePntFinder:MenuStuff:ParmWave
		wave/z Sel=root:DE_RupturePntFinder:MenuStuff:SelWave
	Else
		make/t/n=(4,2) root:DE_RupturePntFinder:MenuStuff:PyParmWave
		wave/t/z Par=root:DE_RupturePntFinder:MenuStuff:PyParmWave
		make/n=(4,2) root:DE_RupturePntFinder:MenuStuff:PySelWave
		wave/z Sel=root:DE_RupturePntFinder:MenuStuff:PySelWave
		
		Par[0][0]={"Smoothing","Decimating","GaussParm","Threshhold"}
		Par[0][1]={"50e-9","1","5","0.1"}
		Sel[][0]=0
		Sel[][1]=2
	endif


end


//Static Function FitHMMButt()

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
					
					
//end

