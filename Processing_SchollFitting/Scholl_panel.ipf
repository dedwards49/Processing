#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma ModuleName = DE_Schollpanel

#include "SchollFitting"
#include "FEC_Fast_Corr"
//#include "DE_Feather"
#include "DE_Filtering"
#include ":\Misc_PanelPrograms\Panel Progs"
#include ":\Misc_PanelPrograms\AsylumNaming"
#include "FastCorrection"


Static Function PMP(ctrlName,popNum,popStr) : PopupMenuControl
	String ctrlName
	Variable popNum
	String popStr

End

Static Function CP(ctrlName,checked) : CheckBoxControl
	String ctrlName
	Variable checked

End

Static Function LBP(ctrlName,row,col,event) : ListBoxControl

	String ctrlName
	Variable row
	Variable col
	Variable event	//1=mouse down, 2=up, 3=dbl click, 4=cell select with mouse or keys
					//5=cell select with shift key, 6=begin edit, 7=end
	wave/t/z FL=root:SchollPanel:FolderList
	wave/t/z LW1=root:SchollPanel:ListWave1
	wave/t/z LW2=root:SchollPanel:ListWave2
	
	DFREF saveDFR = GetDataFolderDFR()	
	if(numpnts(FL)==0)
	return 0
	
	endif
	DFREF Folder=$ FL[row]
	SetDataFolder Folder
		
	controlinfo de_scholl_setvar1
	string search1= S_Value
	controlinfo de_scholl_setvar2
	string search2= S_Value
	string/g CurrFold

	switch (event)
		case -1:
			break
		case 1:
			break
		case 3:
			DE_PanelProgs#ListWaves($Fl[row],Search1,LW1)
			DE_PanelProgs#ListWaves($Fl[row],Search2,LW2)
			SetVariable de_scholl_setvar8,value= _NUM:GuessRamps()
			CurrFold=FL[row]
			
		case 4:	
			DE_PanelProgs#ListWaves($Fl[row],Search1,LW1)
			DE_PanelProgs#ListWaves($Fl[row],Search2,LW2)
			SetVariable de_scholl_setvar8,value= _NUM:GuessRamps()

			CurrFold=FL[row]
	endswitch
	SetDataFolder saveDFR
	return 0
End
//
//



Static Function LBP_1(ctrlName,row,col,event) : ListBoxControl

	String ctrlName
	Variable row
	Variable col
	Variable event	//1=mouse down, 2=up, 3=dbl click, 4=cell select with mouse or keys
	//5=cell select with shift key, 6=begin edit, 7=end
	if(event==10)
		return 0
	endif
	wave/t/z FL=root:SchollPanel:FolderList
	controlinfo de_scholl_list0
	DFREF saveDFR = GetDataFolderDFR()	
	DFREF Folder=$ FL[v_value]
		SetDataFolder Folder

	wave/t/z LW1=root:SchollPanel:ListWave1
	

	make/o/n=0  root:schollpanel:DefV_1, root:schollpanel:ZSnsr_1
	struct ForceWave Name1
	DE_Naming#WavetoStruc(LW1[row],Name1)
	controlinfo de_scholl_setvar5
	string ZsnsrName=DE_Naming#StringCreate(Name1.Name,Name1.VNum,S_Value,Name1.SDirec)
	wave w1=$LW1[row]
	wave w2=$ZsnsrName
	duplicate/o w1 root:schollpanel:DefV_1
	duplicate/o w2 root:schollpanel:Zsnsr_1
	controlinfo de_scholl_setvar3
	Cursor/W=SchollPanel#Plot1/P A,DefV_1,V_Value
	controlinfo de_scholl_setvar4
	Cursor/W=SchollPanel#Plot1/P B,DefV_1,V_Value
	DoUpdate

	SetDataFolder saveDFR
	return 0
End

Static Function LBP_2(ctrlName,row,col,event) : ListBoxControl

	String ctrlName
	Variable row
	Variable col
	Variable event	//1=mouse down, 2=up, 3=dbl click, 4=cell select with mouse or keys
	//5=cell select with shift key, 6=begin edit, 7=end
	if(event==10)
		return 0
	elseif(event==4)
	wave/t/z FL=root:SchollPanel:FolderList
	controlinfo de_scholl_list0
	DFREF saveDFR = GetDataFolderDFR()	
	DFREF Folder=$ FL[v_value]
		SetDataFolder Folder

		wave/t/z LW2=root:SchollPanel:ListWave2
	
		make/o/n=0  root:schollpanel:DefV_2, root:schollpanel:ZSnsr_2

		struct ForceWave Name1
		DE_Naming#WavetoStruc(LW2[row],Name1)
		controlinfo de_scholl_setvar6

		string ZsnsrName=DE_Naming#StringCreate(Name1.Name,Name1.VNum,S_Value,Name1.SDirec)
		string ForceName=DE_Naming#StringCreate(Name1.Name,Name1.VNum,"Force",Name1.SDirec)
		string SepName=DE_Naming#StringCreate(Name1.Name,Name1.VNum,"Sep",Name1.SDirec)


		wave w1=$LW2[row]
		wave w2=$ZsnsrName
		wave w3=$ForceName
		wave w4=$SepName
		wave DrUp=root:Schollpanel:DispRup
		wave DRupX=root:SchollPanel:DispXRup
		if(WaveExists(w3))
			ModifyGraph/W=SchollPanel#Plot2 hideTrace(DispRup)=0
			duplicate/o w3 root:schollpanel:CorForce
			duplicate/o w4 root:schollpanel:CorSep
			DE_SChollpanel#UpdateFilter()
						

		else 
			wave w5=root:schollpanel:CorForce
	
			wave w6=	root:schollpanel:CorSep
			wave w7=	root:schollpanel:DispRup
			wave w8=root:schollpanel:CorForceFil


			w6=Nan
			w7=Nan
			w5=Nan
			w8=Nan

		endif
		UpdatePRHEvents()
		duplicate/o w1 root:schollpanel:DefV_2
		duplicate/o w2 root:schollpanel:Zsnsr_2
	endif
SetDataFolder saveDFR
	return 0
End					
		
Static Function LBP_3(ctrlName,row,col,event) : ListBoxControl

	String ctrlName
	Variable row
	Variable col
	Variable event	//1=mouse down, 2=up, 3=dbl click, 4=cell select with mouse or keys
	//5=cell select with shift key, 6=begin edit, 7=end
	if(event==7)

		wave w1=root:Schollpanel:CorForceFil
		wave w2=root:Schollpanel:CorSepFil
		DE_SChollpanel#UpdateFilter()

	endif
	
	return 0
End			

Static Function LBP_4(ctrlName,row,col,event) : ListBoxControl

	String ctrlName
	Variable row
	Variable col
	Variable event	//1=mouse down, 2=up, 3=dbl click, 4=cell select with mouse or keys
	//5=cell select with shift key, 6=begin edit, 7=end

End		
		
Static Function SVP(ctrlName,varNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable varNum
	String varStr
	String varName
	strswitch(ctrlName)
		case "de_scholl_setvar0":
			DE_PanelProgs#PrintAllFolders(root:SchollPanel:FolderList,varstr)
			break
		case "de_scholl_setvar3":
			Cursor/W=SchollPanel#Plot1/P A,DefV_1,varnum
			doupdate
			break
		case "de_scholl_setvar4":

			Cursor/W=SchollPanel#Plot1/P B,DefV_1,varnum
			doupdate
	
			break
		case "de_scholl_setvar7":
			
			Variable NewNumber=floor(varNum/2)*2+1
			SetVariable de_scholl_setvar7, value= _NUM:NewNumber
			DE_SChollpanel#Updatefilter()
			break
	endswitch
End		


Static Function BP_1(ctrlName) : ButtonControl
	String ctrlName
	controlinfo de_scholl_list1
	variable row1=V_Value
	controlinfo de_scholl_list2
	variable row2=V_Value
	controlinfo de_scholl_setvar3
	variable start=V_Value
	controlinfo de_scholl_setvar4
	variable ends=V_Value
	wave/t/z LW1=root:SchollPanel:ListWave1
	wave/t/z LW2=root:SchollPanel:ListWave2
	struct ForceWave Name1
	DE_Naming#WavetoStruc(LW1[row1],Name1)
	struct ForceWave Name2
	DE_Naming#WavetoStruc(LW2[row2],Name2)
	controlinfo de_scholl_setvar5
	string ZsnsrName1=DE_Naming#StringCreate(Name1.Name,Name1.VNum,S_Value,Name1.SDirec)
	controlinfo de_scholl_setvar6
	string ZsnsrName2=DE_Naming#StringCreate(Name2.Name,Name2.VNum,S_Value,Name2.SDirec)

	
	wave w1=$LW1[row1]
	wave w2=$ZsnsrName1

	wave w3=$LW2[row2]
	wave w4=$ZsnsrName2




	string NewName
	string ForceName
	string SepName
	NewName=DE_Naming#StringCreate(Name2.Name,Name2.VNum,"DeflCor",Name2.SDirec)
	ForceName=DE_Naming#StringCreate(Name2.Name,Name2.VNum,"Force",Name2.SDirec)
	SepName=DE_Naming#StringCreate(Name2.Name,Name2.VNum,"Sep",Name2.SDirec)

	controlinfo de_scholl_popup0

	CorrFEC(w1,w2,w3,w4,start,ends,NewName=NewName,ForceName=ForceName,SepName=SepName,FitType=S_Value,ResName="Resids")
	
	wave w5=$ForceName
	wave w6=$SepName
	wave w7=Resids
	wave NEwDeflWave=$NEwName
	
	variable OffSetGuess= -w2[ends]+nEwDeflWave[ends]-5e-9
	string CurrentNote=note(w5)
	CurrentNote=ReplaceStringbyKey("DE_SchollOffset",CurrentNote,num2str(Offsetguess),":","\r")
	note/K NEwDeflWave, CurrentNote
		note/K w5, CurrentNote
	note/K w6, CurrentNote

	duplicate/o w7 root:schollpanel:FitResiduals
	killwaves w7
	duplicate/o w2 root:schollpanel:FitResidualsX
	duplicate/o w5 root:schollpanel:CorForce
	duplicate/o w6 root:schollpanel:CorSep
	DE_SChollpanel#Updatefilter()


End

//Static Function BP_2(ctrlName) : ButtonControl
//	String ctrlName
//
//	wave/t/z LW1=root:SchollPanel:ListWave1
//	wave/t/z LW2=root:SchollPanel:ListWave2
//	controlinfo de_scholl_list1
//	variable row1=V_Value
//	controlinfo de_scholl_list2
//	variable row2=V_Value
//	controlinfo de_scholl_popup1
//	string Experiment=S_Value
//	
//	struct ForceWave Name1
//	DE_Naming#WavetoStruc(LW1[row1],Name1)
//	variable ExtCurveNumber=Name1.VNum
//
//	//OK, up to here I was trying to be general, at this point...I need to assume that this has the format of a standard glide experiment. 
//	//in truth I need to think of a more generic way to define the events I want to study. 
//	
//	
//	DE_Feather#FitWaveGlide(ExtCurveNumber)
//	DE_Feather#ConvertToSingle(ExtCurveNumber)
//	UpdatePRHEvents()
//
//End


Static Function BP_3(ctrlName) : ButtonControl
	
	String ctrlName
	controlinfo de_scholl_list1
	variable row1=V_Value

	controlinfo de_scholl_setvar3
	variable start=V_Value
	controlinfo de_scholl_setvar4
	variable ends=V_Value
	
	
	wave/t/z LW1=root:SchollPanel:ListWave1
	wave/t/z LW2=root:SchollPanel:ListWave2

	struct ForceWave Name1
	DE_Naming#WavetoStruc(LW1[row1],Name1)


	controlinfo de_scholl_setvar5
	string ZsnsrName1=DE_Naming#StringCreate(Name1.Name,Name1.VNum,S_Value,Name1.SDirec)

	
	wave w1=$LW1[row1]
	wave w2=$ZsnsrName1
	make/o/n=(numpnts(LW2)) root:SchollPanel:SelWave2
	wave SW2=root:SchollPanel:SelWave2
	SW2=0
	NewPanel /W=(187,368,437,831) as "Devin"
	ListBox de_scholl_list3,pos={000,20},size={150,300},proc=DE_Schollpanel#LBP_4,listWave=root:SchollPanel:ListWave2
	ListBox de_scholl_list3,row= 0,selWave= SW2,selRow= 0,mode= 4
	DoWindow/C tmp_Select // Set to an unlikely name

	BatchSelect("tmp_Select",0)
	make/o/n=0 ResultWave
	controlinfo de_scholl_popup0

	DE_FECWiggle#FindWiggleParms(w1,w2, ResultWave,start,ends,FitType=S_Value,ResName="Resids")
	wave w7=Resids

	duplicate/o w7 root:schollpanel:FitResiduals
	duplicate/o w2 root:schollpanel:FitResidualsX
	killwaves w7


	variable n
	
	for(n=0;n<numpnts(SW2);n+=1)
	
		if(SW2[n]==1)
			struct ForceWave Name2

			DE_Naming#WavetoStruc(LW2[n],Name2)
			string NewName
			string ForceName
			string SepName
			controlinfo de_scholl_setvar6
			string ZsnsrName2=DE_Naming#StringCreate(Name2.Name,Name2.VNum,S_Value,Name2.SDirec)
			wave w3=$LW2[n]
			wave w4=$ZsnsrName2
			NewName=DE_Naming#StringCreate(Name2.Name,Name2.VNum,"DeflCor",Name2.SDirec)
			ForceName=DE_Naming#StringCreate(Name2.Name,Name2.VNum,"Force",Name2.SDirec)
			SepName=DE_Naming#StringCreate(Name2.Name,Name2.VNum,"Sep",Name2.SDirec)
			controlinfo de_scholl_popup0

			DE_FECWiggle#SingleProcess(Resultwave,w2,w3,w4,start,ends,NewName=NewName,ForceName=ForceName,SepName=SepName,FitType=S_Value)
		else
		endif
	
	endfor
	killwaves resultwave
End

Static Function SequenceManyFromOne([starttrace,endtrace,FitType])

	variable starttrace,endtrace
	String FitType
		string AllDeflEXtList= wavelist("*Defl_ext",";","")

	variable starting
	if(Paramisdefault(starttrace))
	starttrace=0
	else 
	endif
		if(Paramisdefault(endtrace))
	endtrace=itemsinlist(AllDeflEXtList)
	else 
	endif
	if(endtrace>itemsinlist(AllDeflEXtList))
		endtrace=itemsinlist(AllDeflEXtList)
	endif
	if(Paramisdefault(FitType))
	FitType="LinSin"
	else 
	endif
	
	String ForceWaveList="",SepWaveList=""
	variable n,firstdimdelta,currentdimdelta
	wave DeflExtWave=$stringfromlist(starttrace,AllDeflEXtList)
		wave DeflRetWave=$replacestring("Ext",nameofwave(DeflExtWave),"Ret")
		wave ZExtWave=$replacestring("Defl",nameofwave(DeflExtWave),"ZSnsr")
		wave ZRetWave=$replacestring("Defl",nameofwave(DeflRetWave),"ZSnsr")
		firstdimdelta=dimdelta(DeflRetWave,0)
		dowindow SChollPlot
			if(V_Flag==1)
				killwindow SchollPlot
			endif

			display/N=SchollPlot DeflExtWave vs ZExtWave

			if (DE_Schollpanel#UserCursorAdjust("SchollPlot",0) != 0)
				return -1
			endif
			variable startpnt=pcsr(A,"SchollPlot")
			variable endpnt=pcsr(B,"SchollPlot")
			killwindow SchollPlot
			//variable 	timerRefNum = StartMSTimer


	
	for(n=starttrace;n<endtrace;n+=1)
		wave DeflExtWave=$stringfromlist(n,AllDeflEXtList)
		wave DeflRetWave=$replacestring("Ext",nameofwave(DeflExtWave),"Ret")
		wave ZExtWave=$replacestring("Defl",nameofwave(DeflExtWave),"ZSnsr")
		wave ZRetWave=$replacestring("Defl",nameofwave(DeflRetWave),"ZSnsr")
		
		string NewName=ReplaceString("Defl",nameofwave(DeflExtWave),"DeflCor")
		string ForceName=ReplaceString("Defl",nameofwave(DeflExtWave),"Force")
		string SepName=ReplaceString("Defl",nameofwave(DeflExtWave),"Sep")
		
		if(abs(dimdelta(DeflRetWave,0)-firstdimdelta)>1e-5)
		print/D firstdimdelta
		print/D dimdelta(DeflRetWave,0)
		print/D dimdelta(DeflRetWave,0)-firstdimdelta
		print nameofwave(DeflRetWave)
		return 0
		endif

		CorrFEC(DeflExtWave,ZExtWave,DeflExtWave,ZExtWave,startpnt,endpnt,NewName=NewName,ForceName=ForceName,SepName=SepName,FitType=FitType,ResName="Garbage")
		wave G=Garbage
		killwaves G
		wave w5=$ForceName
		wave w6=$SepName
		wave w7=Resids
		wave NEwDeflWave=$NEwName
	
		variable OffSetGuess= -ZExtWave[endpnt]+nEwDeflWave[endpnt]-5e-9
		print nameofwave(ZExtWave)+":"+num2str(OffSetGuess)
		string CurrentNote=note(w5)
		print ForceName
		CurrentNote=ReplaceStringbyKey("DE_SchollOffset",CurrentNote,num2str(Offsetguess),":","\r")
		note/K NEwDeflWave, CurrentNote
		note/K w5, CurrentNote
		note/K w6, CurrentNote
		
		NewName=ReplaceString("Defl",nameofwave(DeflRetWave),"DeflCor")
		ForceName=ReplaceString("Defl",nameofwave(DeflRetWave),"Force")
		SepName=ReplaceString("Defl",nameofwave(DeflRetWave),"Sep")

		CorrFEC(DeflExtWave,ZExtWave,DeflRetWave,ZRetWave,startpnt,endpnt,NewName=NewName,ForceName=ForceName,SepName=SepName,FitType=FitType,ResName="Garbage")
		wave G=Garbage
		killwaves G
		wave w5=$ForceName
		wave w6=$SepName
		wave w7=Resids
		wave NEwDeflWave=$NEwName
	
		 CurrentNote=note(w5)
		CurrentNote=ReplaceStringbyKey("DE_SchollOffset",CurrentNote,num2str(Offsetguess),":","\r")
		note/K NEwDeflWave, CurrentNote
		note/K w5, CurrentNote
		note/K w6, CurrentNote
	endfor
	
	//variable microSeconds = StopMSTimer(timerRefNum)
	//Print microSeconds/1e6, "seconds Total Run"
	//Print microSeconds/1e6/top, "seconds per iteration"


end

Static Function UserCursorAdjust(graphName,autoAbortSecs)
	String graphName
	Variable autoAbortSecs
	DoWindow/F $graphName // Bring graph to front
	if (V_Flag == 0) // Verify that graph exists
		Abort "UserCursorAdjust: No such graph."
		return -1
	endif
	ShowInfo/W=$graphName
	NewPanel /K=2 /W=(187,368,637,531) as "Pause for Cursor"
	DoWindow/C tmp_PauseforCursor // Set to an unlikely name
	AutoPositionWindow/E/M=1/R=$graphName // Put panel near the graph
	DrawText 21,20,"Adjust the cursors and then"
	DrawText 21,40,"Click Continue."
	Button button0,pos={80,58},size={92,20},title="Continue"
	Button button0,proc=DE_SchollPanel#UserCursorAdjust_ContButtonProc
//	
//	PopupMenu pop0,pos={250,58},size={92,20},title="Garbage"
//	PopupMenu pop0,proc=DE_MultiFEC#PopMenuProc,value= MakeStringList()
//	
//	Button button1,pos={250,88},size={92,20},title="Fix That"
//	Button button1,proc=DE_MultiFEC#UpdateAPoint
//	Button button2,pos={250,110},size={92,20},title="Delete That"
//	Button button2,proc=DE_MultiFEC#DeleteButton
//	Button button3,pos={250,135},size={92,20},title="Add Here"
//	Button button3,proc=DE_MultiFEC#AddButton
	
	Variable didAbort= 0
	if( autoAbortSecs == 0 )
		PauseForUser tmp_PauseforCursor,$graphName
	else
		SetDrawEnv textyjust= 1
		DrawText 162,103,"sec"
		SetVariable sv0,pos={48,97},size={107,15},title="Aborting in "
		SetVariable sv0,limits={-inf,inf,0},value= _NUM:10
		Variable td= 10,newTd
		Variable t0= ticks
		Do
			newTd= autoAbortSecs - round((ticks-t0)/60)
			if( td != newTd )
				td= newTd
				SetVariable sv0,value= _NUM:newTd,win=tmp_PauseforCursor
				if( td <= 10 )
					SetVariable sv0,valueColor= (65535,0,0),win=tmp_PauseforCursor
				endif
			endif
			if( td <= 0 )
				DoWindow/K tmp_PauseforCursor
				didAbort= 1
				break
			endif
			PauseForUser/C tmp_PauseforCursor,$graphName
		while(V_flag)
	endif
	return didAbort
End
Static Function UserCursorAdjust_ContButtonProc(ctrlName) : ButtonControl
	String ctrlName
	DoWindow/K tmp_PauseforCursor // Kill panel
End

Static Function BatchSelect(graphName,autoAbortSecs)
	String graphName
	Variable autoAbortSecs
	DoWindow/F $graphName // Bring graph to front
	if (V_Flag == 0) // Verify that graph exists
		Abort "UserCursorAdjust: No such graph."
		return -1
	endif

	NewPanel /K=2 /W=(187,368,437,531) as "Pause for Cursor"
	DoWindow/C tmp_PauseforCursor // Set to an unlikely name
	AutoPositionWindow/E/M=1/R=$graphName // Put panel near the graph
	DrawText 21,20,"Adjust the cursors and then"
	DrawText 21,40,"Click Continue."
	Button button0,pos={80,58},size={92,20},title="Continue"
	Button button0,proc=DE_Schollpanel#BatchSelect_ContButtonProc
	Variable didAbort= 0
	if( autoAbortSecs == 0 )
		PauseForUser tmp_PauseforCursor,$graphName
	else
		SetDrawEnv textyjust= 1
		DrawText 162,103,"sec"
		SetVariable sv0,pos={48,97},size={107,15},title="Aborting in "
		SetVariable sv0,limits={-inf,inf,0},value= _NUM:10
		Variable td= 10,newTd
		Variable t0= ticks
		Do
			newTd= autoAbortSecs - round((ticks-t0)/60)
			if( td != newTd )
				td= newTd
				SetVariable sv0,value= _NUM:newTd,win=tmp_PauseforCursor
				if( td <= 10 )
					SetVariable sv0,valueColor= (65535,0,0),win=tmp_PauseforCursor
				endif
			endif
			if( td <= 0 )
				DoWindow/K tmp_PauseforCursor
				didAbort= 1
				break
			endif
			PauseForUser/C tmp_PauseforCursor,$graphName
		while(V_flag)
	endif
	return didAbort
End


Static Function BatchSelect_ContButtonProc(ctrlName) : ButtonControl
	String ctrlName
	print "RUNNING"
	DoWindow/K tmp_PauseforCursor // Kill panel
	DoWindow/K tmp_Select // Kill panel

End






Window SchollPanel() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel/N=SchollPanel /W=(10,10,1500,700)
	SetWindow SchollPanel, hook(MyHook)=DE_SChollpanel#WindowHook
	newdatafolder/s/o root:SchollPanel
	make/d/o/n=0 DefV_1,ZSnsr_1,DefV_2,ZSnsr_2,CorForce,CorSep,FitResiduals,FitResidualsX,CorForceFil,CorSepFil,root:SchollPanel:DispRup,root:SchollPanel:DispXRup,root:SchollPanel:ZeroX,root:SchollPanel:ZeroY
	make/o/t/n=0 FolderList,ListWave1,Listwave2


	DE_PanelProgs#PrintAllFolders(root:SchollPanel:FolderList,"Image*")

	Button button0,pos={150,600},size={125,20},proc=DE_SChollpanel#BP_1,title="Fit and Correct BG"
	Button button1,pos={560,0},size={125,16},proc=DE_SChollpanel#BP_2,title="PRH Events"
	Button button2,pos={150,630},size={125,20},proc=DE_SChollpanel#BP_3,title="Batch Process"

	ListBox de_scholl_list0,pos={10,40},size={100,100},proc=DE_SChollpanel#LBP,listWave=root:SchollPanel:FolderList
	ListBox de_scholl_list0,row= 0,mode= 1,selRow= 0
	ListBox de_scholl_list1,pos={120,40},size={150,150},proc=DE_SChollpanel#LBP_1,listWave=root:SchollPanel:ListWave1
	ListBox de_scholl_list1,row= 0,mode= 1,selRow= 0
	
	ListBox de_scholl_list2,pos={280,40},size={150,150},proc=DE_SChollpanel#LBP_2,listWave=root:SchollPanel:ListWave2
	ListBox de_scholl_list2,row= 0,mode= 1,selRow= 0

	SetDrawEnv fillpat= 0, linethick= 3.00;DelayUpdate
	DrawRect -10,-10,440,200
	SetDrawEnv fillpat= 0, linethick= 3.00;DelayUpdate
	DrawRect -10,200,440,800

	PopupMenu de_scholl_popup0,pos={280,207},size={85,20},proc=DE_SChollpanel#PMP,title="Fit Type"
	PopupMenu de_scholl_popup0,mode=1,popvalue="LinSin",value= #"\"LinSin;Poly;Spline;None\""
	
	PopupMenu de_scholl_popup1,pos={760,0},size={125,16},proc=DE_SChollpanel#PMP,title="DataType"
	PopupMenu de_scholl_popup1,mode=1,popvalue="Glide",value= #"\"Glide;MultiRamp;\""
		
	SetVariable de_scholl_setvar0,pos={10,8},size={100,16},proc=DE_SChollpanel#SVP,value= _STR:"Image*",title="Name"
	SetVariable de_scholl_setvar1,pos={120,0},size={150,16},proc=DE_SChollpanel#SVP,value= _STR:"Image*Defl_Ext",title="Y Suff1"
	SetVariable de_scholl_setvar2,pos={280,0},size={150,16},proc=DE_SChollpanel#SVP,value= _STR:"Image*Defl_Ret",title="Y Suffix 2"
	
	SetVariable de_scholl_setvar3,pos={10,210},size={100,20},proc=DE_SChollpanel#SVP,value= _NUM:0,title="BG Start"
	SetVariable de_scholl_setvar4,pos={150,210},size={100,20},proc=DE_SChollpanel#SVP,value=_NUM:100,title="BG End"

	SetVariable de_scholl_setvar5,pos={120,20},size={150,16},proc=DE_SChollpanel#SVP,value= _STR:"Zsnsr",title="X Suff1"
	SetVariable de_scholl_setvar6,pos={280,20},size={150,16},proc=DE_SChollpanel#SVP,value= _STR:"Zsnsr",title="X Suffix 2"
	
	SetVariable de_scholl_setvar7,pos={460,0},size={100,16},proc=DE_SChollpanel#SVP,title="Filtering"
	SetVariable de_scholl_setvar7,limits={-inf,inf,2},value= _NUM:11
	
	SetVariable de_scholl_setvar8,pos={960,0},size={125,16},proc=DE_SChollpanel#SVP,title="Ramps"
	SetVariable de_scholl_setvar8,limits={-inf,inf,2},value= _NUM:3
	
	Display/W=(10,230,420,440)/HOST=SchollPanel/N=Plot1
	Display/W=(450,20,1050,300)/HOST=SchollPanel/N=Plot2
	Display/W=(10,450,420,600)/HOST=SchollPanel/N=Plot3

	appendtograph/W=SchollPanel#Plot3 root:Schollpanel:FitResiduals vs root:Schollpanel:FitResidualsX

	appendtograph/W=SchollPanel#Plot1 root:Schollpanel:DefV_1 vs root:Schollpanel:Zsnsr_1
	appendtograph/W=SchollPanel#Plot1 root:Schollpanel:DefV_2 vs root:Schollpanel:Zsnsr_2
	appendtograph/W=SchollPanel#Plot2 root:Schollpanel:CorForce vs root:Schollpanel:CorSep
	appendtograph/W=SchollPanel#Plot2 root:Schollpanel:CorForceFil vs root:Schollpanel:CorSepFil
	appendtograph/W=SchollPanel#Plot2 root:Schollpanel:ZeroY vs root:Schollpanel:ZeroX
	ModifyGraph/W=SchollPanel#Plot2 rgb(CorForce)=(63232,45824,45824)
	ModifyGraph/W=SchollPanel#Plot2 rgb(CorForceFil)=(58368,6656,7168)
	appendtograph/W=SchollPanel#Plot2 root:Schollpanel:DispRup vs root:SchollPanel:DispXRup
	ModifyGraph/W=SchollPanel#Plot2 mode(DispRup)=3,marker(DispRup)=19,rgb(DispRup)=(58,126,184)
	
EndMacro
////---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


Static Function WindowHook(s)
STRUCT WMWinHookStruct &s
GetWindow $s.winName activeSW
String activeSubwindow = S_value
if (CmpStr(activeSubwindow,"SchollPanel#Plot1") != 0)
return 0
endif
Variable hookResult = 0
switch(s.eventCode)
case 7: // Deactivate
SetVariable de_scholl_setvar3,value= _NUM:str2num(stringbykey("POINT",CsrInfo(A),":",";"))
SetVariable de_scholl_setvar4,value= _NUM:str2num(stringbykey("POINT",CsrInfo(B),":",";"))
break
// And so on . . .
endswitch
return hookResult // 0 if nothing done, else 1
End
////---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


////ListWaves is built to find ALL the waves in a folder of the IPF with the format "ImageXXXXDefl". It then checks/makes
////a separation wave, then pops the ImageXXXX part of the string into the list ListName, which is then called by the 
////ListBox in the panel.
////

Static Function GuessRamps()
	wave/t/z LW2=root:SchollPanel:ListWave2
	variable n=0
	string CurrentName
		if(numpnts(LW2)==0)
	return 0
	endif
	make/free/n=(numpnts(LW2)) Test


	for(n=0;n<numpnts(LW2);n+=1)
		struct ForceWave Name2
		DE_Naming#WavetoStruc(LW2[n],Name2)	
		Test[n]=Name2.VNum
	endfor
	wavestats/Q Test
	return (v_max-v_min-1)

end

Static function Update(w1,w2)
wave w1, w2

	wave Ruptures
	variable m
	wave/t S=root:schollpanel:ParmName
	DE_SchollFitting#Scholl(w1,w2,str2num(S[0][1]),str2num(S[1][1]),str2num(S[2][1]),str2num(S[3][1]))

end

Static function Parms(ForceWave,SepWave,FilSepWave,decimation)
	wave ForceWave,SepWave,FilSepWave
	variable decimation

	wave w1=root:Peaks
	duplicate/o w1 root:SchollPanel:DispRup,  root:SchollPanel:Loading, root:SchollPanel:DispXRup
	wave w2=root:SchollPanel:DispRup
	wave w3=root:SchollPanel:Loading
	wave w5=root:SchollPanel:DispXRup
	wave w6=root:SchollPanel:ZeroX
	wave w7=root:SchollPanel:ZeroY

	
	variable Offset= str2num(stringbykey("NumPtsPerSec",note(forcewave),":","\r"))*.009/str2num(stringbykey("RetractVelocity",note(forcewave),":","\r"))
	w2=imag(DE_SchollFitting#LineFits(ForceWave,SepWave,FilSepWave,decimation,w1))
	w5=SepWave(real(DE_SchollFitting#LineFits(ForceWave,SepWave,FilSepWave,decimation,w1)))
	
	w3=(DE_SchollFitting#Loading(ForceWave,SepWave,FilSepWave,decimation,w1))
	//
	variable StartTime=real(DE_SchollFitting#Loading(ForceWave,SepWave,FilSepWave,decimation,w1[numpnts(w1)-1]))
	variable StartBaseline=x2pnt(ForceWave,StartTime)+Offset
	wavestats/q/r=[StartBaseline,StartBaseline+2*Offset] ForceWave
	duplicate/o/r=[StartBaseline,StartBaseline+2*Offset] SepWave w6
	duplicate/o w6,w7
	w7=v_avg
	
	
end

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Static function UpdateFilter()
	wave w1=root:Schollpanel:CorForceFil
	wave w2=root:Schollpanel:CorForce
	wave w3=root:Schollpanel:CorSep
	wave w4=root:Schollpanel:CorSepFil
	controlinfo de_scholl_setvar7
	DE_Filtering#FilterForceSep(w2,w3,w1,w4,"SVG",V_Value)
end

Static Function UpdatePRHEvents()

	wave/t/z LW2=root:SchollPanel:ListWave2
	controlinfo de_scholl_list1
	variable row1=V_Value
	controlinfo de_scholl_list2
	variable row2=V_Value
	struct ForceWave Name1
		
	DE_Naming#WavetoStruc(LW2[row2],Name1)
	controlinfo de_scholl_setvar6

	string ForceName=DE_Naming#StringCreate(Name1.Name,Name1.VNum,"Force",Name1.SDirec)
	string SepName=DE_Naming#StringCreate(Name1.Name,Name1.VNum,"Sep",Name1.SDirec)

	wave w3=$ForceName
	wave w4=$SepName
	wave DrUp=root:Schollpanel:DispRup
	wave DRupX=root:SchollPanel:DispXRup
	//
	if(WaveExists(w3))

		if(itemsinlist(stringbykey("PRHEvents",note(w3),":","\r"))!=0)
			make/free/n=(itemsinlist(stringbykey("PRHEvents",note(w3),":","\r"))) RupTest,RupXTest
			String EventList=StringbyKey("PRHEvents",note(w3),":","\r")
			Make/n=(itemsinlist(EventList))/Free Events
			Events=str2num(stringfromlist(p,EventList,";"))
			duplicate/o events DrUp,DrUpX
			DrUp=w3[events[p]]
			DrUpX=w4[events[p]]
							
		else
						
			DrUp=Nan
			DrUpX=Nan
			
		endif
						
			
	else 

	endif


end

Menu "Scholl"
	//SubMenu "Processing"
	"Open Scholl", SchollPanel()
	"Runn All", DE_SchollPanel#PopUpForMany()

	//end
	
end

Static Function PopUpForMany()
	String FitType="LinSin"
	Variable starting=0,ending=0
	Prompt FitType,"FitType",popup,"LinSin;Poly;Spline"
	Prompt starting,"Starting Wave"
		Prompt ending,"Ending Wave"

	DoPrompt "Do It!",FitType,starting,ending
	if (V_Flag)
		return 0									// user canceled
	endif
	DE_Schollpanel#SequenceManyFromOne(starttrace=starting,endtrace=ending,FitType=FitType)

end