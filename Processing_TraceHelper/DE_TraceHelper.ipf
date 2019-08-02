#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma modulename=DE_TraceHelper
#include <Readback ModifyStr>
#include ":\Misc_PanelPrograms\Panel Progs"
#include ":\Misc_PanelPrograms\AsylumNaming"
Static Function TraceSelector(ctrlName,row,col,event) : ListBoxControl

	String ctrlName
	Variable row
	Variable col
	Variable event	//1=mouse down, 2=up, 3=dbl click, 4=cell select with mouse or keys
					//5=cell select with shift key, 6=begin edit, 7=end
		wave SWave= root:SelWave
	wave/T TNL=root:TraceNameLists
		variable n
		Svar GraphName=root:STrGN
		String list
		wfprintf list, "%s;", TNL
		string X1= ListMatch(list,"*WLC*",";") 	
		string WLCTrace=stringFromList(0,X1)
		ModifyGraph/W=$GraphName hideTrace=1
				ModifyGraph/W=$GraphName hideTrace($WLCTrace)=0

		for(n=0;n<numpnts(SWave);n+=1)
			if(SWave[n]==1)
				ModifyGraph/W=$GraphName hideTrace($TNL[n])=0
			endif
		endfor
end
Static Function SelectTraceWindow(graphName)
	String graphName
	string/G root:STrGN=graphname
	DoWindow/F $graphName // Bring graph to front
	if (V_Flag == 0) // Verify that graph exists
		Abort "UserCursorAdjust: No such graph."
		return -1
	endif
	NewPanel/N=Devin /W=(187,368,437,831)
	string ListWavesonPlot=traceNameList(graphName,";",1)
	variable IgorVer=(NumberByKey("IGORVERS", IgorInfo(0)))
	
	make/T/o/n=0 root:TraceNameLists
	wave/T TNL=root:TraceNameLists
	DE_ListToTextWave(ListWavesonPlot, TNL,";")

	wave/T TNL=root:TraceNameLists
	make/o/n=(numpnts(TNL)) root:SelWave
	wave Swave=root:SelWave
	Swave=1
	ListBox de_trace win=Devin, pos={000,20},size={150,300},proc=DE_TraceHelper#TraceSelector,listWave=TNL
	ListBox de_trace win=Devin,selWave=Swave,row= 0,selRow= 0,mode= 4
//	DoWindow/C tmp_Select // Set to an unlikely name
//AutoPositionWindow/E/M=1/R=$graphName

//
//	NewPanel /K=2 /W=(187,368,437,531) as "Pause for Cursor"
//	DoWindow/C tmp_PauseforCursor // Set to an unlikely name
//	AutoPositionWindow/E/M=1/R=$graphName // Put panel near the graph
//	DrawText 21,20,"Adjust the cursors and then"
//	DrawText 21,40,"Click Continue."
//	Button button0,pos={80,58},size={92,20},title="Continue"
//	Button button0,proc=DE_Schollpanel#BatchSelect_ContButtonProc
//	Variable didAbort= 0
//	if( autoAbortSecs == 0 )
//		PauseForUser tmp_PauseforCursor,$graphName
//	else
//		SetDrawEnv textyjust= 1
//		DrawText 162,103,"sec"
//		SetVariable sv0,pos={48,97},size={107,15},title="Aborting in "
//		SetVariable sv0,limits={-inf,inf,0},value= _NUM:10
//		Variable td= 10,newTd
//		Variable t0= ticks
//		Do
//			newTd= autoAbortSecs - round((ticks-t0)/60)
//			if( td != newTd )
//				td= newTd
//				SetVariable sv0,value= _NUM:newTd,win=tmp_PauseforCursor
//				if( td <= 10 )
//					SetVariable sv0,valueColor= (65535,0,0),win=tmp_PauseforCursor
//				endif
//			endif
//			if( td <= 0 )
//				DoWindow/K tmp_PauseforCursor
//				didAbort= 1
//				break
//			endif
//			PauseForUser/C tmp_PauseforCursor,$graphName
//		while(V_flag)
//	endif
//	return didAbort
End
Static Function PlotABunch(firstnum,lastnum,Basename,DirectionList)
	variable firstnum,lastnum
	String Basename,DirectionList
	variable n,m
	display
	for(n=firstnum;n<lastnum;n+=1)
	for(m=0;m<itemsinlist(DirectionList,";");m+=1)
		wave w1=$DE_Naming#StringCreate(Basename,n,"Force",stringfromlist(m,DirectionList,";"))
		wave w2=$DE_Naming#StringCreate(Basename,n,"Sep",stringfromlist(m,DirectionList,";"))
		if(waveexists(w1)&&waveexists(w1))
			appendtograph w1 vs w2
		endif
		//

	endfor
	endfor


end
Function FilteringAllonPlot()

	variable n=0
	string list=tracenamelist("",";",1)
	string CurrentForce,CurrentSep,Info
	variable xoff,yoff,xmult,ymult
	for(n=0;n< itemsinlist(list);n+=1)
//	for(n=0;n<1;n+=1)
		CurrentForce=stringfromlist(n,list)
		CurrentSep=Replacestring("Force",CurrentForce,"Sep")
		wave w1=$CurrentForce
		Info=traceinfo("",CurrentForce,0)
		// stringfromlist(0,stringbykey("offset(X)",Info,"=",";"),",")
		xoff= GetNumFromModifyStr(info,"offset","{",0)
		yoff= GetNumFromModifyStr(info,"offset","{",1)
		xmult= GetNumFromModifyStr(info,"muloffset","{",0)
		ymult= GetNumFromModifyStr(info,"muloffset","{",1)

		if(waveexists($CurrentSep)==1)
		print currentsep
			wave w2=$CurrentSep

			make/o/n=0 $(CurrentForce+"_Sm"),$(CurrentSep+"_Sm")
			wave w3=$(CurrentForce+"_Sm")
			wave w4=$(CurrentSep+"_Sm")
			DE_Filtering#FilterForceSep(w1,w2,w3,w4,"SVG",25)
			appendtograph w3 vs w4
			ModifyGraph offset($nameofwave(w3))={xoff,yoff}
						ModifyGraph muloffset($nameofwave(w3))={xmult,ymult}

		else
		endif
	endfor


end
Static Function BatchSelect_ContButtonProc(ctrlName) : ButtonControl
	String ctrlName
	print "RUNNING"
	DoWindow/K tmp_PauseforCursor // Kill panel
	DoWindow/K tmp_Select // Kill panel

End

Static Function DE_ListToTextWave(ListWavesonPlot, Waveout,Divider)

	string ListWavesonPlot,Divider
	wave Waveout
	
	make/T/free/n=(itemsinlist(ListWavesonPlot,Divider)) Test
	
	variable n
	for(n=0;n<itemsinlist(ListWavesonPlot,Divider);n+=1)
		Test[n]=StringFromList(n, ListWavesonPlot,Divider)
		
	endfor
	duplicate/o Test Waveout
end
Menu "TraceHelper"
	//SubMenu "Processing"
	"Open TraceHelp", DE_TraceHelper#SelectTraceWindow(WinName(0,3))
	"Filter All on Plot", DE_TraceHelper#FilteringAllonPlot()

	//end
	
end