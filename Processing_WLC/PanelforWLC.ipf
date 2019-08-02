#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include ":MultiWLCFit"
#include ":\Misc_PanelPrograms\Panel Progs"
#include ":\Misc_PanelPrograms\AsylumNaming"

Function DE_WLC_SVP_1(ctrlName,varNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable varNum
	String varStr
	String varName
	DE_PanelProgs#PrintAllFolders(root:WLCFit:FolderList,varstr)

End

Function DE_WLC_SVP_2(ctrlName,varNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable varNum
	String varStr
	String varName

End

Function DE_WLC_SVP_3(ctrlName,varNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable varNum
	String varStr
	String varName
	DE_WLC_UpdateSens()
End

Function DE_WLC_SVP_4(ctrlName,varNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable varNum
	String varStr
	String varName
	DE_WLC_UpdateSens()
End

Function DE_WLC_SVP_5(ctrlName,varNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable varNum
	String varStr
	String varName

End
Function DE_WLC_SVP_6(ctrlName,varNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable varNum
	String varStr
	String varName
	controlinfo de_wlc_setvar5
	DE_WLC_RemakeCrs(v_value)

	DE_WLC_UpdateParms()
	DE_WLC_RedoColors()
End


Function de_wlc_LBP(ctrlName,row,col,event) : ListBoxControl

	String ctrlName
	Variable row
	Variable col
	Variable event	//1=mouse down, 2=up, 3=dbl click, 4=cell select with mouse or keys
					//5=cell select with shift key, 6=begin edit, 7=end
	wave/t/z FolderList=root:WLCFit:FolderList
	string/g CurrFold
	wave LW= root:WLCFit:ListWave

	controlinfo de_wlc_setvar1 
	switch (event)
		case -1:
			break
		case 1:
			break
		case 3:
			DE_PanelProgs#ListWaves($FolderList[row],S_Value,LW)
			CurrFold=FolderList[row]
		case 4:	
			DE_PanelProgs#ListWaves($FolderList[row],S_Value,LW)
			CurrFold=FolderList[row]
	endswitch
	
	return 0
End


Function de_wlc_LBP_1(ctrlName,row,col,event) : ListBoxControl

	String ctrlName
	Variable row
	Variable col
	Variable event	//1=mouse down, 2=up, 3=dbl click, 4=cell select with mouse or keys
	//5=cell select with shift key, 6=begin edit, 7=end

	
	
	

	if(event==10)
		return 0
	endif
	
	wave/t/z LW=root:WLCFit:ListWave 
	wave/t/z FL= root:WLCFit:FolderList


	String saveDF = GetDataFolder(1)
	

	make/o/n=0  root:WLCFit:X_View, root:WLCFit:Y_View
	controlinfo de_wlc_list0
	SetDataFolder  FL[v_value]
	struct ForceWave Name1
	DE_Naming#WavetoStruc(LW[row],Name1)
	controlinfo de_WLC_setvar4
	string ZsnsrName=DE_Naming#StringCreate(Name1.Name,Name1.VNum,S_Value,Name1.SDirec)
	wave w1=$LW[row]
	nvar LIn=root:WLCFIt:INvols
	Lin=str2num(stringbykey("Invols",note(w1),":","\r"))*1e9
	nvar LSp=root:WLCFIt:Spring
	LSp=str2num(stringbykey("SpringConstant",note(w1),":","\r"))*1e3
	wave w2=$ZsnsrName
	print ZsnsrName
	SetVariable de_wlc_setvar3 value= _num:LSp
	SetVariable de_wlc_setvar2 value= _num:LIn
	duplicate/o w1 root:WLCFit:Y_View
	duplicate/o w2 root:WLCFit:X_View
	wave w3=root:WLCFit:X_View
	wave w4=root:WLCFit:Y_View
	string notes=note(w1)
	variable Xoff,yoff
	if(cmpstr(StringByKey("DE_XOff", notes ,":" ,"\r"),"")==0)
		Xoff=0
	else
		
	
		Xoff=str2num( StringByKey("DE_XOff", notes ,":" ,"\r"))
	endif

	if(cmpstr(StringByKey("DE_YOff", notes ,":" ,"\r"),"")==0)
		yoff=0
	else
		
	
		yoff=str2num( StringByKey("DE_YOff", notes ,":" ,"\r"))
	endif
	w3-=xoff
	w4-=yoff
	DoUpdate
	SetDataFolder saveDF
	return 0

End

Function de_wlc_LBP_2(ctrlName,row,col,event) : ListBoxControl
	String ctrlName
	Variable row
	Variable col
	Variable event	//1=mouse down, 2=up, 3=dbl click, 4=cell select with mouse or keys
					//5=cell select with shift key, 6=begin edit, 7=end
	DE_WLC_RedoColors()
	return 0
End

Function de_wlc_LBP_3(ctrlName,row,col,event) : ListBoxControl
	String ctrlName
	Variable row
	Variable col
	Variable event	//1=mouse down, 2=up, 3=dbl click, 4=cell select with mouse or keys
					//5=cell select with shift key, 6=begin edit, 7=end
	return 0
End

Function de_wlc_LBP_4(ctrlName,row,col,event) : ListBoxControl
	String ctrlName
	Variable row
	Variable col
	Variable event	//1=mouse down, 2=up, 3=dbl click, 4=cell select with mouse or keys
					//5=cell select with shift key, 6=begin edit, 7=end
	return 0
End



Function de_wlc_PMP_1(ctrlName,popNum,popStr) : PopupMenuControl
	String ctrlName
	Variable popNum
	String popStr
	DE_WLC_RemakeCrs(str2num(popstr))
	DE_WLC_UpdateParms()
	DE_WLC_RedoColors()
End


Function de_wlc_PMP_2(ctrlName,popNum,popStr) : PopupMenuControl
	String ctrlName
	Variable popNum
	String popStr
	wave/z/T w1=root:WLCFit:Crs
	if(cmpstr(csrinfo(a,"ImageViewer#autoplot"),"")&&cmpstr(csrinfo(b,"ImageViewer#autoplot"),""))
	variable q=str2num(popStr)
	
	w1[q-1][0]=num2str(pcsr(a,"ImageViewer#autoplot"))
	w1[q-1][1]=num2str(pcsr(b,"ImageViewer#autoplot"))

		endif
			DE_WLC_RedoColors()
End

Function de_wlc_PMP_3(ctrlName,popNum,popStr) : PopupMenuControl
	String ctrlName
	Variable popNum
	String popStr

DE_WLC_UpdateParms()

End


Function de_wlc_BP_1(ctrlName) : ButtonControl
	String ctrlName
	
	wave/T LW= root:WLCFit:ListWave
	string/g CurrFold
	variable m
	String D_Str,Z_Str

	wave Force=root:WLCFit:Y_View
	wave Sep=root:WLCFit:X_View

	wave/T Parms=root:WLCFit:ParmName
	wave/T Curs=root:WLCFit:Crs
	variable Lper
	variable Temp
//	variable n=dimsize(Curs,0)
//
	controlinfo de_wlc_popup2
	variable LStep
	switch(V_Value)
	
		case 1:
			Lper=str2num(Parms[0][1])*1e-9
			Temp=str2num(Parms[3][1])
			variable L0start=str2num(Parms[1][1])*1e-9
			LStep=str2num(Parms[2][1])*1e-9
			DE_MultiWLC#DeltaLC(Force,Sep,L0start,LStep,Lper,Temp,Curs,Folder=root:WLCFit:)
			break
	
		case 2:
			Lper=str2num(Parms[0][1])*1e-9
			Temp=str2num(Parms[1][1])
			DE_MultiWLC#ArbLcs(Force,Sep,Lper,Temp,Curs,Parms,Folder=root:WLCFit:)
			break
		case 3:
			Lper=str2num(Parms[0][1])*1e-9
			Temp=str2num(Parms[3][1])
			variable LFinal=str2num(Parms[1][1])*1e-9
			LStep=str2num(Parms[2][1])*1e-9
			DE_MultiWLC#FinalLC(Force,Sep,LFinal,LStep,Lper,Temp,Curs,Folder=root:WLCFit:)
			break
	endswitch
	wave w3=root:WLCFit:ComSep
	wave w4=root:WLCFit:ResForce

	if(strsearch(TraceNameList("ImageViewer#Autoplot","",1),"ResForce",0)==-1)
		wave w3=root:WLCFit:ComSep
		wave w4=root:WLCFit:ResForce
		
		appendtograph/W=ImageViewer#Autoplot w4 vs w3[][0]
		ModifyGraph/W=ImageViewer#Autoplot rgb(ResForce)=(0,0,65280)
		ModifyGraph/W=ImageViewer#Autoplot lsize(ResForce)=2

	else
	
	endif
	
End

Function de_wlc_BP_2(ctrlName) : ButtonControl
	String ctrlName
	
	wave/T LW= root:WLCFit:ListWave
	string/g CurrFold
	variable m
	String D_Str,Z_Str

	wave Force=root:WLCFit:Y_View
	wave Sep=root:WLCFit:X_View

	wave/T Parms=root:WLCFit:ParmName
	wave/T Curs=root:WLCFit:Crs
	variable Lper
	variable Temp
//	variable n=dimsize(Curs,0)
//
	controlinfo de_wlc_popup2
	variable LStep
	switch(V_Value)
	
		case 1:
			Lper=str2num(Parms[0][1])*1e-9
			Temp=str2num(Parms[3][1])
			variable L0start=str2num(Parms[1][1])*1e-9
			LStep=str2num(Parms[2][1])*1e-9
			DE_MultiWLC#DeltaLC(Force,Sep,L0start,LStep,Lper,Temp,Curs,Folder=root:WLCFit:)
			break
	
		case 2:
			Lper=str2num(Parms[0][1])*1e-9
			Temp=str2num(Parms[1][1])
			DE_MultiWLC#ArbLcs(Force,Sep,Lper,Temp,Curs,Parms,Folder=root:WLCFit:)
			break
		case 3:
			Lper=str2num(Parms[0][1])*1e-9
			Temp=str2num(Parms[3][1])
			variable LFinal=str2num(Parms[1][1])*1e-9
			LStep=str2num(Parms[2][1])*1e-9
			DE_MultiWLC#FinalLC(Force,Sep,LFinal,LStep,Lper,Temp,Curs,Folder=root:WLCFit:)
			break
	endswitch

	wave w3=root:WLCFit:ComSep
	wave w5=root:WLCFit:ComForce

	wave w4=root:WLCFit:ResForce
	
	wave/t w6=root:WLCFit:FitResName
	wave w7=root:WLCFit:FitResSel


	if(strsearch(TraceNameList("ImageViewer#Autoplot","",1),"ResForce",0)==-1)

		appendtograph/W=ImageViewer#Autoplot w4 vs w3 
		ModifyGraph rgb(ResForce)=(0,0,65280)
		ModifyGraph lsize(ResForce)=2


	else
	
	endif
	string FitHold=DE_WLC_MakeCoefsAndHolds()

	wave coefs=root:WLCFIt:W_coef


	controlinfo de_wlc_popup2
	switch(V_Value)
	
		case 1:
			FuncFit/Q/H=FitHold/NTHR=0 DE_Fit_MWC_DLC coefs  w5 /X=w3 /D =w4
			//redimension/n=
			w6[0][1]=num2str(coefs[0])
			note/k w4
			note w4 "Persistance Length="+num2str(coefs[0])
			w6[1][1]=num2str(coefs[1]*1e9)
			note w4 "Starting Contour Length="+num2str(coefs[1])
			w6[2][1]=num2str(coefs[2]*1e9)
			note w4 "Delta LC="+num2str(coefs[2])
			w6[3][1]=num2str(coefs[3])
			note w4 "Temp="+num2str(coefs[3])

			break
		
		case 2:	
			FuncFit/Q/H=FitHold/NTHR=0 DE_Fit_MWC_LCS coefs  w5 /X=w3 /D=w4
			variable q=numpnts(coefs)
			variable b
			note/k w4
			for(b=0;b<q;b+=1)
			w6[b][1]=num2str(coefs[b])
			note w4, w6[b][0]+" "+w6[b][1]
		
			endfor
			break
		case 3:
		print fithold
			FuncFit/Q/H=FitHold/NTHR=0 DE_Fit_MWC_LCF coefs  w5 /X=w3 /D =w4
			note/k w4
			//redimension/n=
			w6[0][1]=num2str(coefs[0])
			note w4 "Persistance Length="+num2str(coefs[0])

			w6[1][1]=num2str(coefs[1]*1e9)
			note w4 "Starting Contour Length="+num2str(coefs[1])
			w6[2][1]=num2str(coefs[2]*1e9)
			note w4 "Delta LC="+num2str(coefs[2])
			w6[3][1]=num2str(coefs[3])
			note w4 "Temp="+num2str(coefs[3])
			break
	endswitch
//	

end
	
	
Function de_wlc_BP_3(ctrlName) : ButtonControl
	String ctrlName
	
	wave/T w1=root:WLCFit:FitResName
	wave/T w2=root:WLCFit:ParmName
	
	w2[][1]=w1[p][1]

end	


Function de_wlc_BP_4(ctrlName) : ButtonControl
	String ctrlName
	
	wave/t/z LW=root:WLCFit:ListWave 
	wave/t/z FL= root:WLCFit:FolderList

	controlinfo de_wlc_list1
	variable row=V_value
	String saveDF = GetDataFolder(1)
	

	controlinfo de_wlc_list0
	SetDataFolder  FL[v_value]
	wave w1=$LW[row]
	//struct ForceWave Name1
	DE_WLC#CalculateDeltaLCs(w1)
	//DE_Naming#WavetoStruc(LW[row],Name1)
	SetDataFolder saveDF

end	


Window panel0() : Panel
	PauseUpdate; Silent 1		// building window...

	NewPanel/N=ImageViewer /W=(10,10,1500,700)
	newdatafolder/o root:WLCFIT
	make/o/t/n=(0,2) root:WLCFit:Crs
	make/o/n=(0,2) root:WLCFit:CrsSel
	make/o/t/n=(4,3) root:WLCFit:ParmName
	make/o/n=(4,3) root:WLCFit:ParmSel
	make/o/t/n=(4,2) root:WLCFit:FitResName
	make/o/n=(4,2) root:WLCFit:FitResSel
	root:WLCFit:ParmName[0][0]={"Persistance Length (nm)","Starting Contour Length (nm)","Contour Length Step (nm)","Temperature(K)"}
	root:WLCFit:ParmName[0][1]={".4","15","18","300"}
	root:WLCFit:ParmName[0][2]={"0","0","0","1"}
	root:WLCFit:ParmSel[][1]=3
	root:WLCFit:FitResName[0][0]={"Persistance Length (nm)","Starting Contour Length (nm)","Contour Length Step (nm)","Temperature(K)"}
	root:WLCFit:FitResName[0][1]={"XX","XX","XX","XXX"}
		
	Button de_wlc_button0,pos={1120,20},size={125,20},proc=de_wlc_BP_1,title="Plot Parms"
	Button de_wlc_button1,pos={1120,50},size={125,20},proc=de_wlc_BP_2,title="Fit!"
	Button de_wlc_button2,pos={1200,350},size={125,20},proc=de_wlc_BP_3,title="Copy Parameter!"
	Button de_wlc_button4,pos={1200,390},size={125,20},proc=de_wlc_BP_4,title="PrintShit!"

	make/o/t/n=0 root:WLCFit:FolderList,root:WLCFit:Listwave
	SetVariable de_wlc_setvar0,pos={10,0},size={150,16},proc=DE_WLC_SVP_1,value= _STR:"Image*",title="Search String"
	SetVariable de_wlc_setvar1,pos={160,0},size={150,16},proc=DE_WLC_SVP_2,value= _STR:"Image*Force*",title="Search String"
	SetVariable de_wlc_setvar4,pos={350,0},size={150,16},proc=DE_WLC_SVP_5,value= _STR:"Sep",title="Search String 2"


	DE_PanelProgs#PrintAllFolders(root:WLCFit:FolderList,"Image*")
	root:WLCFit:ParmSel[][2]=3
	
	ListBox de_wlc_list0,pos={25,35},size={150,150},proc=de_wlc_LBP,listWave=root:WLCFit:FolderList
	ListBox de_wlc_list0,row= 0,mode= 1,selRow= 0
	
	ListBox de_wlc_list1,pos={250,35},size={150,150},proc=de_wlc_LBP_1,listWave=root:WLCFit:ListWave
	ListBox de_wlc_list1,row= 0,mode= 1,selRow= 0
	
	ListBox de_wlc_list2,pos={550,2},size={150,150},proc=de_wlc_LBP_2
	ListBox de_wlc_list2,listWave=root:WLCFit:Crs,selWave=root:WLCFit:CrsSel,editStyle= 2,userColumnResize= 1
	
	ListBox de_wlc_list3,pos={750,2},size={350,150},proc=de_wlc_LBP_3
	ListBox de_wlc_list3,listWave=root:WLCFit:ParmName,selWave=root:WLCFit:ParmSel,editStyle= 2,userColumnResize= 1,widths={35,20,10}
	
	ListBox de_wlc_list4,pos={1110,200},size={350,150},proc=de_wlc_LBP_4
	ListBox de_wlc_list4,listWave=root:WLCFit:FitResName,selWave=root:WLCFit:FitResSel,editStyle= 2,userColumnResize= 1,widths={150,150}

//	PopupMenu de_wlc_popup0,pos={570,160},size={218,21},proc=de_wlc_PMP_1,title="Number of Fits"
//	PopupMenu de_wlc_popup0,mode=1,popvalue="---",value= #"\"1;2;3;4;5\""
	SetVariable de_wlc_setvar5,pos={570,160},size={218,21},title="Num Curves",limits={1,inf,1},proc=DE_WLC_SVP_6,value= _num:5

	PopupMenu de_wlc_popup1,pos={1150,500},size={118,21},proc=de_wlc_PMP_2,title="Assign To"
	PopupMenu de_wlc_popup1,mode=1,popvalue="---",value=DE_WLC_CrsNum()
	
	PopupMenu de_wlc_popup2,pos={1150,450},size={118,21},proc=de_wlc_PMP_3,title="Fit Type"
	PopupMenu de_wlc_popup2,mode=1,popvalue="---",value=#"\"Fixed P, Delta Lc;Fixed P Free Lc;Fixed Lp, Fixed Final LC\""
	
	variable/g root:WLCFit:Invols=Nan
	TitleBox de_Wlc_title0 title="Invols (nm/V)",pos={435,25}
	ValDisplay de_Wlc_valdisp0 pos={409,50},size={35,20},value=root:WLCFit:Invols
	SetVariable de_wlc_setvar2,pos={470,50},size={50,20},proc=DE_WLC_SVP_3,value= _num:Nan

	variable/g root:WLCFit:Spring=Nan
	TitleBox de_Wlc_title1 title="Spring Constant (pN/nm)",pos={410,95}
	ValDisplay de_Wlc_valdisp1 pos={409,120},size={35,20},value=root:WLCFit:Spring
	
	SetVariable de_wlc_setvar3,pos={470,120},size={50,20},proc=DE_WLC_SVP_4,value= _num:Nan

	
	make/o/n=(0) root:WLCFit:X_View,root:WLCFit:Y_View,root:WLCFit:Color,root:WLCFit:ColorColor
	Display/W=(20,200,1100,600)/HOST=#/N=AutoPlot
	appendtograph root:WLCFit:Y_View vs root:WLCFit:X_View
	appendtograph root:WLCFit:Color vs root:WLCFit:X_View
	//ModifyGraph lsize(Color)=2,zmrkSize(Color)={Defl_nm_Hold,*,*,1,10};DelayUpdate
	ModifyGraph lsize(Color)=2;DelayUpdate

	ModifyGraph zColor(Color)={:WLCFIT:ColorColor,0,6,Rainbow,0}
	SetWindow ImageViewer, hook(MyHook)=DE_WLC_Plot_Hook // Install window hook
	ShowInfo/W=ImageViewer
EndMacro






Function DE_WLC_Plot_Hook(s)
	STRUCT WMWinHookStruct &s
	Variable hookResult = 0 // 0 if we do not handle event, 1 if we handle it.
	
	GetWindow $s.winName activeSW
	String activeSubwindow = S_value
	if (CmpStr(activeSubwindow,"ImageViewer#AutoPlot") != 0)
		return 0
	endif
	string items,tracename
	variable xoff,yoff,tracepnt
	switch(s.eventcode)
		case 3: // Keyboard event
			if(s.eventmod==9)
				if(cmpstr(TraceFromPixel(s.mouseloc.h, s.mouseloc.v, ""),"")==0)
					items= "X Offset; Y Offset; X-Y Offset" // 2nd is divider, 3rd is checked
					//variable Xoff=xwave[TracePnt]
					//variable Yoff=ywave[tracepnt]
					TraceName=stringbykey("Trace",TraceFromPixel(s.mouseloc.h, s.mouseloc.v, ""),":",";")
					TracePnt=str2num(stringbykey("HITPOINT",TraceFromPixel(s.mouseloc.h, s.mouseloc.v, ""),":",";"))
					wave Xwave=XWaveRefFromTrace("",TraceName)
					wave Ywave=TraceNameToWaveRef("ImageViewer#AutoPlot",TraceName)
					Xoff=AxisValFromPixel("ImageViewer#AutoPlot","bottom",s.mouseloc.h)
					Yoff=AxisValFromPixel("ImageViewer#AutoPlot","left",s.mouseloc.v)
					PopupContextualMenu items
					switch( V_Flag )
		
						case 1:
							XWave-=Xoff
							DE_WLC_UpdateOffsets(XOff,0)
							break	
						case 2:
							YWave-=Yoff
														DE_WLC_UpdateOffsets(0,Yoff)

							break
						case 3:
							XWave-=Xoff
							YWave-=Yoff

														DE_WLC_UpdateOffsets(XOff,Yoff)

							break
						break
					endswitch		
				else
					
					//print stringbykey("Trace",TraceFromPixel(s.mouseloc.h, s.mouseloc.v, ""),":",";")
					//
					TraceName=stringbykey("Trace",TraceFromPixel(s.mouseloc.h, s.mouseloc.v, ""),":",";")
					TracePnt=str2num(stringbykey("HITPOINT",TraceFromPixel(s.mouseloc.h, s.mouseloc.v, ""),":",";"))
					wave Xwave=XWaveRefFromTrace("",TraceName)
					wave Ywave=TraceNameToWaveRef("ImageViewer#AutoPlot",TraceName)
					items= "A;B;\M1-;X Offset; Y Offset; X-Y Offset" // 2nd is divider, 3rd is checked
					//variable Xoff=xwave[TracePnt]
					//variable Yoff=ywave[tracepnt]
					Xoff=AxisValFromPixel("ImageViewer#AutoPlot","bottom",s.mouseloc.h)
					Yoff=AxisValFromPixel("ImageViewer#AutoPlot","left",s.mouseloc.v)
					PopupContextualMenu items
					switch( V_Flag )
						case 1:
							Execute "Cursor/P/W=ImageViewer#AutoPlot A "+TraceName+" "+num2str( TracePnt)
							break;
						case 2:
							Execute "Cursor/P/W=ImageViewer#AutoPlot B "+TraceName+" "+num2str( TracePnt)
							break
						case 4:
							XWave-=Xoff
														DE_WLC_UpdateOffsets(XOff,0)

							break	
						case 5:
							YWave-=Yoff
														DE_WLC_UpdateOffsets(0,Yoff)

							break
						case 6:
							XWave-=Xoff
							YWave-=Yoff
														DE_WLC_UpdateOffsets(XOff,Yoff)


							break
							break
					endswitch				
				endif
				//String checked= "\\M0:!" + num2char(18) + ":" // checkmark code
				
			endif
			break
	endswitch
	return hookResult // If non-zero, we handled event and Igor will ignore it.
End




function AutoPlotForce(Folder,Number)
DFREF  Folder
variable number
variable Invols, Spring,Zsens
string D_Str,Z_Str
//
sprintf D_Str, "Image%04gForce_Ret",number	
sprintf Z_Str, "Image%04gSep_Ret",number

wave w3=Folder:$D_Str
wave w4=Folder:$Z_Str

appendtograph/W=ImageViewer#Autoplot w3 vs w4
wavestats/q w3
SetAxis/W=ImageViewer#Autoplot left v_min,v_max

end


function DE_WLC_RemakeCrs(num)
	variable num
	newdatafolder/o root:WLCFit
	if(!WaveExists(root:WLCFit:Crs))
		make/t/n=(num,2) root:WLCFit:Crs
		wave/t w1=root:WLCFit:Crs
		w1="0"
	else
	
		wave/t w1=root:WLCFit:Crs
		variable n=dimsize(root:WLCFit:Crs,0)
		n=min(num,n)
		duplicate/t/o root:WLCFit:Crs root:WLCFit:Crs_Sav
		wave/t w2=root:WLCFit:Crs_Sav

		make/t/o/n=(num,2) root:WLCFit:Crs
		w1="0"

		variable i
		for(i=0; i<n;i+=1)
			w1[i][0]=w2[i][0]
			w1[i][1]=w2[i][1]

		endfor
		killwaves root:WLCFit:Crs_Sav
	endif
	make/o/n=(num,2) 	root:WLCFit:CrsSel
	wave w3=root:WLCFit:CrsSel
	w3=2

end

function/S  DE_WLC_CrsNum()
	controlinfo/w=ImageViewer de_wlc_setvar5
	variable n
	string Res=""
	for(n=1;n<=v_value;n+=1)
	Res+=num2str(n)+";"
	endfor
return Res
end




function  DE_WLC_RedoColors()
	wave DWave=TraceNameToWaveRef("ImageViewer#AutoPlot",StringFromList(0,ListMatch(TraceNameList("ImageViewer#Autoplot",";",1),"Y_*")))
	duplicate/o DWave root:WLCFit:Color, root:WLCFit:ColorColor
	wave w1= root:WLCFit:Color
	wave w2= root:WLCFit:ColorColor
	wave/t w3=root:WLCFit:Crs
	w1=NaN
	w2=NaN
	variable n=dimsize(w3,0)
	variable i=0
	for(i=0;i<n;i+=1)
		//w1[str2num(w3[i][0]),str2num(w3[i][1])]=DWave[p]
		if(str2num(w3[i][0])==str2num(w3[i][1]))
		else
			w1[str2num(w3[i][0]),str2num(w3[i][1])]=10e-12
			w2[str2num(w3[i][0]),str2num(w3[i][1])]=i+1
		endif
	endfor

end

function  DE_WLC_UpdateParms()

	wave/z/T w1=root:WLCFit:ParmName
	wave/z w2=root:WLCFit:ParmSel
	wave/z/t w3=root:WLCFit:FitResName
	wave/z w4=root:WLCFit:FitResSel

	controlinfo de_wlc_popup2

	switch (V_Value)
		case 1:
	
			redimension/n=(4,3) w1,w3
			redimension/n=(4,3) w2
			redimension/E=1/n=(4,2) w3
			redimension/n=(4,2) w4
			w2=0
			w2[][0]=0
			w2[][1]=3

			w2[][2]=3
			w4=0
			w1[0][0]={"Persistance Length (nm)","Starting Contour Length (nm)","Contour Length Step (nm)","Temperature(K)"}
			w1[0][1]={".4","15","18","300"}
			w1[0][2]={"0","0","0","1"}
			w3[0][0]={"Persistance Length (nm)","Starting Contour Length (nm)","Contour Length Step (nm)","Temperature(K)"}
			w3[0][1]={"XX","XX","XX","XXX"}
			break
	
		case 2 :
			controlinfo DE_WLC_SetVar5
			redimension/E=1/n=(V_Value+2,3) w1
			redimension/n=(V_Value+2,3) w2
			redimension/E=1/n=(V_Value+2,2) w3
			redimension/n=(V_Value+2,2) w4
			w2[][0]=0
			w2[][1]=3
			w2[][2]=3
			w4=0
			//controlinfo de_wlc_popup0
			controlinfo DE_WLC_SetVar5
			w1[0][0]={"Persistance Length (nm)","Temperature(K)"}
			w1[0][1]={".4","300"}
			w1[0][2]={"0","1"}
			w3[0][0]={"Persistance Length (nm)","Temperature(K)"}
			w3[0][1]={"XX","XX"}
			variable n
			for(n=2;n<=v_value+1;n+=1)
				w1[n][0]={"LC"+num2str(n-1)}
				w1[n][1]={num2str(50+20*(n-2))}
				w1[n][2]={"0"}
				w3[n][0]={"LC"+num2str(n-1)}
				w3[n][1]={"XX"}
			
			endfor
			break
			
	case 3:
	
			redimension/n=(4,3) w1,w3
			redimension/n=(4,3) w2
			redimension/E=1/n=(4,2) w3
			redimension/n=(4,2) w4
			w2=0
			w2[][0]=0
			w2[][1]=3

			w2[][2]=3
			w4=0
			w1[0][0]={"Persistance Length (nm)","Final Contour Length (nm)","Contour Length Step (nm)","Temperature(K)"}
			w1[0][1]={".4","150","18","300"}
			w1[0][2]={"0","0","0","1"}
			w3[0][0]={"Persistance Length (nm)","Final Contour Length (nm)","Contour Length Step (nm)","Temperature(K)"}
			w3[0][1]={"XX","XXX","XX","XXX"}
			break
	
	
	endswitch

end

function/S 	 DE_WLC_MakeCoefsAndHolds()
	wave/T Parms=root:WLCFit:ParmName
	wave/T Curs=root:WLCFit:Crs
	variable n=dimsize(Parms,0)
	string FitHold=""

	controlinfo de_wlc_popup2
	switch(V_Value)
		case 1:
			Make/o/D/N=4 root:WLCFIt:W_coef
			wave coefs=root:WLCFIt:W_coef
			if(dimsize(curs,0)==1)
				FitHold=Parms[0][2]+Parms[1][2]+"1"+Parms[3][2]
			else
				FitHold=Parms[0][2]+Parms[1][2]+Parms[2][2]+Parms[3][2]

			endif
	
			coefs= {str2num(Parms[0][1])*1e-9,str2num(Parms[1][1])*1e-9,str2num(Parms[2][1])*1e-9,str2num(Parms[3][1])}
			break
		case 2:
			Make/O/D/N=(n) root:WLCFIt:W_coef
			wave coefs=root:WLCFIt:W_coef
			variable i
			for(i=0;i<n;i+=1)
				FitHold+=Parms[i][2]
				if(i==1)
			
					coefs[i]=str2num(Parms[i][1])
				else
					coefs[i]=str2num(Parms[i][1])*1e-9
				endif


			endfor
		
			break
		case 3:
			Make/o/D/N=4 root:WLCFIt:W_coef
			wave coefs=root:WLCFIt:W_coef
			if(dimsize(curs,0)==1)
				FitHold=Parms[0][2]+Parms[1][2]+"1"+Parms[3][2]
			else
				FitHold=Parms[0][2]+Parms[1][2]+Parms[2][2]+Parms[3][2]

			endif
	
			coefs= {str2num(Parms[0][1])*1e-9,str2num(Parms[1][1])*1e-9,str2num(Parms[2][1])*1e-9,str2num(Parms[3][1])}
			break
	endswitch
	return FitHold
end

function DE_CorrParms(YIn,XIn,Invols,Spring,ForceOut,SepOut)
	wave Yin, Xin, ForceOut,SepOut
	variable Invols, Spring

	variable origInvols,origspring
	string FindS
	origInvols=str2num(stringbykey("Invols",note(Yin),":","\r"))
	origSpring=str2num(stringbykey("SpringConstant",note(Yin),":","\r"))
	
	if(Invols==0)
		Invols=origInvols
	endif
	if(Spring==0)
		Spring=origspring
	endif
	
	variable m
	sscanf nameofwave(Yin), "Image%04g*",m
	FindS=nameofwave(Yin)[9,13]
	struct ForceWave Name1
	DE_Naming#WavetoStruc(nameofwave(Yin),Name1)
	strswitch(Name1.SMeas)
		case "DeflV":
			duplicate/o Yin Hold_Def, Hold_For,ODefl
			FastOp ODefl=(origInvols)*Yin
			FastOp Hold_Def=(Invols)*Yin
			FastOp Hold_For=(Invols*Spring)*Yin

			break

		case "Defl":

			duplicate/o Yin Hold_Def, Hold_For,ODefl
			FastOp Hold_Def=(Invols/origInvols)*Yin
			FastOp Hold_For=(Invols/origInvols*Spring)*Yin

			break

		case "DeflCor":

			duplicate/o Yin Hold_Def, Hold_For,ODefl
			FastOp Hold_Def=(Invols/origInvols)*Yin
			FastOp Hold_For=(Invols/origInvols*Spring)*Yin

			break

		case "Force":
			duplicate/o Yin Hold_Def, Hold_For,ODefl
			FastOp ODefl=(1/origSpring)*ODefl

			FastOp Hold_Def=(Invols/origInvols/origSpring)*Yin
			FastOp Hold_For=(Invols/origInvols/origSpring*Spring)*Yin

			break
		Default:
			Print "Don't recognize Y Wave Type"
			
			break


	endswitch


	struct ForceWave Name2
	DE_Naming#WavetoStruc(nameofwave(Xin),Name2)
	
	strswitch(Name2.SMeas)
		case "Sep":
			duplicate/o Xin Hold_Zs,Hold_Sep
			FastOp Hold_Zs=Xin-ODefl
			FastOp Hold_Sep=Hold_Zs+Hold_Def
			
			break

		case "ZSns":

			duplicate/o Xin Hold_Zs,Hold_Sep
			FastOp Hold_Sep=Hold_Zs+Hold_Def

			break
		Default:
			Print "Don't recognize X Wave Type"
			
			break



	endswitch

	duplicate/o Hold_For ForceOut
	duplicate/o Hold_Sep SepOut
	killwaves/z Hold_Zs,Hold_Sep, Hold_Def, Hold_For
end

function DE_WLC_UpdateSens()
	wave/t/z LW=root:WLCFit:ListWave 
	controlinfo de_wlc_list1
	variable row=V_value
	controlinfo de_WLC_setvar2
	variable IV=v_value*1e-9
	controlinfo de_WLC_setvar3
	variable SC=v_value*1e-3
	variable m
	
	wave w1=$LW[row]
	
	struct ForceWave Name1
	DE_Naming#WavetoStruc(LW[row],Name1)
	controlinfo de_wlc_setvar4
	string ZsnsrName=DE_Naming#StringCreate(Name1.Name,Name1.VNum,S_Value,Name1.SDirec)
	wave w1=$LW[row]
	wave w2=$ZsnsrName
	wave w3=root:WLCFit:X_View
	wave w4=root:WLCFit:Y_View
	DE_CorrParms(w1,w2,IV,SC,w4,w3)


	string notes=note(w1)
	variable Xoff,yoff
	if(cmpstr(StringByKey("DE_XOff", notes ,":" ,"\r"),"")==0)
	Xoff=0
	else
		
	
	Xoff=str2num( StringByKey("DE_XOff", notes ,":" ,"\r"))
endif

	if(cmpstr(StringByKey("DE_YOff", notes ,":" ,"\r"),"")==0)
	yoff=0
	else
		
	
	yoff=str2num( StringByKey("DE_YOff", notes ,":" ,"\r"))
endif
	w3-=xoff
	w4-=yoff
end


function DE_WLC_UpdateOffsets(xo,yo)
variable xo,yo


	
	controlinfo de_wlc_list1
	variable row=v_value

	wave/t/z LW=root:WLCFit:ListWave 

	struct ForceWave Name1
	DE_Naming#WavetoStruc(LW[row],Name1)
	controlinfo de_WLC_setvar4

	string ZsnsrName=DE_Naming#StringCreate(Name1.Name,Name1.VNum,S_Value,Name1.SDirec)

	wave w1=$LW[row]
	wave w2=$ZsnsrName
variable xoffCurr,yoffCurr
	string notes=note(w1)
	if(cmpstr(StringByKey("DE_XOff", notes ,":" ,"\r"),"")==0)
	xoffCurr=0
	else
		
	
	xoffCurr=str2num( StringByKey("DE_XOff", notes ,":" ,"\r"))
endif

	if(cmpstr(StringByKey("DE_YOff", notes ,":" ,"\r"),"")==0)
	yoffCurr=0
	else
		
	
	yoffCurr=str2num( StringByKey("DE_YOff", notes ,":" ,"\r"))
endif

variable  xoffnew, yoffnew

xoffnew=xo+xoffCurr
yoffnew=yo+yoffCurr

DE_PanelProgs#ReplaceNote(w1,num2str(xoffnew),"DE_Xoff")
DE_PanelProgs#ReplaceNote(w1,num2str(yoffnew),"DE_yoff")

DE_PanelProgs#ReplaceNote(w2,num2str(xoffnew),"DE_Xoff")
DE_PanelProgs#ReplaceNote(w2,num2str(yoffnew),"DE_yoff")

end
