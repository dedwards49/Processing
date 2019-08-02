#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include "C:\HanLabDNPBackUP\Users\Devin\Documents\Software\Boulder\Devin'sIgor\SchollFitting\Scholl_panel"

#include "C:\HanLabDNPBackUP\Users\Devin\Documents\Software\Boulder\Devin'sIgor\PanelProgs\Panel Progs"
//
//Function DE_Rup_CBC(ctrlName,checked) : CheckBoxControl
//	String ctrlName
//	Variable checked
//
//	
//	DoWindow FastWindow
//	if(V_flag==0)
//		return 0
//	endif
//	
//	if(checked==1)
//		string saveDF = GetDataFolder(1)
//		controlinfo de_rup_popup1
//		SetDataFolder s_value
//
//		controlinfo de_rup_popup0
//		wave w1=$S_value
//		variable off=str2num(stringbykey("Invols",note(w1),":","\r"))
//
//		SetDataFolder saveDF
//	
//		ModifyGraph muloffset($S_value)={0,off}
//	else
//		controlinfo de_rup_popup0
//
//		ModifyGraph muloffset($S_value)={0,0}
//	endif
//
//End

Function DE_rup_SVP(ctrlName,varNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable varNum
	String varStr
	String varName

	
	
End

Function/S DE_rup_PMP_1()

	string list=DE_PrintAllFolders_String("*")
	return list
End


Function/S DE_rup_PMP_2()
	String saveDF
	saveDF = GetDataFolder(1)
	controlinfo de_rup_popup1
	
	
	SetDataFolder s_value

	String list = WaveList("*", ";", "")
	SetDataFolder saveDF
	
	return list
End

//Function/S DE_rup_PMP_3()
//	String saveDF
//	saveDF = GetDataFolder(1)
//	controlinfo de_rup_popup1
//	SetDataFolder s_value
//	String list = WaveList("*", ";", "")
//	SetDataFolder saveDF
//	
//	return list
//End
//
//Function/S DE_rup_PMP_4()
//	String saveDF
//	saveDF = GetDataFolder(1)
//	controlinfo de_rup_popup1
//	SetDataFolder s_value
//	String list = WaveList("*", ";", "")
//	SetDataFolder saveDF
//	
//	return list
//End

Function DE_rup_BP(ctrlName) : ButtonControl
	String ctrlName
	String saveDF
	saveDF = GetDataFolder(1)
	controlinfo de_rup_popup1
	SetDataFolder s_value
	controlinfo de_rup_popup0
	string Name1=s_value

	controlinfo de_rup_popup2

//	if(cmpstr(Name1,"X")==0||cmpstr(Name2,"X")==0)
//		print "Waves no Selected"
//		return 0
//	endif 
	Dowindow RuptureWindow
	
	if(v_flag==1)
		killwindow RuptureWindow
	else

	endif
	
	display/N=RuptureWindow $name1
	ModifyGraph offset($name1)={0,0};
	ModifyGraph muloffset($name1)={0,-1}
	SetWindow RuptureWindow, hook(MyHook)=DE_Rup_Plot_Hook // Install window hook

	//appendtograph $name2
	//controlinfo/w=Panel1 de_rup_setvar0
	//ModifyGraph/W=FastWindow offset($Name1)={v_value,0}
	
	//controlinfo/w=Panel1 de_rup_check0
//	if(v_value==1)
//		controlinfo/w=Panel1 de_rup_popup1
//		SetDataFolder s_value
//		wave w1=$name1
//		controlinfo/w=Panel1 de_rup_popup0
//		variable off=str2num(stringbykey("Invols",note(w1),":","\r"))
//		print "HI"
//		print Off
//		SetDataFolder saveDF
//	
//		ModifyGraph muloffset($S_value)={0,off}
//	else
//
//		controlinfo/w=Panel1 de_rup_popup0
//		print s_value
//		ModifyGraph muloffset($S_value)={0,0}
//	endif

	
End

Function DE_rup_BP_2(ctrlName) : ButtonControl
	String ctrlName
	String saveDF
	saveDF = GetDataFolder(1)
	controlinfo de_rup_popup1
	SetDataFolder s_value
		controlinfo de_rup_popup0

	wave w1=$S_value
		controlinfo de_rup_setvar0
variable Offset=v_value
		variable Rupture
	if(cmpstr(CsrInfo(A,"RuptureWindow"),"")==0&&cmpstr(CsrInfo(B,"RuptureWindow"),"")==0)

		print "Nothing?"
	elseif(cmpstr(CsrInfo(B,"RuptureWindow"),"")==0)
		CurveFit/Q/NTHR=0 line  w1[pcsr(A,"RuptureWindow")-Offset,pcsr(A,"RuptureWindow")] /D 
	wave w_coef
	Rupture= w_coef[0]+w_coef[1]*xcsr(A,"RuptureWindow")
	elseif(cmpstr(CsrInfo(A,"RuptureWindow"),"")==0)
		CurveFit/Q/NTHR=0 line  w1[pcsr(B,"RuptureWindow")-Offset,pcsr(B,"RuptureWindow")] /D 
	wave w_coef
	Rupture= w_coef[0]+w_coef[1]*xcsr(B,"RuptureWindow")
	else
		CurveFit/Q/NTHR=0 line  w1[pcsr(A,"RuptureWindow"),pcsr(B,"RuptureWindow")] /D 
	wave w_coef
	if(xcsr(A,"RuptureWindow")>pcsr(B,"RuptureWindow"))
		Rupture= w_coef[0]+w_coef[1]*xcsr(A,"RuptureWindow")

	else
	
		Rupture= w_coef[0]+w_coef[1]*xcsr(B,"RuptureWindow")
endif
	Rupture= w_coef[0]+w_coef[1]*xcsr(B,"RuptureWindow")
	endif
	//wavestats/q/r=[numpnts(w1),numpnts(w1)-100] w1
	//print v_avg
	string fitname="fit_"+nameofwave(w1)
		ModifyGraph/W=RuptureWindow muloffset($fitname)={0,-1},rgb($fitname)=(0,0,0)
	print Rupture*-1
	print w_coef[1]*-1
	killwaves w_coef





	SetDataFolder saveDF
	

	
End

Window DE_RupProc() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(150,77,450,277)
	
	PopupMenu de_rup_popup1,pos={21,2},size={120,21}
	PopupMenu de_rup_popup1,mode=1,popvalue="root:",value= #"DE_rup_PMP_1()"
	PopupMenu de_rup_popup0,pos={150,2},size={150,21}
	PopupMenu de_rup_popup0,mode=1,popvalue="ForceWave",value= #"DE_rup_PMP_2()"
	
//	PopupMenu de_rup_popup2,pos={80,50},size={141,21}
	//PopupMenu de_rup_popup2,mode=3,popvalue="X",value= #"DE_rup_PMP_3()"
//	PopupMenu de_rup_popup3,pos={20,130},size={141,21}
//	PopupMenu de_rup_popup3,mode=3,popvalue="X",value= #"DE_rup_PMP_4()"

	Button de_rup_button0,pos={121,100},size={50,20},proc=DE_Rup_BP,title="Plot"
		Button de_rup_button1,pos={151,170},size={50,20},proc=DE_Rup_BP_2,title="FitLine"

	SetVariable de_rup_setvar0,pos={181,100},size={50,16},proc=DE_Rup_SVP,value= _NUM:0
	//CheckBox de_rup_check0,pos={241,100},size={40,14},proc=DE_Rup_CBC,value= 0
EndMacro


Function DE_Rup_Plot_Hook(s)
	STRUCT WMWinHookStruct &s
	Variable hookResult = 0 // 0 if we do not handle event, 1 if we handle it.
	

	string items,tracename
	variable xoff,yoff,tracepnt
	switch(s.eventcode)
		case 3: // Keyboard event
			if(s.eventmod==9)
				if(cmpstr(TraceFromPixel(s.mouseloc.h, s.mouseloc.v, ""),"")==0)
					
				else
					
					//print stringbykey("Trace",TraceFromPixel(s.mouseloc.h, s.mouseloc.v, ""),":",";")
					//
					TraceName=stringbykey("Trace",TraceFromPixel(s.mouseloc.h, s.mouseloc.v, ""),":",";")
					TracePnt=str2num(stringbykey("HITPOINT",TraceFromPixel(s.mouseloc.h, s.mouseloc.v, ""),":",";"))
					wave Xwave=XWaveRefFromTrace("",TraceName)
					wave Ywave=TraceNameToWaveRef("RuptureWindow",TraceName)
					items= "A;B" // 2nd is divider, 3rd is checked
					//variable Xoff=xwave[TracePnt]
					//variable Yoff=ywave[tracepnt]

					PopupContextualMenu items
					switch( V_Flag )
						case 1:
							Execute "Cursor/P/W=RuptureWindow A "+TraceName+" "+num2str( TracePnt)
							break;
						case 2:
							Execute "Cursor/P/W=RuptureWindow B "+TraceName+" "+num2str( TracePnt)
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