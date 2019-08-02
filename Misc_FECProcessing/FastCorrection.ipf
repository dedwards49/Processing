#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma modulename=FastCorrection
//#include "Scholl_Panel"
#include "FEC_Fast_Corr"
#include ":\Misc_PanelPrograms\Panel Progs"
#include ":\Misc_PanelPrograms\AsylumNaming"
Function DE_FastX_CBC(ctrlName,checked) : CheckBoxControl
	String ctrlName
	Variable checked

	
	DoWindow FastWindow
	if(V_flag==0)
		return 0
	endif
	
	if(checked==1)
		string saveDF = GetDataFolder(1)
		controlinfo de_fastx_popup1
		SetDataFolder s_value

		controlinfo de_fastx_popup0
		wave w1=$S_value
		variable off=str2num(stringbykey("Invols",note(w1),":","\r"))
		SetDataFolder saveDF
	
		ModifyGraph muloffset($S_value)={0,off}
	else
		controlinfo de_fastx_popup0

		ModifyGraph muloffset($S_value)={0,0}
	endif

End

Function DE_FastX_SVP(ctrlName,varNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable varNum
	String varStr
	String varName
	DoWindow FastWindow

	if(V_flag==1)
		controlinfo de_fastx_popup0
		ModifyGraph/W=FastWindow offset($S_Value)={varnum,0}
	else
	
	endif
End

Function/S DE_FastX_PMP_1()

	string list=DE_PanelProgs#PrintAllFolders_String("*")
	return list
End


Function/S DE_FastX_PMP_2()
	String saveDF
	saveDF = GetDataFolder(1)
	controlinfo de_fastx_popup1
	
	
	SetDataFolder s_value

	String list = WaveList("*", ";", "")
	SetDataFolder saveDF
	
	return list
End

Function/S DE_FastX_PMP_3()
	String saveDF
	saveDF = GetDataFolder(1)
	controlinfo de_fastx_popup1
	SetDataFolder s_value
	String list = WaveList("*", ";", "")
	SetDataFolder saveDF
	
	return list
End

Function/S DE_FastX_PMP_4()
	String saveDF
	saveDF = GetDataFolder(1)
	controlinfo de_fastx_popup1
	SetDataFolder s_value
	String list = WaveList("*", ";", "")
	SetDataFolder saveDF
	
	return list
End

Function DE_Fastx_BP(ctrlName) : ButtonControl
	String ctrlName
	String saveDF
	saveDF = GetDataFolder(1)
	controlinfo de_fastx_popup1
	SetDataFolder s_value
	controlinfo de_fastx_popup0
	string Name1=s_value

	controlinfo de_fastx_popup2
	string Name2=s_value

	if(cmpstr(Name1,"X")==0||cmpstr(Name2,"X")==0)
		print "Waves no Selected"
		return 0
	endif 
	Dowindow FastWindow
	
	if(v_flag==1)
		killwindow FastWindow
	else

	endif
	
	display/N=FastWindow $name1
	appendtograph $name2
	controlinfo/w=Panel1 de_fastx_setvar0
	ModifyGraph/W=FastWindow offset($Name1)={v_value,0}
	
	controlinfo/w=Panel1 de_fastx_check0
	if(v_value==1)
		controlinfo/w=Panel1 de_fastx_popup1
		SetDataFolder s_value
		wave w1=$name1
		controlinfo/w=Panel1 de_fastx_popup0
		variable off=str2num(stringbykey("Invols",note(w1),":","\r"))

		SetDataFolder saveDF
	
		ModifyGraph muloffset($S_value)={0,off}
	else

		controlinfo/w=Panel1 de_fastx_popup0
		print s_value
		ModifyGraph muloffset($S_value)={0,0}
	endif

	
End

Function DE_Fastx_BP_2(ctrlName) : ButtonControl
	String ctrlName
	String saveDF
	saveDF = GetDataFolder(1)
	controlinfo de_fastx_popup1
	SetDataFolder s_value
		controlinfo de_fastx_popup0
		wave w1=$s_value
	controlinfo de_fastx_popup3
	wave w2=$s_value
		controlinfo de_fastx_setvar0
	variable offset=v_value
	GenFastX(w2,w1,-offset)
	SetDataFolder saveDF
	

	
End

Window Panel1() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(150,77,450,277)
	ShowTools/A
	PopupMenu de_fastx_popup0,pos={85,2},size={129,21}
	PopupMenu de_fastx_popup0,mode=1,popvalue="X",value= #"DE_FastX_PMP_2()"
	PopupMenu de_fastx_popup1,pos={21,2},size={50,21}
	PopupMenu de_fastx_popup1,mode=1,popvalue="root:",value= #"DE_FastX_PMP_1()"
	PopupMenu de_fastx_popup2,pos={80,50},size={141,21}
	PopupMenu de_fastx_popup2,mode=3,popvalue="X",value= #"DE_FastX_PMP_3()"
		PopupMenu de_fastx_popup3,pos={20,130},size={141,21}
	PopupMenu de_fastx_popup3,mode=3,popvalue="X",value= #"DE_FastX_PMP_4()"

	Button de_fastx_button0,pos={121,100},size={50,20},proc=DE_Fastx_BP,title="Plot"
		Button de_fastx_button1,pos={151,170},size={50,20},proc=DE_Fastx_BP_2,title="MakeX"

	SetVariable de_fastx_setvar0,pos={181,100},size={50,16},proc=DE_FastX_SVP,value= _NUM:0
	CheckBox de_fastx_check0,pos={241,100},size={40,14},proc=DE_FastX_CBC,value= 0
EndMacro



static function Naming(w1)
	wave w1
	wave w2=FastCaptureData
	note/k w2 note(w1)
	struct ForceWave Name1
	DE_Naming#WavetoStruc(nameofwave(w1),Name1)
	string Newname="Fast"+Name1.SNum+"DeflV"+"_"+Name1.SDirec
	rename w2 $NewName  
end