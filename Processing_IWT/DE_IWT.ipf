#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma ModuleName = DE_IWT
#include "C:\Users\dedwards\src_prh\IgorUtil\PythonApplications\InverseWeierstrass\InverseWeierstrass"
#include ":\Misc_PanelPrograms\AsylumNaming"
#include ":\Misc_PanelPrograms\Panel Progs"
#include "DE_Filtering"
Static Function RemovePauses(Forcein,SepIn,ForceOut,Sepout)
	wave Forcein,SepIn,ForceOut,Sepout

	string PauseLocs= stringbykey("DE_PauseLoc",note(ForceIn),":","\r")
	print PauseLocs
	variable n=0
	variable maxpnts=numpnts(ForceIn)
	variable delstart
	variable delnumber
	duplicate/free Forcein TempForceOut
	duplicate/free Sepin TempSepOut

	do
		if(str2num(stringfromlist(n,PauseLocs))>maxpnts)
			PauseLocs=removefromlist((stringfromlist(n,PauseLocs)),PauseLocs)

		else
			n+=1
		endif
	while(n<itemsinlist(PauseLocs))
	for(n=itemsinlist(PauseLocs)-2;n>=0;n-=2)
		delstart=str2num(stringfromlist(n,PauseLocs))
		delnumber=str2num(stringfromlist(n+1,PauseLocs))- str2num(stringfromlist(n,PauseLocs))
		deletepoints delstart,delnumber, TempForceOut,TempSepOut
	endfor
	duplicate/o TempSepOut SepOut
	duplicate/o TempForceOut ForceOut
end

Static Function OutportForce(Finwave,Sinwave)
	wave Finwave,Sinwave
	duplicate/o FInwave '0000Force'
	duplicate/o Sinwave '0000Sep'
	wave ForceWave='0000Force'
	wave SepWave='0000Sep'
	RemovePauses(Finwave,Sinwave,ForceWave,SepWave)
	
	variable K=FindGoodK(ForceWave)
	
	if(cmpstr("",stringbykey("K",note(ForceWave),"=","\r"))==0)
		note/K ForceWave, replacestringbykey("K",note(ForceWave),num2str(K),"=","\r")
	endif


	if(cmpstr("",stringbykey("K",note(SepWave),"=","\r"))==0)
				note/K ForceWave, replacestringbykey("K",note(ForceWave),num2str(K),"=","\r")

	endif
	
	
	if(cmpstr("",stringbykey("TriggerTime",note(ForceWave),":","\r"))==0)
		if(cmpstr("",stringbykey("TriggerTime1",note(ForceWave),":","\r"))!=0)
			note/K ForceWave, replacestringbykey("TriggerTime",note(ForceWave),stringbykey("TriggerTime1",note(ForceWave),":","\r"),":","\r")
		endif
	endif


	if(cmpstr("",stringbykey("TriggerTime",note(SepWave),":","\r"))==0)
		if(cmpstr("",stringbykey("TriggerTime1",note(SepWave),":","\r"))!=0)
			note/K SepWave, replacestringbykey("TriggerTime",note(SepWave),stringbykey("TriggerTime1",note(SepWave),":","\r"),":","\r")
		endif
	endif
	display/N=TMP_D Forcewave vs SepWave 
	String Path="D:\Data\InverseWeierstrass\Hold.pxp"
	SaveGraphCopy/o as Path
	KillWindow TMP_D
	//killwaves ForceWave SepWave
	
end

Static Function RunIWTonOutput()

	String Location ="D:\Data\InverseWeierstrass\Hold.pxp"
	Main("C:/Users/dedwards/src_prh/",Input_file=Location)

end

Static Function FindGoodK(WaveIn)
	wave wavein
	variable returnK
	if(cmpstr("",stringbykey("Spring Constant",note(wavein),":","\r"))!=0)
	returnK=str2num(stringbykey("Spring Constant",note(wavein),":","\r"))
elseif(cmpstr("",stringbykey("SpringConstant",note(wavein),":","\r"))!=0)
	returnK=str2num(stringbykey("SpringConstant",note(wavein),":","\r"))
	elseif(cmpstr("",stringbykey("K",note(wavein),":","\r"))!=0)
	returnK=str2num(stringbykey("K",note(wavein),":","\r"))
		elseif(cmpstr("",stringbykey("k",note(wavein),":","\r"))!=0)
	returnK=str2num(stringbykey("k",note(wavein),":","\r"))
			elseif(cmpstr("",stringbykey("K",note(wavein),"=","\r"))!=0)
	returnK=str2num(stringbykey("K",note(wavein),"=","\r"))
			elseif(cmpstr("",stringbykey("k",note(wavein),"=","\r"))!=0)
	returnK=str2num(stringbykey("k",note(wavein),"=","\r"))
endif
return returnK
end


Static StrConstant DEF_INPUT_REL_TO_BASE =  "IgorUtil/PythonApplications/InverseWeierstrass/Example/input.pxp"

Static Function Main_Windows()
	// Runs a simple IWT on patrick's windows setup
	Main("C:/Users/dedwards/src_prh/")
End Function 

Static Function Main_Mac()
	// Runs a simple IWT on patrick's mac setup 
	Main("/Users/patrickheenan/src_prh/")
End Function
Static Function AppendToWaveWOverlap(StartingWave,WavetoAdd)

	wave StartingWave,WavetoAdd
	
	duplicate/o startingwave test
	insertpoints (numpnts(startingwave)-1), (numpnts(wavetoAdd)-1), Test

	Test[numpnts(startingwave)-1,]=WavetoAdd[p-numpnts(startingwave)+1]

end

Static Function/S CombineWavesforRamp(Number,loops)
	variable number,loops
	variable i

	newdatafolder/o HoldingWaves
	DFREF Folderto=HoldingWaves
	//Initialization withfirst retaction
	String Name="ph6p2"
	string Source=DE_Naming#StringCreate(Name,number+1,"Force","Ret")
	string Target=DE_Naming#StringCreate(Name,number+1,"Force","All")
	wave w1=$Source
	
	duplicate/o w1 $Target
	movewave w1 Folderto
	wave w2=$Target
	String NoteString=note(w2)
	String DirString="Ret;"
	String IDString=num2str(number+1)+"_Ret;"

	String IndStr=num2str(numpnts(w2)-1)+";"
	NoteString=ReplaceStringbyKey("DE_Dir",NoteString,DirString,":","\r")
	NoteString=ReplaceStringbyKey("DE_Ind",NoteString,IndStr,":","\r")
	NoteString=ReplaceStringbyKey("DE_ID",NoteString,IDString,":","\r")

	Note/K  W2,NoteString
	
	
	Source=DE_Naming#StringCreate(Name,number+1,"Force","Ext")
	wave w1=$Source
	AppendToWaveWOverlap(w2,w1)
		movewave w1 Folderto

	wave Test
	duplicate/o Test w2
	NoteString=note(w2)
	DirString+="Ext;"
	IDString+=num2str(number+1)+"_Ext;"
	IndStr+=num2istr(numpnts(w2)-1)+";"
	 
	NoteString=ReplaceStringbyKey("DE_Dir",NoteString,DirString,":","\r")
	NoteString=ReplaceStringbyKey("DE_Ind",NoteString,IndStr,":","\r")
	NoteString=ReplaceStringbyKey("DE_ID",NoteString,IDString,":","\r")

	Note/K  W2,NoteString



	For(i=0;i<(loops-1);i+=1)
	
		Source=DE_Naming#StringCreate(Name,(number+i+2),"Force","Ret")
		wave w1=$Source
		AppendToWaveWOverlap(w2,w1)
			movewave w1 Folderto

		duplicate/o Test w2
		DirString+="Ret;"
		IDString+=num2str(number+i+2)+"_Ret;"

		IndStr+=num2istr(numpnts(w2)-1)+";"
		NoteString=ReplaceStringbyKey("DE_Dir",NoteString,DirString,":","\r")
		NoteString=ReplaceStringbyKey("DE_Ind",NoteString,IndStr,":","\r")
		NoteString=ReplaceStringbyKey("DE_ID",NoteString,IDString,":","\r")

		Note/K  W2,NoteString
		Source=DE_Naming#StringCreate(Name,(number+i+2),"Force","Ext")
		wave w1=$Source
		AppendToWaveWOverlap(w2,w1)
			movewave w1 Folderto

		duplicate/o Test w2
		DirString+="App;"
		IDString+=num2str(number+i+2)+"_Ext;"

		IndStr+=num2istr(numpnts(w2)-1)+";"
		print/D num2istr(numpnts(w2)-1)
		NoteString=ReplaceStringbyKey("DE_Dir",NoteString,DirString,":","\r")
		NoteString=ReplaceStringbyKey("DE_Ind",NoteString,IndStr,":","\r")
		NoteString=ReplaceStringbyKey("DE_ID",NoteString,IDString,":","\r")

		Note/K  W2,NoteString

		//Source=DE_Naming#StringCreate("Image",(number+i+1),"Force","Towd")
		//wave w1=$Source
		//AppendToWaveWOverlap(w2,w1)
		//duplicate/o Test w2

	endfor
//		deletepoints numpnts(w2)-1,1, w2

	Source=DE_Naming#StringCreate(Name,number+1,"Sep","Ret")
	Target=DE_Naming#StringCreate(Name,number+1,"Sep","All")
	wave w1=$Source
	duplicate/o w1 $Target
		movewave w1 Folderto

	wave w2=$Target
	
	//Now we add the first retraction wave
	Source=DE_Naming#StringCreate(Name,number+1,"Sep","Ext")
	wave w1=$Source
	AppendToWaveWOverlap(w2,w1)
		movewave w1 Folderto

	wave Test
	duplicate/o Test w2

	//OK, now we start piling all this shit ON
	For(i=0;i<(loops-1);i+=1)
	
		Source=DE_Naming#StringCreate(Name,(number+i+2),"Sep","Ret")
		wave w1=$Source
		AppendToWaveWOverlap(w2,w1)
			movewave w1 Folderto

		duplicate/o Test w2
		Source=DE_Naming#StringCreate(Name,(number+i+2),"Sep","Ext")
		wave w1=$Source
		AppendToWaveWOverlap(w2,w1)
			movewave w1 Folderto

		duplicate/o Test w2

	endfor
	//deletepoints numpnts(w2)-1,1, w2
	DE_IWT#DumpWavesToFolder(Folderto)
	return nameofwave(w2)
end

Static Function DumpWavestoFolder(TargetFolder)
	dfref TargetFolder
	String list = WaveList("*_Ret", ";", "")
	variable n
	for(n=0;n<itemsinlist(list);n+=1)
		MoveWave $stringfromlist(n,list) TargetFolder
	
	endfor

	
	list = WaveList("*_Ext", ";", "")
	for(n=0;n<itemsinlist(list);n+=1)
		MoveWave $stringfromlist(n,list) TargetFolder
	
	endfor
		 list = WaveList("*_Towd", ";", "")
	for(n=0;n<itemsinlist(list);n+=1)
		MoveWave $stringfromlist(n,list) TargetFolder
	
	endfor
		 list = WaveList("*_Away", ";", "")
	for(n=0;n<itemsinlist(list);n+=1)
		MoveWave $stringfromlist(n,list) TargetFolder
	
	endfor
end



Static Function Main(base,[input_file,Pull])
	// // This function shows how to use the IWT code
	// Args:
	//		base: the folder where the Research Git repository lives 
	//		input_file: the pxp to load. If not present, defaults to 
	//		<base>DEF_INPUT_REL_TO_BASE
	String base,input_file,Pull
	if (ParamIsDefault(input_file))
		input_file  = base +DEF_INPUT_REL_TO_BASE
	EndIf
	if (ParamIsDefault(Pull))
		Pull  = "No"
	EndIf
	//KillWaves /A/Z
	// IWT options
	
	//
	Struct InverseWeierstrassOptions opt

	If(cmpstr(Pull,"Yes")==0)
		if(WaveExists(root:DE_IWT:MenuStuff:ParmWave)==1)
		wave/T parmWave=root:DE_IWT:MenuStuff:ParmWave
		opt.number_of_pairs = str2num(ParmWave[0][1])
		opt.number_of_bins = str2num(ParmWave[1][1])
		opt.z_0 = str2num(ParmWave[2][1])*1e-9
		opt.velocity_m_per_s = str2num(ParmWave[3][1])*1e-6
		opt.kbT = str2num(ParmWave[4][1])
		opt.f_one_half_N = str2num(ParmWave[5][1])*1e-12
		opt.flip_forces = str2num(ParmWave[6][1])
	
	
		else
		opt.number_of_pairs = 10
		opt.number_of_bins = 150
		opt.z_0 = 20e-9
		opt.velocity_m_per_s = 400e-9
		opt.kbT = 4.1e-21
		opt.f_one_half_N = 6.5e-12
		opt.flip_forces = 1
		endif
	
	else
	
		opt.number_of_pairs = 10
		opt.number_of_bins = 150
		opt.z_0 = 20e-9
		opt.velocity_m_per_s = 400e-9
		opt.kbT = 4.1e-21
		opt.f_one_half_N = 6.5e-12
		opt.flip_forces = 1
		
	
	endif
opt.meta.path_to_input_file = input_file
		opt.meta.path_to_research_directory = base
		
	// Make the output waves
	Struct InverseWeierstrassOutput output
	Make /O/N=0, output.molecular_extension_meters
	Make /O/N=0, output.energy_landscape_joules 
	Make /O/N=0, output.tilted_energy_landscape_joules
	// Execte the command
	ModInverseWeierstrass#inverse_weierstrass(opt,output)
	// Make a fun plot wooo
	note/K output.tilted_energy_landscape_joules, num2str(opt.f_one_half_N )
	duplicate/o output.tilted_energy_landscape_joules LS_KB
	LS_KB/=4.11e-21

End Function


Static Function ButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	Strswitch(ba.ctrlname)
	string saveDF
		case "de_IWT_button0":
			switch( ba.eventCode )
				case 2: // mouse up
					saveDF = GetDataFolder(1)
					controlinfo de_IWT_popup0
					SetDataFolder s_value

					controlinfo de_IWT_popup1
					wave w1=$S_value
					controlinfo de_IWT_popup2
					wave w2=$S_value
					
					controlinfo de_IWT_setvar2
					if(V_value==1)					
					OutportForce(w1,w2)
					else
						duplicate/o w1 Testw1
						duplicate/o w2 Testw2
						DE_Filtering#FilterForceSep(w1,w2,Testw1,Testw2,"SVG",V_Value)
						OutportForce(Testw1,Testw2)
						killwaves Testw1,Testw2
					endif
					String Location ="D:\Data\InverseWeierstrass\Hold.pxp"
					Main("C:/Users/dedwards/src_prh/",Input_file=Location,Pull="Yes")
					
						struct ForceWave StartName
						DE_Naming#WavetoStruc(nameofwave(w1),StartName)
						string BaseName=startname.Name+startname.SNum
						wave w3=LS_KB
						wave w4=molecular_extension_meters
						duplicate/o w3 $(Basename+"_LS")
						duplicate/o w4 $(Basename+"_MX")
							doWindow IWTWindow
					if(V_flag==1)
					killwindow IWTWindow
					else
					endif
						Display/N=IWTWindow $(Basename+"_LS") vs $(Basename+"_MX")

						wave energy_landscape_joules,tilted_energy_landscape_joules
						killwaves w3,w4, energy_landscape_joules,tilted_energy_landscape_joules
					SetDataFolder saveDF

					break
				case -1: // control being killed
					break
			endswitch
			break
	
		case "de_IWT_button1":
			switch( ba.eventCode )
				case 2: // mouse up
							saveDF = GetDataFolder(1)
					controlinfo de_IWT_popup0
					SetDataFolder s_value
					controlinfo de_IWT_setvar0
					variable v1=V_value
					controlinfo de_IWT_setvar1
					variable v2=V_value
					string sepwavename
					sepwavename=CombineWavesforRamp(v1,v2)
					wave SepTest=$sepwavename
					SetDataFolder saveDF

				case -1: // control being killed
					break
			endswitch
			break
		case "de_IWT_button2":
			switch( ba.eventCode )
				case 2: // mouse up
				
					saveDF = GetDataFolder(1)
					controlinfo de_IWT_popup0
					SetDataFolder s_value

					controlinfo de_IWT_popup1
					wave w1=$S_value
					controlinfo de_IWT_popup2
					wave w2=$S_value
					display w2
					print itemsinlist(stringbykey("DE_Ind",note(w1),":","\r"))/2
		
					SetDataFolder saveDF
				case -1: // control being killed
					break
			endswitch
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

Static Function PopMenuProc(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa

	switch( pa.eventCode )
		case 2: // mouse up
			Variable popNum = pa.popNum
			String popStr = pa.popStr
			break
		case -1: // control being killed
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
		case "de_IWT_setvar2":
			
			Variable NewNumber=floor(varNum/2)*2+1
			SetVariable de_IWT_setvar2, value= _NUM:NewNumber
		break
	endswitch
End		

Macro IWTPanel() : Panel

	PauseUpdate; Silent 1		// building window...
	NewPanel/N=IWTPanel /W=(697,267,1361,653)
	NewDataFolder/o root:DE_IWT
	NewDataFolder/o root:DE_IWT:MenuStuff

	DE_IWT#UpdateParmWave()
	Button de_IWT_button0,pos={250,110},size={50,20},proc=DE_IWT#ButtonProc,title="IWT!"
	Button de_IWT_button2,pos={250,150},size={50,20},proc=DE_IWT#ButtonProc,title="Numbers"

	PopupMenu de_IWT_popup0,pos={250,2},size={129,21}
	PopupMenu de_IWT_popup0,mode=1,popvalue="X",value= #"DE_IWT#ListFolders()"
	PopupMenu de_IWT_popup1,pos={250,40},size={129,21}
	PopupMenu de_IWT_popup1,mode=1,popvalue="X",value= #"DE_IWT#ListWaves(\"*force*\")"
	PopupMenu de_IWT_popup2,pos={250,80},size={129,21}
	PopupMenu de_IWT_popup2,mode=1,popvalue="X",value= #"DE_IWT#ListWaves(\"*sep*\")"
	
	ListBox list0,pos={400,2},size={175,150},proc=DE_IWT#ListBoxProc,listWave=root:DE_IWT:MenuStuff:ParmWave
	ListBox list0,selWave=root:DE_IWT:MenuStuff:SelWave,editStyle= 2,userColumnResize= 1,widths={70,40,70,40}
	
	
	SetVariable de_IWT_setvar0,pos={2,2},size={150,16},proc=DE_IWT#SVP,title="Initial Ramp Number"
	SetVariable de_IWT_setvar0,limits={0,inf,1},value= _NUM:0
	SetVariable de_IWT_setvar1,pos={2,25},size={150,16},proc=DE_IWT#SVP,title="Number of Ramps"
	SetVariable de_IWT_setvar1,limits={0,inf,1},value= _NUM:10
	
	SetVariable de_IWT_setvar2,pos={250,140},size={100,16},proc=DE_IWT#SVP,title="Filtering"
	SetVariable de_IWT_setvar2,limits={-inf,inf,2},value= _NUM:1
	
	Button de_IWT_button1,pos={2,50},size={80,20},proc=DE_IWT#ButtonProc,title="Stack Curves"

	
	SetDrawEnv fillpat= 0, linethick= 3.00;DelayUpdate
	DrawRect -10,-10,160,160
	SetDrawEnv fillpat= 0, linethick= 3.00;DelayUpdate
	DrawRect 160,-10,600,160
	
	
EndMacro



Static Function UpdateParmWave()
	if(exists("root:DE_IWT:MenuStuff:ParmWave")==1)
		wave/t/z Par=root:DE_IWT:MenuStuff:ParmWave
		wave/z Sel=root:DE_IWT:MenuStuff:SelWave
	Else
		make/t/n=(7,2) root:DE_IWT:MenuStuff:ParmWave
		wave/t/z Par=root:DE_IWT:MenuStuff:ParmWave
		make/n=(7,2) root:DE_IWT:MenuStuff:SelWave
		wave/z Sel=root:DE_IWT:MenuStuff:SelWave
		
		Par[0][0]={"Number of Pairs","Bins","Z0 (nm)","Velocity (um/s)","KbT","Fhalf (pN)","Flip"}
		Par[0][1]={"10","100","20","0.4","4.1e-21","0","1"}
		Sel[][0]=0
		Sel[][1]=2
	endif

	
end


Static Function/S ListWaves(SearchStr)
string SearchStr
	String saveDF
	saveDF = GetDataFolder(1)
	controlinfo de_IWT_popup0
	SetDataFolder s_value
	String list = WaveList(SearchStr, ";", "")
	SetDataFolder saveDF
	return list

end

Static Function/S ListFolders()

	string list=DE_PanelProgs#PrintAllFolders_String("*")
	return list
End

Menu "IWT"
	//SubMenu "Processing"
	"Open IWT", IWTPanel()


	//end
	
end
