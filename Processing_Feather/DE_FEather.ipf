#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma ModuleName = DE_Feather

#include "C:\Users\dedwards\src_prh\IgorUtil\PythonApplications\FEATHER\Example\MainFeather"
#include ":\Misc_PanelPrograms\AsylumNaming"
Static Function OutportForce(ForceWave,SepWave)
	wave ForceWave,SepWave

	duplicate/o ForceWave $(replaceString("Force",nameofwave(ForceWave),"Time"))
	wave TimeWave=$(replaceString("Force",nameofwave(ForceWave),"Time"))
	TimeWave=pnt2x(ForceWave,p)
	display/N=TMP_D Forcewave vs SepWave 
	Appendtograph TimeWave
	String Path="D:\Data\Feather\Hold.pxp"
	SaveGraphCopy/o as Path
	KillWindow TMP_D

end


Static function CombineWaves(ExtW,TW,RetW,AW)

	wave ExtW,TW,RetW,AW

	make/o/n=(numpnts(ExtW)+numpnts(TW)+numpnts(RetW)+numpnts(AW)-3) Test
	Test=0
	Test[0,numpnts(ExtW)-2]=ExtW[p]
	Test[numpnts(ExtW)-1,numpnts(ExtW)+numpnts(TW)-3]=TW[p-(numpnts(ExtW)-1)]
	Test[numpnts(ExtW)+numpnts(TW)-2,numpnts(ExtW)+numpnts(TW)+numpnts(RetW)-4]=RetW[p-(numpnts(ExtW)+numpnts(TW)-2)]
	Test[numpnts(ExtW)+numpnts(TW)+numpnts(RetW)-3,numpnts(ExtW)+numpnts(TW)+numpnts(RetW)+numpnts(AW)-5]=AW[p-(numpnts(ExtW)+numpnts(TW)+numpnts(RetW)-2)]

	print numpnts(ExtW)+numpnts(TW)+numpnts(RetW)-3
	SetScale/P x 0,dimdelta(ExtW,0),"", Test
end


Static function SimpleCombine(ExtW,RetW)

	wave ExtW,RetW

	make/o/n=(numpnts(ExtW)+numpnts(RetW)-1) Test
	Test=0
	Test[0,numpnts(ExtW)-2]=ExtW[p]
	Test[numpnts(ExtW)-1,]=RetW[p-(numpnts(ExtW)-1)]
	SetScale/P x 0,dimdelta(ExtW,0),"", Test
end



Static Function AppendToWaveWOverlap(StartingWave,WavetoAdd)

	wave StartingWave,WavetoAdd
	duplicate/o startingwave test
	insertpoints (numpnts(startingwave)-1), (numpnts(wavetoAdd)-1), Test

	Test[numpnts(startingwave)-1,]=WavetoAdd[p-numpnts(startingwave)+1]

end

Static Function ConvertToSingle(Number,[eventwave,outwave])
	variable Number
	wave eventwave,outwave
		
	string ForceTarget=DE_Naming#StringCreate("Image",number,"Force","All")
	string SepTarget=DE_Naming#StringCreate("Image",number,"Sep","All")
	wave ForceWave=$ForceTarget
	wave SepWave=$SepTarget
	
	if(ParamIsDefault(eventwave))

		String EventList=StringbyKey("PRHEvents",note(ForceWave),":","\r")
		if(itemsinlist(	EventList,":")==0)
			print "No events in wave note, provide event list"
			return 0
		endif	
		Make/n=(itemsinlist(EventList,":"))/Free Events
		Events=str2num(stringfromlist(p,EventList,":"))
		
	else
		duplicate/free eventwave Events
	endif
	
	string SForRet=DE_Naming#StringCreate("Image",number+1,"Force","Ret")
	string SSepRet=DE_Naming#StringCreate("Image",number+1,"Sep","Ret")
	wave ForceRet=$SForRet
	wave SepRet=$SSepRet
	variable n=0
	make/free/n=(numpnts(Events)) ModEvents
	for(n=0;n<numpnts(Events);n+=1)
	
		FindValue/T=1e-17/V=(ForceWave[Events[n]]) ForceRet
		ModEvents[n]=v_value
	endfor
	//duplicate/o ModEvents $("Events_"+num2str(Number)+"_Ret")
	
	duplicate/o ModEvents outwave
	print ModEvents[0]
	EventsIntoNoteRet(Number+1,ModEvents)
end

Static Function FitWaveGlide(Number)
	Variable Number
	string ForceTarget=DE_Naming#StringCreate("Image",number,"Force","All")
	string SepTarget=DE_Naming#StringCreate("Image",number,"Sep","All")
	
	GenerateNewWaveGlide(Number)
	wave ForceWave=$ForceTarget
	wave SepWave=$SepTarget
	OutportForce(ForceWave,SepWave)
	RunFeatheronOutput()
	wave event_starts
	EventsIntoNote(Number,event_starts)
end

Static Function FitWaveRamps(Number,loops)
	Variable Number,loops
	
	string ForceTarget=DE_Naming#StringCreate("Dilute_6p2_100",number,"Force","All")
	string SepTarget=DE_Naming#StringCreate("Dilute_6p2_100",number,"Sep","All")
	
	GenerateNewWaveMultiRamp(Number,loops)
	
	wave ForceWave=$ForceTarget
	wave SepWave=$SepTarget
	
	OutportForce(ForceWave,SepWave)
	RunFeatheronOutput()
	wave event_starts
	EventsIntoNote(Number,event_starts)
	EventsIntoNoteMultiRamp(Number,event_starts)
end

Static Function RuptureLocs(Number,[eventwave])
	variable Number
	wave eventwave
	
	string ForceTarget=DE_Naming#StringCreate("Image",number,"Force","All")
	string SepTarget=DE_Naming#StringCreate("Image",number,"Sep","All")
	wave ForceWave=$ForceTarget
	wave SepWave=$SepTarget
	
	if(ParamIsDefault(eventwave))

		String EventList=StringbyKey("PRHEvents",note(ForceWave),":","\r")
		if(itemsinlist(	EventList,":")==0)
			print "No events in wave note, provide event list"
			return 0
		endif	
		Make/n=(itemsinlist(EventList,":"))/Free Events
		Events=str2num(stringfromlist(p,EventList,":"))
		
	else
		duplicate/free eventwave Events
	endif
	
	duplicate/free events RupForce,RupSep
	RupForce=ForceWave[events[p]]
	RupSep=SepWave[events[p]]
	
	duplicate/o RupForce $("RupForce_"+num2str(Number))
	duplicate/o RupSep $("RupSep_"+num2str(Number))
	
end

Static Function LoadingRate(Number,[eventwave])
	variable Number
	wave eventwave
	
	string ForceTarget=DE_Naming#StringCreate("Image",number,"Force","All")
	string SepTarget=DE_Naming#StringCreate("Image",number,"Sep","All")
	wave ForceWave=$ForceTarget
	wave SepWave=$SepTarget
	
	if(ParamIsDefault(eventwave))
		String EventList=StringbyKey("PRHEvents",note(ForceWave),":","\r")

		if(itemsinlist(	EventList,":")==0)
			print "No events in wave note, provide event list"
			return 0
		endif	
		
		Make/n=(itemsinlist(EventList,":"))/Free Events
		Events=str2num(stringfromlist(p,EventList,":"))
		
	else
		duplicate/free eventwave Events
	endif
	duplicate/free Events Loading
	variable n
	variable backwardspoints=100e-3/dimdelta(ForceWave,0)
	print backwardspoints
	for(n=0;n<numpnts(Events);n+=1)
	duplicate/free/r=[events[n]-backwardspoints,events[n]] ForceWave ForFit
	CurveFit/Q/NTHR=0 line, ForFit  
	wave w_coef
	Loading[n]=w_coef[1]
	endfor
	
	duplicate/o Loading $("Loading_"+num2str(Number))


end


Static Function EventsIntoNote(Number,event_starts)

	Variable Number
	wave event_starts
	if(WaveExists( event_starts)==0)
		Print "No Event Waves"
		return 0
	endif
	
	wave event_starts
	
	variable n=0
	string EventStr=""
	for(n=0;n<numpnts(event_starts);n+=1)
		EventStr+=(num2str(event_starts[n])+":")
		
	
	endfor
	string ForceTarget=DE_Naming#StringCreate("Image",number,"Force","All")
	string SepTarget=DE_Naming#StringCreate("Image",number,"Sep","All")
	wave ForceWave=$ForceTarget
	wave SepWave=$SepTarget
	
		
	String New=ReplaceStringbyKey("PRHEvents",note(ForceWave),EventStr,":","\r")
	note/K ForceWave New
	New=ReplaceStringbyKey("PRHEvents",note(SepWave),EventStr,":","\r")
	note/K SepWave New

end

Static Function EventsIntoNoteRet(Number,event_starts)

	Variable Number
	wave event_starts
	if(WaveExists( event_starts)==0)
		Print "No Event Waves"
		return 0
	endif
	
	
	variable n=0
	string EventStr=""
	for(n=0;n<numpnts(event_starts);n+=1)
		EventStr+=(num2str(event_starts[n])+":")
		
	
	endfor
	string ForceTarget=DE_Naming#StringCreate("Image",number,"Force","Ret")
	string SepTarget=DE_Naming#StringCreate("Image",number,"Sep","Ret")
	wave ForceWave=$ForceTarget
	wave SepWave=$SepTarget
		
	String New=ReplaceStringbyKey("PRHEvents",note(ForceWave),EventStr,":","\r")
	note/K ForceWave New
	New=ReplaceStringbyKey("PRHEvents",note(SepWave),EventStr,":","\r")
	note/K SepWave New

end


Static Function EventsIntoNoteMultiRamp(Number,event_starts)
	variable Number
	wave event_starts
	string ForceTarget=DE_Naming#StringCreate("Image",number,"Force","All")
	string SepTarget=DE_Naming#StringCreate("Image",number,"Sep","All")
	wave ForceWave=$ForceTarget
	wave SepWave=$SepTarget
	if(WaveExists( event_starts)==0)
		Print "No Event Waves"
		return 0
	endif
	String ID,Direc,New,replace
	variable num,ShiftedLoc

	variable n=0
	string EventStr=""
	for(n=0;n<numpnts(event_starts);n+=1)
		ID=PullWaveName(ForceWave,event_starts[n])
		ShiftedLoc=PullEvenSlot(ForceWave,event_starts[n])
		num=str2num(ID[0,3])
		Direc=ID[5,100]
		wave FW=$DE_Naming#StringCreate("Image",Num,"Force",Direc)
		wave SW=$DE_Naming#StringCreate("Image",Num,"Sep",Direc)
		Replace=AddListItem(num2str(ShiftedLoc),StringbyKey("PRHEvents",note(FW),":","\r"),";")
		New=ReplaceStringbyKey("PRHEvents",note(FW),Replace,":","\r")
		note/K FW New
		New=ReplaceStringbyKey("PRHEvents",note(SW),Replace,":","\r")
		note/K SW New
	endfor
end

Static Function/S PullWaveName(FullWave,PntNumber)

	Wave FulLWave
	variable PntNumber
	variable n=0
	variable previous=0
	variable current=0
	do
		if(cmpstr(stringfromlist(n,stringbykey("DE_Ind",note(FulLWave),":","\r")),"")==0)
			break
		endif

		if(n==0)
			current=str2num(stringfromlist(n,stringbykey("DE_Ind",note(FulLWave),":","\r")))
			previous=0

		else
		
			current=str2num(stringfromlist(n,stringbykey("DE_Ind",note(FulLWave),":","\r")))
			previous=str2num(stringfromlist(max(0,(n-1)),stringbykey("DE_Ind",note(FulLWave),":","\r")))
		endif

		if(PntNumber>=previous&&PntNumber<current)
			return stringfromlist(n,stringbykey("DE_ID",note(FulLWave),":","\r"))

		else
		endif

		n+=1
	while (n>-1)
end

Static Function PullEvenSlot(FullWave,PntNumber)

	Wave FulLWave
	variable PntNumber
	variable n=0
	variable previous=0
	variable current=0
	do
		if(cmpstr(stringfromlist(n,stringbykey("DE_Ind",note(FulLWave),":","\r")),"")==0)
			break
		endif

		if(n==0)
			current=str2num(stringfromlist(n,stringbykey("DE_Ind",note(FulLWave),":","\r")))
			previous=0

		else
		
			current=str2num(stringfromlist(n,stringbykey("DE_Ind",note(FulLWave),":","\r")))
			previous=str2num(stringfromlist(max(0,(n-1)),stringbykey("DE_Ind",note(FulLWave),":","\r")))
		endif

		if(PntNumber>=previous&&PntNumber<current)
			return (PntNumber-previous)

		else
		endif

		n+=1
	while (n>-1)
end



Static Function MakeTIme(number)
variable number
	string ForceTarget=DE_Naming#StringCreate("Image",number,"Force","All")
	string SepTarget=DE_Naming#StringCreate("Image",number,"Sep","All")
	wave ForceWave=$ForceTarget
	wave SepWave=$SepTarget	
			String EventList=StringbyKey("PRHEvents",note(ForceWave),":","\r")

	if(itemsinlist(	EventList,":")==0)
			print "No events in wave note, provide event list"
			return 0
		endif	
		Make/o/n=(itemsinlist(EventList,":")) Events,Stuff
		Events=pnt2x(ForceWave,str2num(stringfromlist(p,EventList,":")))
		Stuff=ForceWave(Events)
	
end


Static Function GenerateNewWaveGlide(Number)
	variable number
	string Source=DE_Naming#StringCreate("Image",number,"Force","Ext")
	string Target=DE_Naming#StringCreate("Image",number,"Force","All")
	wave w1=$Source
	duplicate/o w1 $Target
	wave w2=$Target
	Source=DE_Naming#StringCreate("Image",number,"Force","Towd")
	wave w1=$Source

	AppendToWaveWOverlap(w2,w1)
	wave Test
	duplicate/o Test w2
	Source=DE_Naming#StringCreate("Image",number,"Force","Ret")
	wave w1=$Source
	AppendToWaveWOverlap(w2,w1)
	duplicate/o Test w2
	Source=DE_Naming#StringCreate("Image",number+1,"Force","Ret")
	wave w1=$Source
	AppendToWaveWOverlap(w2,w1)
	duplicate/o Test w2

	Source=DE_Naming#StringCreate("Image",number,"Sep","Ext")
	Target=DE_Naming#StringCreate("Image",number,"Sep","All")
	wave w1=$Source
	duplicate/o w1 $Target

	wave w2=$Target
	Source=DE_Naming#StringCreate("Image",number,"Sep","Towd")
	wave w1=$Source
	AppendToWaveWOverlap(w2,w1)
	wave Test
	duplicate/o Test w2
	Source=DE_Naming#StringCreate("Image",number,"Sep","Ret")
	wave w1=$Source
	AppendToWaveWOverlap(w2,w1)
	duplicate/o Test w2
	Source=DE_Naming#StringCreate("Image",number+1,"Sep","Ret")
	wave w1=$Source
	AppendToWaveWOverlap(w2,w1)
	duplicate/o Test w2

end


Static Function GenerateNewWaveJustRamps(Number,loops)
	variable number,loops
	variable i
	
	//Initialization withfirst retaction
	string Source=DE_Naming#StringCreate("Image",number+1,"Force","Ret")
	string Target=DE_Naming#StringCreate("Image",number+1,"Force","All")
	wave w1=$Source
	duplicate/o w1 $Target
	wave w2=$Target
	String NoteString=note(w2)
	String DirString="Ret;"
	String IDString=num2str(number+1)+"_Ret;"

	String IndStr=num2str(numpnts(w2)-1)+";"
	NoteString=ReplaceStringbyKey("DE_Dir",NoteString,DirString,":","\r")
	NoteString=ReplaceStringbyKey("DE_Ind",NoteString,IndStr,":","\r")
	NoteString=ReplaceStringbyKey("DE_ID",NoteString,IDString,":","\r")

	Note/K  W2,NoteString
	
	
	Source=DE_Naming#StringCreate("Image",number+1,"Force","Ext")
	wave w1=$Source
	AppendToWaveWOverlap(w2,w1)
	wave Test
	duplicate/o Test w2
	NoteString=note(w2)
	DirString+="Ext;"
	IDString+=num2str(number+1)+"_Ext;"
	IndStr+=num2str(numpnts(w2)-1)+";"
	 
	NoteString=ReplaceStringbyKey("DE_Dir",NoteString,DirString,":","\r")
	NoteString=ReplaceStringbyKey("DE_Ind",NoteString,IndStr,":","\r")
	NoteString=ReplaceStringbyKey("DE_ID",NoteString,IDString,":","\r")

	Note/K  W2,NoteString
	

	For(i=0;i<loops;i+=1)
	
		Source=DE_Naming#StringCreate("Image",(number+i+2),"Force","Ret")
		wave w1=$Source
		AppendToWaveWOverlap(w2,w1)
		duplicate/o Test w2
		DirString+="Ret;"
		IDString+=num2str(number+i+2)+"_Ret;"

		IndStr+=num2str(numpnts(w2)-1)+";"
		NoteString=ReplaceStringbyKey("DE_Dir",NoteString,DirString,":","\r")
		NoteString=ReplaceStringbyKey("DE_Ind",NoteString,IndStr,":","\r")
		NoteString=ReplaceStringbyKey("DE_ID",NoteString,IDString,":","\r")

		Note/K  W2,NoteString
		Source=DE_Naming#StringCreate("Image",(number+i+2),"Force","Ext")
		wave w1=$Source
		AppendToWaveWOverlap(w2,w1)
		duplicate/o Test w2
		DirString+="App;"
		IDString+=num2str(number+i+2)+"_Ext;"

		IndStr+=num2str(numpnts(w2)-1)+";"
		NoteString=ReplaceStringbyKey("DE_Dir",NoteString,DirString,":","\r")
		NoteString=ReplaceStringbyKey("DE_Ind",NoteString,IndStr,":","\r")
		NoteString=ReplaceStringbyKey("DE_ID",NoteString,IDString,":","\r")

		Note/K  W2,NoteString

		//Source=DE_Naming#StringCreate("Image",(number+i+1),"Force","Towd")
		//wave w1=$Source
		//AppendToWaveWOverlap(w2,w1)
		//duplicate/o Test w2

	endfor


	Source=DE_Naming#StringCreate("Image",number+1,"Sep","Ret")
	Target=DE_Naming#StringCreate("Image",number+1,"Sep","All")
	wave w1=$Source
	duplicate/o w1 $Target
	wave w2=$Target
	
	//Now we add the first retraction wave
	Source=DE_Naming#StringCreate("Image",number+1,"Sep","Ext")
	wave w1=$Source
	AppendToWaveWOverlap(w2,w1)
	wave Test
	duplicate/o Test w2

	//OK, now we start piling all this shit ON
	For(i=0;i<loops;i+=1)
	
		Source=DE_Naming#StringCreate("Image",(number+i+2),"Sep","Ret")
		wave w1=$Source
		AppendToWaveWOverlap(w2,w1)
		duplicate/o Test w2
		Source=DE_Naming#StringCreate("Image",(number+i+2),"Sep","Ext")
		wave w1=$Source
		AppendToWaveWOverlap(w2,w1)
		duplicate/o Test w2

	endfor

	
end


Static Function GenerateNewWaveMultiRamp(Number,loops)
	variable number,loops
	variable i
	
	//Initialization with the extension wave.
	string Source=DE_Naming#StringCreate("Dilute_6p2_100",number,"Force","Ext")
	string Target=DE_Naming#StringCreate("Dilute_6p2_100",number,"Force","All")
	wave w1=$Source
	duplicate/o w1 $Target
	wave w2=$Target

	String NoteString=note(w2)
	String DirString="App;"
	String IDString=num2str(number)+"_Ext;"

	String IndStr=num2str(numpnts(w2)-1)+";"
	NoteString=ReplaceStringbyKey("DE_Dir",NoteString,DirString,":","\r")
	NoteString=ReplaceStringbyKey("DE_Ind",NoteString,IndStr,":","\r")
	NoteString=ReplaceStringbyKey("DE_ID",NoteString,IDString,":","\r")

	Note/K  W2,NoteString
	
	
	//Now we add the first retraction wave
	Source=DE_Naming#StringCreate("Dilute_6p2_100",number,"Force","Towd")
	wave w1=$Source
	AppendToWaveWOverlap(w2,w1)
	wave Test
	duplicate/o Test w2
	NoteString=note(w2)
	DirString+="Surf;"
	IDString+=num2str(number)+"_Towd;"
	IndStr+=num2str(numpnts(w2)-1)+";"
	 
	NoteString=ReplaceStringbyKey("DE_Dir",NoteString,DirString,":","\r")
	NoteString=ReplaceStringbyKey("DE_Ind",NoteString,IndStr,":","\r")
	NoteString=ReplaceStringbyKey("DE_ID",NoteString,IDString,":","\r")

	Note/K  W2,NoteString
	
	
	
	Source=DE_Naming#StringCreate("Dilute_6p2_100",number,"Force","Ret")
	wave w1=$Source
	AppendToWaveWOverlap(w2,w1)

	wave Test
	duplicate/o Test w2
	DirString+="IRet;"
	IDString+=num2str(number)+"_Ret;"

	IndStr+=num2str(numpnts(w2)-1)+";"
	NoteString=ReplaceStringbyKey("DE_Dir",NoteString,DirString,":","\r")
	NoteString=ReplaceStringbyKey("DE_Ind",NoteString,IndStr,":","\r")
	NoteString=ReplaceStringbyKey("DE_ID",NoteString,IDString,":","\r")

	Note/K  W2,NoteString
	//OK, now we start piling all this shit ON
	For(i=0;i<loops;i+=1)
	
		Source=DE_Naming#StringCreate("Dilute_6p2_100",(number+i+1),"Force","Ret")
		wave w1=$Source
		AppendToWaveWOverlap(w2,w1)
		duplicate/o Test w2
		DirString+="Ret;"
		IDString+=num2str(number+i+1)+"_Ret;"

		IndStr+=num2str(numpnts(w2)-1)+";"
		NoteString=ReplaceStringbyKey("DE_Dir",NoteString,DirString,":","\r")
		NoteString=ReplaceStringbyKey("DE_Ind",NoteString,IndStr,":","\r")
		NoteString=ReplaceStringbyKey("DE_ID",NoteString,IDString,":","\r")

		Note/K  W2,NoteString
		Source=DE_Naming#StringCreate("Dilute_6p2_100",(number+i+1),"Force","Ext")
		wave w1=$Source
		AppendToWaveWOverlap(w2,w1)
		duplicate/o Test w2
		DirString+="App;"
		IDString+=num2str(number+i+1)+"_Ext;"

		IndStr+=num2str(numpnts(w2)-1)+";"
		NoteString=ReplaceStringbyKey("DE_Dir",NoteString,DirString,":","\r")
		NoteString=ReplaceStringbyKey("DE_Ind",NoteString,IndStr,":","\r")
		NoteString=ReplaceStringbyKey("DE_ID",NoteString,IDString,":","\r")

		Note/K  W2,NoteString

		//Source=DE_Naming#StringCreate("Image",(number+i+1),"Force","Towd")
		//wave w1=$Source
		//AppendToWaveWOverlap(w2,w1)
		//duplicate/o Test w2

	endfor
	//	
	//	
	//
	Source=DE_Naming#StringCreate("Dilute_6p2_100",(number+loops+1),"Force","Towd")
	wave w1=$Source
	AppendToWaveWOverlap(w2,w1)
	duplicate/o Test w2
	DirString+="Final;"
	IDString+=num2str(number+loops+1)+"_Towd;"

	IndStr+=num2str(numpnts(w2)-1)+";"
	NoteString=ReplaceStringbyKey("DE_Dir",NoteString,DirString,":","\r")
	NoteString=ReplaceStringbyKey("DE_Ind",NoteString,IndStr,":","\r")
	NoteString=ReplaceStringbyKey("DE_ID",NoteString,IDString,":","\r")

	Note/K  W2,NoteString

	Source=DE_Naming#StringCreate("Dilute_6p2_100",number,"Sep","Ext")
	Target=DE_Naming#StringCreate("Dilute_6p2_100",number,"Sep","All")
	wave w1=$Source
	duplicate/o w1 $Target
	wave w2=$Target
	
	//Now we add the first retraction wave
	Source=DE_Naming#StringCreate("Dilute_6p2_100",number,"Sep","Towd")
	wave w1=$Source
	AppendToWaveWOverlap(w2,w1)
	wave Test
	duplicate/o Test w2
	Source=DE_Naming#StringCreate("Dilute_6p2_100",number,"Sep","Ret")
	wave w1=$Source
	AppendToWaveWOverlap(w2,w1)
	wave Test
	duplicate/o Test w2
	
	//OK, now we start piling all this shit ON
	For(i=0;i<loops;i+=1)
	
		Source=DE_Naming#StringCreate("Dilute_6p2_100",(number+i+1),"Sep","Ret")
		wave w1=$Source
		AppendToWaveWOverlap(w2,w1)
		duplicate/o Test w2
		Source=DE_Naming#StringCreate("Dilute_6p2_100",(number+i+1),"Sep","Ext")
		wave w1=$Source
		AppendToWaveWOverlap(w2,w1)
		duplicate/o Test w2
		//Source=DE_Naming#StringCreate("Image",(number+i+1),"Sep","Towd")
		//wave w1=$Source
		//AppendToWaveWOverlap(w2,w1)
		//duplicate/o Test w2

	endfor
	//	
	//	
	//
	Source=DE_Naming#StringCreate("Dilute_6p2_100",(number+loops+1),"Sep","Towd")
	wave w1=$Source
	AppendToWaveWOverlap(w2,w1)
	duplicate/o Test w2
	
end

Static Function RunFeatheronOutput()

	String Location = "D:\Data\Feather\Hold.pxp"
	ModMainFEATHER#Main("C:/Users/dedwards/src_prh/",Input_file=Location)

end


Static Function ClearNotes(StartNum,Loops)
	variable startnum,Loops

	variable i=0
	string Source
	for(i=0;i<=Loops;i+=1)
		Source=DE_Naming#StringCreate("Image",(StartNum+i),"Force","All")
		if(Waveexists($Source))
			wave w1= $Source
			note/k w1, replacestringbykey("PRHEvents",note(w1),"",":","\r")
		endif
		Source=DE_Naming#StringCreate("Image",(StartNum+i),"Force","Ext")
		if(Waveexists($Source))
			wave w1= $Source
			note/k w1, replacestringbykey("PRHEvents",note(w1),"",":","\r")
		endif
		Source=DE_Naming#StringCreate("Image",(StartNum+i),"Force","Ret")
		if(Waveexists($Source))
			wave w1= $Source
			note/k w1, replacestringbykey("PRHEvents",note(w1),"",":","\r")

		endif
		Source=DE_Naming#StringCreate("Image",(StartNum+i),"Force","Towd")
		if(Waveexists($Source))
			wave w1= $Source
			note/k w1, replacestringbykey("PRHEvents",note(w1),"",":","\r")
		endif
		Source=DE_Naming#StringCreate("Image",(StartNum+i),"Force","Away")
	
		if(Waveexists($Source))
			wave w1= $Source
			note/k w1, replacestringbykey("PRHEvents",note(w1),"",":","\r")
		endif
		
		Source=DE_Naming#StringCreate("Image",(StartNum+i),"Sep","All")
		if(Waveexists($Source))
			wave w1= $Source
			note/k w1, replacestringbykey("PRHEvents",note(w1),"",":","\r")
		endif
		Source=DE_Naming#StringCreate("Image",(StartNum+i),"Sep","Ext")
		if(Waveexists($Source))
			wave w1= $Source
			note/k w1, replacestringbykey("PRHEvents",note(w1),"",":","\r")
		endif
		Source=DE_Naming#StringCreate("Image",(StartNum+i),"Sep","Ret")
		if(Waveexists($Source))
			wave w1= $Source
			note/k w1, replacestringbykey("PRHEvents",note(w1),"",":","\r")
		endif
		Source=DE_Naming#StringCreate("Image",(StartNum+i),"Sep","Towd")
		if(Waveexists($Source))
			wave w1= $Source
			note/k w1, replacestringbykey("PRHEvents",note(w1),"",":","\r")
		endif
		Source=DE_Naming#StringCreate("Image",(StartNum+i),"Sep","Away")
	
		if(Waveexists($Source))
			wave w1= $Source
			note/k w1, replacestringbykey("PRHEvents",note(w1),"",":","\r")
		endif
	
	endfor

End
