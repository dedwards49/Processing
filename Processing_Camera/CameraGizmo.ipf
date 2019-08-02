#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma ModuleName = CameraGizmo
#include "Camera_SavingProcs"


Static function LOADFile(Folder,basename,num)
	string Folder,basename
	variable num
	string File
	NewPath/o/Q Hold Folder
	File=basename+num2str(num)
	LoadWave/Q/H/P=Hold/N File+".ibw"
	wave w1=$stringfromlist(0,S_waveNames)
	imagetransform decompress w1
	wave w_decompressed
	duplicate/o w_decompressed CurrentImage
	killwaves w_decompressed,w1
	//dowindow Test
	//if( V_flag==0)
	//Display/N=Test;AppendImage CurrentImage
	//endif

end

Window CameraViewer() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(50,50,1200,700)
	Button  button0,pos={550,250},size={150,25},proc= CameraGizmo#ButtonProc
	Button button0 title="Load Data"
	SetVariable PathName proc=CameraGizmo#SetVarProc, pos={550,50}, size={450,50},  value= _STR:""
	SetVariable BaseName proc=CameraGizmo#SetVarProc, pos={550,100}, size={250,50},  value= _STR:"Image"
	CheckBox  check0,pos={850,250}, proc=CameraGizmo#CheckProc

	//	
EndMacro

Static function Check(Folder,Basename)
	string Folder,Basename
	string ALL="No"
	Prompt ALL,"Load all Data?",popup,"Yes;No"
	DoPrompt "Load",All
	if(cmpstr(All,"Yes")==0)
		newdatafolder/o/s Data

		LOADAll(Folder,basename)
		setdatafolder root:

		return 1

	else
		return 0
	endif
end

Static Function MyWindowHook(s)
	STRUCT WMWinHookStruct &s
	
	Variable hookResult = 0	// 0 if we do not handle event, 1 if we handle it.

	switch(s.eventCode)
		case 7:					// Keyboard event
			UpdateThatShit()
			break
	endswitch

	return hookResult	// If non-zero, we handled event and Igor will ignore it.
End

Static function LOADAll(Folder,basename)
	string Folder,basename
	variable n=0
	variable num=0
	string File,WaveNames,fName
	variable loadflag
	NewPath/o/Q Hold Folder
	do
		fName= IndexedFile(Hold,n,".ibw")
		sscanf fName, basename+"%g.ibw", num
		if(cmpstr(fname,"")==0)
			return -1
		endif

		WaveNames="Image"+num2str(num)
		LoadWave/Q/H/P=Hold/N fName
		wave w1=$stringfromlist(0,S_waveNames)
		imagetransform decompress w1
		wave w_decompressed
		duplicate/o w_decompressed $WaveNames
		killwaves w_decompressed,w1

		n+=1
	while(1>0)
	WaveNames="Image"+num2str(0)
	duplicate $WaveNames CurrentImage
 
	////dowindow Test
	////if( V_flag==0)
	////Display/N=Test;AppendImage CurrentImage
	////endif

end

Static function LOADAllJPEG(Folder,basename)
	string Folder,basename
	variable n=0
	variable num=0
	string File,WaveNames,fName
	variable loadflag
	NewPath/o/Q Hold Folder
	
	do
		fName= IndexedFile(Hold,n,".jpg")
		sscanf fName, basename+"%g.jpg", num
		if(cmpstr(fname,"")==0)

			return -1
		endif

		WaveNames="Image"+num2str(num)
		ImageLoad/Q/P=Hold fName
		
		wave w1=$stringfromlist(0,S_waveNames)

		n+=1
	while(1<3)
	WaveNames="Image"+num2str(0)
	duplicate $WaveNames CurrentImage

end


Static function UpdateThatShit()
	Controlinfo BaseName
	string BN=S_value
	Controlinfo PathName
	string S_path=S_value
	variable point=pcsr(A,"CameraViewer#New")
	Controlinfo LoadAll
	variable all=v_value
	string wavenames
	controlinfo check0
	
	if(v_value==1)
		wavenames="Image"+num2str(point)+".jpg"
		setdatafolder root:Data
		duplicate/o $wavenames root:CurrentImage
		setdatafolder root:
	else
		if(all==0)
			LOADFile(S_path,BN,point)

		else
			wavenames="Image"+num2str(point)
			setdatafolder root:Data
			duplicate/o $wavenames root:CurrentImage
			setdatafolder root:
		endif
	endif

end

Static Function SetVarProc(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Static Function ButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			killdatafolder/Z Data
			LoadWave/Q/G/N=Testing
			wave Testing0,Testing1,Testing2
			SetVariable PathName proc=SetVarProc, pos={550,50}, size={450,50},  value= _STR:S_path
			if(cmpstr(TraceNameList("CameraViewer#New", ";",1 ),"")==0)
				display/W=(50,50,500,300)/Host=CameraViewer/N=New
				SetWindow CameraViewer, hook(MyHook) = CameraGizmo#MyWindowHook
				Appendtograph/W=CameraViewer#New Testing1 vs Testing0
				AppendToGraph/W=CameraViewer#New/L=L1/B=B1 Testing2 vs Testing0
				ModifyGraph axisEnab(left)={0,0.45}
				ModifyGraph axisEnab(l1)={0.55,1}
				ModifyGraph freePos(L1)=0,freePos(B1)=-110
				ModifyGraph mode=4,marker=19,rgb(Testing1)=(58368,6656,7168);DelayUpdate
				ModifyGraph rgb(Testing2)=(14848,32256,47104)
				ModifyGraph lblPos(L1)=50;DelayUpdate
				Label L1 "\\f01Sum (V)"
				ModifyGraph noLabel(B1)=2
				Label bottom "\\f01Zsensor (µm)"
				ModifyGraph prescaleExp(bottom)=6
				Label left "\\f01Deflection (V)"
				ModifyGraph fSize=11,axisOnTop=1,standoff=0,font="Arial"
			endif
			//
			//
			//
			//	
			make/o/n=(0,0,0) CurrentImage
			controlinfo check0
			variable JPG=V_Value
			Controlinfo BaseName

			if(JPG==1)
	
				SetVariable LoadAll proc=SetVarProc, pos={550,150}, size={100,100},  value= _num:1
				newdatafolder/o/s Data
				LOADAllJPEG(S_path,S_Value)
				setdatafolder root:

	
			else
				variable all=Check(S_path,S_Value)
				SetVariable LoadAll proc=SetVarProc, pos={550,150}, size={100,100},  value= _num:all
				if(all==0)
					//
					LOADFile(S_path,S_Value,0)
				else
					//	
				endif
			endif
			display/W=(50,350,500,800)/Host=CameraViewer/N=Image
			Cursor/P/W=CameraViewer#New  A Testing1 0
			//				UpdateThatShit()

			//
			Appendimage CurrentImage
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function CheckProc(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End