#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma ModuleName = SimpleCameraGizmo


Window SimpleCameraViewer() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(50,50,1200,700)
	Button  button0,pos={550,250},size={150,25},proc= SimpleCameraGizmo#ButtonProc
	Button button0 title="Load Data"
	SetVariable PathName proc=SimpleCameraGizmo#SetVarProc, pos={550,50}, size={450,50},  value= _STR:""
	SetVariable BaseName proc=SimpleCameraGizmo#SetVarProc, pos={550,100}, size={250,50},  value= _STR:"Image"

	//	
EndMacro

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


Static function LoadJPEG(Folder,basename,num)
	string Folder,basename
	variable num
	string File
	NewPath/o/Q Hold Folder
	File=basename+num2str(num)+".jpg"
	ImageLoad/Q/P=Hold File
	wave w1=$stringfromlist(0,S_waveNames)
	duplicate/o w1 Currentimage
	killwaves w1
	
end


Static function UpdateThatShit()
	Controlinfo BaseName
	string BN=S_value
	Controlinfo PathName
	string S_path=S_value
	variable point=pcsr(A,"SimpleCameraViewer#New")
	Controlinfo LoadAll
	variable all=v_value
	string wavenames
	controlinfo check0
	LoadJPEG(S_path,BN,point)

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
			if(cmpstr(TraceNameList("SimpleCameraViewer#New", ";",1 ),"")==0)
				display/W=(50,50,500,300)/Host=SimpleCameraViewer/N=New
				SetWindow SimpleCameraViewer, hook(MyHook) = SimpleCameraGizmo#MyWindowHook
				Appendtograph/W=SimpleCameraViewer#New Testing1 vs Testing0
				AppendToGraph/W=SimpleCameraViewer#New/L=L1/B=B1 Testing2 vs Testing0
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
			LOADJPEG(S_path,S_Value,0)

			
			display/W=(50,350,500,800)/Host=SimpleCameraViewer/N=Image
			Cursor/P/W=SimpleCameraViewer#New  A Testing1 0
			//				UpdateThatShit()

			//
			Appendimage CurrentImage
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

