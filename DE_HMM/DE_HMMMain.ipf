// Use modern global access method, strict compilation
#pragma rtGlobals=3	
#pragma ModuleName = DE_HMMPython
#include "C:\Users\dedwards\src_prh\IgorUtil\IgorCode\Util\IoUtil"
#include "C:\Users\dedwards\src_prh\IgorUtil\IgorCode\Util\PlotUtil"
#include "C:\Users\dedwards\src_prh\IgorUtil\IgorCode\Util:OperatingSystemUtil"
#include "DE_HMM"

//#include "C:\Users\dedwards\src_prh\IgorUtil\PythonApplications\FEATHER\Example\MainFeather"

Static StrConstant path_to_HMM_folder = "D:\Devin\Python\HMM\\"
//Static StrConstant path_to_feather_folder = "Research/Perkins/Projects/PythonCommandLine/FEATHER/"
Static StrConstant HMM_file = "HMM_DE.py"

Static Function OutportTrace(TraceWave)
	wave TraceWave
	NewPath/Q/O Test "D:Devin:Python:HMM:"
	Save/O/J/P=Test TraceWave as "HMMData.txt"
	Killpath Test

end

Static Function RunHMMOutput(OptionsWave)
	wave OptionsWave
	String Location = "D:\Devin\Python\HMM\HMMData.txt"
	
	//HMMMain("C:/Users/dedwards/src_prh/",Location,OptionsWave=OptionsWave)
	HMMMain("C:/Users/dedwards/src_prh/",Location,OptionsWave=OptionsWave)

end


//
//Static StrConstant DEF_INPUT_REL_TO_BASE =  "IgorUtil/PythonApplications/FEATHER/Example/feather_example.pxp"
//
//Static Function Main_Windows()
//	// Runs a simple IWT on patrick's windows setup
//	ModMainFEATHER#Main("C:/Users/pahe3165/src_prh/")
//End Function 
//
//Static Function Main_Mac()
//	// Runs a simple IWT on patrick's mac setup 
//	ModMainFEATHER#Main("/Users/patrickheenan/src_prh/")
//End Function
//
Static Function HMMMain(base,input_file,[OptionsWave])
	// // This function shows how to use the IWT code
	// Args:
	//		base: the folder where the Research Git repository lives 
	//		input_file: the pxp to load. If not present, defaults to 
	//		<base>DEF_INPUT_REL_TO_BASE
	String base,input_file
	wave OptionsWave

//	//KillWaves /A/Z
//	ModPlotUtil#KillAllGraphs()
//	// IWT options
	Struct HMMOptions opt
//
//	// Add the meta information
	if (ParamIsDefault(OptionsWave))
//		opt.srate=1234
//		opt.tau=2
	else
		opt.srate=OptionsWave[0]
		opt.SState1=OptionsWave[1]
		opt.SState2=OptionsWave[2]

	EndIf
	
	Struct HMMOutput output
	//Make /O/N=0, output.event_starts
	opt.meta.path_to_input_file = input_file
	opt.meta.path_to_research_directory = base
	HmmRun(opt,output)
	LoadTheHMMOutput(opt)

End Function

Static Function HmmRun(user_options,output)
	// Function that calls the IWT python code
	//
	// Args:
	// 		options : instance of the InverseWeierstrassOptions struct. 
	//		output: instance of InverseWeierstrassOutput. Note that the waves 
	//		should already be allocated
	// Returns:
	//		Nothing, but sets output appropriately. 
	////
	Struct HMMOutput & output
	Struct HMMOptions & user_options
	// Manually set the main path and output file
	user_options.meta.path_to_main = full_path_to_HMM_main(user_options)
	user_options.meta.path_to_output_file = output_file_name(user_options)
	// make a local copy of user_options, since we have to mess with paths (ugh)
	// and we want to adhere to principle of least astonishment
	Struct HMMOptions options 
	options = user_options
	ModOperatingSystemUtil#get_updated_options(options.meta)
	// Run the python code 
	String PythonCommand = DE_HMMPython#python_command(options)	
	ModOperatingSystemUtil#execute_python(PythonCommand)
	//Make /O/FREE/Wave wave_tmp = {output.event_starts}
	//ModOperatingSystemUtil#get_output_waves(wave_tmp,options.meta,skip_lines=2)
End Function

Static Function /S python_command(opt)
//	// Function that turns an options structure into a python command string
//	//
//	// Args:
//	//		options: the InverseWeierstrassOptions structure, initialized (see inverse_weierstrass_options function)
//	// Returns:
//	//		string to command, OS-specific
	Struct HMMOptions & opt
	String PythonCommand
	String python_str  = ModOperatingSystemUtil#python_binary_string()
	String FolderPath = DE_HMMPython#full_path_to_HMM_folder(opt)
	String FullPath = DE_HMMPython#full_path_to_HMM_main(opt)
//	// Get just the python portion of the command
	String Output
	sprintf Output,"%s %s ",python_str,FullPath
	ModOperatingSystemUtil#append_argument(Output,"StartingRate",num2str(opt.srate))
	ModOperatingSystemUtil#append_argument(Output,"StartingState1",num2str(opt.SState1))
	ModOperatingSystemUtil#append_argument(Output,"StartingState2",num2str(opt.SState2))
//	ModOperatingSystemUtil#append_argument(Output,"dwell_time",num2str(opt.dwell_time))
//	ModOperatingSystemUtil#append_argument(Output,"trigger_time",num2str(opt.trigger_time))
	String output_file = output_file_name(opt)
	String input_file = opt.meta.path_to_input_file
	// Windows is a special flower and needs its paths adjusted
	if (running_windows())
		output_file = ModOperatingSystemUtil#sanitize_path_for_windows(output_file)
		input_file = ModOperatingSystemUtil#sanitize_path_for_windows(input_file)
	endif
	ModOperatingSystemUtil#append_argument(Output,"file_input",input_file)
	ModOperatingSystemUtil#append_argument(Output,"file_output",output_file,AddSpace=0)
	return Output
End Function

//Static StrConstant path_to_feather_folder = "D:/PatApp/AppFeather/AppPython/"
////Static StrConstant path_to_feather_folder = "Research/Perkins/Projects/PythonCommandLine/FEATHER/"
//Static StrConstant feather_file = "main_feather.py"
//
Structure HMMOutput
	Wave event_starts
EndStructure

Structure HMMOptions
	// Structure describing all the options the feather takes
	//
	//	Args: 
	//		threshold: the probability (between 0 and 1, log spacing is a good idea)
	//		tau: the number of points to use for spline fitting XXX currently not supported
	//
	//		path_to_research_directory: absolute path to directory one above Research (parent of 
	//		Research/Perkins/Projects/ ...)
	//
	//		path_to_input_file: where the pxp you want to analyze is. Should have a single wave 
	// 		like <x><d>_Sep and a single wave like <x><d>_Force. <x> can be anything,
	//  		<d> should be a 4-digit identifier (e.g. "Image0001_Sep" would be OK)
	variable srate
	Variable SState1
	Variable SState2
	Struct RuntimeMetaInfo meta

EndStructure
//

Static Function /S full_path_to_HMM_folder(options)
	// Function that gives the full path to the inverse weierstrass python folder
	//
	// Args:
	//		options: the FeatherOptions structure, initialized (see inverse_weierstrass_options function)
	// Returns:
	//		string to full path
	Struct HMMOptions & options
	return path_to_HMM_folder
	//return options.meta.path_to_research_directory + path_to_feather_folder
End Function
//
Static Function /S output_file_name(options)
//	//
//	// Returns:
//	//		the output file to read from
	Struct HMMOptions & options
	String FolderPath = DE_HMMPython#full_path_to_HMM_folder(options)
	return FolderPath//+ "events.csv"
End Function
//

//
Static Function /S full_path_to_HMM_main(options)
	// Function that gives the full path to the inverse weierstrass python folder
	//
	// Args:
	//		options: the InverseWeierstrassOptions structure, initialized (see inverse_weierstrass_options function)
	// Returns:
	//		string, full path to feather python main
	Struct HMMOptions & options
	return full_path_to_HMM_folder(options) + HMM_file
End Function
//
Static Function LoadTheHMMOutput(options)
	
	Struct HMMOptions & options
	String Folder=options.meta.path_to_output_file
	String TransMatrixOut=Folder+"Trans.txt"
	String  FitOut=Folder+"Fit.txt"
	String StateOut=Folder+"State.txt"
		String MeanOut=Folder+"Means.txt"

	LoadWave/O/Q/G/M/D/N=TransitionMatrix TransMatrixOut
	LoadWave/o/Q/G/N=FitWave FitOut
	LoadWave/o/Q/G/N=StateWave StateOut
		LoadWave/o/Q/G/N=MeanWave MeanOut

end

Static Function ButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			string saveDF 
			strswitch(ba.ctrlName)
				case "de_HMMPy_Button0":
					
					
					saveDF = GetDataFolder(1)
					controlinfo de_HMMPy_popup0
					SetDataFolder s_value
					controlinfo de_HMMPy_popup1
					wave w1=$S_value
					wave/T parmWave=root:DE_HMMPython:MenuStuff:ParmWave

					if(waveexists($(nameofwave(w1)+"_drift"))==0)
						DetermineDrift()
					else
					endif
					wave DriftWave=$(S_value+"_drift")
					variable filtering=str2num(parmWave[6][1])

					RunTheHMM(w1,DriftWave,filtering)


					duplicate/o $"FitWave0" $(S_value+"_Fit")
					duplicate/o $"StateWave0" $(S_value+"_States")
					duplicate/o $"TransitionMatrix0" $(S_value+"_Trans")
					duplicate/o $"MeanWave0" $(S_value+"_Mean")
					wave MeanWave= $(S_value+"_Mean")
					//print MeanWave
					duplicate/o $"Test1" $(S_value+"_Perfect")
					DE_HMM#PlotandProcessHMM($(S_value+"_States"),w1,$(S_value+"_Perfect"),S_value,MeanWave[0]*1e-9,MeanWave[1]*1e-9)
//
					killwaves $"FitWave0", $"StateWave0",$"TransitionMatrix0",$"test1"
					SetDataFolder saveDF
					break
				
				case "de_HMMPy_Button1":
					saveDF = GetDataFolder(1)
					controlinfo de_HMMPy_popup0
					SetDataFolder s_value
					controlinfo de_HMMPy_popup1
					wave w1=$S_value
					DetermineDrift()
					break
				case "de_HMMPy_Button2":
					
					
					saveDF = GetDataFolder(1)
					controlinfo de_HMMPy_popup0
					SetDataFolder s_value
					controlinfo de_HMMPy_popup1
					wave w1=$S_value
					wave/T parmWave=root:DE_HMMPython:MenuStuff:ParmWave

					if(waveexists($(nameofwave(w1)+"_drift"))==0)
						DetermineDrift()
					else
					endif
					wave DriftWave=$(S_value+"_drift")
					make/o/n=0 TestSeries_TVD,TestSeries_SVG
					  TestSeries_TVD[0]= {5e-10,1e-09,2e-09,5e-09,6e-09,7e-09,8e-09,9e-09,1e-08,1.2e-08,1.5e-08,2e-08,2.5e-08,3e-08}
 					 TestSeries_SVG[0]= {25,51,75,101,151,201,251,301,401,501,601,701,801,901,1001,1251,1501,1751,2001,2251,2501,2751,3001,3251,3501,3751,4001}
					make/o/c/n=0 TestResult_TVD,TestResult_SVG
					RunASeries(w1,driftwave,TestSeries_TVD,TestResult_TVD)
					RunASeries(w1,driftwave,TestSeries_SVG,TestResult_SVG)

					killwaves $"FitWave0", $"StateWave0",$"TransitionMatrix0",$"test1"
					SetDataFolder saveDF
					break
			endswitch
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Static Function Process(StateWave)

	wave StateWave
	wave RawWave=$Stringfromlist(0,nameofwave(StateWave),"_")+"_"+Stringfromlist(1,nameofwave(StateWave),"_")
	string S_value=replacestring("_States",nameofwave(StateWave),"")
	wave PerfectWave=$replacestring("States",nameofwave(StateWave),"Perfect")
	wave MeanWave=$replacestring("States",nameofwave(StateWave),"Mean")

	DE_HMM#PlotandProcessHMM(StateWave,RawWave,$(S_value+"_Perfect"),S_value,MeanWave[0]*1e-9,MeanWave[1]*1e-9)

end
Function RunASeries(w1,driftwave,filteringwave,ResWave)
	wave w1,driftwave,filteringwave
	wave/c ResWave
	string wavenamestring=nameofwave(w1)
	string Basename
	variable n=0
	variable filtering
	make/c/free/n=(numpnts(filteringwave)) FreeResults
	for(n=0;n<numpnts(filteringwave);n+=1)
	
		filtering=filteringwave[n]
		RunTheHMM(w1,DriftWave,filtering)
		wave trans=$"TransitionMatrix0"
		variable toprate=trans[0][1]/dimdelta($"Test1",0)
		variable botrate=trans[1][0]/dimdelta($"Test1",0)
		FreeResults[n]=cmplx(1/toprate,1/botrate)

		Basename=wavenamestring+"_"+num2str(filtering)
		duplicate/o $"FitWave0" $(Basename+"_Fit")
		duplicate/o $"StateWave0" $(Basename+"_States")
		duplicate/o $"TransitionMatrix0" $(Basename+"_Trans")
		duplicate/o $"MeanWave0" $(Basename+"_Mean")
		wave MeanWave= $(Basename+"_Mean")
	//	print MeanWave
		duplicate/o $"Test1" $(Basename+"_Perfect")
	endfor
	duplicate/o/c FreeResults ResWave

end

Static Function RunTheHMM(w1,DriftWave,filtering)
	Wave w1,DriftWave
	variable filtering
	String S_Value=nameofwave(w1)
	Interpolate2/T=2/N=(numpnts(w1))/Y=$(nameofwave(DriftWave)+"_INT") DriftWave
	wave DriftWaveInt=$(nameofwave(w1)+"_drift_INT")
	duplicate/o w1 $(S_value+"_BGC")
	wave BGC=$(S_value+"_BGC")
	BGC-=DriftWaveInt
	//killwaves DriftWaveInt
					
	duplicate/o BGC Test1
	if(filtering>=1)
		Smooth/S=2 Filtering,Test1
	else
		DE_Filtering#TVD1D_denoise(BGC,filtering,Test1)
	endif
					
	OutportTrace(Test1)
	make/free/n=0 options
	options={str2num(stringbykey("State1Rate",note(DriftWave),":","\r")),1e9*str2num(stringbykey("State1Loc",note(DriftWave),":","\r")),1e9*str2num(stringbykey("State2Loc",note(DriftWave),":","\r"))}
	RunHMMOutput(options)
	SetScale/P x pnt2x(Test1,0),dimdelta(Test1,0),"", $"FitWave0", $"StateWave0"
	wave trans=$"TransitionMatrix0"
	variable toprate=trans[0][1]/dimdelta(Test1,0)
	variable botrate=trans[1][0]/dimdelta(Test1,0)
	print filtering
	print 1/toprate
	print 1/botrate
				
//	duplicate/o $"FitWave0" $(S_value+"_Fit")
//	duplicate/o $"StateWave0" $(S_value+"_States")
//	duplicate/o $"TransitionMatrix0" $(S_value+"_Trans")
//	duplicate/o $"MeanWave0" $(S_value+"_Mean")
//	wave MeanWave= $(S_value+"_Mean")
//	duplicate/o Test1 $(S_value+"_Perfect")
//	DE_HMM#PlotandProcessHMM($(S_value+"_States"),w1,$(S_value+"_Perfect"),S_value,MeanWave[0]*1e-9,MeanWave[1]*1e-9)
//
//	killwaves $"FitWave0", $"StateWave0",$"TransitionMatrix0"
//	SetDataFolder saveDF


end

Static Function DetermineDrift()
	string saveDF = GetDataFolder(1)
	controlinfo de_HMMPy_popup0
	SetDataFolder s_value
	controlinfo de_HMMPy_popup1
	wave w1=$S_value
	duplicate/free w1 Test,Test1	
	wave/T parmWave=root:DE_HMMPython:MenuStuff:ParmWave
	variable filtering=str2num(parmWave[6][1])

	if(filtering>=1)
		Smooth/S=2 Filtering,Test1
	else
		DE_Filtering#TVD1D_denoise(Test,filtering,Test1)
	endif
			
	Test1*=1e9
	setscale/P x dimoffset(w1,0), dimdelta(w1,0), "s", Test1//ensuring scaling of input and output wave are the same
	Resample/DOWN=(str2num(parmWave[7][1]))/N=1/WINF=None Test1
	DE_HMM#DriftMarkovFitter( Test1, str2num(parmWave[0][1]), str2num(parmWave[1][1]), dimdelta(w1,0),str2num(parmWave[2][1]),str2num(parmWave[3][1]), str2num(parmWave[4][1]), str2num(parmWave[5][1]))
	wave HidMar0,HidMar1,HidMar2,HidMar3,HidMar4,HidMarParms0
	variable state0=HidMarParms0[1][0]/1e12
	variable state1=state0+HidMarParms0[1][1]/1e12
	HidMar1/=1e9
	//HidMar3/=1e12
	Test1/=1e9
	//state0+=mean(HidMar3)
	//state1+=mean(HidMar3)
						
	duplicate/o HidMar1 $(S_value+"_Driftfit")
	//duplicate/o HidMar2 $(S_value+"_st")
	duplicate/o HidMar3 $(S_value+"_drift")
	WAVE Driftwave=$(S_value+"_drift")
	Driftwave/=1e9
	//wave States=$(S_value+"_st")
	variable State1Loc=(HidMarParms0[1][0])*1e-9
	variable State2Loc=((HidMarParms0[1][0])+(HidMarParms0[1][1]))*1e-9
	variable State1Rate=(HidMarParms0[1][dimsize(HidMarParms0,1)-3])*dimdelta(Test1,0)
	variable State2Rate=(HidMarParms0[1][dimsize(HidMarParms0,1)-2])*dimdelta(Test1,0)
	string NoteString=""
	NoteString=ReplaceStringbyKey("State1Loc",NoteString,num2str(State1Loc),":","\r")
	NoteString=ReplaceStringbyKey("State2Loc",NoteString,num2str(State2Loc),":","\r")
	NoteString=ReplaceStringbyKey("State1Rate",NoteString,num2str(State1Rate),":","\r")
	NoteString=ReplaceStringbyKey("State2Rate",NoteString,num2str(State2Rate),":","\r")
	note/K Driftwave NoteString
	setscale/P x dimoffset(test1,0), dimdelta(test1,0), "s",  $(S_value+"_drift"),$(S_value+"_DriftFit")//ensuring scaling of input and output wave are the same
	duplicate/o Test1 $(S_value+"_HmmDriftFit")
	dowindow DriftFit
	if(V_Flag==1)
		killwindow DriftFit
	endif

	display/N=DriftFit $(S_value+"_HmmDriftFit")
	Appendtograph/W=DriftFit $(S_value+"_DriftFit")
	ModifyGraph/W=DriftFit rgb($(S_value+"_DriftFit"))=(0,0,0), lsize($(S_value+"_DriftFit"))=1.5

	if (DE_HMMPython#UserCursorAdjust("DriftFit",0) != 0)
		return -1
	endif
			
			
	killwaves HidMarParms0,HidMar0,HidMar1,HidMar2,HidMar3,HidMar4
	SetDataFolder saveDF
	
end

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
	DrawText 21,40,"Click Continue."
	Button button0,pos={80,58},size={92,20},title="Continue"
	Button button0,proc=DE_HMMPython#UserCursorAdjust_ContButtonProc

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

Window PythonHMMPanel() : Panel

	PauseUpdate; Silent 1		// building window...
	NewPanel/N=PythonHMMPanel /W=(697,267,1361,653)
	NewDataFolder/o root:DE_HMMPython
	NewDataFolder/o root:DE_HMMPython:MenuStuff

	DE_HMMPython#UpdateParmWave()
	Button de_HMMPy_button0,pos={250,110},size={150,20},proc=DE_HMMPython#ButtonProc,title="HMM!"
	Button de_HMMPy_button1,pos={250,140},size={150,20},proc=DE_HMMPython#ButtonProc,title="Recalc Drift!"
	Button de_HMMPy_button2,pos={250,170},size={150,20},proc=DE_HMMPython#ButtonProc,title="Run a Series!"

	PopupMenu de_HMMPy_popup0,pos={250,2},size={129,21}
	PopupMenu de_HMMPy_popup0,mode=1,popvalue="X",value= #"DE_PanelProgs#ListFolders()"

	PopupMenu de_HMMPy_popup1,pos={250,40},size={129,21}
	PopupMenu de_HMMPy_popup1,mode=1,popvalue="X",value= #"DE_HMMPython#ListWaves()"
	
	ListBox DE_HMMPy_list0,pos={400,2},size={175,150},proc=DE_HMMPython#ListBoxProc,listWave=root:DE_HMMPython:MenuStuff:ParmWave
	ListBox DE_HMMPy_list0,selWave=root:DE_HMMPython:MenuStuff:SelWave,editStyle= 2,userColumnResize= 1,widths={70,40,70,40}

EndMacro

Static Function/S ListWaves()

	String saveDF
	saveDF = GetDataFolder(1)
	controlinfo de_HMMPy_popup0
	SetDataFolder s_value
	String list = WaveList("*", ";", "")
	SetDataFolder saveDF
	return list

end

Static Function UpdateParmWave()
	if(exists("root:DE_HMMPython:MenuStuff:ParmWave")==1)
		wave/t/z Par=root:DE_HMMPython:MenuStuff:ParmWave
		wave/z Sel=root:DE_HMMPython:MenuStuff:SelWave
	Else
		make/t/n=(8,2) root:DE_HMMPython:MenuStuff:ParmWave
		wave/t/z Par=root:DE_HMMPython:MenuStuff:ParmWave
		make/n=(8,2) root:DE_HMMPython:MenuStuff:SelWave
		wave/z Sel=root:DE_HMMPython:MenuStuff:SelWave
		
		Par[0][0]={"Number of States","Number of Modes","Drift Bound (nm)","Sd. Deviation (nm)","Transition Bound","Iterations","Smoothing","Decimation"}
		Par[0][1]={"2","4",".5",".2",".5","3","10e-9","1"}
		Sel[][0]=0
		Sel[][1]=2
	endif


end

Menu "Equilibrium"
	//SubMenu "Processing"
	"Open Python HMM", PythonHMMPanel()

end
	