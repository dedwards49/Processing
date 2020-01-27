#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
//#include "C:\Users\dedwards\src_prh\IgorUtil\PythonApplications\FEATHER\Example\MainFeather"
#include "C:\Devin\Documents\Software\AppFEATHER\AppIgor\Example\MainFeather"
#pragma modulename=DE_NewFeather
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
	killwaves Timewave

end

Static Function RunFeatheronOutput(OptionsWave)
	wave OptionsWave
	String Location = "C:\Data\Feather\Hold.pxp"
	
	///ModMainFEATHER#Main(base="C:/Users/dedwards/src_prh/",Input_file=Location,OptionsWave=OptionsWave)

	ModMainFEATHER#Main(base="C:/Devin/Documents/Software/AppFEATHER/",Input_file=Location,OptionsWave=OptionsWave)
end

Static Function RunFeatheronOutputFolder(OptionsWave)
	wave OptionsWave
	String LocationFolder = "C:\Data\Feather"
	
	///ModMainFEATHER#Main(base="C:/Users/dedwards/src_prh/",Input_file=Location,OptionsWave=OptionsWave)

	RunMultiFeather("C:/Devin/Documents/Software/AppFEATHER/",LocationFolder,OptionsWave)
end

Static Function MakeSingleWaves(SearchString)
	STRING SearchString
	string AllForceRet= wavelist(SearchString,";","")
	String ForceWaveList="",SepWaveList=""
	variable n
	string ForceWave, SepWave
	for(n=0;n<itemsinlist(AllForceRet);n+=1)
		//for(n=0;n<itemsinlist(AllFOrceRet);n+=1)
		wave ForceRetWave=$stringfromlist(n,AllForceRet)
		wave ForceExtWave=$replacestring("Ret",nameofwave(ForceRetWave),"Ext")
		wave SepRetWave=$replacestring("Force",nameofwave(ForceRetWave),"Sep")
		wave SepExtWave=$replacestring("Force",nameofwave(ForceExtWave),"Sep")
		
		make/free/n=0 ForceAll,SepAll
		Concatenate/o {ForceExtWave,ForceRetWave},ForceAll
		Concatenate/o {SepExtWave,SepRetWave},SepAll
		duplicate/o ForceAll $replacestring("Force_Ret",nameofwave(ForceRetWave),"Force")
		duplicate/o SepAll $replacestring("Force_ret",nameofwave(ForceRetWave),"Sep")

	endfor
	killwaves ForceAll,SepAll
end



Static Function SaveOutAllWaves(SearchString)
	String SearchString
	string AllForce= wavelist(SearchString,";","")
	CleanOutFolder("C:Data:Feather:")

	variable n,tot=itemsinlist(AllForce)
	display/N=TMP_D
	for(n=0;n<tot;n+=1)
		Wave ForceWave=$stringfromlist(n,AllForce)
		if(waveexists($ReplaceString("Force",nameofwave(ForceWave),"Sep"))==0)
		else
		Wave SepWave=$ReplaceString("Force",nameofwave(ForceWave),"Sep")
		duplicate/o ForceWave $(replaceString("Force",nameofwave(ForceWave),"Time"))
		wave TimeWave=$(replaceString("Force",nameofwave(ForceWave),"Time"))
		Appendtograph/W=TMP_D TimeWave
		Appendtograph/W=TMP_D ForceWave vs SepWave
		TimeWave=pnt2x(ForceWave,p)
		endif
	endfor
	
	


	String Path="C:\Data\Feather\Hold.pxp"
	SaveGraphCopy/o as Path
	KillWindow TMP_D
	for(n=0;n<tot;n+=1)
		Wave ForceWave=$stringfromlist(n,AllForce)
		if(waveexists($ReplaceString("Force",nameofwave(ForceWave),"Sep"))==0)
		else
		wave TimeWave=$(replaceString("Force",nameofwave(ForceWave),"Time"))

		TimeWave=pnt2x(ForceWave,p)
		killwaves timewave
		endif
	endfor
end

Static Function LoadTheWaves(SearchString)
	String SearchString
	string AllForce= wavelist(SearchString,";","")
	NewPath/o/Q DataPath "C:Data:Feather:"
	string AllFiles= IndexedFile(DataPath, -1, ".txt")
	string CurrentFile,CurrentName
	variable n,tot=itemsInList(AllFiles)

	for(n=0;n<tot;n+=1)
		CurrentFile=stringfromlist(n,AllFiles)
		CurrentName=CurrentFile[0,strlen(CurrentFile)-5]
		if(FindListItem(Currentname,AllForce)!=-1)
			CurrentName=replaceString("Force",CurrentName,"Starts")
			Wave ForceWave=$replacestring("Starts",CurrentName,"Force_ret")
			Wave SepWave=$replacestring("Starts",CurrentName,"Sep_ret")

			LoadWave/Q/G/P=DataPath/N=GARBAGE/L={0, 1, 0, 0, 0 } CurrentFile
			
			wave garbage0
			if(Garbage0[0]==0)
				if(numpnts(Garbage0)==1)
					Garbage0[0]=0
				else
					deletepoints 0,1, Garbage0
				endif
			elseif(Garbage0[0]<0)
				String NoteIn=replacestringbykey("DE_FeatherZero",note(ForceWave),num2str(Garbage0[0]),":","\r")
				note/K ForceWave NoteIn
				note/K SepWave NoteIn
				if(numpnts(Garbage0)==1)
					Garbage0[0]=0
				else
					deletepoints 0,1, Garbage0
				endif
			endif
			
			duplicate/o GARBAGE0 $CurrentName
		else
		endif
	endfor
	killpath Datapath


end

Static Function CleanOutFolder(PathString)
	String PathString

	NewPath/o/Q DataPath PathString
	string AllFiles= IndexedFile(DataPath, -1,"????")
	string CurrentFile,CurrentName
	variable n,tot=itemsInList(AllFiles)

	for(n=0;n<tot;n+=1)
		CurrentFile=stringfromlist(n,AllFiles)
		//print CurrentFile
		DeleteFile/P=DataPath CurrentFile
	endfor
	killpath Datapath
end

Static Function InitMultiFeather(user_options)
	// Function that calls the IWT python code
	//
	// Args:
	// 		options : instance of the InverseWeierstrassOptions struct. 
	//		output: instance of InverseWeierstrassOutput. Note that the waves 
	//		should already be allocated
	// Returns:
	//		Nothing, but sets output appropriately. 
	//
	Struct FeatherOptions & user_options
	// Manually set the main path and output file
	user_options.meta.path_to_main = full_path_to_feather_Devin(user_options)
	user_options.meta.path_to_output_file = user_options.meta.path_to_input_file
	// make a local copy of user_options, since we have to mess with paths (ugh)
	// and we want to adhere to principle of least astonishment
	Struct FeatherOptions options 
	options = user_options
	//ModOperatingSystemUtil#get_updated_options(options.meta)

	// Run the python code 
	String PythonCommand = Devin_python_command(options)	
	ModOperatingSystemUtil#execute_python(PythonCommand,options.meta)
	//Make /O/FREE/Wave wave_tmp = {output.event_starts}
	//ModOperatingSystemUtil#get_output_waves(wave_tmp,options.meta,skip_lines=2)
End Function

Static Function /S full_path_to_feather_Devin(options)
	// Function that gives the full path to the inverse weierstrass python folder
	//
	// Args:
	//		options: the InverseWeierstrassOptions structure, initialized (see inverse_weierstrass_options function)
	// Returns:
	//		string, full path to feather python main
	Struct FeatherOptions & options
	return Full_path_to_PythonIgor_Feather(options) + "main_devin.py"
End Function

Static Function /S Full_path_to_PythonIgor_Feather(options)
	// Function that gives the full path to the inverse weierstrass python folder
	//
	// Args:
	//		options: the InverseWeierstrassOptions structure, initialized (see inverse_weierstrass_options function)
	// Returns:
	//		string, full path to feather python main
	Struct FeatherOptions & options
	return  "C:/Devin/Documents/Software/Library/Python/Feather/"
End Function


Static Function /S Devin_python_command(opt)
	// Function that turns an options structure into a python command string
	//
	// Args:
	//		options: the InverseWeierstrassOptions structure, initialized (see inverse_weierstrass_options function)
	// Returns:
	//		string to command, OS-specific
	Struct FeatherOptions & opt
	String PythonCommand
	String python_str  = ModOperatingSystemUtil#python_binary_string(opt.meta)
	String FolderPath = ModFeather#full_path_to_feather_folder(opt)
	String FullPath = full_path_to_feather_Devin(opt)
//	// Get just the python portion of the command
	String Output
	sprintf Output,"%s %s ",python_str,FullPath
	ModOperatingSystemUtil#append_argument(Output,"threshold",num2str(opt.threshold))
	ModOperatingSystemUtil#append_argument(Output,"tau",num2str(opt.tau))
	String output_file = opt.meta.path_to_output_file
	String input_folder = opt.meta.path_to_input_file
//	// Windows is a special flower and needs its paths adjusted
//	if (running_windows())
//		output_file = ModOperatingSystemUtil#sanitize_path_for_windows(output_file)
//		input_file = ModOperatingSystemUtil#sanitize_path_for_windows(input_file)
//	endif
	ModOperatingSystemUtil#append_argument(Output,"input_directory",input_folder)
	ModOperatingSystemUtil#append_argument(Output,"folder_output",output_file,AddSpace=0)
	return Output
End Function


Static Function RunMultiFeather(base,input_file,OptionsWave)
	// // This function shows how to use the IWT code
	// Args:
	//		base: the folder where the Research Git repository lives 
	//		input_file: the pxp to load. If not present, defaults to 
	//		<base>DEF_INPUT_REL_TO_BASE
	String base,input_file
	wave OptionsWave
//	if (ParamIsDefault(base))
//		base = ModIoUtil#pwd_igor_path(DEF_PATH_NAME,n_up_relative=3)
//	EndIf
//	if (ParamIsDefault(input_file))
//		input_file  = base +DEF_INPUT_REL_TO_BASE
//	EndIf

//	
	Struct FeatherOptions opt
//
//	
//	if (ParamIsDefault(OptionsWave))
//		opt.threshold = 1e-3
//		opt.tau = 1e-3
//		opt.trigger_time = 1.05
//		opt.dwell_time = .01
//		opt.spring_constant = 5.1e-3
//	else
		opt.threshold = OptionsWave[0]
		opt.tau = OptionsWave[1]
//		opt.trigger_time = OptionsWave[2]
//		opt.dwell_time = OptionsWave[3]
//		opt.spring_constant = OptionsWave[4]
//	EndIf
//	// add the file information
	opt.meta.path_to_input_file = input_file
	opt.meta.path_to_research_directory = base
//	// Make the output waves
//	Struct FeatherOutput output
//	Make /O/N=0, output.event_starts
//	// Execte the command
	InitMultiFeather(opt)
//	// Make a fun plot wooo
//	//LoadData /Q/O/R (ModOperatingSystemUtil#sanitize_path(input_file))
//	//Wave Y =  $("Image0994Force")
//	//Wave X =  $("Image0994Sep")
//	//Display Y vs X
//	Variable n_events = DimSize(output.event_starts,0)
//	Variable i
//	//for (i=0; i<n_events; i+=1)
//	//	Variable event_idx = output.event_starts[i]
//	//	ModPlotUtil#axvline(X[event_idx])
//	//endfor
//	//Edit output.event_starts as "Predicted Event Indices in Wave"
//	
	
	
End Function