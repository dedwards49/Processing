#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma rtGlobals=3	

#include "c:\Users\dedwards\src_prh\IgorUtil\PythonApplications\InverseBoltzmann\Boltzmann"
#include "c:\Users\dedwards\src_prh\IgorUtil\IgorCode\Util\PlotUtil"
#pragma ModuleName = DE_InvBolt

Static StrConstant DEF_INPUT_REL_TO_BASE =  "IgorUtil/PythonApplications/InverseBoltzmann/Example/Experiment.pxp"

#include ":\Misc_PanelPrograms\Panel Progs"
#include ":\Misc_PanelPrograms\AsylumNaming"


Static Function OutportForce(Sinwave)
	wave Sinwave
	duplicate/o Sinwave 'Image0000Sep'
	wave SepWave='Image0000Sep'

	display/N=TMP_D SepWave 
	String Path="D:\Data\InverseBoltz\Data.pxp"
	SaveGraphCopy/o as Path
	KillWindow TMP_D
	killwaves SepWave

end

Static Function Main_Windows()
	// Runs a simple IWT on patrick's windows setup
	DE_InvBolt#Main("C:/Users/dedwards/src_prh/")
End Function 

Static Function Main_Mac()
	// Runs a simple IWT on patrick's mac setup 
	DE_InvBolt#Main("/Users/patrickheenan/src_prh/")
End Function

Static Function Main(base,[input_file,Pull])
	// // This function shows how to use the IWT code
	// Args:
	//		base: the folder where the Research Git repository lives 
	//		input_file: the pxp to load. If not present, defaults to 
	//		<base>DEF_INPUT_REL_TO_BASE
	String base,input_file,pull
	if (ParamIsDefault(input_file))
		input_file  = base +DEF_INPUT_REL_TO_BASE
	EndIf
	if (ParamIsDefault(Pull))
		Pull  = "No"
	EndIf
	//KillWaves /A/Z
	//ModPlotUtil#KillAllGraphs()
	// IWT options
	Struct BoltzmannOptions opt
	
	
	If(cmpstr(Pull,"Yes")==0)
		if(WaveExists(root:DE_IB:MenuStuff:ParmWave)==1)
			wave/T parmWave=root:DE_IB:MenuStuff:ParmWave

			opt.number_of_bins =str2num(ParmWave[0][1])
			opt.interpolation_factor = str2num(ParmWave[1][1])
      			opt.smart_interpolation =  str2num(ParmWave[2][1])
   		      opt.gaussian_stdev  =  str2num(ParmWave[3][1])*1e-9
      			opt.output_interpolated = str2num(ParmWave[4][1])
      			opt.n_iters = str2num(ParmWave[5][1])

		else
			opt.number_of_bins =500
			opt.interpolation_factor = 1
			opt.smart_interpolation = 1
			// Note: this stdev is from simulated data; 20nm would be
			// a huge (and awful) gaussian / point-spread function
			opt.gaussian_stdev  = .002e-8
			opt.output_interpolated = 1
			opt.n_iters = 1000
		endif
	
	else
	
		opt.number_of_bins =500
		opt.interpolation_factor = 1
		opt.smart_interpolation = 1
		// Note: this stdev is from simulated data; 20nm would be
		// a huge (and awful) gaussian / point-spread function
		opt.gaussian_stdev  = .002e-8
		opt.output_interpolated = 1
		opt.n_iters = 1000
		
	
	endif
	

      	// add the file information
	opt.meta.path_to_input_file = input_file
	opt.meta.path_to_research_directory = base
	// Make the output waves
	Struct BoltzmannOutput output
	Make /O/N=0,  output.extension_bins,output.distribution,output.distribution_deconvolved
	// Execte the command
	ModBoltzmann#inverse_boltzmann(opt,output)
	// plot the distributions and energy landscapes, units of kT
	Variable kT =  1
	Make /O/N=(DimSize(output.distribution,0)) landscape,landscape_deconvolved
	landscape[] = -ln(output.distribution[p]) * kT
	landscape_deconvolved[] = -ln(output.distribution_deconvolved[p]) * kT
	dowindow UtilDisp0
	if(V_flag==1)
	else
	ModPlotUtil#figure(hide=0)
	ModPlotUtil#subplot(2,1,1)
	ModPlotUtil#plot(output.distribution,mX=output.extension_bins)
	ModPlotUtil#plot(output.distribution_deconvolved,mX=output.extension_bins,color="r",linestyle="--")
	ModPlotUtil#pLegend(labelStr="Measured,Deconvolved")
	ModPlotUtil#xlabel("Extension (m)")
	ModPlotUtil#ylabel("Probability (1/m)")
	ModPlotUtil#subplot(2,1,2)
	ModPlotUtil#plot(landscape,mX=output.extension_bins)
	ModPlotUtil#plot(landscape_deconvolved,mX=output.extension_bins,color="r",linestyle="--")
	ModPlotUtil#xlabel("Extension (m)")
	ModPlotUtil#ylabel("Energy (kT)")
	endif

End Function


Static Function ButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
				string saveDF = GetDataFolder(1)
				controlinfo de_IB_popup0
				SetDataFolder s_value

				controlinfo de_IB_popup1
				wave w1=$S_value

				OutportForce(w1)

				Main("C:/Users/dedwards/src_prh/",input_file="D:\Data\InverseBoltz\Data.pxp",Pull="Yes")
				SetDataFolder saveDF

			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

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

Window InverseBoltzmannPanel() : Panel

	PauseUpdate; Silent 1		// building window...
	NewPanel/N=IBPanel /W=(697,267,1361,653)
	NewDataFolder/o root:DE_IB
	NewDataFolder/o root:DE_IB:MenuStuff

	DE_InvBolt#UpdateParmWave()
	Button de_IB_button0,pos={250,110},size={150,20},proc=DE_InvBolt#ButtonProc,title="Inverse that Boltz!"
	PopupMenu de_IB_popup0,pos={250,2},size={129,21}
	PopupMenu de_IB_popup0,mode=1,popvalue="X",value= #"DE_InvBolt#ListFolders()"
	PopupMenu de_IB_popup1,pos={250,40},size={129,21}
	PopupMenu de_IB_popup1,mode=1,popvalue="X",value= #"DE_InvBolt#ListWaves()"

	ListBox DE_IB_list0,pos={400,2},size={175,150},proc=DE_InvBolt#ListBoxProc,listWave=root:DE_IB:MenuStuff:ParmWave
	ListBox DE_IB_list0,selWave=root:DE_IB:MenuStuff:SelWave,editStyle= 2,userColumnResize= 1,widths={70,40,70,40}
	
	
//	SetVariable de_IWT_setvar0,pos={2,2},size={150,16},proc=DE_IWT#SVP,title="Initial Ramp Number"
//	SetVariable de_IWT_setvar0,limits={0,inf,1},value= _NUM:0
//	SetVariable de_IWT_setvar1,pos={2,25},size={150,16},proc=DE_IWT#SVP,title="Number of Ramps"
//	SetVariable de_IWT_setvar1,limits={0,inf,1},value= _NUM:10
//	
//	SetVariable de_IWT_setvar2,pos={250,140},size={100,16},proc=DE_IWT#SVP,title="Filtering"
//	SetVariable de_IWT_setvar2,limits={-inf,inf,2},value= _NUM:1
//	
//	Button de_IWT_button1,pos={2,50},size={80,20},proc=DE_IWT#ButtonProc,title="Stack Curves"
//
//	
//	SetDrawEnv fillpat= 0, linethick= 3.00;DelayUpdate
//	DrawRect -10,-10,160,160
//	SetDrawEnv fillpat= 0, linethick= 3.00;DelayUpdate
//	DrawRect 160,-10,600,160

EndMacro

Static Function/S ListWaves()

	String saveDF
	saveDF = GetDataFolder(1)
	controlinfo de_IB_popup0
	SetDataFolder s_value
	String list = WaveList("*", ";", "")
	SetDataFolder saveDF
	return list

end

Static Function/S ListFolders()

	string list=DE_PanelProgs#PrintAllFolders_String("*")
	return list
End

Static Function UpdateParmWave()
	if(exists("root:DE_IB:MenuStuff:ParmWave")==1)
		wave/t/z Par=root:DE_IB:MenuStuff:ParmWave
		wave/z Sel=root:DE_IB:MenuStuff:SelWave
	Else
		make/t/n=(6,2) root:DE_IB:MenuStuff:ParmWave
		wave/t/z Par=root:DE_IB:MenuStuff:ParmWave
		make/n=(6,2) root:DE_IB:MenuStuff:SelWave
		wave/z Sel=root:DE_IB:MenuStuff:SelWave
		
		Par[0][0]={"Number of Bins","Interpolation","Interp Factor","Sd. Deviation (nm)","Output Interp","Iterations"}
		Par[0][1]={"500","1","1",".2","1","300"}
		Sel[][0]=0
		Sel[][1]=2
	endif


end

Menu "Equilibrium"
	//SubMenu "Processing"
	"Open IB Decon", InverseBoltzmannPanel()


	//end
	
end