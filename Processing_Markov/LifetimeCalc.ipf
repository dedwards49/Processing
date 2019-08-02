#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma ModuleName = DE_HMM

#include ":\Misc_PanelPrograms\Panel Progs"
#include ":\Misc_PanelPrograms\AsylumNaming"

Static Function DriftMarkovFitter( UseWave, stateCount, modeCount, timeStep, driftBound, sigmaBound, transitionBound, iterationCount, [RAM, Threads])//Variables demanded by MarkovFit Code
	Wave UseWave//Input wave
	Variable stateCount, modeCount, timeStep, driftBound, sigmaBound, iterationCount, RAM, Threads
	Variable TransitionBound
	killwaves /z HidMar0, HidMar1, HidMar2, HidMar3, HidMar4, usable//Getting rid of generated waves to generate new ones
	RAM = paramIsDefault(RAM) ? 4:RAM
	Threads = paramIsDefault(Threads) ? 1000:Threads
	killwaves /z Used
	duplicate /o UseWave Used
	if(timeStep==0)
		timestep = 1.0
	endif
	Variable hold
	if(iterationCount==0)
		iterationCount = 4
	endif
	if(RAM == 0)
		RAM = 4
	endif
	if(modeCount ==0)
		Variable i
		for(i=0;i<numpnts(Used); i+=1)
			Used[i] += -driftBound*i
		endfor
	endif
	String InfoPass = "java -Xmx" + num2str(RAM) +"g -jar C:\MarkovFitter\DriftMarkov2.jar C:\MarkovFitter\UseWave.txt " + num2str(stateCount)+" 0 "//infopass exists to hold the command sent to DOS
	InfoPass = InfoPass + num2str(modeCount)+" "+num2str(timeStep)+" "+num2str(driftBound)+" "+num2str(sigmaBound)+" "+num2str(transitionBound)+" "+num2str(iterationCount)+" "+num2str(Threads)
	Save/J/W Used as "C:\MarkovFitter\UseWave.txt"//saving the wave that was given to  proper location
	print(InfoPass)//gives view of command line in case anything is wrong
	executescripttext InfoPass//sendng command to command line
	LoadWave/A=HidMar/J/D/W/K=0 "C:MarkovFitter:DriftMarkovOut.txt"//getting waves from location jar tosses them to(waves have base name HidMar
	//Display UseWave//displaying wave given
	variable Temp
	duplicate/o $"HidMar1" usable//while wave1 is created through this code it cannot regonize it so it must be duplicated
	Temp =dimoffset(UseWave,0)
	setscale/P x dimoffset(UseWave,0), dimdelta(UseWave,0), "s", usable//ensuring scaling of input and output wave are the same
	if(modeCount ==0)
		for(i=0;i<numpnts(UseWave);i+=1)
			usable[i] += driftBound*i
		endfor
	endif
	//AppendToGraph usable//putting on same graph
	//ModifyGraph rgb(usable)=(0,0,65280)//changing color so both waves are visible
	//display $"HidMar2"//displaying simple jump wave
	executescripttext "java -jar C:\MarkovFitter\GetRidOfUseWave.jar"//Eliminates file created earlier to prevent problems on future runs
end

Static Function ButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			string saveDF = GetDataFolder(1)
			controlinfo de_HMM_popup0
			SetDataFolder s_value
			controlinfo de_HMM_popup1
			wave w1=$S_value
			duplicate/o w1 Test
			wave test
			Test*=1e9
			setscale/P x dimoffset(w1,0), dimdelta(w1,0), "s", Test//ensuring scaling of input and output wave are the same
			DriftMarkovFitter( Test, 2, 3, dimdelta(w1,0), .5, .3, .2, 4)
			//killwaves test
			wave HidMar0,HidMar1,HidMar2,HidMar3,HidMar4
			
			HidMar1/=1e9
			duplicate/o/o HidMar1 $(S_value+"_fit")
			duplicate/o HidMar2 $(S_value+"_st")
			duplicate/o HidMar3 $(S_value+"_dr")
			setscale/P x dimoffset(w1,0), dimdelta(w1,0), "s", $(S_value+"_fit"),$(S_value+"_st"), $(S_value+"_dr")//ensuring scaling of input and output wave are the same
			DoWindow HMMPlot
			if(V_flag==0)
				display/n=HMMPlot w1
				Appendtograph/w=HMMPlot $(S_value+"_fit")
				ModifyGraph rgb($nameofwave(w1))=(58368,6656,7168);DelayUpdate
				ModifyGraph rgb($(S_value+"_fit"))=(0,0,0)
			else
				killwindow HMMPlot
				display/n=HMMPlot w1
				Appendtograph/w=HMMPlot $(S_value+"_fit")
				ModifyGraph rgb($nameofwave(w1))=(58368,6656,7168);DelayUpdate
				ModifyGraph rgb($(S_value+"_fit"))=(0,0,0)
			endif
			Killwaves HidMar0,HidMar1,HidMar2,HidMar3,HidMar4
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

Window HMMPanel() : Panel

	PauseUpdate; Silent 1		// building window...
	NewPanel/N=HMMPanel /W=(697,267,1361,653)
	NewDataFolder/o root:DE_HMM
	NewDataFolder/o root:DE_HMM:MenuStuff

	DE_HMM#UpdateParmWave()
	Button de_HMM_button0,pos={250,110},size={150,20},proc=DE_HMM#ButtonProc,title="HMM!"
	PopupMenu de_HMM_popup0,pos={250,2},size={129,21}
	PopupMenu de_HMM_popup0,mode=1,popvalue="X",value= #"DE_HMM#ListFolders()"
	PopupMenu de_HMM_popup1,pos={250,40},size={129,21}
	PopupMenu de_HMM_popup1,mode=1,popvalue="X",value= #"DE_HMM#ListWaves()"

	ListBox DE_HMM_list0,pos={400,2},size={175,150},proc=DE_HMM#ListBoxProc,listWave=root:DE_HMM:MenuStuff:ParmWave
	ListBox DE_HMM_list0,selWave=root:DE_HMM:MenuStuff:SelWave,editStyle= 2,userColumnResize= 1,widths={70,40,70,40}

EndMacro

Static Function/S ListWaves()

	String saveDF
	saveDF = GetDataFolder(1)
	controlinfo de_HMM_popup0
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
	if(exists("root:DE_HMM:MenuStuff:ParmWave")==1)
		wave/t/z Par=root:DE_HMM:MenuStuff:ParmWave
		wave/z Sel=root:DE_HMM:MenuStuff:SelWave
	Else
		make/t/n=(6,2) root:DE_HMM:MenuStuff:ParmWave
		wave/t/z Par=root:DE_HMM:MenuStuff:ParmWave
		make/n=(6,2) root:DE_HMM:MenuStuff:SelWave
		wave/z Sel=root:DE_HMM:MenuStuff:SelWave
		
		Par[0][0]={"Number of Bins","Interpolation","Interp Factor","Sd. Deviation (nm)","Output Interp","Iterations"}
		Par[0][1]={"500","1","1",".2","1","300"}
		Sel[][0]=0
		Sel[][1]=2
	endif


end

Menu "Equilibrium"
	//SubMenu "Processing"
	"Open HMM", HMMPanel()


	//end
	
end
//
//
//function SplitLifetimes(States,UpperLT,LowerLT)
//wave States,UpperLT,LowerLT
//
//variable n=1
//variable un=0,lun=0
//variable Current
//if(States[0]==1)
//Current=1
//elseif(States[0]==0)
//Current=-1
//else
//return -1
//endif
//make/o/n=0 $nameofwave(UpperLT),$nameofwave(LowerLT)
//variable last=pnt2x(States,0)
//for(n=1;n<numpnts(States);n+=1)
//if(States[n]!=States[n-1])
//
//
//if(Current==1)
//InsertPoints un,1,UpperLT
//UpperLT[un]=(pnt2x(States,n)-last)
//last=pnt2x(States,n)
//un+=1
//else
//InsertPoints un,1,LowerLT
//
//LowerLT[lun]=(pnt2x(States,n)-last)
//lun+=1
//last=pnt2x(States,n)
//endif
//Current*=-1
//endif
//
//endfor
//
//end