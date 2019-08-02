#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma ModuleName = DE_HMM
#include "DE_Filtering"

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
	LoadWave/M/A=HidMarParms/J/D/K=1 "C:MarkovFitter:DriftMarkovProperties.txt"
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
	killwaves usable, used
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
			duplicate/free w1 Test,Test1	
			wave/T parmWave=root:DE_HMM:MenuStuff:ParmWave
			variable filtering=str2num(parmWave[6][1])
			if(filtering>=1)
	
				Smooth/S=2 Filtering,Test1

			else
				DE_Filtering#TVD1D_denoise(Test,filtering,Test1)
	
			endif
			
			Test1*=1e12
			setscale/P x dimoffset(w1,0), dimdelta(w1,0), "s", Test1//ensuring scaling of input and output wave are the same
			Resample/DOWN=(str2num(parmWave[7][1]))/N=1/WINF=None Test1
			DriftMarkovFitter( Test1, str2num(parmWave[0][1]), str2num(parmWave[1][1]), dimdelta(w1,0),str2num(parmWave[2][1]),str2num(parmWave[3][1]), str2num(parmWave[4][1]), str2num(parmWave[5][1]))
			wave HidMar0,HidMar1,HidMar2,HidMar3,HidMar4,HidMarParms0
			variable state0=HidMarParms0[1][0]/1e12
			variable state1=state0+HidMarParms0[1][1]/1e12
			HidMar1/=1e12
			HidMar3/=1e12
			Test1/=1e12
			state0+=mean(HidMar3)
			state1+=mean(HidMar3)
			
			duplicate/o HidMar1 $(S_value+"_fit")
			duplicate/o HidMar2 $(S_value+"_st")
			duplicate/o HidMar3 $(S_value+"_dr")

			wave States=$(S_value+"_st")
			setscale/P x dimoffset(test1,0), dimdelta(test1,0), "s", $(S_value+"_fit"),$(S_value+"_st"), $(S_value+"_dr")//ensuring scaling of input and output wave are the same
			duplicate/o Test1 $(S_value+"_Hmm")
			wave HmmIN=$((S_value+"_Hmm")[0,30])
			PlotandProcessHMM($(S_value+"_st"),w1,HmmIN,S_value,state0,state1)

			Killwaves HidMar0,HidMar1,HidMar2,HidMar3,HidMar4,HidMarParms0
			
			SetDataFolder saveDF

			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Static Function PlotandProcessHMM(StateWave,RawForceWave,ForceWave,NameString,state1,state2)
	wave statewave,ForceWave,RawForceWave

	variable state1,state2

	string NameString
	
	make/o/n=0 $(NameString+"_HMM_ULT"),$(NameString+"_HMM_LLT")
	wave fitwave=$(NameString+"_fit")
	wave UnfoldedLifetime=$(NameString+"_HMM_ULT")
	
	wave FoldedLifetime=$(NameString+"_HMM_LLT")
	CalcLifetimes(statewave,UnfoldedLifetime,FoldedLifetime)

	make/free/n=0 OutUHist,OutLHist
	variable totallength=dimdelta(statewave,0)*dimsize(statewave,0)
	ReturnStateLifetimes(statewave,totallength/2000,OutUHist,OutLHist)
	make/o/n=40 $(nameofwave(ForceWave)+"_Hist")
	wave FHist=$(nameofwave(ForceWave)+"_Hist")
	Histogram/C/B=1 ForceWave,FHist

	wavestats/q/R=[0,numpnts(FHist)/2] FHist
	variable P1=state1
	variable H1=v_max
	wavestats/q/R=[numpnts(FHist)/2,numpnts(FHist)-1] FHist
	variable P2=state2
	variable H2=v_max
	make/D/o/n=7 w_coefs
	w_coefs[0]={0,H1,H2,P1,P2,.2e-9,.2e-9}
	FuncFit/Q/W=2/H="1000000"/N/NTHR=0 DGauss w_coefs  FHist/D
	//wave WFIT=$(("fit_"+nameofwave(FHist))[0,30])
	wave WFIT=$(("fit_"+nameofwave(FHist)))
	variable totaltransition=numpnts(FoldedLifetime)+numpnts(unFoldedLifetime)
	variable FoldedLTAvg=mean(FoldedLifetime)
	variable UnfoldedLTAvg= mean(UnFoldedLifetime)
	make/o/n=(20) $(nameofwave(RawForceWave)+"_Hmm_FH")
	wave FLTHist=$(nameofwave(RawForceWave)+"_Hmm_FH")
	Histogram/C/B={1e-3,(dimdelta(RawForceWave,0)*numpnts(RawForceWave)/totaltransition/3),20} FoldedLifetime FLTHist
	make/o/n=(20) $(nameofwave(RawForceWave)+"_HMM_UH")
	wave UFLTHist=$(nameofwave(RawForceWave)+"_HMM_UH")
	Histogram/C/B={1e-3,(dimdelta(RawForceWave,0)*numpnts(RawForceWave)/totaltransition)/3,20} UnFoldedLifetime UFLTHist
	make/D/free/n=3 w_coefL1
	w_coefL1[0]={0,100,1e-2}
	CurveFit/W=2/Q/G/H="100"/NTHR=0 exp_XOffset kwCWave=w_coefL1 ,FLTHist /D
	wave FLFit= $(("fit_"+nameofwave(FLTHist))[0,30])	
	make/D/free/n=3 w_coefL2
	w_coefL2[0]={0,100,1e-2}
	CurveFit/W=2/Q/G/H="100"/NTHR=0 exp_XOffset kwCWave=w_coefL2,UFLTHist /D 
	wave UFLFit= $(("fit_"+nameofwave(UFLTHist))[0,30])	
	make/o/n=(3,2) $(nameofwave(RawForceWave)+"_HMM_S1"),$(nameofwave(RawForceWave)+"_HMM_S2")
	wave State1Pos=$(nameofwave(RawForceWave)+"_HMM_S1")
	wave State2pos=$(nameofwave(RawForceWave)+"_HMM_S2")
	State1Pos[0][0]=State1-.05e-9
	State1Pos[0][1]=0
	State1Pos[1][0]=State1
	State1Pos[1][1]=wavemax(FHist)
	State1Pos[2][0]=State1+.05e-9
	State1Pos[2][1]=0
	State2Pos[0][0]=State2-.05e-9
	State2Pos[0][1]=0
	State2Pos[1][0]=State2
	State2Pos[1][1]=wavemax(FHist)
	State2Pos[2][0]=State2+.05e-9
	State2Pos[2][1]=0

	DE_HMM#MakeNicePlot(RawForceWave,ForceWave,fitwave,FHist,WFIT,FLTHist,FLFit,UFLTHist,UFLFit,State1Pos,State2Pos)

	TextBox/N=Populations/X=35/Y=1/C/N=text0/F=0 num2str(round(100*w_coefs[1]/(w_coefs[1]+w_coefs[2])))+"%"
	string lifetime1S,lifetime2S,lifetime1AS,lifetime2AS
	sprintf lifetime1S, "%0.2f",-1e3*w_coefL1[2]*ln(1/2)
	sprintf lifetime2S, "%0.2f",-1e3*w_coefL2[2]*ln(1/2)
	sprintf lifetime1AS, "%0.2f",1e3*FoldedLTAvg
	sprintf lifetime2AS, "%0.2f",1e3*UnfoldedLTAvg		
	TextBox/N=Lifetimes/X=0/Y=1/C/N=text0/F=0 "\\K(19712,44800,18944)Folded: "+lifetime1AS+"("+lifetime1S+") ms\r\\K(14848,32256,47104)UnFolded: "+lifetime2AS+"("+lifetime2S+") ms"
	wave W_Coef,W_Sigma,w_coefs,W_fitConstants
	killwaves W_Coef,W_Sigma,w_coefs,W_fitConstants				
					
	//killwaves FHist
end


Static Function/C ReturnStateLifetimes(HMM_FIt,HistStepO,OutUHist,OutLHist)
	Wave HMM_Fit,OutUHist,OutLHist
	variable HistStepO
	make/free/n=0 UOut,LOut
	CalcLifetimes(HMM_FIT,UOut,LOut)
	variable HistStepL=HistStepO
	variable HistStepU=HistStepO
	variable HisStepC
	make/free/n=40 UOutHist,LOutHist

	
//	do
		HisStepC=HistStepL
		make/free/n=40 LOutHist
		Histogram/C/B={0,HistStepL,40} LOut,LOutHist
		make/o/n=3 W_coef
		W_Coef={0,LoutHist[0],1/HistStepL}
		CurveFit/Q/W=2/N/H="100"/NTHR=0 exp  LOutHist
//		if(1/W_coef[2]>20*HistStepL)
//			HistStepL*=5
		
//		elseif(1/W_coef[2]<HistStepL*4)
//			HistStepL/=3
		
//		endif
	//while(1/W_coef[2]>20*HisStepC||1/W_coef[2]<HisStepC*4)
	
	variable/C Result=Cmplx(1/w_coef[2],0)

//	do
			HisStepC=HistStepU
		make/free/n=40 UOutHist
		Histogram/C/B={0,HistStepU,40} UOut,UOutHist
		W_Coef={0,UoutHist[0],1/HistStepU}
		CurveFit/Q/W=2/N/H="100"/NTHR=0 exp  UOutHist
//		print 1/W_coef[2]
//		if(1/W_coef[2]>20*HistStepU)
//			HistStepU*=5
//		
//		elseif(1/W_coef[2]<HistStepU*4)
//			HistStepU/=3
//
		
//		endif
//	while(1/W_coef[2]>20*HisStepC||1/W_coef[2]<HisStepC*4)
	Result+=cmplx(0,1/w_coef[2])
	duplicate/o UOutHist OutUHist
	duplicate/o LOutHist OutLHist
	return result
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
	PopupMenu de_HMM_popup0,mode=1,popvalue="X",value= #"DE_PanelProgs#ListFolders()"
	

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


//Function DriftMarkovFitter( UseWave, stateCount, modeCount, timeStep, driftBound, sigmaBound, transitionBound, iterationCount, [RAM, Threads])
Static Function UpdateParmWave()
	if(exists("root:DE_HMM:MenuStuff:ParmWave")==1)
		wave/t/z Par=root:DE_HMM:MenuStuff:ParmWave
		wave/z Sel=root:DE_HMM:MenuStuff:SelWave
	Else
		make/t/n=(8,2) root:DE_HMM:MenuStuff:ParmWave
		wave/t/z Par=root:DE_HMM:MenuStuff:ParmWave
		make/n=(8,2) root:DE_HMM:MenuStuff:SelWave
		wave/z Sel=root:DE_HMM:MenuStuff:SelWave
		
		Par[0][0]={"Number of States","Number of Modes","Drift Bound (nm)","Sd. Deviation (nm)","Transition Bound","Iterations","Smoothing","Decimation"}
		Par[0][1]={"2","4",".5",".2",".5","3","10e-9","1"}
		Sel[][0]=0
		Sel[][1]=2
	endif


end

Menu "Equilibrium"
	//SubMenu "Processing"
	"Open HMM", HMMPanel()


	//end
	
end
Static Function MakeNicePlot(RawForce,SmForce,ForceFit,HistWave,HistFit,LT1,LT1Fit,LT2,LT2Fit,state1,state2)
	wave RawForce,SmForce,HistWave,HistFit,LT1,LT2,LT1Fit,LT2Fit,ForceFit, state1,state2
	string WindowName=nameofwave(RawForce)+"_HMM_Win"
	dowindow $Windowname
	if(V_flag==1)
		killwindow $windowname
	else
	endif
	Display/N=$WindowName RawForce,SmForce,ForceFit

	AppendToGraph/W=$WindowName/B=B1/L=L1/VERT HistWave
	AppendToGraph/W=$WindowName/B=B1/L=L1/VERT HistFit
	AppendToGraph/W=$WindowName/B=B1/L=L1/Vert State1[][1] vs State1[][0]
	AppendToGraph/W=$WindowName/B=B1/L=L1/Vert State2[][1] vs State2[][0]

	AppendToGraph/W=$WindowName/B=B2/L=L2 LT1
	AppendToGraph/W=$WindowName/B=B2/L=L2 LT1Fit

	AppendToGraph/W=$WindowName/B=B2/L=L2 LT2
	AppendToGraph/W=$WindowName/B=B2/L=L2 LT2Fit

	//	//SetAxis L1 2.1312133e-08,2.9592798e-08
	ModifyGraph/W=$WindowName mode($nameofwave(State1))=5,hbFill($nameofwave(State1))=4,rgb($nameofwave(State1))=(0,65280,65280);DelayUpdate
	ModifyGraph/W=$WindowName mode($nameofwave(State2))=5,hbFill($nameofwave(State2))=4,rgb($nameofwave(State2))=(0,65280,65280)

	ModifyGraph/W=$WindowName tick=2,fSize=9,lblPosMode=1,lblPos=42,standoff=0,font="Arial"
	ModifyGraph/W=$WindowName axisEnab(bottom)={0,0.5}
	ModifyGraph/W=$WindowName axisEnab(B1)={0.55,0.7}
	ModifyGraph/W=$WindowName axisEnab(B2)={0.75,1}
	ModifyGraph/W=$WindowName freePos(B1)={0,L1}
	ModifyGraph/W=$WindowName freePos(L1)={0,B1}
	ModifyGraph/W=$WindowName freePos(B2)={0,L2}
	ModifyGraph/W=$WindowName freePos(L2)={0,B2}
	ModifyGraph/W=$WindowName rgb($nameofwave(HistFit))=(0,0,0)
	ModifyGraph/W=$WindowName rgb($(nameofwave(RawForce)))=(65280,48896,48896)
	ModifyGraph/W=$WindowName hideTrace($(nameofwave(RawForce)))=1
	ModifyGraph/W=$WindowName rgb($(nameofwave(ForceFit)))=(0,0,0)
	ModifyGraph/W=$WindowName margin(left)=36,margin(bottom)=29,margin(top)=14,margin(right)=14;DelayUpdate
	ModifyGraph/W=$WindowName mode($nameofwave(LT1))=3,marker($nameofwave(LT1))=16;DelayUpdate
	ModifyGraph/W=$WindowName rgb($nameofwave(LT1))=(19712,44800,18944);DelayUpdate
	ModifyGraph/W=$WindowName useMrkStrokeRGB($nameofwave(LT1))=1,mode($nameofwave(LT2))=3;DelayUpdate
	ModifyGraph/W=$WindowName marker($nameofwave(LT2))=16;DelayUpdate
	ModifyGraph/W=$WindowName rgb($nameofwave(LT2))=(14848,32256,47104);DelayUpdate
	ModifyGraph/W=$WindowName useMrkStrokeRGB($nameofwave(LT2))=1
	ModifyGraph/W=$WindowName width=576,height=144
	ModifyGraph/W=$WindowName noLabel(L1)=2
	ModifyGraph/W=$WindowName tickUnit(left)=1,prescaleExp(left)=9;DelayUpdate
	Label/W=$WindowName left "\\f01Extension (nm)"
	Label/W=$WindowName bottom "\\f01Time (s)"
	ModifyGraph/W=$WindowName mode($nameofwave(HistWave))=3,marker($nameofwave(HistWave))=19;DelayUpdate
	ModifyGraph/W=$WindowName rgb($nameofwave(HistWave))=(58368,6656,7168);DelayUpdate
	ModifyGraph/W=$WindowName useMrkStrokeRGB($nameofwave(HistWave))=1
	ModifyGraph/W=$WindowName lsize($Nameofwave(HistFit))=2
	
	
	ModifyGraph/W=$WindowName lsize($nameofwave(LT1Fit))=1.5;DelayUpdate
	ModifyGraph/W=$WindowName rgb($nameofwave(LT1Fit))=(19712,44800,18944);DelayUpdate
	ModifyGraph/W=$WindowName lsize($nameofwave(LT2Fit))=1.5;DelayUpdate
	ModifyGraph/W=$WindowName rgb($nameofwave(LT2Fit))=(14848,32256,47104)
	
	
	DoUpdate 
	GetAxis/W=$WindowName/Q L1
	SetAxis/W=$WindowName left, v_min,v_max
end

Static Function CalcLifetimes(States,UOut,LOut)

	Wave States,UOut,LOut
	variable n
	variable StartPoint=0
	make/free/n=0 UpperLifetime,LowerLifetime,transitionpnts
	for(n=1;n<numpnts(States);n+=1)
		if(States[n]==States[n-1])
		else
			
			
			if(States[n-1]==1)
			insertpoints 0,1,UpperLifetime
			UpperLifetime[0]=((n-1-StartPoint))
			StartPoint=n
			else
			insertpoints 0,1,LowerLifetime
			LowerLifetime[0]=((n-1-StartPoint))
			StartPoint=n
			endif
		endif

	endfor
	UpperLifetime*=dimdelta(States,0)
	LowerLifetime*=dimdelta(States,0)
//	string UName= nameofwave(states)[0,strlen(nameofwave(States))-3]+"UL"
//	string LName= nameofwave(states)[0,strlen(nameofwave(States))-3]+"LL"
	duplicate/o UpperLifetime UOut
	duplicate/o LowerLifetime LOut
	//duplicate/o transitionpnts Trxs
end
