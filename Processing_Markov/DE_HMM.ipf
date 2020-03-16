#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma ModuleName = DE_HMM
#include "DE_Filtering"
#include "EquilibriumJumpData"
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
	
	print InfoPass
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
			controlinfo/W=HmmPanel de_HMM_popup0
			String BaseFoldStr=S_value
			SetDataFolder BaseFoldStr
			controlinfo/W=HmmPanel de_HMM_popup1
			wave w1=$S_value
			string BaseName=nameofwave(w1)

			duplicate/free w1 Test,Test1	
						setscale/P x dimoffset(w1,0), dimdelta(w1,0), "s", Test,Test1//ensuring scaling of input and output wave are the same

			wave/T parmWave=root:DE_HMM:MenuStuff:ParmWave
			variable filtering=str2num(parmWave[6][1])

//					
//			if(waveexists($(nameofwave(w1)+"_drift"))==0)
//				DetermineDrift()
//			else
//			endif
//			wave DriftWave=$(S_value+"_drift")
	
			if(filtering>=1)
	
				Smooth/S=2 Filtering,Test1
				if(str2num(parmWave[7][1])==1)	
					if(Filtering>4)
						Resample/DOWN=(floor(Filtering/2))/N=1/WINF=None Test1

					endif
					
				else
				
					Resample/DOWN=(str2num(parmWave[7][1]))/N=1/WINF=None Test1
				endif

			else
				DE_Filtering#TVD1D_denoise(Test,filtering,Test1)
				Resample/DOWN=(str2num(parmWave[7][1]))/N=1/WINF=None Test1

			endif
			
			Test1*=1e9
			newdataFolder/o/S $(BaseFoldStr+"JavaHMMFit")
			print "Dim"+num2str(dimdelta(test1,0))
			DriftMarkovFitter(Test1, str2num(parmWave[0][1]), str2num(parmWave[1][1]), 1,str2num(parmWave[2][1]),str2num(parmWave[3][1]), str2num(parmWave[4][1]), str2num(parmWave[5][1]))
			wave HidMar0,HidMar1,HidMar2,HidMar3,HidMar4,HidMarParms0
			variable state0=HidMarParms0[1][0]/1e9
			variable state1=state0+HidMarParms0[1][1]/1e9
			HidMar1/=1e9
			HidMar3/=1e9
			Test1/=1e9
			variable Rate1=HidMarParms0[1][11]/dimdelta(Test1,0)
			variable Rate2=HidMarParms0[1][12]/dimdelta(Test1,0)
			print rate1
			print rate2
			state0+=mean(HidMar3)
			state1+=mean(HidMar3)
			String NewNoteInfo=note(test1)
			NewNoteInfo=replacestringbykey("DE_JavaHMM_UpperRate",NewNoteInfo,num2str(rate1),":","\r")
			NewNoteInfo=replacestringbykey("DE_JavaHMM_LowerRate",NewNoteInfo,num2str(rate2),":","\r")		
			NewNoteInfo=replacestringbykey("DE_JavaHMM_State0",NewNoteInfo,num2str(state0),":","\r")
			NewNoteInfo=replacestringbykey("DE_JavaHMM_State1",NewNoteInfo,num2str(state1),":","\r")
			duplicate/o HidMar1 $(S_value+"_Javafit")
			duplicate/o HidMar2 $(S_value+"_Javast")
			duplicate/o HidMar3 $(S_value+"_Javadr")

			wave States=$(S_value+"_Javast")
			setscale/P x dimoffset(test1,0), dimdelta(test1,0), "s", $(S_value+"_Javafit"),$(S_value+"_Javast"), $(S_value+"_Javadr")//ensuring scaling of input and output wave are the same
			duplicate/o Test1 $(S_value+"_JavaHMM")
			wave HmmIN=$((S_value+"_JavaHMM")[0,30])
		note/k HmmIN,NewNoteInfo
				note/k $(S_value+"_Javafit"),NewNoteInfo
				note/k $(S_value+"_Javadr"),NewNoteInfo
				note/k $(S_value+"_Javast"),NewNoteInfo

			PlotandProcessJavaHMM($(S_value+"_Javast"),w1,HmmIN,BaseName,state0,state1)

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

	wave fitwave=$(NameString+"_fit")

	//Position histograms
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
	
	//Make lifetimes
	make/o/n=0 $(NameString+"_HMM_ULT"),$(NameString+"_HMM_FLT")
	wave UnfoldedLifetime=$(NameString+"_HMM_ULT")
	wave FoldedLifetime=$(NameString+"_HMM_FLT")
	CalcLifetimes(statewave,UnfoldedLifetime,FoldedLifetime)
	variable totaltransition=numpnts(FoldedLifetime)+numpnts(unFoldedLifetime)
	variable FoldedLTAvg=mean(FoldedLifetime)
	variable UnfoldedLTAvg= mean(UnFoldedLifetime)
	String HMMUpperTime=num2str(1/str2num(Stringbykey("DE_HMM_URate",note(ForceWave),":","\r"))*1e3)
	String HMMLowerTime=num2str(1/str2num(Stringbykey("DE_HMM_LRate",note(ForceWave),":","\r"))*1e3)

	//Make Lifetime histograms and FIT
	make/free/n=0 OutUHist,OutLHist,OutUHist,FitOutUHist,FitOutLHist
	variable totallength=dimdelta(statewave,0)*dimsize(statewave,0)
	make/o/n=(20) $(nameofwave(RawForceWave)+"_Hmm_FH"),$(nameofwave(RawForceWave)+"_HMM_UH"),$("fit_"+nameofwave(RawForceWave)+"_Hmm_FH"),$("fit_"+nameofwave(RawForceWave)+"_Hmm_UH")
	wave FLTHist=$(nameofwave(RawForceWave)+"_Hmm_FH")
	wave UFLTHist=$(nameofwave(RawForceWave)+"_HMM_UH")
	wave FLFit= $("fit_"+nameofwave(RawForceWave)+"_Hmm_FH")	
	wave UFLFit= $("fit_"+nameofwave(RawForceWave)+"_Hmm_UH")	

	variable/C Lifetimes=ReturnStateLifetimes(statewave,totallength/2000,UFLTHist,UFLFit,FLTHist,FLFit)


	//print/C Lifetimes

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
	String WindowName=nameofwave(RawForceWave)+"_HMM_Win"
	DE_HMM#MakeNicePlot(RawForceWave,ForceWave,fitwave,FHist,WFIT,FLTHist,FLFit,UFLTHist,UFLFit,State1Pos,State2Pos,WindowName)
	variable Spacing=abs(State2-State1)*1e9
	variable StatePercentage=round(100*(1-sum(StateWave)/numpnts(StateWave)))
	variable fitPerCentage=round(100*w_coefs[1]/(w_coefs[1]+w_coefs[2]))
	TextBox/W=$WindowName/N=Populations/X=40/Y=0/C/N=text0/F=0 num2str(StatePercentage)+"%\r"+num2str(Spacing)+" nm"
	string lifetimeUS,lifetimeFS,lifetimeUAS,lifetimeFAS
	print lifetimes
	sprintf lifetimeUS, "%0.2f",real(lifetimes)*1e3
	sprintf lifetimeFS, "%0.2f",imag(Lifetimes)*1e3
	sprintf lifetimeUAS, "%0.2f",1e3*UnfoldedLTAvg
	sprintf lifetimeFAS, "%0.2f",1e3*FoldedLTAvg	
			
	TextBox/W=$WindowName/N=FLifetimes/X=22/Y=0/C/N=text0/F=0 "\\K(19712,44800,18944)Folded:\r"+HMMUpperTime+" ms (HMM)\r"+lifetimeFAS+" ms (Avg)\r"+lifetimeFS+" ms (Fit)"
	TextBox/W=$WindowName/N=ULifetimes/X=-1/Y=0/C/N=text0/F=0 "\\K(14848,32256,47104)UnFolded:\r"+HMMLowertime+" ms (HMM)\r"+lifetimeUAS+" ms (Avg)\r"+lifetimeUS+" ms (Fit)"
	wave W_Coef,W_Sigma,w_coefs,W_fitConstants
	
	make/o/n=(1,8) Test
	Test[0][0]=StatePercentage
	Test[0][1]=Spacing
	Test[0][2]=str2num(HMMUpperTime)
	Test[0][3]=str2num(HMMLowertime)
	Test[0][4]=str2num(lifetimeFAS)
	Test[0][5]=str2num(lifetimeUAS)
	Test[0][6]=str2num(lifetimeFS)
	Test[0][7]=str2num(lifetimeUS)
	DoWindow CurrRes
	if(V_flag==1)
		killwindow CurrRes
	endif
	edit/W=(10,500,800,625)/N=CurrRes	 Test
	
	killwaves W_Coef,W_Sigma,w_coefs,W_fitConstants				
	
					
	//killwaves FHist
end




Static Function PlotandProcessJavaHMM(StateWave,RawForceWave,ForceWave,NameString,state1,state2)
	wave statewave,ForceWave,RawForceWave

	variable state1,state2

	string NameString
	
	wave fitwave=$(NameString+"_Javafit")

	//Position histograms
	make/o/n=40 $(nameofwave(ForceWave)+"_JavaHist")
	wave FHist=$(nameofwave(ForceWave)+"_JavaHist")

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
	
	//Make lifetimes
	make/o/n=0 $(NameString+"_JavaHMM_ULT"),$(NameString+"_JavaHMM_FLT")
	wave UnfoldedLifetime=$(NameString+"_JavaHMM_ULT")
	wave FoldedLifetime=$(NameString+"_JavaHMM_FLT")
	CalcLifetimes(statewave,UnfoldedLifetime,FoldedLifetime)
	variable totaltransition=numpnts(FoldedLifetime)+numpnts(unFoldedLifetime)
	variable FoldedLTAvg=mean(FoldedLifetime)
	variable UnfoldedLTAvg= mean(UnFoldedLifetime)
	print nameofwave(ForceWave)
	String HMMUpperTime=num2str(1/str2num(Stringbykey("DE_JavaHMM_UpperRate",note(ForceWave),":","\r"))*1e3)
	String HMMLowerTime=num2str(1/str2num(Stringbykey("DE_JavaHMM_LowerRate",note(ForceWave),":","\r"))*1e3)

	//Make Lifetime histograms and FIT
	make/free/n=0 OutUHist,OutLHist,OutUHist,FitOutUHist,FitOutLHist
	variable totallength=dimdelta(statewave,0)*dimsize(statewave,0)
	make/o/n=(20) $(nameofwave(RawForceWave)+"_Hmm_FH"),$(nameofwave(RawForceWave)+"_HMM_UH"),$("fit_"+nameofwave(RawForceWave)+"_Hmm_FH"),$("fit_"+nameofwave(RawForceWave)+"_Hmm_UH")
	wave FLTHist=$(nameofwave(RawForceWave)+"_Hmm_FH")
	wave UFLTHist=$(nameofwave(RawForceWave)+"_HMM_UH")
	wave FLFit= $("fit_"+nameofwave(RawForceWave)+"_Hmm_FH")	
	wave UFLFit= $("fit_"+nameofwave(RawForceWave)+"_Hmm_UH")	

	variable/C Lifetimes=ReturnJavaStateLifetimes(statewave,totallength/2000,UFLTHist,UFLFit,FLTHist,FLFit)


	//print/C Lifetimes

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

	DE_HMM#MakeNicePlot(RawForceWave,ForceWave,fitwave,FHist,WFIT,FLTHist,FLFit,UFLTHist,UFLFit,State1Pos,State2Pos,nameofwave(RawForceWave)+"_JavaHMM_Win")
	TextBox/N=Populations/X=45/Y=0/C/N=text0/F=0 num2str(round(100*w_coefs[1]/(w_coefs[1]+w_coefs[2])))+"%"
	string lifetimeUS,lifetimeFS,lifetimeUAS,lifetimeFAS
	print lifetimes
	sprintf lifetimeUS, "%0.2f",real(lifetimes)*1e3
	sprintf lifetimeFS, "%0.2f",imag(Lifetimes)*1e3
	sprintf lifetimeUAS, "%0.2f",1e3*UnfoldedLTAvg
	sprintf lifetimeFAS, "%0.2f",1e3*FoldedLTAvg	
			
	TextBox/N=FLifetimes/X=22/Y=0/C/N=text0/F=0 "\\K(19712,44800,18944)Folded:\r"+HMMUpperTime+" ms (HMM)\r"+lifetimeFAS+" ms (Avg)\r"+lifetimeFS+" ms (Fit)"
	TextBox/N=ULifetimes/X=-1/Y=0/C/N=text0/F=0 "\\K(14848,32256,47104)UnFolded:\r"+HMMLowertime+" ms (HMM)\r"+lifetimeUAS+" ms (Avg)\r"+lifetimeUS+" ms (Fit)"
	wave W_Coef,W_Sigma,w_coefs,W_fitConstants
	killwaves W_Coef,W_Sigma,w_coefs,W_fitConstants				
					
	//killwaves FHist
end

Static Function/C ReturnStateLifetimes(HMM_FIt,HistStepO,OutUHist,FitOutUHist,OutFHist,FitOutFHist)
	Wave HMM_Fit,OutUHist,OutFHist,FitOutFHist,FitOutUHist
	variable HistStepO
	make/free/n=0 FreeUTime,FreeFTime
	CalcLifetimes(HMM_FIT,FreeUTime,FreeFTime)
	variable HistStepF=HistStepO
	variable HistStepU=HistStepO
	String UpperRate=Stringbykey("DE_HMM_URate",note(HMM_FIT),":","\r")
	if(cmpstr(UpperRate,"")==0)
	else
		HistStepU=1/str2num(Stringbykey("DE_HMM_URate",note(HMM_FIT),":","\r"))/5
		HistStepF=1/str2num(Stringbykey("DE_Hmm_LRate",note(HMM_FIT),":","\r"))/5

	endif
	variable LTBins=10
	make/free/n=(LTBins) FreeFHist,FreeFHistFit
	Histogram/C/B={10e-4,3*HistStepF,LTBins} FreeFTime,FreeFHist
	make/o/n=3 W_coef
	W_Coef={0,wavemax(FreeFHist),1/HistStepF}
	CurveFit/Q/W=2/N/H="100"/NTHR=0 exp  FreeFHist/D=FreeFHistFit
	Setscale/P x, pnt2x(FreeFHist,0), dimdelta(FreeFHist,0),FreeFHistFit
	
	variable/C Result=Cmplx(0,1/w_coef[2])
	make/free/n=(LTBins) FreeUHist,FreeUHistFit
	Histogram/C/B={10e-4,3*HistStepU,LTBins} FreeUTime,FreeUHist
	W_Coef={0,wavemax(FreeFHist),1/HistStepU}
	CurveFit/Q/W=2/N/H="100"/NTHR=0 exp  FreeUHist/D=FreeUHistFit
	Setscale/p x, pnt2x(FreeUHist,0), dimdelta(FreeUHist,0),FreeUHistFit
	Result+=cmplx(1/w_coef[2],0)
	duplicate/o FreeUHist OutUHist
	duplicate/o FreeFHist OutFHist
	duplicate/o FreeUHistFit FitOutUHist
	duplicate/o FreeFHistFit FitOutFHist
	return result
end
Static Function/C ReturnJavaStateLifetimes(HMM_FIt,HistStepO,OutUHist,FitOutUHist,OutFHist,FitOutFHist)
	Wave HMM_Fit,OutUHist,OutFHist,FitOutFHist,FitOutUHist
	variable HistStepO
	make/free/n=0 FreeUTime,FreeFTime
	CalcLifetimes(HMM_FIT,FreeUTime,FreeFTime)
	variable HistStepF=HistStepO
	variable HistStepU=HistStepO
	String UpperRate=Stringbykey("DE_JavaHMM_UpperRate",note(HMM_FIT),":","\r")
	if(cmpstr(UpperRate,"")==0)
	else
		HistStepU=1/str2num(Stringbykey("DE_JavaHMM_UpperRate",note(HMM_FIT),":","\r"))/5
		HistStepF=1/str2num(Stringbykey("DE_JavaHmm_LowerRate",note(HMM_FIT),":","\r"))/5

	endif
	print note(HMM_FIT)
	variable LTBins=10
	make/free/n=(LTBins) FreeFHist,FreeFHistFit
	Histogram/C/B={10e-4,3*HistStepF,LTBins} FreeFTime,FreeFHist
	make/o/n=3 W_coef
	W_Coef={0,wavemax(FreeFHist),1/HistStepF}
	CurveFit/Q/W=2/N/H="100"/NTHR=0 exp  FreeFHist/D=FreeFHistFit
	Setscale/P x, pnt2x(FreeFHist,0), dimdelta(FreeFHist,0),FreeFHistFit
	
	variable/C Result=Cmplx(0,1/w_coef[2])
	make/free/n=(LTBins) FreeUHist,FreeUHistFit
	Histogram/C/B={10e-4,3*HistStepU,LTBins} FreeUTime,FreeUHist
	W_Coef={0,wavemax(FreeFHist),1/HistStepU}
	CurveFit/Q/W=2/N/H="100"/NTHR=0 exp  FreeUHist/D=FreeUHistFit
	Setscale/p x, pnt2x(FreeUHist,0), dimdelta(FreeUHist,0),FreeUHistFit
	Result+=cmplx(1/w_coef[2],0)
	duplicate/o FreeUHist OutUHist
	duplicate/o FreeFHist OutFHist
	duplicate/o FreeUHistFit FitOutUHist
	duplicate/o FreeFHistFit FitOutFHist
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
	String list = WaveList("*Sep_1", ";", "")+WaveList("*Sep_2", ";", "")+WaveList("*Sep_3", ";", "")
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
Static Function MakeNicePlot(RawForce,SmForce,ForceFit,HistWave,HistFit,LT1,LT1Fit,LT2,LT2Fit,state1,state2,WindowName)
	wave RawForce,SmForce,HistWave,HistFit,LT1,LT2,LT1Fit,LT2Fit,ForceFit, state1,state2
	string WindowName
//	string WindowName=nameofwave(RawForce)+"_HMM_Win"
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

	AppendToGraph/W=$WindowName/B=B3/L=L3 LT2
	AppendToGraph/W=$WindowName/B=B3/L=L3 LT2Fit

	//	//SetAxis L1 2.1312133e-08,2.9592798e-08
	ModifyGraph/W=$WindowName mode($nameofwave(State1))=5,hbFill($nameofwave(State1))=4,rgb($nameofwave(State1))=(0,65280,65280);DelayUpdate
	ModifyGraph/W=$WindowName mode($nameofwave(State2))=5,hbFill($nameofwave(State2))=4,rgb($nameofwave(State2))=(0,65280,65280)

	ModifyGraph/W=$WindowName tick=2,fSize=9,lblPosMode=1,lblPos=42,standoff=0,font="Arial"
	ModifyGraph/W=$WindowName axisEnab(bottom)={0,0.45}
	ModifyGraph/W=$WindowName axisEnab(B1)={0.47,0.60}
	ModifyGraph/W=$WindowName axisEnab(B2)={0.63,.8}
		ModifyGraph/W=$WindowName axisEnab(B3)={0.83,1}

	ModifyGraph/W=$WindowName freePos(B1)={0,L1}
	ModifyGraph/W=$WindowName freePos(L1)={0,B1}
	ModifyGraph/W=$WindowName freePos(B2)={0,L2}
	ModifyGraph/W=$WindowName freePos(L2)={0,B2}
	ModifyGraph/W=$WindowName freePos(B3)={0,L3}
	ModifyGraph/W=$WindowName freePos(L3)={0,B3}
	ModifyGraph/W=$WindowName rgb($nameofwave(HistFit))=(0,0,0)
	ModifyGraph/W=$WindowName rgb($(nameofwave(RawForce)))=(65280,48896,48896)
	ModifyGraph/W=$WindowName hideTrace($(nameofwave(RawForce)))=1
	ModifyGraph/W=$WindowName rgb($(nameofwave(ForceFit)))=(0,0,0)
	ModifyGraph/W=$WindowName margin(left)=36,margin(bottom)=29,margin(top)=21,margin(right)=14;DelayUpdate
	ModifyGraph/W=$WindowName mode($nameofwave(LT1))=5,hbfill($nameofwave(LT1))=5;DelayUpdate
	ModifyGraph/W=$WindowName rgb($nameofwave(LT1))=(19712,44800,18944);DelayUpdate
	ModifyGraph/W=$WindowName useMrkStrokeRGB($nameofwave(LT1))=1,mode($nameofwave(LT2))=5;DelayUpdate
	ModifyGraph/W=$WindowName hbfill($nameofwave(LT2))=5;DelayUpdate
	ModifyGraph/W=$WindowName rgb($nameofwave(LT2))=(14848,32256,47104);DelayUpdate
	ModifyGraph/W=$WindowName useMrkStrokeRGB($nameofwave(LT2))=1
	ModifyGraph/W=$WindowName width=800,height=144
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
	variable DownTrans,UpTrans
	for(n=1;n<numpnts(States);n+=1)
		if(States[n]==States[n-1])
		else
			if(States[n-1]==1)
			insertpoints 0,1,UpperLifetime
			UpperLifetime[0]=((n-1-StartPoint))
			StartPoint=n
			DownTrans+=1
			else
			insertpoints 0,1,LowerLifetime
			LowerLifetime[0]=((n-1-StartPoint))
			StartPoint=n
			UpTrans+=1
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
