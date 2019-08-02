#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma ModuleName = DE_EquilJump
#include ":\Misc_PanelPrograms\AsylumNaming"
#include ":\Misc_PanelPrograms\Panel Progs"
#include "DE_Filtering"
#include <Multi-peak fitting 2.0>

function TresEquil(ZsnsrWave,SepWave,SepAddOff,ForceWave,ForceAddOff,BaseOutString)
	wave ZsnsrWave,SepWave,ForceWave
	variable SepAddOff,ForceAddOff
	String BaseOutString
	GenerateRegions(ZsnsrWave)
	wave MidPoints
	variable q,FOffset,SepOffset,length
	string Fdup,Sdup
	string FWaveName=nameofwave(ForceWave)+"_cut"
	string SWaveName=nameofwave(SepWave)+"_cut"
	string ZWaveName=nameofwave(ZsnsrWave)+"_cut"
	
	variable totaltime= pnt2x(ForceWave,numpnts(ForceWave)-1)
	variable stepdown,stepup
	if(totaltime>3&&totaltime<4)
		stepdown=.1
		stepup=.97
	elseif(totaltime>15&&totaltime<16)
		stepdown=.1
		stepup=4.9
	elseif(totaltime>2&&totaltime<3)
		stepdown=.15
		stepup=.45
	elseif(totaltime>9.5&&totaltime<10.5)
		stepdown=.15
		stepup=2.9
	elseif(totaltime>30.5&&totaltime<31.5)
		stepdown=.15
		stepup=9.9
	else
		stepdown=.15
		stepup=totaltime/3-1.5/3
	endif

//	DoWindow EqJump
//	if(V_flag==0)
//		display/n=EqJump
//	else
//		killwindow EqJump
//		Display/n=EqJump
//	endif
	for(q=0;q<numpnts(midpoints);q+=1)
	
		if(mod(q,2)==0)//even
			CutoutRegion(SepWave,MidPoints,q,stepdown)
			wave/z CutS=$SWaveName
			wavestats/q CutS
			SepOffset=(v_avg+SepAddOff)
			
			CutoutRegion(ForceWave,MidPoints,1,stepdown)
			wave/z CutF=$FWaveName
			wavestats/q CutF
			//FOffset=v_avg+ForceAddOff
			FOffset=ForceAddOff
		else //odd
			CutoutRegion(SepWave,MidPoints,q,stepup)
			wave/z CutS=$SWaveName
			CutS-=SepAddOff
			
			length=strlen(nameofwave(SepWave))-1
			Sdup=BaseOutString+"Sep_"+num2str(ceil(q/2))
			//Sdup=nameofwave(SepWave)[0,length]+num2str(ceil(q/2))

			duplicate/o CutS $Sdup

			CutoutRegion(ForceWave,MidPoints,q,stepup)
			wave/z CutF=$FWaveName
			CutF-=ForceAddOff
			CutF-=FOffset

			length=strlen(nameofwave(ForceWave))-1
			Fdup=BaseOutString+"Force_"+num2str(ceil(q/2))
			
			duplicate/o CutF $Fdup
	
		endif
		
		
	
	
	endfor
	variable np=numpnts(midpoints)/2
	killwaves CutS,CutF,MidPoints

	return np
end

Static function GenerateRegions(Wavein)
	wave Wavein
	FilterWaveDown(wavein)
	string filteredname=nameofwave(wavein)+"_dif_samp"
	wave wfiltered=$filteredname
	insertpoints (numpnts(wfiltered)),10, wfiltered
	FindPeaks(wfiltered)
	wave MidPoints
	killwaves wfiltered
end

Static function CutoutRegion(WavetoCut,Midvalues,cutoutpoint,range)
	wave WavetoCut,Midvalues
	variable cutoutpoint,range
	string Resultwave=nameofwave(wavetocut)+"_cut"
	variable Start=real(ApplyRange(Midvalues,range,cutoutpoint))
	if(cutoutpoint==0)
	//start=0
	endif
	variable Ends=imag(ApplyRange(Midvalues,range,cutoutpoint))
	duplicate/o/r=(Start,Ends) WavetoCut $Resultwave
end

function StepCutoutRegion(WavetoCut,Midvalues,cutoutpoint,range)
	wave WavetoCut,Midvalues
	variable cutoutpoint,range
	string Resultwave=nameofwave(wavetocut)+"_cut"
	variable Start=real(ApplyRange(Midvalues,range,cutoutpoint))
	//if(cutoutpoint==0)
	//start=0
	//
	//endif
	variable Ends=imag(ApplyRange(Midvalues,range,cutoutpoint))
	duplicate/o/r=(Start,Ends) WavetoCut $Resultwave
end

function FilterWaveDown(wavein)
	wave wavein
	variable rate=1/dimdelta(wavein,0)/1000
	variable down=rate*2*2
	variable Ns=rate*4*2
	if((round(Ns/2)*2+1)>=40000)
		rate=1/dimdelta(wavein,0)/10000
		down=rate*2*2
		Ns=rate*4*2
	endif
	duplicate/o wavein Test
	Resample/DOWN=(down)/N=(round(Ns/2)*2+1) Test

	
	string difname=nameofwave(wavein)+"_DIF"
	wave difwave=$difname
	Differentiate Test/D=difwave
	string name=difname+"_samp"
	Duplicate/O difwave,$name;
	wave waveout=$name
	


	waveout=abs(waveout)
	killwaves difwave,Test
	return 0
end

function FindPeaks(wavein)
	wave wavein
	string xdata="_calculated_;"
	
	WAVE/Z wx=$xdata
	Variable pBegin=0, pEnd= numpnts(wavein)-1
	Variable maxPeaks=150, minPeakPercent=30
	Variable noiselevel=0
	Variable smoothingFactor=10
	AutoFindPeaksWorker(wavein, wx, pBegin, pEnd, maxPeaks, minPeakPercent, noiseLevel, smoothingFactor)
	killwindow showpeaks
	killwindow table0
	wave WA_PeakCentersX, W_AutoPeakInfo,WA_PeakCentersY
	KillWaves  W_AutoPeakInfo,WA_PeakCentersY

	reorderwave(WA_PeakCentersX)
	duplicate/o WA_PeakCentersX MidPoints
	MidPointFind(WA_PeakCentersX,MidPoints)
	killwaves WA_PeakCentersX
	return 0
end
Static function ReorderWave(w1)
	wave w1
	duplicate/o w1 Hold1
	Hold1=.11
	variable i
	wavestats/q w1
variable ends=v_npnts
	for(i=0;i<(ends);i+=1)
		wavestats/q w1
		Hold1[i]=w1[v_minloc]
		w1[v_minloc]=Nan
	endfor
	w1=Hold1
	killwaves hold1
	
end
function/C ApplyRange(wavein,range,number)
	wave wavein
	variable range,number
	wavestats/q wavein
	return cmplx(wavein[number]-range/2,wavein[number]+range/2)

end

function MidPointFind(inwave,outwave)
	wave inwave,outwave
	variable n
	outwave[n]=inwave[n]/2

	for(n=1;n<numpnts(inwave);n+=1)
		outwave[n]=inwave[n-1]/2+inwave[n]/2
	endfor
end



//
Static Function ButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	Strswitch(ba.ctrlname)
			string saveDF
			
		case "de_EQJ_button0":
			switch( ba.eventCode )
				case 2: // mouse up
					saveDF = GetDataFolder(1)
					controlinfo/W=EqJPanel de_EQJ_popup0
					SetDataFolder s_value
					controlinfo/W=EqJPanel de_EQJ_popup1
					wave w1=$S_value
					wave w2=$ReplaceString("Force", S_Value, "Sep")
					wave w3=$ReplaceString("Force", S_Value, "ZSnsr")
					print nameofwave(w1)

					controlinfo/W=EqJPanel de_EQJ_setvar0
					variable ForceOff =V_value
					controlinfo/W=EqJPanel de_EQJ_setvar1
					variable SepOff=V_value
					controlinfo/W=EqJPanel de_EQJ_setvar2
					variable Filtering=V_value
					
					struct ForceWave Name1
					DE_Naming#WavetoStruc(nameofwave(w1),Name1)
					
	
					string basename=name1.Name+name1.Snum
					print basename
					
					variable number=DE_EQuilJump#TresEquil(w3,w2,SepOff,w1,ForceOff,basename)
					variable n=0
					variable StartingState,FoldedLTAvg,UnfoldedLTAvg
					
					for(n=1;n<=number;n+=1)
						//make/o/n=0 $(basename+"Sep_"+num2str(n)+"_Sm")
						wave w0=$(basename+"Sep_"+num2str(n))
						duplicate/o w0 $(basename+"Sep_"+num2str(n)+"_Sm")
						wave w1=$(basename+"Sep_"+num2str(n)+"_Sm")
						//DE_FILTERING#TVD1D_denoise(w0,Filtering,w1)
						Smooth/S=2 Filtering,w1
						make/o/n=100 $(nameofwave(w1)+"_Hist")
						wave w1Hist=$(nameofwave(w1)+"_Hist")
						wavestats/q w1
						variable SMax=v_max
						variable SMin=v_min
						Histogram/C w1 w1Hist
						Integrate/METH=1 w1Hist/D=w1Hist_Int;DelayUpdate
						w1Hist_Int/=w1Hist_Int[numpnts(w1Hist_Int)-1]
						FindLevel/P/Q w1Hist_Int,0.33
						variable firstthird=v_levelx
							FindLevel/P/Q w1Hist_Int,0.66
						variable lastthird=v_levelx
						killwaves w1Hist_Int
						//wavestats/Q/R=[0,firstthird] w1Hist

						wavestats/q/R=[0,numpnts(w1Hist)/2] w1Hist
						variable P1=v_maxloc
						variable H1=v_max
						//wavestats/Q/R=[lastthird,numpnts(w1Hist)-1] w1Hist

						wavestats/q/R=[numpnts(w1Hist)/2,numpnts(w1Hist)-1] w1Hist
						variable P2=v_maxloc
						variable H2=v_max
						make/D/o/n=7 w_coefs
						w_coefs[0]={0,H1,H2,P1,P2,.2e-9,.2e-9}
						print w_coefs
						FuncFit/ODR=1/W=2/Q/H="1000000"/NTHR=0 DGauss w_coefs  w1Hist /D
						wave WFIT=$("fit_"+nameofwave(w1Hist))
						variable midpoint=(w_coefs[3]+w_coefs[4])/2
						make/free/n=0 OutWave,FoldedLiftime,UnFoldedLiftime
						StartingState=GenerateThreeCuts(w1,w_coefs[3],w_coefs[4],OutWave)
						//StartingState=FindTransitions(w1,midpoint,OutWave)
						duplicate/o OutWave $(nameofwave(w1)+"TransitionX"),$(nameofwave(w1)+"TransitionY")
						wave trx=$(nameofwave(w1)+"TransitionX")
						wave trys=$(nameofwave(w1)+"Transitiony")
						trys=w1(trx)
						CalculateTwoStateLifeTimes(OutWave,FoldedLiftime,UnFoldedLiftime,StartingState)
						FoldedLTAvg=mean(FoldedLiftime)
						UnfoldedLTAvg= mean(UnFoldedLiftime)
						make/o/n=(40) $(nameofwave(w1)+"_FLHist")
						wave FLTHist=$(nameofwave(w1)+"_FLHist")
						Histogram/C/B={1e-5,(dimdelta(w1,0)*numpnts(w1)/numpnts(Outwave)/4),40} FoldedLiftime FLTHist
						make/o/n=(40) $(nameofwave(w1)+"_UFLHist")
						wave UFLTHist=$(nameofwave(w1)+"_UFLHist")
						Histogram/C/B={1e-5,(dimdelta(w1,0)*numpnts(w1)/numpnts(Outwave))/4,40} UnFoldedLiftime UFLTHist
						make/D/o/n=3 w_coefL1
						w_coefL1[0]={0,100,1e-2}
						CurveFit/W=2/Q/G/H="100"/NTHR=0 exp_XOffset kwCWave=w_coefL1 ,FLTHist /D
						wave FLFit= $("fit_"+nameofwave(FLTHist))	
						make/D/o/n=3 w_coefL2
						w_coefL2[0]={0,100,1e-2}
						CurveFit/W=2/Q/G/H="100"/NTHR=0 exp_XOffset kwCWave=w_coefL2,UFLTHist /D 
						wave UFLFit= $("fit_"+nameofwave(UFLTHist))	
						wave WFIT=$("fit_"+nameofwave(w1Hist))				
						MakeNicePlot(w0,w1,w1Hist,WFIT,FLTHist,FLFit,UFLTHist,UFLFit,Trx,trys)
						TextBox/N=Populations/X=35/Y=1/C/N=text0/F=0 num2str(round(100*w_coefs[1]/(w_coefs[1]+w_coefs[2])))+"%"
						string lifetime1S,lifetime2S,lifetime1AS,lifetime2AS
						sprintf lifetime1S, "%0.2f",-1e3*w_coefL1[2]*ln(1/2)
						sprintf lifetime2S, "%0.2f",-1e3*w_coefL2[2]*ln(1/2)
						sprintf lifetime1AS, "%0.2f",1e3*FoldedLTAvg
						sprintf lifetime2AS, "%0.2f",1e3*UnfoldedLTAvg		
						TextBox/N=Lifetimes/X=0/Y=1/C/N=text0/F=0 "\\K(19712,44800,18944)Folded: "+lifetime1AS+"("+lifetime1S+") ms\r\\K(14848,32256,47104)UnFolded: "+lifetime2AS+"("+lifetime2S+") ms"
					endfor
					
					wave W_sigma,W_FindLevels,W_fitConstants
					killwaves w_coefs,w_coefL2,w_coefL1,W_sigma,W_FindLevels,W_fitConstants
					//SetDataFolder saveDF

					break
				case -1: // control being killed
					break
			endswitch
			break

	
	endswitch
	return 0
End

 Function StateandLifetimeswithFiltering(w0,filtering)
	wave w0
	variable filtering
	variable StartingState,FoldedLTAvg,UnfoldedLTAvg
	DFREF startingDFR=GetDataFolderDFR( )
	setdatafolder GetWavesDataFolderDFR(w0)
	duplicate/o w0 $(nameofwave(w0)+"_Sm")
	wave w1=$(nameofwave(w0)+"_Sm")
	if(filtering>=1)
	
		Smooth/S=2 Filtering,w1

	else
		DE_FILTERING#TVD1D_denoise(w0,Filtering,w1)
	
	endif

	make/o/n=100 $(nameofwave(w1)+"_Hist")
	wave w1Hist=$(nameofwave(w1)+"_Hist")
	wavestats/q w1
	variable SMax=v_max
	variable SMin=v_min
	Histogram/C w1 w1Hist
	Integrate/METH=1 w1Hist/D=w1Hist_Int;DelayUpdate
	w1Hist_Int/=w1Hist_Int[numpnts(w1Hist_Int)-1]
	FindLevel/P/Q w1Hist_Int,0.33
	variable firstthird=v_levelx
	FindLevel/P/Q w1Hist_Int,0.66
	variable lastthird=v_levelx
	killwaves w1Hist_Int
	wavestats/Q/R=[0,firstthird] w1Hist

	//wavestats/q/R=[0,numpnts(w1Hist)/2] w1Hist
	variable P1=v_maxloc
	variable H1=v_max
	wavestats/Q/R=[lastthird,numpnts(w1Hist)-1] w1Hist

	//wavestats/q/R=[numpnts(w1Hist)/2,numpnts(w1Hist)-1] w1Hist
	variable P2=v_maxloc
	variable H2=v_max
	make/D/o/n=7 w_coefs
	w_coefs[0]={0,H1,H2,P1,P2,.2e-9,.2e-9}
	print w_coefs
	FuncFit/ODR=1/W=2/Q/H="1000000"/NTHR=0 DGauss w_coefs  w1Hist /D
	wave WFIT=$("fit_"+nameofwave(w1Hist))
	variable midpoint=(w_coefs[3]+w_coefs[4])/2
	make/free/n=0 OutWave,FoldedLiftime,UnFoldedLiftime
	StartingState=GenerateThreeCuts(w1,w_coefs[3],w_coefs[4],OutWave)
	//StartingState=FindTransitions(w1,midpoint,OutWave)
	duplicate/o OutWave $(nameofwave(w1)+"TransitionX"),$(nameofwave(w1)+"TransitionY")
	wave trx=$(nameofwave(w1)+"TransitionX")
	wave trys=$(nameofwave(w1)+"Transitiony")
	trys=w1(trx)
	CalculateTwoStateLifeTimes(OutWave,FoldedLiftime,UnFoldedLiftime,StartingState)
	FoldedLTAvg=mean(FoldedLiftime)
	UnfoldedLTAvg= mean(UnFoldedLiftime)
	make/o/n=(40) $(nameofwave(w1)+"_FLHist")
	wave FLTHist=$(nameofwave(w1)+"_FLHist")
	Histogram/C/B={1e-5,(dimdelta(w1,0)*numpnts(w1)/numpnts(Outwave)/4),40} FoldedLiftime FLTHist
	make/o/n=(40) $(nameofwave(w1)+"_UFLHist")
	wave UFLTHist=$(nameofwave(w1)+"_UFLHist")
	Histogram/C/B={1e-5,(dimdelta(w1,0)*numpnts(w1)/numpnts(Outwave))/4,40} UnFoldedLiftime UFLTHist
	make/D/o/n=3 w_coefL1
	w_coefL1[0]={0,100,1e-2}
	CurveFit/W=2/Q/G/H="100"/NTHR=0 exp_XOffset kwCWave=w_coefL1 ,FLTHist /D
	wave FLFit= $("fit_"+nameofwave(FLTHist))	
	make/D/o/n=3 w_coefL2
	w_coefL2[0]={0,100,1e-2}
	CurveFit/W=2/Q/G/H="100"/NTHR=0 exp_XOffset kwCWave=w_coefL2,UFLTHist /D 
	wave UFLFit= $("fit_"+nameofwave(UFLTHist))	
	wave WFIT=$("fit_"+nameofwave(w1Hist))				
	MakeNicePlot(w0,w1,w1Hist,WFIT,FLTHist,FLFit,UFLTHist,UFLFit,Trx,trys)
	TextBox/N=Populations/X=35/Y=1/C/N=text0/F=0 num2str(round(100*w_coefs[1]/(w_coefs[1]+w_coefs[2])))+"%"
	string lifetime1S,lifetime2S,lifetime1AS,lifetime2AS
	sprintf lifetime1S, "%0.2f",-1e3*w_coefL1[2]*ln(1/2)
	sprintf lifetime2S, "%0.2f",-1e3*w_coefL2[2]*ln(1/2)
	sprintf lifetime1AS, "%0.2f",1e3*FoldedLTAvg
	sprintf lifetime2AS, "%0.2f",1e3*UnfoldedLTAvg		
	TextBox/N=Lifetimes/X=0/Y=1/C/N=text0/F=0 "\\K(19712,44800,18944)Folded: "+lifetime1AS+"("+lifetime1S+") ms\r\\K(14848,32256,47104)UnFolded: "+lifetime2AS+"("+lifetime2S+") ms"
	wave W_sigma,W_FindLevels,W_fitConstants,M_Jacobian
	killwaves w_coefs,w_coefL2,w_coefL1,W_sigma,W_FindLevels,W_fitConstants,M_Jacobian
	SetDataFolder startingDFR

end


Static Function CalculateLifetimesofOne(RawSep,SepWaveIn,FoldedLoc,UnfoldedLoc)

	wave RawSep,SepWaveIn
	variable FoldedLoc,UnfoldedLoc
	
	
	make/o/n=100 $(nameofwave(SepWaveIn)+"_Hist")
	wave w1Hist=$(nameofwave(SepWaveIn)+"_Hist")
	
	Histogram/C SepWaveIn w1Hist

	wavestats/q/R=[0,numpnts(w1Hist)/2] w1Hist
	variable P1=FoldedLoc
	variable H1=v_max
	//wavestats/Q/R=[lastthird,numpnts(w1Hist)-1] w1Hist

	wavestats/q/R=[numpnts(w1Hist)/2,numpnts(w1Hist)-1] w1Hist
	variable P2=UnfoldedLoc
	variable H2=v_max
	make/D/o/n=7 w_coefs
	w_coefs[0]={0,H1,H2,P1,P2,.2e-9,.2e-9}
	//print w_coefs
	FuncFit/ODR=1/W=2/Q/H="1000000"/NTHR=0 DGauss w_coefs  w1Hist /D
	wave WFIT=$("fit_"+nameofwave(w1Hist))
	variable midpoint=(w_coefs[3]+w_coefs[4])/2
	make/free/n=0 OutWave,FoldedLiftime,UnFoldedLiftime
						
						
						
	make/free/n=0 OutWave,FoldedLiftime,UnFoldedLiftime
	variable StartingState=GenerateThreeCuts(SepWaveIn,FoldedLoc,UnfoldedLoc,OutWave)
	//StartingState=FindTransitions(w1,midpoint,OutWave)
	duplicate/o OutWave $(nameofwave(SepWaveIn)+"TransitionX"),$(nameofwave(SepWaveIn)+"TransitionY")
	wave trx=$(nameofwave(SepWaveIn)+"TransitionX")
	wave trys=$(nameofwave(SepWaveIn)+"Transitiony")
	trys=SepWaveIn(trx)
	CalculateTwoStateLifeTimes(OutWave,FoldedLiftime,UnFoldedLiftime,StartingState)
	variable FoldedLTAvg=mean(FoldedLiftime)
	variable UnfoldedLTAvg= mean(UnFoldedLiftime)
	make/o/n=(15) $(nameofwave(SepWaveIn)+"_FLHist")
	wave FLTHist=$(nameofwave(SepWaveIn)+"_FLHist")
	Histogram/C/B={1e-5,(dimdelta(SepWaveIn,0)*numpnts(SepWaveIn)/numpnts(Outwave)/1.5),15} FoldedLiftime FLTHist
	make/o/n=(15) $(nameofwave(SepWaveIn)+"_UFLHist")
	wave UFLTHist=$(nameofwave(SepWaveIn)+"_UFLHist")
	Histogram/C/B={1e-5,(dimdelta(SepWaveIn,0)*numpnts(SepWaveIn)/numpnts(Outwave))/1.5,15} UnFoldedLiftime UFLTHist
	make/D/o/n=3 w_coefL1
	w_coefL1[0]={0,100,1e-2}
	CurveFit/W=2/Q/G/H="100"/NTHR=0 exp_XOffset kwCWave=w_coefL1 ,FLTHist /D
	wave FLFit= $("fit_"+nameofwave(FLTHist))	
	make/D/o/n=3 w_coefL2
	w_coefL2[0]={0,100,1e-2}
	CurveFit/W=2/Q/G/H="100"/NTHR=0 exp_XOffset kwCWave=w_coefL2,UFLTHist /D 
	wave UFLFit= $("fit_"+nameofwave(UFLTHist))	
	wave WFIT=$("fit_"+nameofwave(w1Hist))		
	MakeNicePlot(RawSep,SepWaveIn,w1Hist,WFIT,FLTHist,FLFit,UFLTHist,UFLFit,trx,trys)
TextBox/N=Populations/X=35/Y=1/C/N=text0/F=0 num2str(round(100*w_coefs[1]/(w_coefs[1]+w_coefs[2])))+"%"
						string lifetime1S,lifetime2S,lifetime1AS,lifetime2AS
						sprintf lifetime1S, "%0.2f",-1e3*w_coefL1[2]*ln(1/2)
						sprintf lifetime2S, "%0.2f",-1e3*w_coefL2[2]*ln(1/2)
						sprintf lifetime1AS, "%0.2f",1e3*FoldedLTAvg
						sprintf lifetime2AS, "%0.2f",1e3*UnfoldedLTAvg		
						TextBox/N=Lifetimes/X=0/Y=1/C/N=text0/F=0 "\\K(19712,44800,18944)Folded: "+lifetime1AS+"("+lifetime1S+") ms\r\\K(14848,32256,47104)UnFolded: "+lifetime2AS+"("+lifetime2S+") ms"

	print num2str((w_coefs[1]/(w_coefs[1]+w_coefs[2])))+"\r"+num2str(w_coefs[3])+"\r"+num2str(w_coefs[4])+"\r"+"Folded: "+lifetime1AS+"("+lifetime1S+") ms\rUnFolded: "+lifetime2AS+"("+lifetime2S+") ms"
end
//
	
//						
//						
Static Function CalculateTwoStateLifeTimes(IndexWave,FoldedLT,UnFoldedLT,Starting)
	wave IndexWave,FoldedLT,UnFoldedLT
	variable Starting
	variable n
	variable currlife
	variable odd=0
	make/free/n=0 State1,State2
	for(n=1;n<(numpnts(IndexWave)-1);n+=1)
		currlife=(INdexWave[n]-IndexWave[n-1])
		if(odd==0)
			Insertpoints numpnts(State1),1, State1
			State1[numpnts(State1)-1]=currlife
			odd=1
		else
			Insertpoints numpnts(State2),1 ,State2
			State2[numpnts(State2)-1]=currlife
			odd=0

		endif
	endfor
	if(Starting==1)
		duplicate/o State1 FoldedLT
		duplicate/o State2 UnFoldedLT
	else
		duplicate/o State2 FoldedLT
		duplicate/o State1 UnFoldedLT
	endif

	return 1

end

Static Function GenerateThreeCuts(InWave,Folded,Unfolded,Transitions)
	wave Inwave,Transitions
	variable Folded,Unfolded
	variable center=(Folded+Unfolded)/2
	variable UnFoldedBarrier=Unfolded-(Unfolded-folded)*.25
	variable FoldedBarrier=Unfolded-(Unfolded-folded)*.75
	FindLevels/edge=1/Q InWave,center
	wave w_findLevels
	duplicate/free w_findlevels CEnterLevelsUp
	FindLevels/edge=2/Q InWave,center
	duplicate/free w_findlevels CEnterLevelsDown
	FindLevels/Q InWave,UnFoldedBarrier
	duplicate/free w_findlevels UnfoldedLevels
	FindLevels/Q  InWave,FoldedBarrier
	duplicate/free w_findlevels FoldedLevels


	//Yuck is all crossings, Yuck2tells what kind of crossing it is. We then sort in order of the crossing point
	Concatenate/o/NP {CEnterLevelsUp,CEnterLevelsDown,UnfoldedLevels,FoldedLevels}, YUCK
	CEnterLevelsUp=1
	CEnterLevelsDown=-1
	UnfoldedLevels=2
	FoldedLevels=-2
	Concatenate/o/NP {CEnterLevelsUp,CEnterLevelsDown,UnfoldedLevels,FoldedLevels}, YUCK2
	Sort YUCK, YUCK,YUCK2
	variable n,p
	make/free/n=0 Up,Down
	//for(n=1;n<1000;n+=1)
	for(n=1;n<numpnts(YUCK)-1;n+=1)
		if(YUCK2[n]==1)//If this is a level up
			if(YUCK2[n+1]==2&&Yuck2[n-1]==-2)//If the next point is unfolded and the previous point was folded, count it as a transition
				insertpoints numpnts(Up),1,Up
				Up[numpnts(Up)-1]=Yuck[n]
				//print "PLUS: "+num2str(YUCK[n])
			endif
			if(YUCK2[n+1]==2&&Yuck2[n-1]==-1) //If the next point was unfolded but the previous was a negative crossing, count backwards till we
				//where the last 
				p=1
				do
					if(n==1||p>=n)
					break
					endif
					p+=1	
				
				
				while(Yuck2[n-p]==1||Yuck2[n-p]==-1)
				if(p==n)
				elseif(Yuck2[n-p]==-2)
					insertpoints numpnts(Up),1,Up
					Up[numpnts(Up)-1]=Yuck[n]
				endif
			endif
		elseif(YUCK2[n]==-1)
			if(YUCK2[n+1]==-2&&Yuck2[n-1]==2)
				insertpoints numpnts(Down),1,Down
				Down[numpnts(Down)-1]=Yuck[n]
				//print "Minus: "+num2str(YUCK[n])
			endif
			if(YUCK2[n+1]==-2&&Yuck2[n-1]==1)
				p=1
				do
					if(n==1||p>=n)
					break
					endif
					p+=1	
				
				
				while(Yuck2[n-p]==1||Yuck2[n-p]==-1||p<n)
			
				if(p==n)
				elseif(Yuck2[n-p]==2)
					insertpoints numpnts(Down),1,Down
					Down[numpnts(Down)-1]=Yuck[n]
					//print "Minus: "+num2str(YUCK[n])

				endif
			endif
		endif
	endfor
	Concatenate/o/NP {Up,Down}, YUCK3
	Sort YUCK3, YUCK3

	duplicate/o YUCK3,Transitions
	variable firstpoint=x2pnt(Inwave,Transitions[0])-1
	variable secondpoint=x2pnt(Inwave,Transitions[0])+1
	variable result
	if(Inwave[firstpoint]<Inwave[secondpoint])
		result=0
	else
		result=1
	endif
	killwaves Yuck,Yuck2,Yuck3,W_FindLevels
	return result
end


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


static function FindTransitions(INputWave,Level,OutWave)
	wave InputWave,Outwave
	variable Level
	duplicate/free InputWave SmoothedWave
	FindLevels/Q  SmoothedWave Level
	wave W_FindLevels
	variable firstpoint=x2pnt(SmoothedWave,W_FindLevels[0])-1
	variable secondpoint=x2pnt(SmoothedWave,W_FindLevels[0])+1
	variable result
	if(SmoothedWave[firstpoint]<SmoothedWave[secondpoint])
	result=0
	else
	result=1
	endif
	
	duplicate/o W_FindLevels OutWave
	return result
end

Static Function GuessDGaussParms(Hist,OutCoefs)
	wave Hist,OutCoefs

	


end

Static Function MakeNicePlot(RawForce,SmForce,HistWave,HistFit,LT1,LT1Fit,LT2,LT2Fit,TransX,TransY)
	wave RawForce,SmForce,HistWave,HistFit,LT1,LT2,LT1Fit,LT2Fit,TransX,TransY
	string WindowName=nameofwave(RawForce)+"_Win"
	dowindow $Windowname
	if(V_flag==1)
		killwindow $windowname
	else
	endif
	Display/N=$WindowName RawForce,SmForce
	AppendToGraph/W=$WindowName TransY vs TransX

	AppendToGraph/W=$WindowName/B=B1/L=L1/VERT HistWave
	AppendToGraph/W=$WindowName/B=B1/L=L1/VERT HistFit
	AppendToGraph/W=$WindowName/B=B2/L=L2 LT1
	AppendToGraph/W=$WindowName/B=B2/L=L2 LT1Fit

	AppendToGraph/W=$WindowName/B=B2/L=L2 LT2
	AppendToGraph/W=$WindowName/B=B2/L=L2 LT2Fit

	//	//SetAxis L1 2.1312133e-08,2.9592798e-08
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
	ModifyGraph mode($nameofwave(TransY))=3;
	ModifyGraph rgb($nameofwave(TransY))=(0,0,0)
	ModifyGraph marker($nameofwave(TransY))=19
	
	DoUpdate 
	GetAxis/W=$WindowName/Q L1
	SetAxis/W=$WindowName left, v_min,v_max
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

Static Function SVP(ctrlName,varNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable varNum
	String varStr
	String varName
	strswitch(ctrlName)

	endswitch
End		

Macro EqJumpPanel() : Panel

	PauseUpdate; Silent 1		// building window...
	NewPanel/N=EqJPanel /W=(600,250,900,500)
	NewDataFolder/o root:DE_EQJ
	NewDataFolder/o root:DE_EQJ:MenuStuff

	Button de_EQJ_button0,pos={5,200},size={50,20},proc=DE_EquilJump#ButtonProc,title="Process!"
	PopupMenu de_EQJ_popup0,pos={5,2},size={129,21},title="Folder"
	PopupMenu de_EQJ_popup0,mode=1,popvalue="X",value= #"DE_EquilJump#ListFolders()"
	PopupMenu de_EQJ_popup1,pos={5,40},size={129,21},title="Force"
	PopupMenu de_EQJ_popup1,mode=1,popvalue="X",value= #"DE_EquilJump#ListWaves()"
	//PopupMenu de_EQJ_popup2,pos={5,80},size={129,21},title="Sep"
	//PopupMenu de_EQJ_popup2,mode=1,popvalue="X",value= #"DE_EquilJump#ListWaves()"
	//PopupMenu de_EQJ_popup3,pos={5,110},size={129,21},title="ZSnsr"
	//PopupMenu de_EQJ_popup3,mode=1,popvalue="X",value= #"DE_EquilJump#ListWaves()"
	SetVariable de_EQJ_setvar0,pos={5,140},size={150,16},proc=DE_EquilJump#SVP,title="Force Offset"
	SetVariable de_EQJ_setvar0,limits={-inf,inf,0},value= _NUM:0
	SetVariable de_EQJ_setvar1,pos={5,170},size={150,16},proc=DE_EquilJump#SVP,title="Sep Offset"
	SetVariable de_EQJ_setvar1,limits={-inf,inf,0},value= _NUM:0
	SetVariable de_EQJ_setvar2,pos={175,170},size={150,16},proc=DE_EquilJump#SVP,title="Smoothing"
	SetVariable de_EQJ_setvar2,limits={-inf,inf,0},value= _NUM:5e-9
//	
	
EndMacro
//

Static Function/S ListWaves()

	String saveDF
	saveDF = GetDataFolder(1)
	controlinfo de_EQJ_popup0
	SetDataFolder s_value
	String list = WaveList("*Equil*", ";", "")
	SetDataFolder saveDF
	return list

end

Static Function/S ListFolders()

	string list=DE_PanelProgs#PrintAllFolders_String("*")
	return list
End

Menu "Equilibrium"
	//SubMenu "Processing"
	"Open EquilJump Panel", EqJumpPanel()


	//end
	
end
