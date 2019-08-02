#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma ModuleName=DE_Dudko
#include ":\Misc_PanelPrograms\Panel Progs"
#include ":\Misc_PanelPrograms\AsylumNaming"
#include "DE_CompileHistograms"
#include "DE_Filtering"


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Dudko

Static Function CalcRates(Var,bins,pare)
	variable Var,bins,pare
	make/T/Free/n=0 StateWaves,RefWaves
	DFREF Here1= $stringfromlist(0, DE_PanelProgs#PrintAllFolders_String("*_2States_Adj"))
	DE_PanelProgs#ListWaves(Here1,"*_2States_Adj",StateWaves)
	print StateWaves
	variable n=0,j	, divj=0,USlopehold,UNumHold

	make/free/n=0 JustForce,StateImIn
	for(n=0;n<numpnts(StateWaves);n+=1)
		wave SW=$StateWaves[n]	
		wave RW=$ReplaceString("2States_adj",StateWaves[n]	,"Ref")
		print nameofwave(SW)
		Make/free/n=(dimsize(SW,0)) Forces,TempState
		Forces= -1*SW[p][1]
		TempState= SW[p][2]

		Extract/o Forces, LocalUnFoldingForces, TempState[p]==-1
		Extract/o Forces, LocalFoldingForces, TempState[p]==1
		make/o/n=(bins) $ReplaceString("States",StateWaves[n]	,"UHist"),$ReplaceString("States",StateWaves[n]	,"FHist")
		
		wave w1=$ReplaceString("States",StateWaves[n]	,"UHist")
		wave w2=$ReplaceString("States",StateWaves[n]	,"FHist")
		Histogram/C/B=1 LocalUnFoldingForces w1
		Histogram/C/B=1 LocalFoldingForces w2

		
		duplicate/o w1 $ReplaceString("States",StateWaves[n]	,"USlope"),$ReplaceString("States",StateWaves[n]	,"UNum"),$ReplaceString("States",StateWaves[n]	,"URates")
		duplicate/o w2 $ReplaceString("States",StateWaves[n]	,"FSlope"),$ReplaceString("States",StateWaves[n]	,"FNum"),$ReplaceString("States",StateWaves[n]	,"FRates")
		wave w3=$ReplaceString("States",StateWaves[n]	,"USlope")
		wave w4=$ReplaceString("States",StateWaves[n]	,"FSlope")
		wave w5=$ReplaceString("States",StateWaves[n]	,"UNum")
		wave w6=$ReplaceString("States",StateWaves[n]	,"FNum")
		wave w7=$ReplaceString("States",StateWaves[n]	,"URates")
		wave w8=$ReplaceString("States",StateWaves[n]	,"FRates")

	
		for(j=0;j<numpnts(w1);j+=1)
			make/free/n=(dimsize(RW,0)) AllForces,UnAllSlopes,UnAllNumbers
			AllForces=real(RW[p][0])
			if(Var==1)
			
				UnAllSlopes=real(RW[p][1])
			else
				UnAllSlopes=real(RW[p][5])

			endif
			UnAllNumbers=real(RW[p][3])

			FindLevel/Q/p AllForces pnt2x(w1,j)
			if(numtype(v_levelx)!=0)
				if(pnt2x(w1,j)>=wavemax(AllForces))
					w5[j]=wavemin(UnAllNumbers)
				elseif(pnt2x(w1,j)<wavemax(AllForces))
					w5[j]=wavemax(UnAllNumbers)
				else
					print "WARNING"
				endif
			elseif(numtype(UnAllSlopes[v_levelx])==2)
				w5[j]=UnAllNumbers[v_levelx]

			else
				//print UnAllSlopes[v_levelx]
				w3[j]=UnAllSlopes[v_levelx]
				w5[j]=UnAllNumbers[v_levelx]
			endif

		endfor
		j=0
		for(j=0;j<numpnts(w2);j+=1)
			make/free/n=(dimsize(RW,0)) AllForces,FAllSlopes,FAllNumbers
			AllForces=real(RW[p][0])
			if(Var==1)
			
				FAllSlopes=real(RW[p][2])
			else
				FAllSlopes=real(RW[p][6])

			endif
			FAllNumbers=real(RW[p][4])

			FindLevel/Q/p AllForces pnt2x(w1,j)
			if(numtype(v_levelx)!=0)
				if(pnt2x(w2,j)>=wavemax(AllForces))
					w6[j]=wavemin(FAllNumbers)
				elseif(pnt2x(w2,j)<wavemax(AllForces))
					w6[j]=wavemax(FAllNumbers)
				else
					print "WARNING"
				endif
			elseif(numtype(FAllSlopes[v_levelx])==2)
				w6[j]=FAllNumbers[v_levelx]

			else
				w4[j]=FAllSlopes[v_levelx]
				w6[j]=FAllNumbers[v_levelx]
			endif

		endfor
		w7=w1/dimdelta(w1,0)*w3/w5
		w8=w2/dimdelta(w2,0)*w4/w6
		duplicate/free JustForce TestConcatH
		Concatenate/o/NP {TestConcatH,Forces}, JustForce
		duplicate/o StateImIn StateImInH
		Concatenate/o/NP {StateImInH,TempState}, StateImIn
	endfor
	
	Extract/o JustForce, UnFoldingForces, StateImIn[p]==-1
	Extract/o JustForce, FoldingForces, StateImIn[p]==1
	
	make/o/n=(bins) FHist,UHist
	Histogram/C/B=1 UnFoldingForces UHist
	Histogram/C/B=1 FoldingForces FHist
	duplicate/o UHist USlopeF,UNumF,URate
	duplicate/o FHist FSlopeF,FNumF,FRate
	
	variable i

	//for(j=1;j<2;j+=1)
	for(j=0;j<numpnts(UHist);j+=1)
		divj=0
		for(n=0;n<numpnts(Statewaves);n+=1)
			wave SW=$StateWaves[n]	
			wave RW=$ReplaceString("2States_Adj",StateWaves[n],"Ref")
			make/free/n=(dimsize(RW,0)) AllForces,UnAllSlopes,UnAllNumbers
			AllForces=real(RW[p][0])
			if(Var==1)
			
				UnAllSlopes=real(RW[p][1])
			else
				UnAllSlopes=real(RW[p][5])

			endif
			UnAllNumbers=real(RW[p][3])

			FindLevel/Q/p AllForces pnt2x(UHist,j)
			if(numtype(v_levelx)!=0)
				if(pnt2x(UHist,j)>=wavemax(AllForces))
					UNumHold+=wavemin(UnAllNumbers)
				elseif(pnt2x(UHist,j)<wavemax(AllForces))
					UNumHold+=wavemax(UnAllNumbers)
				else
					print "WARNING"
				endif
			elseif(numtype(UnAllSlopes[v_levelx])==2)
				UNumHold+=UnAllNumbers[v_levelx]

			else
				//print UnAllSlopes[v_levelx]
				USlopehold+=UnAllSlopes[v_levelx]
				UNumHold+=UnAllNumbers[v_levelx]
				divj+=1

			endif

		endfor
		USlopeF[j]=USlopehold/divj
		UNumF[j]=UNumHold
		USlopehold=0
		UNumHold=0
		divj=0
	endfor
	variable FSlopehold,FNumHold,		divi=0
	for(i=0;i<numpnts(FHist);i+=1)
		for(n=0;n<numpnts(Statewaves);n+=1)
			wave SW=$StateWaves[n]	
			wave RW=$ReplaceString("2States_Adj",StateWaves[n],"Ref")
			make/free/n=(dimsize(RW,0)) AllForces,FAllSlopes,FAllNumbers
			//	make/free/n=(12) USlope,Fslope,UNum,FNum
			AllForces=real(RW[p][0])
			if(Var==1)
			
				FAllSlopes=real(RW[p][2])
			else
				FAllSlopes=real(RW[p][6])

			endif			
			FAllNumbers=real(RW[p][4])
			FindLevel/Q/p AllForces pnt2x(FHist,i)
			if(numtype(v_levelx)!=0)
				if(pnt2x(FHist,i)>=wavemax(AllForces))
					FNumHold+=wavemax(FAllNumbers)
				elseif(pnt2x(FHist,i)<+wavemax(AllForces))
					FNumHold+=wavemin(FAllNumbers)
				else
					print "WARNING"
				endif
			elseif(numtype(FAllSlopes[v_levelx])==2)
				FNumHold+=FAllNumbers[v_levelx]

			else
				//print FAllSlopes[v_levelx]
				FSlopehold+=FAllSlopes[v_levelx]
				FNumHold+=FAllNumbers[v_levelx]
				divi+=1

			endif

		endfor
		FSlopeF[i]=FSlopehold/divi
		FNumF[i]=FNumHold
		FSlopehold=0
		FNumHold=0
		divi=0
	endfor
	
		
	URate=UHist/dimdelta(UHist,0)*USlopeF/UNumF
	FRate=FHist/dimdelta(FHist,0)*FSlopeF/FNumF
	if(pare==0)
	
	else
		for(n=0;n<numpnts(Urate);n+=1)
			if(UHist[n]<pare)
			URate[n]=Nan
			endif
			endfor
	for(n=0;n<numpnts(FRate);n+=1)
			if(FHist[n]<pare)
			FRate[n]=Nan
			endif
			endfor
	endif
	
	killwaves JustForce,StateImIN,StateImInH,LocalUnFoldingForces,LocalFoldingForces
end


static Function MakeSingleStateKey(Forcewave,RuptureWaveLeaving,RuptureWaveEntering,CombinedOut,[smoothing])
	wave Forcewave,RuptureWaveLeaving,RuptureWaveEntering,CombinedOut
	variable smoothing
	variable n,points
	make/free/n=(numpnts(RuptureWaveLeaving)+numpnts(RuptureWaveEntering),2) TestCombined
	make/free/n=0 PauseWave,TestWave,Combined,CombinedV
	DE_Dudko#grabIndices(ForceWave,TestWave)
	DE_Dudko#grabPauseIndices(ForceWave,PauseWave)
	make/free/n=0 SurfStop,ExtStop
	
	duplicate/free ForceWave Fw

	if( ParamIsDefault(Smoothing) )
	elseif(Smoothing<5)
		print "Smoothing <5 not allowed"
	elseif(2*round((smoothing)/2)==smoothing)
		Smoothing=2*round((smoothing)/2)+1
		print "Smoothing moved to nearest Odd number "+num2str(smoothing)
		Smooth/S=2 Smoothing, Fw
		Fw*=-1

	else
		
		Smooth/S=2 Smoothing, Fw
		Fw*=-1

	endif
	
	for(n=0;n<numpnts(TestWave);n+=1)
		if (mod(n, 2) == 0)
			points=n/2
			insertpoints points,1, SurfStop
			SurfStop[points]=TestWave[n]
		else
			points=(n-1)/2
			insertpoints points,1, ExtStop
			ExtStop[points]=TestWave[n]
		endif
	endfor
	duplicate/free SurfStop SurfStopV
	SurfStopV=2
	duplicate/free ExtStop ExtStopV
	ExtStopV=-2
	duplicate/free PauseWave PauseWaveV
	PauseWaveV=0
	duplicate/free RuptureWaveLeaving RuptureWaveLeavingV
	RuptureWaveLeavingV=-1
	duplicate/free RuptureWaveEntering RuptureWaveEnteringV
	RuptureWaveEnteringV=1
	Concatenate/NP {PauseWave,SurfStop,ExtStop,RuptureWaveLeaving,RuptureWaveEntering}, Combined
	Concatenate/NP {PauseWaveV,SurfStopV,ExtStopV,RuptureWaveLeavingV,RuptureWaveEnteringV},CombinedV
	Sort Combined,CombinedV
	Sort Combined,Combined
	make/free/n=(numpnts(Combined),5) AboutToFinish
	if(Combined[dimsize(Combined,0)-1]==numpnts(FW))
	Combined[dimsize(Combined,0)-1]-=1
	endif
	AboutToFinish[][0]=Combined[p]
	AboutToFinish[][4]=pnt2x(FW,AboutToFinish[p][0])
	AboutToFinish[][1]=FW[AboutToFinish[p][0]]
	AboutToFinish[][2]=CombinedV[p]
	print/D AboutToFinish[1][0]
	variable States=-1
	for(n=0;n<dimsize(Combined,0);n+=1)
		if(AboutToFinish[n][2]==2)
		States+=1
		endif
		AboutToFinish[n][3]=States
	
	endfor

	note AboutToFinish, nameofwave(Forcewave)

	duplicate/o AboutToFinish CombinedOut
	wave '_free_'
	killwaves '_free_'
End



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Slope Calculations

Static Function/C NewAvgSlope(ForceWaveSM,CombinedWave,Force,State,Smoothed,TimeStep,minspacing)
	wave ForceWaveSM, CombinedWave
	variable Force,State,Smoothed,TimeStep,minspacing
	variable n,FoldedMax,UnfoldedMax,q
	String StateList=""
	Variable StartTime,EndTime,EventBefore,EventAfter
	If(Smoothed==1)
		duplicate/free ForceWaveSM FWAll
	else
		duplicate/free ForceWaveSM FWAll
		Smooth/S=2 101, FWALL
		FWALL*=-1
	endif
	
	make/free/n=(dimsize(CombinedWave,0)) Points,RupForce,Type,Trace
	make/free/n=0 HoldOnToSlopes,Folded,Unfolded
	Points=CombinedWave[p][0]
	RupForce=CombinedWave[p][1]
	Type=CombinedWave[p][2]
	Trace=CombinedWave[p][3]
	variable lastpoint=0,firstone,firstmone
	for(n=0;n<numpnts(points);n+=1)
		//for(n=0;n<5;n+=1)
	
		 
		if(Type[n]==-1)
			if((Points[n]-lastpoint)<=0)

			else
				if(firstmone==0)

					duplicate/free/r=[0,Points[n]] FWAll, Folded
					firstmone=1
				else		
	
					duplicate/free/r=[lastpoint,Points[n]] FWAll, NewSection
					AppendToWaveWOverlapOW(Folded,NewSection,lastpoint)
				endif
			endif
			
			lastpoint=Points[n]+minspacing
		endif	
	
		if(Type[n]==1)
			if((Points[n]-lastpoint)<=0)
			else
				if(firstone==0)
					duplicate/free/r=[0,Points[n]] FWAll, UnFolded
					UnFolded[0,lastpoint]=NaN
					firstone=1
				else		

					duplicate/free/r=[lastpoint,Points[n]] FWAll, NewSection
					AppendToWaveWOverlapOW(UnFolded,NewSection,lastpoint)
				endif
				
			endif
			lastpoint=Points[n]+minspacing
		endif
		
	endfor
	
	if(State==-1)
		FindLevels/M=0.1/Q Folded, Force
	elseif(State==1)
		FindLevels/M=0.1/Q UnFolded, Force

	endif
	wave W_FindLevels

	if(numpnts(W_FINDLEVELS)==0)
	
	else
		for(q=0;q<numpnts(W_FINDLEVELS);q+=1)
			//for(q=7;q<8;q+=1)
		
			if(State==-1)
				StartTime=W_FindLevels[q]-TimeStep/2
				EndTime=W_FindLevels[q]+TimeStep/2
				duplicate/free/r=(StartTime,EndTime) Folded FWFit
			elseif(State==1)
				StartTime=W_FindLevels[q]-TimeStep/2
				EndTime=W_FindLevels[q]+TimeStep/2
				duplicate/free/r=(StartTime,EndTime) UnFolded FWFit
			endif
			CurveFit/Q/NTHR=0 line  FWFit 
			wave w_coef
			if((Endtime-Starttime)<0.001)
			else
				insertpoints 0,1, HoldOnToSlopes
				HoldOnToSlopes[0]=abs(w_coef[1])
			endif
		endfor
	endif

 
	if( numpnts(Holdontoslopes)==0)
		return NaN
	endif
	wavestats/q HoldOnToSlopes
	//print v_Sdev
	wave W_Coef,W_Sigma
	killwaves W_FindLevels,W_Coef,W_Sigma
	return cmplx(v_avg,v_sdev)
	
end


Static Function/C ReturnSeparationSlopes(SepWaveSm,CombinedWave,MinSpacing)
	wave SepWaveSm,CombinedWave
	variable MinSpacing
	make/free/n=(dimsize(CombinedWave,0)) Points,RupForce,Type,Trace
	make/free/n=0 HoldontoFoldedSlope,HoldontoUnFoldedSlope
	make/free/n=0 UnfoldedSep,FoldedSep
	Points=CombinedWave[p][0]
	variable lastpoint=0,firstone,firstmone,NewZero,firstfit

	Type=CombinedWave[p][2]
	Trace=CombinedWave[p][3]
	variable i,n
	for(n=0;n<Trace[numpnts(Trace)-1];n+=1)
		Extract/INDX/FREE Points, LocalSet, Trace==n
		make/free/n=(numpnts(LocalSet)) LocalPoints,LocalType,LocalTrace
		Extract/INDX/FREE Points, LocalSet, Trace==n
		make/free/n=(numpnts(LocalSet)) LocalPoints,LocalType,LocalTrace
		LocalPoints[]=Points[LocalSet[1]+p][0]
		LocalType[]=Type[LocalSet[1]+p][2]
		LocalTrace[]=Trace[LocalSet[1]+p][2]
		for(i=0;i<numpnts(LocalPoints);i+=1)
			if(LocalType[i]==2||LocalType[i]==-2)
				if((LocalPoints[i]-lastpoint)<=0)

				else
					if(LocalType[i]==2)
						if(firstmone==0)
							duplicate/free/r=[lastpoint,LocalPoints[i]] SepWaveSm, FoldedSep
							firstmone=1
						else		
							duplicate/free/r=[lastpoint,LocalPoints[i]] SepWaveSm, NewSectionSep
							AppendToWaveWOverlap(FoldedSep,NewSectionSep)
						endif
					else
						if(firstone==0)

							duplicate/free/r=[lastpoint,LocalPoints[i]] SepWaveSm, unFoldedSep
							firstone=1
						else		
							duplicate/free/r=[lastpoint,LocalPoints[i]] SepWaveSm, NewSectionSep
							AppendToWaveWOverlap(UnFoldedSep,NewSectionSep)
						endif
					endif
			
				endif
				
				if(numpnts(FoldedSep)<30||numpnts(UnFoldedSep)<30)
					print "SKIPPED: "+num2str(n)
					deletepoints 0,1e8,  FoldedSep,UnfoldedSep
					NewZero=Points[n]
					firstmone=0
					firstone=0
				else
					lastpoint=Points[n]+minspacing
					CurveFit/W=2/Q/NTHR=0 line  FoldedSep 
					wave w_Coef
					insertpoints 0,1, HoldontoFoldedSlope,HoldontoUnFoldedSlope
					if(LocalType[i]==-2)
						HoldontoFoldedSlope[0]= w_coef[1]
					else
						HoldontoFoldedSlope[0]= -1*w_coef[1]
					endif
					CurveFit/W=2/Q/NTHR=0 line  UnfoldedSep 
					if(LocalType[i]==-2)
						HoldontoUnFoldedSlope[0]= w_coef[1]
					else
						HoldontoUnFoldedSlope[0]= -1*w_coef[1]
					endif
					firstfit+=1
					deletepoints 0,1e8, FoldedSep,UnfoldedSep
					NewZero=LocalPoints[i]
					lastpoint=LocalPoints[i]
					firstmone=0
					firstone=0
				endif	
			
			elseif(LocalType[i]==-1)
				if((LocalPoints[i]-lastpoint)<=0)
				else
					if(firstmone==0)
						duplicate/free/r=[NewZero,LocalPoints[i]] SepWaveSm, FoldedSep
						firstmone=1
					else		
						duplicate/free/r=[lastpoint,LocalPoints[i]] SepWaveSm, NewSectionSep
						AppendToWaveWOverlap(FoldedSep,NewSectionSep)
					endif
				endif
				lastpoint=LocalPoints[i]+minspacing
			elseif(LocalType[i]==1)
				if((LocalPoints[i]-lastpoint)<=0)

				else
					if(firstone==0)
						duplicate/free/r=[lastpoint,LocalPoints[i]] SepWaveSm, UnFoldedSep
						firstone=1
					else		
						duplicate/free/r=[lastpoint,LocalPoints[i]] SepWaveSm, NewSectionSep
						AppendToWaveWOverlap(UnFoldedSep,NewSectionSep)
					endif
				endif
				lastpoint=LocalPoints[i]+minspacing
			endif
		endfor
	endfor	
	wavestats/Q HoldontoFoldedSlope
	variable/C result=cmplx(v_avg,0)
	wavestats/Q HoldontoUnFoldedSlope
	result+=cmplx(0,v_avg)
	return result

end

Static Function FitFreeWLCSlope(LC,Offset,Force,rate)
	variable LC,Offset,Force,rate
	variable z=DE_WLC#ReturnExtentionatForce(Force+Offset,.4e-9,LC,298)
	return WLCSlope(z,.4e-9,LC,298)*rate
end

Static Function WLCSlopeTimeDomain(WLCWave,Force)
	wave WLCWave
	variable Force
	FindLevels/M=0.05 WLCWave, Force
	wave W_FindLevels
	duplicate/o W_FindLevels Here1,Here2
	Here2=WLCWave(Here1)
	make/o/n=(numpnts(Here1)) BLAG
	variable distance=.1
	variable n
	for(n=0;n<numpnts(Here1);n+=1)
		duplicate/o/R=(Here1[n]-distance,Here1[n]+distance) WLCWave CutBlag
		CurveFit/Q/W=2/NTHR=0 line  CutBlag /D
		wave w_coef
		BLAG[n]=abs(w_coef[1])
		
	endfor
	return Mean(Blag)
end

Static Function MakeWLCFvT(WLCParms,FSep,USep,OutWLC1,OutWLC2)
	wave WLCParms,FSep,USep,,OutWLC1,OutWLC2
	duplicate/free FSep WLC1
	duplicate/free USep WLC2
	WLC1=WLC(FSep-WLCParms[3],.4e-9,WLCParms[0],298)
	WLC2=WLC(USep-WLCParms[3],.4e-9,WLCParms[1],298)
	duplicate/o WLC1 OutWLC1
	duplicate/o WLC2 OutWLC2
end

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Number of Traces in a State

Static Function NumberofTracesinState(ForceWaveSM,CombinedWave,Force,State,Smoothed)
	wave ForceWaveSM, CombinedWave
	variable Force,State,Smoothed
	variable n,FoldedMax,UnfoldedMax,q
	String StateList=""
	
	If(Smoothed==1)
		duplicate/free ForceWaveSM FWAll
	else
		duplicate/free ForceWaveSM FWAll
		Smooth/S=2 101, FWALL
		FWALL*=-1
	endif
	
	make/free/n=(dimsize(CombinedWave,0)) Points,RupForce,Type,Trace
	make/T/free/n=0 HoldOnToStates
	Points=CombinedWave[p][0]
	RupForce=-CombinedWave[p][1]
	Type=CombinedWave[p][2]
	Trace=CombinedWave[p][3]
	//for(n=10;n<11;n+=1)
	for(n=0;n<Trace[numpnts(Trace)-1][2];n+=1)
		Extract/INDX/Free Points, LocalSet, Trace==n
		make/free/n=(numpnts(LocalSet)+1) LocalPoints,LocalType,LocalTrace
		LocalPoints[]=Points[LocalSet[0]+p][0]
		LocalType[]=Type[LocalSet[0]+p][2]
		LocalTrace[]=Trace[LocalSet[0]+p][2]
		FindValue/V=-2/T=.01 LocalType
		duplicate/free/r=[,V_value] LocalPoints UpDir
		duplicate/free/r=[V_value+1,] LocalPoints DownDir
		duplicate/free/r=[UpDir[0],UpDir[numpnts(UpDir)-1]] FWALL FW
		FindLevels/Q FW,Force		
		wave W_FindLevels

		if(numpnts(W_FindLevels)==0)
			wavestats/q FW
			if(v_max<Force)
				StateList+="1"
			elseif(Force<v_min)
				StateList+="-1"

			else
			endif
		else
			for(q=0;q<numpnts(W_FindLevels);q+=1)
				StateList+=num2str(StateAtPoint(x2pnt(FWALL,W_FindLevels[q]),CombinedWave))+";"
			endfor
		endif
		if(FindListItem(num2str(State), StateList)==-1)
		
		else
			InsertPoints 0, 1,HoldOnToStates
			HoldOnToStates[0]=num2str(LocalTrace[0])+"U"
		endif
		StateList=""
		
		
		duplicate/free/r=[UpDir[numpnts(UpDir)-1],DownDir[numpnts(DownDir)-1]] FWALL FW
		
		FindLevels/Q FW,Force		
		wave W_FindLevels

		if(numpnts(W_FindLevels)==0)
			wavestats/Q FW
			if(v_max<Force)
				StateList+="1"
			elseif(Force<v_min)
				StateList+="-1"

			else
			endif
		else

			for(q=0;q<numpnts(W_FindLevels);q+=1)
				StateList+=num2str(StateAtPoint(x2pnt(FWALL,W_FindLevels[q]),CombinedWave))+";"
			endfor
	
			endif
			if(FindListItem(num2str(State), StateList)==-1)
			
		else
			InsertPoints 0, 1,HoldOnToStates
			HoldOnToStates[0]=num2str(LocalTrace[0])+"D"
		endif
		StateList=""
		
	endfor
	killwaves W_FindLevels
	
	return numpnts(HoldOnToStates)
end

Static Function GenerateFakeSeparation(SepIn,States,SepOut)

	wave SepIn,States,SepOut
	duplicate/free SepIn Assembled

	
	variable n,prevcut=0,currcut
	for(n=1;n<dimsize(States,0);n+=1)
		if(States[n][2]==2||States[n][2]==-2)
			currcut=States[n][0]
			duplicate/free/r=[prevcut,currcut] Sepin TestSep
			
			CurveFit/Q/W=2/NTHR=0 line  TestSep/D
			wave w_coef
			
			Assembled[prevcut,min(currcut,numpnts(Assembled)-1)]=w_coef[0]+w_coef[1]*x
			prevcut=States[n][0]
			
		
		
		endif
	
	endfor
	duplicate/o Assembled SepOut

end


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//WLC Fitting

Static Function MakeWLCs(SepWave,CoefficientWave,OutWave1,Outwave2)
	wave SepWave,CoefficientWave,OutWave1,Outwave2
	make/free/n=1000 WLC1,WLC2
	WAVESTATS/Q sEPwAVE
	SetScale/P x v_Min,(v_max-v_min)/999, WLC1,WLC2
	WLC1=WLC(x-CoefficientWave[3],0.4e-9,CoefficientWave[0],298)-CoefficientWave[2]
	WLC2=WLC(x-CoefficientWave[3],0.4e-9,CoefficientWave[1],298)-CoefficientWave[2]
	duplicate/o WLC1 Outwave1
	duplicate/o WLC2 Outwave2


end

Static Function MakeWLCs2(SepWave,CoefficientWave,OutWave1,Outwave2)
	wave SepWave,CoefficientWave,OutWave1,Outwave2
	print CoefficientWave
	make/free/n=1000 WLC1,WLC2
	WAVESTATS/Q sEPwAVE
	SetScale/P x v_Min,(v_max-v_min)/999, WLC1,WLC2
	WLC1=WLC(x-CoefficientWave[5],CoefficientWave[0],CoefficientWave[1],CoefficientWave[3])-CoefficientWave[4]
	WLC2=WLC(x-CoefficientWave[5],CoefficientWave[0],CoefficientWave[2],CoefficientWave[3])-CoefficientWave[4]
	duplicate/o WLC1 Outwave1
	duplicate/o WLC2 Outwave2


end


Static Function TMakeWLCs(SepWave,CoefficientWave,OutWave1,Outwave2)
	wave SepWave,CoefficientWave,OutWave1,Outwave2
	print	CoefficientWave

	make/free/n=1000 WLC1,WLC2
	WAVESTATS/Q sEPwAVE
	SetScale/P x v_Min,(v_max-v_min)/999, WLC1,WLC2
	WLC1=WLC(x-CoefficientWave[5],0.4e-9,CoefficientWave[1],298)-CoefficientWave[4]
	WLC2=WLC(x-CoefficientWave[5],0.4e-9,CoefficientWave[2],298)-CoefficientWave[4]
	duplicate/o WLC1 Outwave1
	duplicate/o WLC2 Outwave2


end

Static Function/C ContourLengthDetermine(ForceWaveSM,SepWaveSm,CombinedWave,MinSpacing,ResOut)
	wave ForceWaveSM, SepWaveSm,CombinedWave,ResOut
	variable MinSpacing
	variable n,FoldedMax,UnfoldedMax,q
	
	make/free/n=(dimsize(CombinedWave,0)) Points,RupForce,Type,Trace
	make/free/n=0 HoldOnToCF,HoldOnToCU,HoldOnToF,HoldOnToOff,Unfolded,UnfoldedSep,Folded,FoldedSep,HoldN
	Points=CombinedWave[p][0]
	RupForce=CombinedWave[p][1]
	Type=CombinedWave[p][2]
	Trace=CombinedWave[p][3]
	
	variable fitstart=SepWaveSm[0]-3e-9,firstfit,FOFF
	variable lastpoint=0,firstone,firstmone,NewZero
	for(n=0;n<numpnts(points);n+=1)
		if(Type[n]==2)
			if(n==0)
			
			elseif(numpnts(Folded)<10||numpnts(UnFolded)<10)
				print "SKIPPED: "+num2str(n)
				deletepoints 0,1e8,  Folded,Unfolded,FoldedSep,UnfoldedSep
				NewZero=Points[n]
				firstmone=0
				firstone=0
			else

				if((Points[n]-lastpoint)<=0)

				else
					duplicate/free/r=[lastpoint,Points[n]] ForceWaveSM, NewSection
					duplicate/free/r=[lastpoint,Points[n]] SepWaveSm, NewSectionSep
					AppendToWaveWOverlap(Folded,NewSection)
					AppendToWaveWOverlap(FoldedSep,NewSectionSep)

			
				endif
				lastpoint=Points[n]+minspacing
				Make/D/N=5/O W_coef
				//FOFF=ReturnOffForceForState(ForceWaveSM,Trace[n-1])
				//print FOFF
				//Unfolded=Unfolded+FOFF
				//folded=folded+FOFF
				if(Firstfit==0)
					W_coef[0] = {.4e-9,100e-9,298,0,fitstart}
					FuncFit/Q/H="10110"/NTHR=0 WLC_FIT W_coef  Unfolded /X=UnFoldedSep
					fitstart=w_coef[4]
					
				else
					W_coef[0] = {.4e-9,100e-9,298,0,fitstart}
					FuncFit/Q/H="10111"/NTHR=0 WLC_FIT W_coef  Unfolded /X=UnFoldedSep
				endif
				
				InsertPoints 0,1, HoldOnToCU,HoldOnToCF,HoldN,HoldOnToF,HoldOnToOff

				HoldOnToF[0]=w_coef[3]
				HoldN[0]=n
				HoldOnToCU[0]=w_coef[1]
				Make/D/N=5/O W_coef
				W_coef[0] = {.4e-9,100e-9,298,w_coef[3],w_coef[4]}
				FuncFit/Q/H="10111"/NTHR=0 WLC_FIT W_coef  folded /X=FoldedSep
				HoldOnToCF[0]=w_coef[1]
				HoldOnToOff[0]=w_coef[4]
				//if(firstfit==0)
				//return 0
				//endif
				firstfit+=1
				deletepoints 0,1e8, Folded,Unfolded,FoldedSep,UnfoldedSep
				NewZero=Points[n]
				firstmone=0
				firstone=0
			endif
		
		endif

		if(Type[n]==-1)
			if((Points[n]-lastpoint)<=0)

			else
				if(firstmone==0)

					duplicate/free/r=[NewZero,Points[n]] ForceWaveSM, Folded
					duplicate/free/r=[NewZero,Points[n]] SepWaveSm, FoldedSep

					firstmone=1
				else		
	
					duplicate/free/r=[lastpoint,Points[n]] ForceWaveSM, NewSection
					//NewSection-=ForceShift
					duplicate/free/r=[lastpoint,Points[n]] SepWaveSm, NewSectionSep

					AppendToWaveWOverlap(Folded,NewSection)
					AppendToWaveWOverlap(FoldedSep,NewSectionSep)

				endif
			endif
			
			lastpoint=Points[n]+minspacing
		endif	
	
		if(Type[n]==1)
			if((Points[n]-lastpoint)<=0)
			else
				if(firstone==0)
					duplicate/free/r=[lastpoint,Points[n]] ForceWaveSM, UnFolded
					//UnFolded[0,lastpoint-NewZero]=NaN
					duplicate/free/r=[lastpoint,Points[n]] SepWaveSm, UnFoldedSep
					//UnFoldedSep[0,lastpoint-NewZero]=NaN
					firstone=1
				else		

					duplicate/free/r=[lastpoint,Points[n]] ForceWaveSM, NewSection

					duplicate/free/r=[lastpoint,Points[n]] SepWaveSm, NewSectionSep

					AppendToWaveWOverlap(UnFolded,NewSection)
					AppendToWaveWOverlap(UnFoldedSep,NewSectionSep)

				endif
				
			endif
			lastpoint=Points[n]+minspacing
		endif
		
	endfor
	make/o/n=(dimsize(HoldOntoCF,0),4) CombinedWLCParms
	
	make/free/n=4 HereYouAre
	wavestats/Q HoldOnToCF
	CombinedWLCParms[][0]=HoldOnToCF[p]
	HereYouAre[0]=v_avg

	wavestats/Q HoldOnToCU
	CombinedWLCParms[][1]=HoldOnToCU[p]

	HereYouAre[1]=v_avg
	wavestats/Q HoldOnToF
	CombinedWLCParms[][2]=HoldOnToF[p]

	HereYouAre[2]=v_avg
	wave w_coef,W_Sigma 

	wavestats/Q HoldOnToOff
	CombinedWLCParms[][3]=HoldOnToOff[p]

	HereYouAre[3]=v_avg

	killwaves w_coef,W_Sigma 
	duplicate/o HereYouAre ResOut
end

Static Function/C ContourLengthDetermineCombined(ForceWaveSM,SepWaveSm,CombinedWave,MinSpacing,ResOut,FitFOldedFirst,[CopyWavesOUt])
	wave ForceWaveSM, SepWaveSm,CombinedWave,ResOut
	variable MinSpacing,CopyWavesOUt,FitFOldedFirst
	variable n,FoldedMax,UnfoldedMax,q
	
	make/free/n=(dimsize(CombinedWave,0)) Points,RupForce,Type,Trace
	make/free/n=0 HoldOnToCF,HoldOnToCU,HoldOnToF,HoldOnToOff,Unfolded,UnfoldedSep,Folded,FoldedSep,HoldN
	Points=CombinedWave[p][0]
	RupForce=CombinedWave[p][1]
	Type=CombinedWave[p][2]
	Trace=CombinedWave[p][3]
	
	variable fitstart=SepWaveSm[0]-10e-9,firstfit,FOFF
	//fitstart=6.12917e-07
	variable lastpoint=0,firstone,firstmone,NewZero
	variable LastNType
	//for(n=1;n<50;n+=1)

	for(n=1;n<numpnts(points);n+=1)
		
	
	
		if(Type[n]==-1)
			if((Points[n]-lastpoint)<=0)
				lastpoint=Points[n]+minspacing
			else
				if(firstmone==0)

					duplicate/free/r=[0,Points[n]] ForceWaveSM, Folded
					duplicate/free/r=[0,Points[n]] SepWaveSm, FoldedSep
					LastNType=Type[n]

					firstmone=1
				else		
	
					duplicate/free/r=[lastpoint,Points[n]] ForceWaveSM, NewSection
					//NewSection-=ForceShift
					duplicate/free/r=[lastpoint,Points[n]] SepWaveSm, NewSectionSep
					LastNType=Type[n]

					AppendToWaveWOverlap(Folded,NewSection)
					AppendToWaveWOverlap(FoldedSep,NewSectionSep)

				endif
			endif
			
			lastpoint=Points[n]+minspacing
		endif	
	
		if(Type[n]==1)
			if((Points[n]-lastpoint)<=0)
				lastpoint=Points[n]+minspacing


			else
				if(firstone==0)
					duplicate/free/r=[lastpoint,Points[n]] ForceWaveSM, UnFolded
					//UnFolded[0,lastpoint-NewZero]=NaN
					duplicate/free/r=[lastpoint,Points[n]] SepWaveSm, UnFoldedSep
					//UnFoldedSep[0,lastpoint-NewZero]=NaN
					firstone=1
					LastNType=Type[n]

				else		

					duplicate/free/r=[lastpoint,Points[n]] ForceWaveSM, NewSection

					duplicate/free/r=[lastpoint,Points[n]] SepWaveSm, NewSectionSep
					LastNType=Type[n]

					AppendToWaveWOverlap(UnFolded,NewSection)
					AppendToWaveWOverlap(UnFoldedSep,NewSectionSep)

				endif

			endif
			lastpoint=Points[n]+minspacing
		endif

	endfor
	if((Points[n-1]-lastpoint)<=0)

	else
		duplicate/free/r=[lastpoint,Points[n-1]] ForceWaveSM, NewSection
		duplicate/free/r=[lastpoint,Points[n-1]] SepWaveSm, NewSectionSep
		AppendToWaveWOverlap(Folded,NewSection)
		AppendToWaveWOverlap(FoldedSep,NewSectionSep)
	endif
	Make/D/N=5/O W_coef
	//FOFF=ReturnOffForceForState(ForceWaveSM,Trace[n-1])
	//print FOFF
	//Unfolded=Unfolded+FOFF
	//folded=folded+FOFF


	W_coef[0] = {.4e-9,100e-9,298,0,fitstart}
	Make/free/T/N=2 T_Constraints
	//T_Constraints[0] = {"K4<"+num2str(fitstart+3e-9),"K4>"+num2str(fitstart-30e-9),"K3<"+num2str(1e-12),"K3>"+num2str(-1e-12)}
	T_Constraints[0] = {"K4<"+num2str(fitstart+10e-9),"K4>"+num2str(fitstart-30e-9)}
	if(FitFOldedFirst==1)
		FuncFit/Q/H="10100"/NTHR=0 WLC_FIT W_coef  folded /X=FoldedSep/C=T_Constraints 
		print w_coef

		fitstart=w_coef[4]
				
		InsertPoints 0,1, HoldOnToCU,HoldOnToCF,HoldN,HoldOnToF,HoldOnToOff

		HoldOnToF[0]=w_coef[3]
		HoldN[0]=n
		HoldOnToCF[0]=w_coef[1]
		W_coef[0] = {.4e-9,w_coef[1]+40e-9,298,w_coef[3],w_coef[4]}
		FuncFit/Q/H="10111"/NTHR=0 WLC_FIT W_coef  unfolded /X=unFoldedSep
		HoldOnToCU[0]=w_coef[1]
		HoldOnToOff[0]=w_coef[4]
		print w_coef

	else
		FuncFit/Q/H="10100"/NTHR=0 WLC_FIT W_coef  unfolded /X=unFoldedSep/C=T_Constraints 
		print w_coef

		fitstart=w_coef[4]
				
		InsertPoints 0,1, HoldOnToCU,HoldOnToCF,HoldN,HoldOnToF,HoldOnToOff

		HoldOnToF[0]=w_coef[3]
		HoldN[0]=n
		HoldOnToCU[0]=w_coef[1]
		W_coef[0] = {.4e-9,w_coef[1]-10e-9,298,w_coef[3],w_coef[4]}
		FuncFit/Q/H="10111"/NTHR=0 WLC_FIT W_coef  folded /X=FoldedSep
		HoldOnToCF[0]=w_coef[1]
		HoldOnToOff[0]=w_coef[4]
		print w_coef
	endif

			
		
	//
	make/free/n=4 HereYouAre
	wavestats/Q HoldOnToCU
	HereYouAre[0]=v_avg

	wavestats/Q HoldOnToCF
	HereYouAre[1]=v_avg
	wavestats/Q HoldOnToF
	HereYouAre[2]=v_avg
	wave w_coef,W_Sigma 

	wavestats/Q HoldOnToOff
	HereYouAre[3]=v_avg

	killwaves w_coef,W_Sigma 
	if(ParamisDefault(CopyWavesOUt)||CopyWavesOUt==0)
	
	else
		duplicate/o FoldedSep FoldedSepa
		duplicate/o Folded	Foldeda
		duplicate/o UnfoldedSep UnfoldedSepa
		duplicate/o Unfolded Unfoldeda
	
	endif
	duplicate/o HereYouAre ResOut
end


Static Function/C MaxOverlapwrtWLC(ForceWaveSM,SepWaveSm,CombinedWave,MinSpacing,WLCParms,FoldedorUn,SepShift,Shift,[CopyWavesOUt])
	wave ForceWaveSM, SepWaveSm,CombinedWave,Shift,WLCParms
	variable MinSpacing,FoldedorUn,CopyWavesOUt,SepShift
	variable n,FoldedMax,UnfoldedMax,q
	
	make/free/n=(dimsize(CombinedWave,0)) Points,RupForce,Type,Trace
	make/free/n=0 Unfolded,UnfoldedSep,Folded,FoldedSep,HoldN
	Points=CombinedWave[p][0]
	RupForce=CombinedWave[p][1]
	Type=CombinedWave[p][2]
	Trace=CombinedWave[p][3]
	String InitialHolds
	if(SepShift==1)
	InitialHolds="11100"
	else
	InitialHolds="11101"
	endif
	//	variable fitstart=SepWaveSm[0]-3e-9,firstfit,FOFF
	//	fitstart=6.12917e-07
	variable lastpoint=0,firstone,firstmone,NewZero
	variable LastNType
	//for(n=1;n<50;n+=1)

	for(n=1;n<numpnts(points);n+=1)
		
	
	
		if(Type[n]==-1)
			if((Points[n]-lastpoint)<=0)
				lastpoint=Points[n]+minspacing
			else
				if(firstmone==0)

					duplicate/free/r=[0,Points[n]] ForceWaveSM, Folded
					duplicate/free/r=[0,Points[n]] SepWaveSm, FoldedSep
					LastNType=Type[n]

					firstmone=1
				else		
	
					duplicate/free/r=[lastpoint,Points[n]] ForceWaveSM, NewSection
					//NewSection-=ForceShift
					duplicate/free/r=[lastpoint,Points[n]] SepWaveSm, NewSectionSep
					LastNType=Type[n]

					AppendToWaveWOverlap(Folded,NewSection)
					AppendToWaveWOverlap(FoldedSep,NewSectionSep)

				endif
			endif
			
			lastpoint=Points[n]+minspacing
		endif	
	
		if(Type[n]==1)
			if((Points[n]-lastpoint)<=0)
				lastpoint=Points[n]+minspacing


			else
				if(firstone==0)
					duplicate/free/r=[lastpoint,Points[n]] ForceWaveSM, UnFolded
					//UnFolded[0,lastpoint-NewZero]=NaN
					duplicate/free/r=[lastpoint,Points[n]] SepWaveSm, UnFoldedSep
					//UnFoldedSep[0,lastpoint-NewZero]=NaN
					firstone=1
					LastNType=Type[n]

				else		

					duplicate/free/r=[lastpoint,Points[n]] ForceWaveSM, NewSection

					duplicate/free/r=[lastpoint,Points[n]] SepWaveSm, NewSectionSep
					LastNType=Type[n]

					AppendToWaveWOverlap(UnFolded,NewSection)
					AppendToWaveWOverlap(UnFoldedSep,NewSectionSep)

				endif

			endif
			lastpoint=Points[n]+minspacing
		endif

	endfor
	if((Points[n-1]-lastpoint)<=0)

	else
		duplicate/free/r=[lastpoint,Points[n-1]] ForceWaveSM, NewSection
		duplicate/free/r=[lastpoint,Points[n-1]] SepWaveSm, NewSectionSep
		AppendToWaveWOverlap(Folded,NewSection)
		AppendToWaveWOverlap(FoldedSep,NewSectionSep)
	endif
	Make/D/N=5/O W_coef
	//FOFF=ReturnOffForceForState(ForceWaveSM,Trace[n-1])
	//print FOFF
	//Unfolded=Unfolded+FOFF
	//folded=folded+FOFF


	W_coef= WLCParms
	make/free/n=6 HereYouAre

	Make/free/T/N=2 T_Constraints
	//T_Constraints[0] = {"K4<"+num2str(fitstart+3e-9),"K4>"+num2str(fitstart-30e-9),"K3<"+num2str(1e-12),"K3>"+num2str(-1e-12)}
	//T_Constraints[0] = {"K4<"+num2str(W_coef[4]+1e-9),"K4>"+num2str(W_coef[4]-1e-9)}
	
	if(FoldedorUn==0)
	FuncFit/Q/H=InitialHolds/NTHR=0 WLC_FIT W_coef  unfolded /X=unFoldedSep///C=T_Constraints 
	HereYouAre[0]=w_coef[0]
	HereYouAre[2]=w_coef[1]
	HereYouAre[3]=w_coef[2]
	//Make/D/N=5/O W_coef
	W_coef[0] = {.4e-9,w_coef[1]*2,298,w_coef[3],w_coef[4]}
//	FuncFit/Q/H="10111"/NTHR=0 WLC_FIT W_coef  folded /X=FoldedSep
	HereYouAre[1]=w_coef[1]
	HereYouAre[4]=w_coef[3]
	HereYouAre[5]=w_coef[4]
	else
	FuncFit/Q/H=InitialHolds/NTHR=0 WLC_FIT W_coef  folded /X=FoldedSep///C=T_Constraints 
	HereYouAre[0]=w_coef[0]
	HereYouAre[1]=w_coef[1]
	HereYouAre[3]=w_coef[2]
	//Make/D/N=5/O W_coef
	W_coef[0] = {.4e-9,w_coef[1]*2,298,w_coef[3],w_coef[4]}
//	FuncFit/Q/H="10111"/NTHR=0 WLC_FIT W_coef  unfolded /X=unFoldedSep
	HereYouAre[2]=w_coef[1]
	HereYouAre[4]=w_coef[3]
	HereYouAre[5]=w_coef[4]
	endif

	wave w_coef,W_Sigma 
	killwaves w_coef,W_Sigma 
	if(ParamisDefault(CopyWavesOUt)||CopyWavesOUt==0)
	
	else
		duplicate/o FoldedSep FoldedSepa
		duplicate/o Folded	Foldeda
		duplicate/o UnfoldedSep UnfoldedSepa
		duplicate/o Unfolded Unfoldeda
	
	endif

	duplicate/o HereYouAre Shift

end


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Dudko Panel Stuff
Static Function PopMenuProc(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa

	switch( pa.eventCode )
		case 2: // mouse up

			break
		case -1: // control being killed
			break
	endswitch

	return 0
End
Static Function ButtonProc_2(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	string saveDF
	variable FOffset, Soffset

	switch( ba.eventCode )
		case 2: // mouse up
			saveDF = GetDataFolder(1)
			wave/T parmWave=root:DE_Dudko:MenuStuff:ParmWave


			controlinfo de_Dudko_popup0
			SetDataFolder s_value
			controlinfo de_Dudko_popup1
			wave ForceWave=$S_value
			wave SepWave=$ReplaceString("Force",S_value,"Sep")
			controlinfo de_Dudko_popup2
			wave UpPoints=$S_value
			controlinfo de_Dudko_popup3
			wave DownPoints=$S_value
			duplicate/o ForceWave root:De_DUDKO:ForceWaveS
			duplicate/o SepWave root:De_DUDKO:SepWaveS
			wave ForceWaveS=root:De_DUDKO:ForceWaveS
			wave SepWaveS=root:De_DUDKO:SepWaveS
			make/o/n=0 root:De_DUDKO:ForceWaveS_SH,root:De_DUDKO:ForceWave_SH
			wave ForceWaveS_SH=root:De_DUDKO:ForceWaveS_SH
			wave ForceWave_SH= root:De_DUDKO:ForceWave_SH
			DE_Filtering#FilterForceSep(ForceWave,SepWave,ForceWaveS,SepWaveS,"TVD",str2num(parmWave[0][1]))
			ForceWaveS*=-1;
			make/o/n=0 root:DE_Dudko:States
			wave States=root:DE_Dudko:States

			DE_DUDKO#MakeSingleStateKey(ForceWaveS,UpPoints,DownPoints,States)
			print DE_OverlapRamps#AddForceOffsetstoForceWave(ForceWaveS,SepWaveS,States)
			note/k ForceWave note(ForceWaveS)

			DE_OverlapRamps#OffsetEachStepinForceWave(ForceWaveS,ForceWaveS_SH,States)
			DE_OverlapRamps#OffsetEachStepinForceWave(ForceWave,ForceWave_SH,States)

			DE_OverlapRamps#AddShiftToStates(ForceWaveS,States)
			struct ForceWave StartName
			DE_Naming#WavetoStruc(nameofwave(ForceWave),StartName)
			string BaseName=startname.Name+startname.SNum
			//DE_DUdko#ConstructHistograms(States,Reference,10,NameBase=BaseName)
			duplicate/o States $(Basename+"_States")
			wave w1=$(Basename+"_Ref")
			duplicate/o ForceWaveS $(Basename+"_Sm")
			duplicate/o ForceWaveS_SH $(Basename+"_Sm_SH")
			duplicate/o SepWaveS $(Basename+"Sep_Sm")
			
			duplicate/o ForceWave_SH $(Basename+"_SH")

			SetDataFolder saveDF
				
		
				
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Static Function ButtonProc_1(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	string saveDF
	variable FOffset, Soffset
	switch( ba.eventCode )
		case 2: // mouse up
			saveDF = GetDataFolder(1)
			string Here1
			dowindow AllShift
			if(V_flag==1)
				killwindow AllShift
				display/n=AllShift

			else
				display/n=AllShift
			endif
			wave/T parmWave=root:DE_Dudko:MenuStuff:ParmWave
			make/T/Free/n=0 AvailableForceWaves
			controlinfo/W=Dudko_Panel de_Dudko_popup0
			SetDataFolder s_value
			Here1=GetDataFolder(1)
			DE_PanelProgs#ListWaves($Here1,"*Force_Adj",AvailableForceWaves)
			variable n=0
			for(n=0;n<numpnts(AvailableForceWaves);n+=1)
				
				print "At: "+num2str(n)+"  TOTAL: "+num2str(numpnts(AvailableForceWaves))
				wave ForceWave=$AvailableForceWaves[n]	
				print nameofwave(ForceWave)
				wave SepWave=$ReplaceString("Force",AvailableForceWaves[n],"Sep")
				wave UpPoints=$ReplaceString("Force",AvailableForceWaves[n],"RupPntU")
				wave DownPoints=$ReplaceString("Force",AvailableForceWaves[n],"RupPntD")
				duplicate/o ForceWave root:De_DUDKO:ForceWaveS
				duplicate/o SepWave root:De_DUDKO:SepWaveS
				wave ForceWaveS=root:De_DUDKO:ForceWaveS
				wave SepWaveS=root:De_DUDKO:SepWaveS
				make/o/n=0 root:De_DUDKO:ForceWaveS_SH
				wave ForceWaveS_SH=root:De_DUDKO:ForceWaveS_SH
			
				variable pointstoignore=floor(str2num(parmWave[2][1])/str2num(parmWave[1][1])/dimdelta(ForceWaveS,0))
				variable FitTime=str2num(parmWave[3][1])/str2num(parmWave[1][1])
				//DE_Filtering#TVD1D_denoise(ForceWave,str2num(parmWave[0][1]),ForceWaveS)
				DE_Filtering#FilterForceSep(ForceWave,SepWave,ForceWaveS,SepWaveS,"TVD",str2num(parmWave[0][1]))
				ForceWaveS*=-1;
				make/o/n=0 root:DE_Dudko:States
				wave States=root:DE_Dudko:States

				DE_DUDKO#MakeSingleStateKey(ForceWaveS,UpPoints,DownPoints,States)
				DE_OverlapRamps#AddForceOffsetstoForceWave(ForceWaveS,SepWaveS,States)
				DE_OverlapRamps#OffsetEachStepinForceWave(ForceWaveS,ForceWaveS_SH,States)
				DE_OverlapRamps#AddShiftToStates(ForceWaveS,States)
				make/o/n=(60,7) root:DE_DUDKO:Reference	
				wave Reference=root:DE_DUDKO:Reference
			
				Reference=0;
			
				Reference[][0]=1.5*DE_DUDKO#MaxRupForce(States)/59*p
			
				make/free/n=3 Results
				//ContourLengthDetermine(ForceWaveS_SH,SepWaveS,States,pointstoignore,Results)
				variable/C slopes=DE_DUDKO#ReturnSeparationSlopes(SepWaveS,States,pointstoignore)
				print/C slopes
				DE_DUDKO#ContourLengthDetermineCombined(ForceWaveS_SH,SepWaveS,States,pointstoignore,Results,1)
				make/free/n=0 WLC1,WLC2
				MakeWLCs(SepWaveS,Results,WLC1,WLC2)
				variable LC,Offset
				LC=Results[0]
				offSet=Results[2]
				print results
				//Reference[][1]=ReturnSlopesFromWLC(LC,Offset,real(Reference[p][0]),str2num(parmWave[1][1]))
				Reference[][1]=FitFreeWLCSlope(LC,Offset,real(Reference[p][0]),real(slopes))
				LC=Results[1]
				//Reference[][2]=ReturnSlopesFromWLC(LC,Offset,real(Reference[p][0]),str2num(parmWave[1][1]))
				Reference[][2]=FitFreeWLCSlope(LC,Offset,real(Reference[p][0]),imag(slopes))
			
				Reference[][5]=real(NewAvgSlope(ForceWaveS_SH,States,real(Reference[p][0]),-1,1,FitTime,pointstoignore))
				Reference[][6]=real(NewAvgSlope(ForceWaveS_SH,States,real(Reference[p][0]),1,1,FitTime,pointstoignore))


				Reference[][3]=real(DE_DUDKO#NumberofTracesinState(ForceWaveS_SH,States,real(Reference[p][0]),-1,1))	
				Reference[][4]=real(DE_DUDKO#NumberofTracesinState(ForceWaveS_SH,States,real(Reference[p][0]),1,1))
				//			
				struct ForceWave StartName
				DE_Naming#WavetoStruc(nameofwave(ForceWave),StartName)
				string BaseName=startname.Name+startname.SNum
				//DE_DUdko#ConstructHistograms(States,Reference,10,NameBase=BaseName)
				duplicate/o States $(Basename+"_States")
				duplicate/o Reference $(Basename+"_Ref")
				wave w1=$(Basename+"_Ref")
				duplicate/o ForceWaveS $(Basename+"_Sm")
				duplicate/o ForceWaveS_SH $(Basename+"_Sm_SH")
				duplicate/o SepWaveS $(Basename+"Sep_Sm")
				duplicate/o WLC1 $(Basename+"WLC1")
				duplicate/o WLC2 $(Basename+"WLC2")
				wave w2= $(Basename+"_Sm_SH")
				wave w3= $(Basename+"Sep_Sm")
				wave w4= $(Basename+"WLC1")
				wave w5=$(Basename+"WLC2")
				Display w1[][1] vs w1[][0]	
				appendtograph w1[][2] vs w1[][0]	
				appendtograph w1[][5] vs w1[][0]	
				appendtograph w1[][6] vs w1[][0]	
				ModifyGraph rgb($(nameofwave(w1)))=(58368,6656,7168);DelayUpdate
				ModifyGraph lsize($(nameofwave(w1)))=2,lsize($(nameofwave(w1)+"#1"))=2
				ModifyGraph rgb($(nameofwave(w1)+"#1"))=(14848,32256,47104);DelayUpdate
			
				ModifyGraph mode($(nameofwave(w1)+"#2"))=3,marker($(nameofwave(w1)+"#2"))=8;DelayUpdate
				ModifyGraph rgb($(nameofwave(w1)+"#2"))=(58368,6656,7168);DelayUpdate
				ModifyGraph mode($(nameofwave(w1)+"#3"))=3,marker($(nameofwave(w1)+"#3"))=8;DelayUpdate
				ModifyGraph rgb($(nameofwave(w1)+"#3"))=(14848,32256,47104)
				display w2 vs w3
				appendtograph w4
				appendtograph w5
				ModifyGraph rgb($nameofwave(w4))=(0,0,0),rgb($nameofwave(w5))=(0,0,0)
				SetAxis left 0,5e-11
				
				Appendtograph/w=AllShift w2 vs w3
				SetDataFolder saveDF
				//				
			endfor
				
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End
Static Function ButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	string saveDF
	variable FOffset, Soffset
//
	switch( ba.eventCode )
		case 2: // mouse up
			saveDF = GetDataFolder(1)
			wave/T parmWave=root:DE_Dudko:MenuStuff:ParmWave
//
//			//find the right waves
			controlinfo de_Dudko_popup0
			SetDataFolder s_value
			controlinfo de_Dudko_popup1
			wave ForceWave=$S_value
			wave SepWave=$ReplaceString("Force",S_value,"Sep")
			controlinfo de_Dudko_popup2
			wave States=$S_value
			controlinfo de_Dudko_popup3
			wave Results=$S_value
//
//			//make from additional waves
			duplicate/o ForceWave root:De_DUDKO:ForceWave_SH
			duplicate/o ForceWave root:De_DUDKO:ForceWaveS
			duplicate/o SepWave root:De_DUDKO:SepWaveS
			wave ForceWaveS=root:De_DUDKO:ForceWaveS
			//wave ForceWaveSH=root:De_DUDKO:ForceWave_SH
			wave SepWaveS=root:De_DUDKO:SepWaveS
			make/o/n=0 root:De_DUDKO:ForceWaveS_SH
			wave ForceWaveS_SH=root:De_DUDKO:ForceWaveS_SH	
//					
//			//define some variables
			variable pointstoignore=floor(str2num(parmWave[2][1])/str2num(parmWave[1][1])/dimdelta(ForceWaveS,0))
			variable FitTime=str2num(parmWave[3][1])/str2num(parmWave[1][1])
//			
//			//filter and correct waves
			DE_Filtering#FilterForceSep(ForceWave,SepWave,ForceWaveS,SepWaveS,"TVD",str2num(parmWave[0][1]))
			ForceWaveS*=-1;
			ForceWave*=-1
//
//			//make the state key
//			make/o/n=0 root:DE_Dudko:States
//			wave States=root:DE_Dudko:States
//			DE_DUDKO#MakeSingleStateKey(ForceWaveS,UpPoints,DownPoints,States)
//			
//			//handle offsetting both the smoothed wave and the unsmoothed wave, also corrects our estimates of the rupture forces
//			String Offsets=DE_OverlapRamps#AddForceOffsetstoForceWave(ForceWaveS,SepWaveS,States)
//			note/K ForceWaveSH,ReplaceStringByKey("DE_FOff", note(ForceWaveSH), Offsets,":","\r" )
//			duplicate/free ForceWaveSH HoldSH
//			DE_OverlapRamps#OffsetEachStepinForceWave(HoldSH,ForceWaveSH,States)
//			DE_OverlapRamps#OffsetEachStepinForceWave(ForceWaveS,ForceWaveS_SH,States)
//			DE_OverlapRamps#AddShiftToStates(ForceWaveS,States)
//			
//			//Makes our reference wave
			make/o/n=(60,7) root:DE_DUDKO:Reference	
			wave Reference=root:DE_DUDKO:Reference
			Reference=0;
			Reference[][0]=-1.5*DE_DUDKO#MaxRupForce(States)/59*p
//			
//			make/free/n=3 Results
//			//ContourLengthDetermine(ForceWaveS_SH,SepWaveS,States,pointstoignore,Results)
			variable/C slopes=DE_DUDKO#ReturnSeparationSlopes(SepWaveS,States,pointstoignore)
			print/c slopes
//			DE_DUDKO#ContourLengthDetermineCombined(ForceWaveS_SH,SepWaveS,States,pointstoignore,Results)
			make/o/n=0 WLC1a,WLC2a
			MakeWLCs(SepWaveS,Results,WLC1a,WLC2a)
			variable LC,Offset
			print Results
			LC=Results[0]
			offSet=Results[2]
//			print results
//			//Reference[][1]=ReturnSlopesFromWLC(LC,Offset,real(Reference[p][0]),str2num(parmWave[1][1]))
			Reference[][1]=FitFreeWLCSlope(LC,Offset,real(Reference[p][0]),real(slopes))
			LC=Results[1]
			//Reference[][2]=ReturnSlopesFromWLC(LC,Offset,real(Reference[p][0]),str2num(parmWave[1][1]))
			Reference[][2]=FitFreeWLCSlope(LC,Offset,real(Reference[p][0]),imag(slopes))
//			//Reference[][5]=real(NewAvgSlope(ForceWaveS_SH,States,real(Reference[p][0]),-1,1,FitTime,pointstoignore))
//			//Reference[][6]=real(NewAvgSlope(ForceWaveS_SH,States,real(Reference[p][0]),1,1,FitTime,pointstoignore))
//
//
			Reference[][3]=real(DE_DUDKO#NumberofTracesinState(ForceWaveS,States,real(Reference[p][0]),-1,1))	
			Reference[][4]=real(DE_DUDKO#NumberofTracesinState(ForceWaveS,States,real(Reference[p][0]),1,1))
//			//			
			struct ForceWave StartName
			DE_Naming#WavetoStruc(nameofwave(ForceWave),StartName)
			string BaseName=startname.Name+startname.SNum
//			//DE_DUdko#ConstructHistograms(States,Reference,10,NameBase=BaseName)
print BaseName
//			duplicate/o States $(Basename+"_States")
			duplicate/o Reference $(Basename+"_Ref")
//			wave w1=$(Basename+"_Ref")
			duplicate/o ForceWaveS $(Basename+"_Sm")
						ForceWave*=-1

//			duplicate/o ForceWaveSH $(Basename+"_Sh")
//
//			duplicate/o ForceWaveS_SH $(Basename+"_Sm_SH")
//			duplicate/o SepWaveS $(Basename+"Sep_Sm")
//			duplicate/o WLC1 $(Basename+"WLC1")
//			duplicate/o WLC2 $(Basename+"WLC2")
//			duplicate/o results $(Basename+"WLCParms")
//			wave w2= $(Basename+"_Sm_SH")
//			wave w3= $(Basename+"Sep_Sm")
//			wave w4= $(Basename+"WLC1")
//			wave w5=$(Basename+"WLC2")
//			Display w1[][1] vs w1[][0]	
//
//			appendtograph w1[][2] vs w1[][0]	
//			appendtograph w1[][5] vs w1[][0]	
//			appendtograph w1[][6] vs w1[][0]	
//
//			ModifyGraph rgb($(nameofwave(w1)))=(58368,6656,7168);DelayUpdate
//			ModifyGraph lsize($(nameofwave(w1)))=2,lsize($(nameofwave(w1)+"#1"))=2
//			ModifyGraph rgb($(nameofwave(w1)+"#1"))=(14848,32256,47104);DelayUpdate
//			
//			ModifyGraph mode($(nameofwave(w1)+"#2"))=3,marker($(nameofwave(w1)+"#2"))=8;DelayUpdate
//			ModifyGraph rgb($(nameofwave(w1)+"#2"))=(58368,6656,7168);DelayUpdate
//			ModifyGraph mode($(nameofwave(w1)+"#3"))=3,marker($(nameofwave(w1)+"#3"))=8;DelayUpdate
//			ModifyGraph rgb($(nameofwave(w1)+"#3"))=(14848,32256,47104)
//						
//			display w2 vs w3
//			appendtograph w4
//			appendtograph w5
//			ModifyGraph rgb($nameofwave(w4))=(0,0,0),rgb($nameofwave(w5))=(0,0,0)
//			SetAxis left 0,5e-11
//			SetDataFolder saveDF
//				
//			break
//		case -1: // control being killed
//			break
	endswitch
//
	return 0
end

//Static Function ButtonProc(ba) : ButtonControl
//	STRUCT WMButtonAction &ba
//	string saveDF
//	variable FOffset, Soffset
//
//	switch( ba.eventCode )
//		case 2: // mouse up
//			saveDF = GetDataFolder(1)
//			wave/T parmWave=root:DE_Dudko:MenuStuff:ParmWave
//
//			//find the right waves
//			controlinfo de_Dudko_popup0
//			SetDataFolder s_value
//			controlinfo de_Dudko_popup1
//			wave ForceWave=$S_value
//			wave SepWave=$ReplaceString("Force",S_value,"Sep")
//			controlinfo de_Dudko_popup2
//			wave UpPoints=$S_value
//			controlinfo de_Dudko_popup3
//			wave DownPoints=$S_value
//
//			//make from additional waves
//			duplicate/o ForceWave root:De_DUDKO:ForceWave_SH
//			duplicate/o ForceWave root:De_DUDKO:ForceWaveS
//			duplicate/o SepWave root:De_DUDKO:SepWaveS
//			wave ForceWaveS=root:De_DUDKO:ForceWaveS
//			wave ForceWaveSH=root:De_DUDKO:ForceWave_SH
//			wave SepWaveS=root:De_DUDKO:SepWaveS
//			make/o/n=0 root:De_DUDKO:ForceWaveS_SH
//			wave ForceWaveS_SH=root:De_DUDKO:ForceWaveS_SH	
//					
//			//define some variables
//			variable pointstoignore=floor(str2num(parmWave[2][1])/str2num(parmWave[1][1])/dimdelta(ForceWaveS,0))
//			variable FitTime=str2num(parmWave[3][1])/str2num(parmWave[1][1])
//			
//			//filter and correct waves
//			DE_Filtering#FilterForceSep(ForceWave,SepWave,ForceWaveS,SepWaveS,"TVD",str2num(parmWave[0][1]))
//			ForceWaveS*=-1;
//			ForceWaveSH*=-1
//
//			//make the state key
//			make/o/n=0 root:DE_Dudko:States
//			wave States=root:DE_Dudko:States
//			DE_DUDKO#MakeSingleStateKey(ForceWaveS,UpPoints,DownPoints,States)
//			
//			//handle offsetting both the smoothed wave and the unsmoothed wave, also corrects our estimates of the rupture forces
//			String Offsets=DE_OverlapRamps#AddForceOffsetstoForceWave(ForceWaveS,SepWaveS,States)
//			note/K ForceWaveSH,ReplaceStringByKey("DE_FOff", note(ForceWaveSH), Offsets,":","\r" )
//			duplicate/free ForceWaveSH HoldSH
//			DE_OverlapRamps#OffsetEachStepinForceWave(HoldSH,ForceWaveSH,States)
//			DE_OverlapRamps#OffsetEachStepinForceWave(ForceWaveS,ForceWaveS_SH,States)
//			DE_OverlapRamps#AddShiftToStates(ForceWaveS,States)
//			
//			//Makes our reference wave
//			make/o/n=(60,7) root:DE_DUDKO:Reference	
//			wave Reference=root:DE_DUDKO:Reference
//			Reference=0;
//			Reference[][0]=1.5*DE_DUDKO#MaxRupForce(States)/59*p
//			
//			make/free/n=3 Results
//			//ContourLengthDetermine(ForceWaveS_SH,SepWaveS,States,pointstoignore,Results)
//			variable/C slopes=DE_DUDKO#ReturnSeparationSlopes(SepWaveS,States,pointstoignore)
//			print/c slopes
//			DE_DUDKO#ContourLengthDetermineCombined(ForceWaveS_SH,SepWaveS,States,pointstoignore,Results)
//			make/free/n=0 WLC1,WLC2
//			MakeWLCs(SepWaveS,Results,WLC1,WLC2)
//			variable LC,Offset
//			LC=Results[0]
//			offSet=Results[2]
//			print results
//			//Reference[][1]=ReturnSlopesFromWLC(LC,Offset,real(Reference[p][0]),str2num(parmWave[1][1]))
//			Reference[][1]=FitFreeWLCSlope(LC,Offset,real(Reference[p][0]),real(slopes))
//			LC=Results[1]
//			//Reference[][2]=ReturnSlopesFromWLC(LC,Offset,real(Reference[p][0]),str2num(parmWave[1][1]))
//			Reference[][2]=FitFreeWLCSlope(LC,Offset,real(Reference[p][0]),imag(slopes))
//			//Reference[][5]=real(NewAvgSlope(ForceWaveS_SH,States,real(Reference[p][0]),-1,1,FitTime,pointstoignore))
//			//Reference[][6]=real(NewAvgSlope(ForceWaveS_SH,States,real(Reference[p][0]),1,1,FitTime,pointstoignore))
//
//
//			Reference[][3]=real(DE_DUDKO#NumberofTracesinState(ForceWaveS_SH,States,real(Reference[p][0]),-1,1))	
//			Reference[][4]=real(DE_DUDKO#NumberofTracesinState(ForceWaveS_SH,States,real(Reference[p][0]),1,1))
//			//			
//			struct ForceWave StartName
//			DE_Naming#WavetoStruc(nameofwave(ForceWave),StartName)
//			string BaseName=startname.Name+startname.SNum
//			//DE_DUdko#ConstructHistograms(States,Reference,10,NameBase=BaseName)
//			duplicate/o States $(Basename+"_States")
//			duplicate/o Reference $(Basename+"_Ref")
//			wave w1=$(Basename+"_Ref")
//			duplicate/o ForceWaveS $(Basename+"_Sm")
//			duplicate/o ForceWaveSH $(Basename+"_Sh")
//
//			duplicate/o ForceWaveS_SH $(Basename+"_Sm_SH")
//			duplicate/o SepWaveS $(Basename+"Sep_Sm")
//			duplicate/o WLC1 $(Basename+"WLC1")
//			duplicate/o WLC2 $(Basename+"WLC2")
//			duplicate/o results $(Basename+"WLCParms")
//			wave w2= $(Basename+"_Sm_SH")
//			wave w3= $(Basename+"Sep_Sm")
//			wave w4= $(Basename+"WLC1")
//			wave w5=$(Basename+"WLC2")
//			Display w1[][1] vs w1[][0]	
//
//			appendtograph w1[][2] vs w1[][0]	
//			appendtograph w1[][5] vs w1[][0]	
//			appendtograph w1[][6] vs w1[][0]	
//
//			ModifyGraph rgb($(nameofwave(w1)))=(58368,6656,7168);DelayUpdate
//			ModifyGraph lsize($(nameofwave(w1)))=2,lsize($(nameofwave(w1)+"#1"))=2
//			ModifyGraph rgb($(nameofwave(w1)+"#1"))=(14848,32256,47104);DelayUpdate
//			
//			ModifyGraph mode($(nameofwave(w1)+"#2"))=3,marker($(nameofwave(w1)+"#2"))=8;DelayUpdate
//			ModifyGraph rgb($(nameofwave(w1)+"#2"))=(58368,6656,7168);DelayUpdate
//			ModifyGraph mode($(nameofwave(w1)+"#3"))=3,marker($(nameofwave(w1)+"#3"))=8;DelayUpdate
//			ModifyGraph rgb($(nameofwave(w1)+"#3"))=(14848,32256,47104)
//						
//			display w2 vs w3
//			appendtograph w4
//			appendtograph w5
//			ModifyGraph rgb($nameofwave(w4))=(0,0,0),rgb($nameofwave(w5))=(0,0,0)
//			SetAxis left 0,5e-11
//			SetDataFolder saveDF
//				
//			break
//		case -1: // control being killed
//			break
//	endswitch
//
//	return 0
//End

Window Dudko_Panel() : Panel

	PauseUpdate; Silent 1		// building window...
	NewPanel/N=Dudko_Panel /W=(0,0,300,450)
	NewDataFolder/o root:DE_Dudko
	NewDataFolder/o root:DE_Dudko:MenuStuff

	DE_Dudko#UpdateParmWave()
	Button de_RupRamp_button0,pos={75,130},size={150,20},proc=DE_Dudko#ButtonProc,title="GO!"
		Button de_RupRamp_button2,pos={75,155},size={150,20},proc=DE_Dudko#ButtonProc_2,title="GO!"

	Button de_RupRamp_button1,pos={75,400},size={150,20},proc=DE_Dudko#ButtonProc_1,title="All!"
	PopupMenu de_Dudko_popup0,pos={75,2},size={129,21},Title="Folder"
	PopupMenu de_Dudko_popup0,mode=1,popvalue="X",value= #"DE_PanelProgs#ListFolders()"
	PopupMenu de_Dudko_popup1,pos={75,40},size={129,21},Title="ForceWave"
	PopupMenu de_Dudko_popup1,mode=1,popvalue="X",value= #"DE_dudko#ListWaves(\"de_Dudko_popup0\")"
	PopupMenu de_Dudko_popup2,pos={75,70},size={129,21},Title="Up"
	PopupMenu de_Dudko_popup2,mode=1,popvalue="X",value= #"DE_Dudko#ListWaves(\"de_Dudko_popup0\")"
	PopupMenu de_Dudko_popup3,pos={75,100},size={129,21},Title="Down"
	PopupMenu de_Dudko_popup3,mode=1,popvalue="X",value= #"DE_Dudko#ListWaves(\"de_Dudko_popup0\")"
	//ListBox DE_Dudko_list0,pos={50,200},size={175,150},proc=DE_Dudko#ListBoxProc,listWave=root:DE_Dudko:MenuStuff:ParmWave
	ListBox DE_Dudko_list0,pos={50,200},size={175,150},listWave=root:DE_Dudko:MenuStuff:ParmWave
	ListBox DE_Dudko_list0,selWave=root:DE_dudko:MenuStuff:SelWave,editStyle= 2,userColumnResize= 1,widths={70,40,70,40}

EndMacro

Static Function/S ListWaves(ControlStr)
	string ControlStr
	String saveDF

	saveDF = GetDataFolder(1)
	controlinfo $ControlStr
	SetDataFolder s_value
	String list = WaveList("*", ";", "")
	SetDataFolder saveDF
	return list

end

Static Function UpdateParmWave()
	if(exists("root:DE_Dudko:MenuStuff:ParmWave")==1)
		wave/t/z Par=root:DE_Dudko:MenuStuff:ParmWave
		wave/z Sel=root:DE_Dudko:MenuStuff:SelWave
	Else
		make/t/n=(4,2) root:DE_Dudko:MenuStuff:ParmWave
		wave/t/z Par=root:DE_Dudko:MenuStuff:ParmWave
		make/n=(4,2) root:DE_Dudko:MenuStuff:SelWave
		wave/z Sel=root:DE_Dudko:MenuStuff:SelWave
		
		Par[0][0]={"Smoothing","Pulling Rate","Distance to Ignore","Distance to Fit"}
		Par[0][1]={"10e-9","50e-9","3e-9","3e-9"}
		Sel[][0]=0
		Sel[][1]=2
	endif

end




 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 //Visualization
 
 Static Function OverlayTwoCurves(ForceIn,Sepin,n1,n2)
	wave Forcein,SepIn
	variable n1,n2
	
	make/o/n=0 FDisp1,FDisp2,SDisp1,SDisp2
	
	ExtractForceCurve(ForceIn,SepIn,n1,FDisp1,SDisp1)
	ExtractForceCurve(ForceIn,SepIn,n2,FDisp2,SDisp2)
	SetScale/P x 0,dimdelta(FDisp1,0),"s", FDisp1
	SetScale/P x 0,dimdelta(FDisp2,0),"s", FDisp2

	SetScale/P x 0,dimdelta(sDisp1,0),"s", SDisp1
	SetScale/P x 0,dimdelta(sDisp2,0),"s", SDisp2

	variable offset1= str2num(stringfromlist(n1,stringbykey("DE_FOff",note(ForceIN),":","\r")))
	variable offset2= str2num(stringfromlist(n2,stringbykey("DE_FOff",note(ForceIN),":","\r")))
	FDisp1+=offset1
	FDIsp2+=offset2
		dowindow DEDISP

	if(v_flag==1)
	else 
		Display/n=DEDISP FDisp1 vs SDisp1
		Appendtograph/w=DEDISP FDisp2 vs SDIsp2
		ModifyGraph rgb(FDisp1)=(58368,6656,7168),rgb(FDisp2)=(14848,32256,47104)
	endif
end

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Offset Forces
//
//
//Static Function ExtractForceCurveWithGaps(ForceIn,SepIn,StateWave,num,CutWave,ForceOut,SepOut)
//	Wave Forcein,ForceOut,StateWave,SepIn,SepOut,CutWave
//	variable num
//
//	variable startret,mid,endext
//	if(num==0)
//		startRet=0
//		Mid= str2num(StringfromList(0,stringbykey("DE_Ind",note(Forcein),":","\r")	))
//		EndExt=str2num(StringfromList(1,stringbykey("DE_Ind",note(Forcein),":","\r")	))
//	else
//		startRet=str2num(StringfromList(2*num-1,stringbykey("DE_Ind",note(Forcein),":","\r")	))
//		Mid= str2num(StringfromList(2*num,stringbykey("DE_Ind",note(Forcein),":","\r")	))
//		EndExt=str2num(StringfromList(2*num+1,stringbykey("DE_Ind",note(Forcein),":","\r")	))
//	endif
//	
//	duplicate/free/r=[startRet,EndExt] ForceIn FO
//	duplicate/free/r=[startRet,EndExt] SepIn SO
//	SetScale/P x 0,dimdelta(ForceIn,0),"s", FO
//	SetScale/P x 0,dimdelta(ForceIn,0),"s", SO
//
//	variable start1= CutWave[0]
//	variable end1= CutWave[1]
//	variable start2= CutWave[2]
//	variable end2= CutWave[3]
//
//	deletepoints start2,(end2-start2), FO
//	deletepoints start1,(end1-start1), FO
//	deletepoints start2,(end2-start2), SO
//	deletepoints start1,(end1-start1), SO
//	duplicate/o FO ForceOut
//	duplicate/o SO SepOut
//
//end
//
//Static Function ReturnGaps(ForceIn,SepIn,StateWave,num,Results)
//	Wave Forcein,SepIn,StateWave,Results
//	variable num
//
//	variable startret,mid,endext
//	if(num==0)
//		startRet=0
//		Mid= str2num(StringfromList(0,stringbykey("DE_Ind",note(Forcein),":","\r")	))
//		EndExt=str2num(StringfromList(1,stringbykey("DE_Ind",note(Forcein),":","\r")	))
//	else
//		startRet=str2num(StringfromList(2*num-1,stringbykey("DE_Ind",note(Forcein),":","\r")	))
//		Mid= str2num(StringfromList(2*num,stringbykey("DE_Ind",note(Forcein),":","\r")	))
//		EndExt=str2num(StringfromList(2*num+1,stringbykey("DE_Ind",note(Forcein),":","\r")	))
//	endif
//	
//	duplicate/free/r=[startRet,EndExt] ForceIn FO
//	duplicate/free/r=[startRet,EndExt] SepIn SO
//	make/free/n=0 ReturnPoints
//	ExtractPointsFromThisState(StateWave,num,ReturnPoints)
//	FindLevel/q ReturnPoints, Mid
//	variable start1= ReturnPoints[1]-250
//	start1-=startRet
//	variable end1= ReturnPoints[v_levelx-1]+250
//	end1-=startRet
//
//	variable start2= ReturnPoints[v_levelx+1]	-250
//	start2-=startRet
//
//	variable end2= ReturnPoints[numpnts(returnpoints)-1]+250
//	end2-=startRet
//
//	make/free/n=4 Garbage
//	Garbage={start1,end1,start2,end2}
//	duplicate/o Garbage Results
//
//end
//
//Static Function ReturnForcesWidestGaps(ForceWave1,SepWave1,StateWave1,num1,num2,OutForce1,OutSep1,OutForce2,OutSep2)
//	wave ForceWave1,SepWave1,StateWave1,OutForce1,OutSep1,OutForce2,OutSep2
//	variable num1,num2
//	make/free/N=0 Points1,Points2,RF1,RF2,RS1,RS2
//	make/free/n=4 FinalPoints
//	ReturnGaps(ForceWave1,SepWave1,StateWave1,num1,Points1)
//	ReturnGaps(ForceWave1,SepWave1,StateWave1,num2,Points2)
//	FinalPoints[0]=min(Points1[0],Points2[0])
//	FinalPoints[1]=max(Points1[1],Points2[1])
//	FinalPoints[2]=min(Points1[2],Points2[2])
//	FinalPoints[3]=max(Points1[3],Points2[3])
//	ExtractForceCurveWithGaps(ForceWave1,SepWave1,StateWave1,num1,FinalPoints,RF1,RS1)
//	ExtractForceCurveWithGaps(ForceWave1,SepWave1,StateWave1,num2,FinalPoints,RF2,RS2)
//	duplicate/o RF1 OutForce1
//	duplicate/o RS1 OutSep1
//	duplicate/o RF2 OutForce2
//	duplicate/o RS2 OutSep2
//
//end
//
//
//
//Static Function CalculateMaxIndexforForceWave(ForceWave)
//	wave ForceWave
//	
//	
//	string FullInd=Stringbykey("DE_Ind",note(ForceWave),":","\r")
//	variable numberofpulls=itemsinlist(FullINd,";")/2
//	return numberofpulls
//end
//
//
//
//Static Function CalculateForceOffsetIndex(ForceWave,SepWave,Statewave,index1,index2)
//	wave ForceWave,SepWave,Statewave
//	variable index1,index2
//	make/free/n=0 F1,F2,S1,S2;ReturnForcesWidestGaps(ForceWave,SepWave,Statewave,index1,index2,F1,S1,F2,S2)
//	make/o/n=2 W_Coef
//	W_Coef={0,1}
//	CurveFit/Q/H="01"/NTHR=0 line  kwCWave=W_Coef F1 /X=F2 /D 
//	return W_Coef[0]
//end
//
//Static Function ExtractPointsFromThisState(StateWave,Index,ReturnPoints)
//	wave StateWave,ReturnPoints
//	variable Index
//	
//	make/free/n=(dimsize(StateWave,0)) Points,Type,States
//	Points=StateWave[p][0]
//	Type=StateWave[p][2]
//	States=StateWave[p][3]
//	Extract/Free/Indx States, HiDevin, States==index
//	make/free/n=(numpnts(HiDevin)) Results
//	Results=Points[HiDevin]
//	duplicate/o Results ReturnPoints
//
//end
//
//Static Function OffsetEachStepinForceWave(ForceIn,ForceOut,States)
//	wave ForceIn,ForceOut,States
//	
//	variable n,prevcut=0,currcut,FOFF
//	duplicate/o ForceIn ForceOut
//	for(n=1;n<dimsize(States,0);n+=1)
//		if(States[n][2]==2)
//			currcut=States[n][0]
//			FOFF=ReturnOffForceForState(ForceIn,States[n][3]-1)
//			ForceOut[prevcut,currcut]+=FOFF
//			prevcut=States[n][0]
//		endif
//	
//	endfor
//end
//
//Static Function/S AddForceOffsetstoForceWave(ForceWave,SepWave,StateWave)
//	wave ForceWave,SepWave,StateWave
//	variable numsteps=CalculateMaxIndexforForceWave(ForceWave)
//	variable n=0
//	string Offsets="0;"
//	for(n=1;n<numsteps;n+=1)
//		Offsets+=num2str(CalculateForceOffsetIndex(ForceWave,SepWave,StateWave,0,n))+";"
//	
//	
//	endfor
//	note/K ForceWave,ReplaceStringByKey("DE_FOff", note(ForceWave), Offsets,":","\r" )
//	note/K SepWave,ReplaceStringByKey("DE_FOff", note(ForceWave), Offsets,":","\r" )
//	return (Offsets)
//end
//
//Static Function ReturnOffForceForState(ForceWave,index)
//
//		wave ForceWave
//		variable index
//		return str2num(stringfromlist(index,stringbykey("DE_FOff",note(ForceWave),":","\r")))
//end
//
//Static Function ReturnForcewithOffset(ForceWave,StateWave,Point)
//	wave ForceWave,StateWave	
//	variable Point
//	
//	make/free/n=(dimsize(StateWave,0)) Points
//	Points=StateWave[p][0]
//	FindLevel/Q Points Point
//	variable index=StateWave[ceil(V_levelx)][3]
//	return (ForceWave[Point]+ReturnOffForceForState(ForceWave,index))
//end
//
//Static Function AddShiftToStates(ForceWave,StateWave)
//	wave ForceWave,StateWave
//	
//	make/free/n=(dimsize(StateWave,0),6) AboutToFinish
//	AboutToFinish[][0]=StateWave[p][0]
//	AboutToFinish[][1]=StateWave[p][1]
//	AboutToFinish[][2]=StateWave[p][2]
//	AboutToFinish[][3]=StateWave[p][3]
//	AboutToFinish[][4]=StateWave[p][4]
//	AboutToFinish[][5]=StateWave[p][1]+ReturnOffForceForState(ForceWave,StateWave[p][3])
//
//	duplicate/o AbouttoFinish StateWave
//end
//
//

//Misc

Static Function ExtractForceCurve(ForceIn,SepIn,num,ForceOut,SepOut)
	Wave Forcein,ForceOut,SepIn,SepOut
	variable num

	variable startret,mid,endext
	if(num==0)
		startRet=0
		Mid= str2num(StringfromList(0,stringbykey("DE_Ind",note(Forcein),":","\r")	))
		EndExt=str2num(StringfromList(1,stringbykey("DE_Ind",note(Forcein),":","\r")	))
	else
		startRet=str2num(StringfromList(2*num-1,stringbykey("DE_Ind",note(Forcein),":","\r")	))
		Mid= str2num(StringfromList(2*num,stringbykey("DE_Ind",note(Forcein),":","\r")	))
		EndExt=str2num(StringfromList(2*num+1,stringbykey("DE_Ind",note(Forcein),":","\r")	))
	endif

	duplicate/o/r=[startRet,EndExt] ForceIn ForceOut
	duplicate/o/r=[startRet,EndExt] SepIn SepOut

end

Static Function AppendToWaveWOverlapOW(StartingWave,WavetoAdd,PointtoStart)

	wave StartingWave,WavetoAdd
	variable PointtoStart

	variable a=numpnts(Startingwave)
	variable b=numpnts(wavetoadd)
	duplicate/free startingwave test
	insertpoints (numpnts(startingwave)-1), (PointtoStart-(numpnts(startingwave)-1)), Test

	Test[(numpnts(startingwave)-1),]=NaN
	insertpoints PointtoStart, (numpnts(wavetoAdd))-1, Test
	//	Test[PointtoStart,]=1

	Test[PointtoStart,]=WavetoAdd[p-PointtoStart]

	duplicate/o Test StartingWave
end

Static Function AppendToWaveWOverlap(StartingWave,WavetoAdd)

	wave StartingWave,WavetoAdd
	
	duplicate/free startingwave test
	insertpoints (numpnts(startingwave)-1), (numpnts(wavetoAdd)-1), Test

	Test[numpnts(startingwave)-1,]=WavetoAdd[p-numpnts(startingwave)+1]
	duplicate/o Test StartingWave

end


static function grabIndices(ForceWave,SaveWave)
	wave ForceWave,SaveWave
	string Combined="0;"+stringbykey("DE_Ind",note(ForceWave),":","\r")
	make/o/n=0 $nameofwave(SaveWave)
	variable n
	for(n=0;n<itemsinlist(Combined);n+=1)
		if(str2num(StringFromList(n, Combined))>numpnts(ForceWave))
		return 0
		endif

		InsertPoints n,1,Savewave
		Savewave[n]=str2num(StringFromList(n, Combined))
	endfor

End
static function grabPauseIndices(ForceWave,SaveWave)
	wave ForceWave,SaveWave
	string Combined="0;"+stringbykey("DE_PauseLoc",note(ForceWave),":","\r")
	make/free/n=0 Freesave
	
	variable n
	for(n=1;n<itemsinlist(Combined);n+=2)
		if(str2num(StringFromList(n, Combined))>numpnts(ForceWave))
		return 0
		endif

		InsertPoints n,1,Freesave
		Freesave[(n-1)/2]=str2num(StringFromList(n, Combined))
	endfor
	duplicate/o Freesave SaveWave
End

Static Function StateAtPoint(Point,CombinedWave)
	variable Point
	Wave CombinedWave
	
	make/free/n=(dimsize(CombinedWave,0)) Locations,States
	
	Locations=CombinedWave[p][0]
	if(Point==0)
		return -1
	endif
	
	FindLevel/Q Locations,Point
	if(floor(V_levelX)==0)
		return -1
	endif
	
	variable before=floor(V_LevelX)
	if(CombinedWave[before][2]==-2||CombinedWave[before][2]==2)
		before-=1
	endif
		
	if(CombinedWave[before][2]==-1)
		return 1
	else
		return -1
	endif

end


Static Function MaxRupForce(SingleKey)

	Wave SingleKey
	make/free/n=(dimsize(Singlekey,0)) Force, Type
	Force=Singlekey[p][1]
	Type=Singlekey[p][2]
	
	Extract/Free Force, ForceOut, Type[p]==-1
	waveStats/q ForceOut
	return v_max

end

Static Function/C SeparateintoStates(ForceWaveSM,SepWaveSm,CombinedWave,MinSpacing,FForce,FSep,UForce,USep)
	wave ForceWaveSM, SepWaveSm,CombinedWave,FSep,USep,FForce,UForce
	variable MinSpacing
	variable n,FoldedMax,UnfoldedMax,q
	
	make/free/n=(dimsize(CombinedWave,0)) Points,RupForce,Type,Trace
	make/free/n=0 HoldOnToCF,HoldOnToCU,HoldOnToF,HoldOnToOff,Unfolded,UnfoldedSep,Folded,FoldedSep,HoldN
	Points=CombinedWave[p][0]
	RupForce=CombinedWave[p][1]
	Type=CombinedWave[p][2]
	Trace=CombinedWave[p][3]
	
	variable fitstart=SepWaveSm[0]-3e-9,firstfit,FOFF
	variable lastpoint=0,firstone,firstmone,NewZero
	for(n=0;n<numpnts(points);n+=1)
		if(Type[n]==-1)


			if((Points[n]-lastpoint)<=0)

			else
				if(firstmone==0)

					duplicate/free/r=[NewZero,Points[n]] ForceWaveSM, Folded
					duplicate/free/r=[NewZero,Points[n]] SepWaveSm, FoldedSep

					firstmone=1
				else		
	
					duplicate/free/r=[lastpoint,Points[n]] ForceWaveSM, NewSection
					//NewSection-=ForceShift
					duplicate/free/r=[lastpoint,Points[n]] SepWaveSm, NewSectionSep

					AppendToWaveWOverlapOW(Folded,NewSection,lastpoint)
					AppendToWaveWOverlapOW(FoldedSep,NewSectionSep,lastpoint)

				endif
			endif
			
			lastpoint=Points[n]+minspacing
		endif	
	
		if(Type[n]==1)

			if((Points[n]-lastpoint)<=0)
			else
				if(firstone==0)
					//duplicate/free/r=[lastpoint,Points[n]] ForceWaveSM, UnFolded
					//UnFolded[0,lastpoint-NewZero]=NaN
					//duplicate/free/r=[lastpoint,Points[n]] SepWaveSm, UnFoldedSep
					//UnFoldedSep[0,lastpoint-NewZero]=NaN
					duplicate/free/r=[0,Points[n]] ForceWaveSM, UnFolded
					UnFolded[0,lastpoint]=NaN
					duplicate/free/r=[0,Points[n]] SepWaveSm, UnFoldedSep
					UnFoldedSep[0,lastpoint]=NaN
					firstone=1
				else		

					duplicate/free/r=[lastpoint,Points[n]] ForceWaveSM, NewSection

					duplicate/free/r=[lastpoint,Points[n]] SepWaveSm, NewSectionSep

					AppendToWaveWOverlapOW(UnFolded,NewSection,lastpoint)
					AppendToWaveWOverlapOW(UnFoldedSep,NewSectionSep,lastpoint)

				endif
				
			endif
			lastpoint=Points[n]+minspacing
		endif
		
	endfor
	if((Points[n-1]-lastpoint)<=0)

	else
		duplicate/free/r=[lastpoint,Points[n-1]] ForceWaveSM, NewSection
		duplicate/free/r=[lastpoint,Points[n-1]] SepWaveSm, NewSectionSep
		AppendToWaveWOverlap(Folded,NewSection)
		AppendToWaveWOverlap(FoldedSep,NewSectionSep)
	endif
//	FuncFit/Q/H="10111"/NTHR=0 WLC_FIT W_coef  folded /X=FoldedSep
	duplicate/o FoldedSep FSep
	duplicate/o UnFoldedSep USep
	duplicate/o Folded FForce
	duplicate/o UnFolded UForce
end

Menu "Dudko"
	"Open Dudko", dudko_panel()
end
//Static Function MinFoldForce(SingleKey)
//
//	Wave SingleKey
//	make/free/n=(dimsize(Singlekey,0)) Force, Type
//	Force=Singlekey[p][1]
//	Type=Singlekey[p][2]
//	
//	Extract/Free Force, ForceOut, Type[p]==1
//	waveStats/q ForceOut
//	return v_min
//
//end
//
//Static Function MakeMultiWLCFvT(AllWLCParms,SepWave,States,OutWLC1,OutWLC2)
//	wave AllWLCParms,SepWave,OutWLC1,OutWLC2,States
//	duplicate/free SepWave WLC1,WLC2
//	//WLC1=WLC(SepWave-WLCParms[3],.4e-9,WLCParms[0],298)
//	//WLC2=WLC(SepWave-WLCParms[3],.4e-9,WLCParms[1],298)
//	make/free/n=4 WLCParms
//	variable n,currcut,prevcut=0,i
//	for(n=1;n<dimsize(States,0);n+=1)
//		if(States[n][2]==2)
//		print i
//		WLCParms=AllWLCParms[i][p]
//			currcut=States[n][0]
//			//FOFF=ReturnOffForceForState(ForceIn,States[n][3]-1)
//			WLC1[prevcut,currcut]=WLC(SepWave-WLCParms[3],.4e-9,WLCParms[0],298)-WLCParms[2]
//			WLC2[prevcut,currcut]=WLC(SepWave-WLCParms[3],.4e-9,WLCParms[1],298)	-WLCParms[2]	
//			prevcut=States[n][0]
//			i+=1
//		
//		endif
//	
//	endfor
//		
//	duplicate/o WLC1 OutWLC1
//	duplicate/o WLC2 OutWLC2
//
//end
//Static Function FindSlopeNearRupture(ForceWave,CombinedStates,[smoothing])
//	wave ForceWave,CombinedStates
//	
//	variable Smoothing
//	duplicate/free ForceWave Fw
//	FW*=-1
//
//	if( ParamIsDefault(Smoothing) )
//	elseif(Smoothing<5)
//		print "Smoothing <5 not allowed"
//	elseif(2*round((smoothing)/2)==smoothing)
//		Smoothing=2*round((smoothing)/2)+1
//		print "Smoothing moved to nearest Odd number "+num2str(smoothing)
//		Smooth/S=2 Smoothing, Fw
//	else
//		Smooth/S=2 Smoothing, Fw
//
//	endif
//	
//	
//	variable n,StartIndex,EndIndex,i,j
//	make/o/n=(0,4) UpResults,DownResults
//	variable preferedsize=3000
//	For(n=0;n<dimsize(CombinedStates,0);n+=1)
//		if(CombinedStates[n][2]==1)
//			insertpoints i,1, UpResults
//			UpResults[i][0]=CombinedStates[n][0]
//			UpResults[i][1]=FW[CombinedStates[n][0]]
//
//			EndIndex=CombinedStates[n][0]
//			StartIndex=EndIndex-preferedsize
//		
//			if(StartIndex<(CombinedStates[n-1][0]+50))
//				startIndex=(CombinedStates[n-1][0]+50)
//			endif
//			duplicate/free/r=[startIndex,Endindex] FW, FWFit
//			CurveFit/Q/NTHR=0 line  FWFit /D 
//			wave w_coef
//			UpResults[i][2]=w_coef[1]
//			UpResults[i][3]=numpnts(FWFIT)
//		
//			i+=1
//		
//		elseif(CombinedStates[n][2]==-1)
//			insertpoints i,1, DownResults
//			DownResults[i][0]=CombinedStates[n][0]
//			DownResults[i][1]=FW[CombinedStates[n][0]]
//
//			EndIndex=CombinedStates[n][0]
//			StartIndex=EndIndex-preferedsize
//		
//			if(StartIndex<(CombinedStates[n-1][0]+50))
//				startIndex=(CombinedStates[n-1][0]+50)
//			endif
//			duplicate/free/r=[startIndex,Endindex] FW, FWFit
//			CurveFit/Q/NTHR=0 line  FWFit /D 
//			wave w_coef
//			DownResults[i][2]=w_coef[1]
//			DownResults[i][3]=numpnts(FWFIT)
//
//			j+=1
//		endif
//	endfor
//end

//
//Static Function NewSlopeatForce(ForceWave,Force,CombinedKey,DesiredStates,Direction)
//	wave ForceWave,CombinedKey	
//	variable Force,DesiredStates,Direction
//	make/o/n=0 InvestigateWaveStart,ResultWave,InvestigateWave
//	AssembleAllStatesNaN(Forcewave,CombinedKey,DesiredStates,InvestigateWaveStart)
//
//	Extract /o InvestigateWaveStart, InvestigateWaveStart, numtype(InvestigateWaveStart)==0
//	setscale/P x pnt2x(Forcewave,0), dimdelta(Forcewave,0), InvestigateWaveStart
//
//	variable smoothing=2001
//	do
//	duplicate/o InvestigateWaveStart InvestigateWave
//	Smooth/E=0/S=2 smoothing, InvestigateWave
//	
//	if(Direction==1)
//		FindLevels/edge=1/P/Q InvestigateWave,Force
//
//	
//	elseif(Direction==0)
//		FindLevels/edge=2/P/Q InvestigateWave,Force
//
//	else
//	endif
//	
//	wave W_FindLevels
//	
//	smoothing=2*round(smoothing/4)+1
//	
//	
//	while(numpnts(W_FindLevels)==0&&smoothing>51)
//	duplicate/o W_FindLevels W_Values
//	W_Values=InvestigateWave[W_FindLevels]
//	variable n,Center,PointShift,Pointstart,PointEnd
//	for(n=0;n<1;n+=1)
//	//for(n=0;n<numpnts(w_FindLevels);n+=1)
//		Center=W_FindLevels[n]
//		PointShift=500
//		if((Center-PointShift)<0)
//			Pointstart=0
//		else
//			PointStart=Center-PointShift
//		endif
//		if((Center)>numpnts(InvestigateWave))
//			Pointend=numpnts(InvestigateWave)
//		else
//			PointEnd=Center
//		endif
//
//		duplicate/o/r=[PointStart,PointEnd] InvestigateWave,GarbagetoLookat,GarbagetoLookatSm
//		Smooth/S=2 11, GarbagetoLookatSm
//		CurveFit/Q/NTHR=0 line  GarbagetoLookatSm /D 	
//		wave w_Coef
//		insertpoints n,1,ResultWave
//		ResultWave[n]=w_coef[1]
//	endfor
//	
//	if(numpnts(ResultWave)==0)
//		print "Hey, we didn't find any crossings" 
//		return 0
//	endif
//		
//		
//	wavestats/Q ResultWave
//		
//	for(n=0;n<numpnts(ResultWave);n+=1)
//
//		if(abs(ResultWave[n])<abs(v_avg)/20)
//			deletepoints n,1, ResultWave
//			wavestats/Q ResultWave
//			n-=1
//		endif
//	endfor
//	wavestats/Q ResultWave
//
//	for(n=0;n<numpnts(ResultWave);n+=1)
//
//		if(sign(ResultWave[n])!=sign(v_avg))
//			deletepoints n,1, ResultWave
//			wavestats/Q ResultWave
//			n-=1
//		endif
//	endfor
//	wavestats/Q ResultWave
//
//	if( abs(v_sdev/v_avg)>.25)
//		print "Warning, sdev ratio= "+num2str(abs(v_sdev/v_avg))
//	endif
//	return v_avg
//end
//
//Static Function FindSlopeAtForce(ForceWave,Force,RuptureWaveLeaving,RuptureWaveEntering,[num])
//	Wave ForceWave,RuptureWaveLeaving,RuptureWaveEntering
//	Variable Force
//	variable num
//	variable startn,endn, n,index,m
//	
//	make/free/n=0 ResultWave,TestWave
//
//	DE_Dudko#grabIndices(ForceWave,TestWave)
//	
//	if(ParamisDefault(num)==1)
//		startn=0
//		endn=numpnts(RuptureWaveLeaving)
//	else
//		startn=num
//		endn=num+1
//	endif
//
//	for(n=startn;n<endn;n+=1)
//		index=RuptureWaveLeaving[n]
//		variable PreviousIndex=FindPreceedingFromOtherList(index,RuptureWaveEntering,TestWave)
//		duplicate/free/r=[PreviousINdex,Index] ForceWave,GarbagetoLookat,GarbagetoLookatSm
//		if(numpnts(GarbagetoLookat)<2100)
//	
//		else 
//			Smooth/S=2 501, GarbagetoLookatSm
//			FindLevels/P/Q GarbagetoLookatSm,Force
//	
//			wave W_FindLevels
//
//			if(V_levelsFound==1)
//				variable Center=W_FindLevels[0]
//				variable Pointstart,PointEnd,PointShift=500
//				if((Center-PointShift)<0)
//					Pointstart=0
//				else
//					PointStart=Center-PointShift
//				endif
//				if((Center+PointShift)>numpnts(GarbagetoLookatsm))
//					Pointend=numpnts(GarbagetoLookatsm)
//				else
//					PointEnd=Center+PointShift
//				endif
//	
//				duplicate/free/r=[PreviousINdex,Index] ForceWave,GarbagetoLookat,GarbagetoLookatSm
//				Smooth/S=2 11, GarbagetoLookatSm
//
//				duplicate/o/r=[PointStart,PointEnd] GarbagetoLookatsm,GarbagetoLookata
//
//				CurveFit/Q/NTHR=0 line  GarbagetoLookata /D 	
//				wave w_Coef
//				insertpoints m,1,ResultWave
//				ResultWave[m]=w_coef[1]
//				m+=1
//			endif
//		endif	
//
//	endfor
//	if(numpnts(ResultWave)==0)
//		print "Hey, we didn't find any crossings" 
//		return 0
//	endif
//	wavestats/Q ResultWave
//	
//	if( abs(v_sdev/v_avg)>.25)
//		print "Warning, sdev ratio= "+num2str(abs(v_sdev/v_avg))
//	endif
//	return v_avg
//end

//static function ExtractRuptureForces(ForceWave,RupWave1,RupWave2,RupOut1,RupOut2,[Smoothing])
//
//	wave ForceWave,RupWave1,RupWave2,RupOut1,RupOut2
//	variable Smoothing
//	duplicate/free ForceWave Fw
//
//	if( ParamIsDefault(Smoothing) )
//	elseif(Smoothing<5)
//		print "Smoothing <5 not allowed"
//	elseif(2*round((smoothing)/2)==smoothing)
//		Smoothing=2*round((smoothing)/2)+1
//		print "Smoothing moved to nearest Odd number "+num2str(smoothing)
//		Smooth/S=2 Smoothing, Fw
//	else
//		Smooth/S=2 Smoothing, Fw
//
//	endif
//	Fw*=-1
//	duplicate/free RupWave1 RupForce1
//	duplicate/free RupWave2 RupForce2
//	RupForce1=Fw[RupWave1[p]]
//	RupForce2=Fw[RupWave2[p]]
//
//	make/free/n=(numpnts(RupWave1),2) Test1
//	make/free/n=(numpnts(RupWave2),2) Test2
//
//	Test1[][0]=RupWave1
//	Test1[][1]=Fw[RupWave1[p]]
//	Test2[][0]=RupWave2
//	Test2[][1]=Fw[RupWave2[p]]
//
//	duplicate/o Test1 RupOut1
//	duplicate/o Test2 RupOut2
//	variable n=10
////	make/o/n=(n) $(nameofwave(RupOut1)+"_Hist"),$(nameofwave(RupOut2)+"_Hist")
////	wave Hist1= $(nameofwave(RupOut1)+"_Hist")
////	wave Hist2=$(nameofwave(RupOut2)+"_Hist")
//
//
//
////	Histogram/C/B=1 RupForce1,Hist1
////	Histogram/C/B=1 RupForce2,Hist2
////	Hist1/=dimdelta(Hist1,0)
////	Hist2/=dimdelta(Hist2,0)
//
//end


//Static Function FindPreceedingFromOtherList(Point,PreceedingWave,EdgeWave)
//	variable Point
//	wave PreceedingWave,EdgeWave
//	variable v1,v2
//
//	FindLevel/Q PreceedingWave,Point
//	if(numtype(v_levelx)==2)
//		v1= 0
//	else
//		v1= PreceedingWave[floor(V_levelx)]
//	endif
//	FindLevel/Q EdgeWave,Point
//	if(numtype(v_levelx)==2)
//		v2= 0
//	else
//		v2= EdgeWave[floor(V_levelx)]
//	endif
//
//	return max(v1,v2)
//
//end
//
//

//
//Static Function AssembleAllStatesNan(Forcewave,CombinedKey,DesiredStates,AllDesiredStates)
//	
//	wave Forcewave,CombinedKey,AllDesiredStates
//	variable DesiredStates
//	
//	variable n
//	variable previousindex=0,deletelength
//	duplicate/o ForceWave AllDesiredStates
//	for(n=0;n<dimsize(CombinedKey,0);n+=1)
//		if((CombinedKey[n][1]==-2)||(CombinedKey[n][1]==2))
//		elseif(CombinedKey[n][1]==DesiredStates)
//
//		//	AllDesiredStates[PreviousIndex,CombinedKey[n][0]]=NaN
//			previousindex=CombinedKey[n][0]
//		else
//			AllDesiredStates[PreviousIndex,CombinedKey[n][0]+100]=NaN
//			previousindex=CombinedKey[n][0]+100
//		endif
//	
//	endfor
//end


//Static Function MakeTwoStateForceMask(Forcewave,CombinedKey,Mask)
//	wave Forcewave,CombinedKey,Mask
//	variable n
//	variable previousindex=0
//	duplicate/o ForceWave Mask
//	for(n=0;n<dimsize(CombinedKey,0);n+=1)
//		if(CombinedKey[n][1]==-1)
//			Mask[PreviousIndex,CombinedKey[n][0]]=-1
//			previousindex=CombinedKey[n][0]
//		endif
//		
//		if(CombinedKey[n][1]==1)
//			Mask[PreviousIndex,CombinedKey[n][0]]=1
//			previousindex=CombinedKey[n][0]
//		endif
//	
//	endfor
//end
//Static Function/C NewCLSlope(ForceWaveSM,SepWaveSm,CombinedWave,Force,State,Smoothed,rate)
//	wave ForceWaveSM, SepWaveSm,CombinedWave
//	variable Force,State,Smoothed,rate
//	variable LC,OffSet
//	make/free/n=3 Results
//	ContourLengthDetermine(ForceWaveSM,SepWaveSm,CombinedWave,Smoothed,Results)
//	if(State==-1)
//		LC=Results[0]
//		
//	else
//		LC=Results[1]
//	endif
//	OffSet=Results[2]
//	return ReturnSlopesFromWLC(LC,Offset,Force,rate)
//end
//Static Function ReturnSlopesFromWLC(LC,Offset,Force,rate)
//	variable LC,Offset,Force,rate
//	make/free/n=1000 WLCX,WLCF
//	WLCX=2e-7/999*x
//	WLCF=WLC(WLCX,0.4e-9,LC,298)-Offset
//	FindLevel/Q WLCF(Force) 
//	variable StartPnt,EndPnt
//	if(numtype(V_levelX)==2)
//	return NaN
//	
//	
//	else
//	StartPnt=max(v_levelx-10,0)
//	EndPnt=min(v_levelx+10,numpnts(WLCF)-1)
//	return (WLCF[EndPnt]-WLCF[StartPnt])/(WLCX[EndPnt]-WLCX[StartPnt])*rate
//	endif
//end


//Static Function Generate(ForceWave,UpPoints,DownPoints)
//	wave ForceWave,UpPoints,DownPoints
//	
//	duplicate/o ForceWave SmoothForce
//	string SmoothingType="TVD"
//	variable SmoothingValue=50e-12
//	StrSwitch (SmoothingType)
//	case "TVD":
//		DE_Filtering#TVD1D_denoise(ForceWave,SmoothingValue,SmoothForce)
//	break
//	
//	case "SVG":
//		Smooth/S=2 SmoothingValue,SmoothForce;
//
//	break
//	
//	endswitch
//	SmoothForce*=-1 
//	
//	make/o/n=0 $(nameofwave(ForceWave)+"_State")
//	wave StateWave=$(nameofwave(ForceWave)+"_State")
//	DE_DUDKO#MakeSingleStateKey(SmoothForce,UpPoints,DownPoints,StateWave)
//	make/o/n=(40,5) $(nameofwave(ForceWave)+"_Ref")
//	wave RefWave=$(nameofwave(ForceWave)+"_Ref")	
//	RefWave=0
//	RefWave[][0]=1.1*DE_DUDKO#MaxRupForce(StateWave)/39*p
//	RefWave[][1]=real(DE_DUDKO#AverageSlopeAtForce(SmoothForce,StateWave,(RefWave[p][0]),-1,1))	
//	RefWave[][2]=real(DE_DUDKO#AverageSlopeAtForce(SmoothForce,StateWave,(RefWave[p][0]),1,1))	
//	
//	RefWave[][3]=real(DE_DUDKO#NumberofTracesinState(SmoothForce,StateWave,(RefWave[p][0]),-1,1))	
//	RefWave[][4]=real(DE_DUDKO#NumberofTracesinState(SmoothForce,StateWave,(RefWave[p][0]),1,1))
//	
//
//
//	
//	DE_DUdko#ConstructHistograms(StateWave,RefWave,10)
//end

//Static Function/C AverageSlopeAtForce(ForceWaveSM,CombinedWave,Force,State,Smoothed)
//	wave ForceWaveSM, CombinedWave
//	variable Force,State,Smoothed
//	variable n,FoldedMax,UnfoldedMax,q
//	String StateList=""
//	Variable TimeStep=0.05
//	Variable StartTime,EndTime,EventBefore,EventAfter
//	If(Smoothed==1)
//		duplicate/free ForceWaveSM FWAll
//	else
//		duplicate/free ForceWaveSM FWAll
//		Smooth/S=2 101, FWALL
//		FWALL*=-1
//	endif
//	
//	make/o/n=(dimsize(CombinedWave,0)) Points,RupForce,Type,Trace
//	make/o/n=0 HoldOnToSlopes
//	Points=CombinedWave[p][0]
//	RupForce=CombinedWave[p][1]
//	Type=CombinedWave[p][2]
//	Trace=CombinedWave[p][3]
//	//for(n=1;n<2;n+=1)
//	for(n=0;n<Trace[numpnts(Trace)-1][2];n+=1)
//		Extract/INDX/FREE Points, LocalSet, Trace==n
//		make/o/n=(numpnts(LocalSet)+1) LocalPoints,LocalType,LocalTrace
//		LocalPoints[]=Points[LocalSet[0]+p][0]
//		LocalType[]=Type[LocalSet[0]+p][2]
//		LocalTrace[]=Trace[LocalSet[0]+p][2]
//		FindValue/V=-2/T=.01 LocalType
//		duplicate/free/r=[,V_value] LocalPoints UpDir
//		duplicate/free/r=[V_value+1,] LocalPoints DownDir
//		duplicate/o/r=[UpDir[0],UpDir[numpnts(UpDir)-1]] FWALL FW
//		FindLevels/edge=1/Q FW,Force		
//		wave W_FindLevels
//
//		if(numpnts(W_FindLevels)==0)
//			print "No Levels: "+num2str(Force)
//		else
//			duplicate/FREE W_FindLevels W_FindValues
//			W_FindValues=FW(W_FindLevels)
//			//for(q=1;q<2;q+=1)
//			for(q=0;q<numpnts(W_FindLevels);q+=1)
//				if(State==StateAtPoint(x2pnt(FWALL,W_FindLevels[q]),CombinedWave))
//					EventBefore=EventBeforeTime(W_FindLevels[q],LocalPoints,LocalType,FWAll)
//					EventAfter=EventAfterTime(W_FindLevels[q],LocalPoints,LocalType,FWAll)
//					if(State==-1)
//						StartTime=W_FindLevels[q]-TimeStep/2
//						EndTime=W_FindLevels[q]+TimeStep/2
//						
//						if(StartTime<(EventBefore+0.001))
//							EndTime+=((EventBefore+0.001)-StartTime)
//					
//							StartTime=(EventBefore+0.001)
//						endif
//						
//						if(EndTime>(EventAfter-0.001))
//							EndTime=EventAfter-0.001
//
//						endif
//						duplicate/o/r=(StartTime,EndTime) FW FWFit
//					
//					elseif(State==1)
//					
//						StartTime=W_FindLevels[q]-TimeStep/2
//						EndTime=W_FindLevels[q]+TimeStep/2
//						
//						if(StartTime<(EventBefore+0.001))
//							EndTime+=((EventBefore+0.001)-StartTime)
//							StartTime=(EventBefore+0.001)
//						endif
//						
//						if(EndTime>(EventAfter-0.001))
//							EndTime=EventAfter-0.001
//
//						endif
//						duplicate/o/r=(StartTime,EndTime) FW FWFit
//					endif
//					CurveFit/Q/NTHR=0 line  FWFit 
//					wave w_coef
//					if((Endtime-Starttime)<0.001)
//					else
//						insertpoints 0,1, HoldOnToSlopes
//						HoldOnToSlopes[0]=w_coef[1]
//					endif
//				endif
//
//			endfor
//		endif
//		
//
//				
//		duplicate/o/r=[UpDir[numpnts(UpDir)-1],DownDir[numpnts(DownDir)-1]] FWALL FW
//				
//		FindLevels/edge=2/Q FW,Force		
//		wave W_FindLevels
//				
//		if(numpnts(W_FindLevels)==0)
//		else
//			duplicate/o W_FindLevels W_FindValues
//			W_FindValues=FW(W_FindLevels)
//			//for(q=1;q<2;q+=1)
//			for(q=0;q<numpnts(W_FindLevels);q+=1)
//				if(State==StateAtPoint(x2pnt(FWALL,W_FindLevels[q]),CombinedWave))
//					EventBefore=EventBeforeTime(W_FindLevels[q],LocalPoints,LocalType,FWAll)
//					EventAfter=EventAfterTime(W_FindLevels[q],LocalPoints,LocalType,FWAll)
//					if(State==-1)
//						StartTime=W_FindLevels[q]-TimeStep/2
//						EndTime=W_FindLevels[q]+TimeStep/2
//								
//						if(StartTime<(EventBefore+0.001))
//							EndTime+=((EventBefore+0.001)-StartTime)
//							StartTime=(EventBefore+0.001)
//						endif
//							
//						if(EndTime>(EventAfter-0.001))
//		
//							EndTime=EventAfter-0.001
//		
//						endif
//		
//						duplicate/o/r=(StartTime,EndTime) FW FWFit
//					elseif(State==1)
//						StartTime=W_FindLevels[q]-TimeStep/2
//						EndTime=W_FindLevels[q]+TimeStep/2
//							
//						if(StartTime<(EventBefore+0.001))
//							EndTime+=((EventBefore+0.001)-StartTime)
//							StartTime=(EventBefore+0.001)
//						endif
//								
//						if(EndTime>(EventAfter-0.001))
//							EndTime=EventAfter-0.001
//						endif
//						duplicate/free/r=(StartTime,EndTime) FW FWFit
//					endif
//					CurveFit/Q/NTHR=0 line  FWFit
//					wave w_coef
//					if((Endtime-Starttime)<0.001)
//					else
//						insertpoints 0,1, HoldOnToSlopes
//						HoldOnToSlopes[0]=-1*w_coef[1]
//					endif
//				endif
//			endfor
//		endif
//					
//	
//	endfor
//	wave w_sigma
//	killwaves W_FindLevels,w_coef,w_sigma
//	 
//	 
//	if( numpnts(Holdontoslopes)==0)
//		return NaN
//	endif
//	wavestats/q HoldOnToSlopes
//	//print v_Sdev
//	return cmplx(v_avg,v_sdev)
//	
//end
//
//Static Function EventAfterTime(CurrTime,Points,Type,ForceWave)
//	Variable CurrTime
//	wave Points,Type,ForceWave
//	variable n
//	for(n=0;n<numpnts(Points);n+=1)
//		if(CurrTime<pnt2x(Forcewave,Points[n]))
//			return pnt2x(Forcewave,Points[(n)])
//		endif
//	endfor
//end
//
//Static Function EventBeforeTime(CurrTime,Points,Type,ForceWave)
//	Variable CurrTime
//	wave Points,Type,ForceWave
//	variable n
//	for(n=0;n<numpnts(Points);n+=1)
//		if(CurrTime<pnt2x(Forcewave,Points[n]))
//			return pnt2x(Forcewave,Points[(n-1)])
//		 
//		endif
//	endfor
//end
//
//Static Function ConstructHistograms(StateWave,RefWave,bins,[NameBase])
//	wave StateWave,RefWave
//	variable bins
//	string NameBase
//	variable n,Location
//	make/o/n=(bins) FHist,UHist
//	make/free/n=(bins) NumberofSlopes,TestSlopes
//
//	make/free/n=(dimsize(StateWave,0)) JustForce,StateImIn
//
//	JustForce=StateWave[p][1]
//	StateImIn=StateWave[p][2]
//
//	//make/free/n=0 UnfoldingPoints,FoldingPoints
//
//	Extract/free JustForce, UnFoldingForces, StateImIn[p]==-1
//	Extract/free JustForce, FoldingForces, StateImIn[p]==1
//
//
//	Histogram/C/B=1 UnFoldingForces UHist
//	Histogram/C/B=1 FoldingForces FHist
//
//	duplicate/o FHist FSlope, FNum
//	duplicate/o UHist USlope, UNum
//	
//	make/free/n=(dimsize(RefWave,0)) AllForces,UnAllSlopes,UnAllNumbers,FAllSlopes,FAllNumbers
//	
//	AllForces=real(RefWave[p][0])
//	UnAllSlopes=real(RefWave[p][1])
//	UnAllNumbers=real(RefWave[p][3])
//	FAllSlopes=real(RefWave[p][2])
//	FAllNumbers=real(RefWave[p][4])
//	for(n=0;n<numpnts(UHist);n+=1)
//		FindLevel/Q/p AllForces pnt2x(UHist,n)
//		if(numtype(v_levelx)!=0)
//		print pnt2x(UHist,n)
//		USlope[n]=Nan
//		UNum[n]=Nan
//		else
//		USlope[n]=UnAllSlopes[v_levelx]
//		UNum[n]=UnAllNumbers[v_levelx]
//
//		endif
//	endfor
//	
//		for(n=0;n<numpnts(FHist);n+=1)
//		FindLevel/Q/p AllForces pnt2x(FHist,n)
//		if(numtype(v_levelx)!=0)
//				print pnt2x(UHist,n)
//
//		FSlope[n]=Nan
//		FNum[n]=Nan
//		else
//		FSlope[n]=FAllSlopes[v_levelx]
//		FNum[n]=FAllNumbers[v_levelx]
//		endif
//		
//	endfor
//	
//	Duplicate/o UHist URate
//	Duplicate/o FHist FRate
//	URate=UHist/dimdelta(UHist,0)*USlope/UNum
//	FRate=FHist/dimdelta(FHist,0)*FSlope/FNum
//	
//	if(ParamisDefault(NameBase))
//	
//	else
//		duplicate/o URate $(Namebase+"_URate")
//		duplicate/o FRate $(Namebase+"_FRate")
//
//		duplicate/o UHist $(Namebase+"_UHist")
//		duplicate/o FHIST $(Namebase+"_FHist")
//		
//		duplicate/o USlope $(Namebase+"_USlope")
//		duplicate/o FSlope $(Namebase+"_FSlope")		
//		
//		duplicate/o UNum $(Namebase+"_UNum")
//		duplicate/o FNum $(Namebase+"_FNum")
//		
//		killwaves URate,FRate,UHist,FHIST,USlope,FSlope,UNum,FNum
//	endif
//end

//Static Function CalculateBestForceOffset(ForceWave1,ForceWave2,[startindex,endindex])
//	wave ForceWave1,ForceWave2
//	variable startindex,endindex
//
//	if(paramisdefault(startindex)||paramisdefault(endindex))
//		startindex=0
//		endindex=numpnts(ForceWave1)-1
//	endif
//	variable K1 = 1;
//	CurveFit/Q/H="01"/NTHR=0 line  ForceWave1[startindex,endindex] /X=ForceWave2  
//	wave w_coef
//	return w_coef[0]
//end

//
//Static Function WhatBinDoesThisGoIN(HistogramWave,Force)
//	wave HistogramWave
//	variable Force
//
//	if(Force>(pnt2x(Histogramwave,numpnts(Histogramwave)-1)+dimdelta(Histogramwave,0)*(.5)))
//		print "Force too high"
//		return 0
//	
//	elseif(Force<(pnt2x(Histogramwave,0)-dimdelta(Histogramwave,0)*(.5)))
//		print "Force too low"
//		return 0
//	endif
//
//	make/free/n=(numpnts( HistogramWave)+1) XValues
//	Xvalues=pnt2x(Histogramwave,0)+dimdelta(Histogramwave,0)*(p-.5)
//	FindLevel/Q  Xvalues,Force
//	if(numtype(v_levelx)!=0)
//		print Force
//		print "HUH"
//	endif
//	if(floor(v_levelx)>9)
//	return (floor(v_levelx))-1
//	endif
//	return (floor(v_levelx))
//end