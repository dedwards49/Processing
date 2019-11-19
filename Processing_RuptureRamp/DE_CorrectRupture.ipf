#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma ModuleName=DE_CorrRup
Static Function FixPicks(UpWave,ForceWave,SmoothFW,WLC,SepWave)
	wave Upwave,ForceWave,SepWave,WLC,SmoothFW
	wave DownWave=$ReplaceString("PntU", nameofwave(Upwave), "PntD")
	variable n
	make/D/o/n=(numpnts(DownWave),6) DownCrap
	make/D/o/n=(numpnts(Upwave),6) UpCrap

	variable prepnt
	variable nextpnt
	//	if(Direction==-1)
	//		wave OtherPnt=$replacestring("PntD",nameofwave(pntWave),"PntU")

	//	elseif(Direction==1)
	//		wave OtherPnt=$replacestring("PntU",nameofwave(pntWave),"PntD")

	//		
	//	endif


	make/o/n=0 Test

	for(n=0;n<numpnts(DownWave);n+=1)
		prepnt=Upwave[n]
		if(n==(numpnts(DownWave)-1))
			nextpnt=numpnts(ForceWave)
		else
			nextpnt=Upwave[n+1]

		endif

		//DE_CorrRup#CalculateBestCrossingGuess(ForceWave,SepWave,WLC,DownWave[n],prepnt,nextpnt,-1,Test)
		if((DE_CorrRup#CalculateCrossingfromWLC(ForceWave,SepWave,WLC,101,DownWave[n],prepnt,nextpnt,-1,Test))==-1)
			print "WLC Down: "+num2str(n)
			DownCrap[n][0,2]=NaN
		else
			DownCrap[n][0,2]=Test[q]
		endif
		
		

		if((DE_CorrRup#CalculateCrossingLinear(ForceWave,SepWave,WLC,101,DownWave[n],prepnt,nextpnt,-1,Test))==-1)
			print "Line Down: "+num2str(n)
		
			DownCrap[n][3,5]=Nan

		else
			DownCrap[n][3,5]=Test[q-3]

		endif

	endfor
	
	for(n=0;n<numpnts(Upwave);n+=1)
		if(n==0)
			prepnt=0
			nextpnt=DownWave[n]
		else
			prepnt=DownWave[n-1]
			nextpnt=DownWave[n]
		endif
		if((DE_CorrRup#CalculateCrossingfromWLC(ForceWave,SepWave,WLC,101,Upwave[n],prepnt,nextpnt,1,Test))==-1)
			print "WLC UP: "+num2str(n)
			UpCrap[n][0,2]=Nan

		else
			UpCrap[n][0,2]=Test[q]

		endif


		if((DE_CorrRup#CalculateCrossingLinear(ForceWave,SepWave,WLC,101,Upwave[n],prepnt,nextpnt,1,Test))==-1)
			print "Line UP: "+num2str(n)
			UpCrap[n][3,5]=NaN
		else
			UpCrap[n][3,5]=Test[q-3]

		endif


	endfor
	wavestats/Q Test
	if(V_numNans>0)
		print "We missed something"
	endif
	make/o/D/n=(numpnts(Upwave)) $(nameofwave(UpWave)+"_Mod")
	wave w1=$(nameofwave(UpWave)+"_Mod")
	w1=UpCrap[p][5]
	//w1=(UpCrap[p][3]-dimoffset(ForceWave,0))/dimdelta(ForceWave,0)
	//w1=x2pnt(ForceWave,UpCrap[p][0])
	make/o/D/n=(numpnts(DownWave)) $(nameofwave(DownWave)+"_Mod")
	wave w2=$(nameofwave(DownWave)+"_Mod")
	w2=(DownCrap[p][5])
	//w2=(DownCrap[p][3]-dimoffset(ForceWave,0))/dimdelta(ForceWave,0)
	wavestats/Q w1
	if(V_numNans>0)
		print "We missed something"
	endif
	wavestats/Q w2
	if(V_numNans>0)
		print "We missed something"
	endif
	//w2=x2pnt(ForceWave,DownCrap[p][0])


	
end
//Can remove prevpoint

Static Function CalculateCrossingfromWLC(ForceWave,SepWave,WLCParms,FilterPoints,EstimatePoint,PrevPoint,NextPoint,TransType,Waveout,[CopyWavesOUt])
	wave ForceWave,SepWave,WLCParms,Waveout
	variable EstimatePoint,TransType,PrevPoint,NextPoint, CopyWavesOUt,FilterPoints
	variable startingscan,endingscan

 	variable/C Result
		variable edge //Depending on whether it's folding or unfolding we pick a direction to find levels
	if(TransType==1)
		edge=2
	elseif(transtype==-1)
		edge=1
	else
		return 0
	endif
	
	variable currenttime=pnt2x(ForceWave,EstimatePoint)
	variable nexttime=pnt2x(Forcewave,NextPoint)
	variable Prevtime=pnt2x(Forcewave,PrevPoint)

	//variable FirstWLC_backtime=0//2e-4
	variable FirstWLC_backtime=2e-4

	variable FirstWLC_Fwdtime=min(1e-3,(nexttime-currenttime)*.3)  //In the first WLC, look ahead EITHER 2 ms, or 30% to the next transition
	//variable SecondWLC_backtime=-0e-4
	variable SecondWLC_backtime=-2e-4

	variable SecondWLC_Fwdtime=min(2e-3,(nexttime-currenttime)*.6) //In the second WLC look ahead EITHER 2 ms, or 60% to the next transition
	variable Smooth_start=currenttime-1e-2
	variable Smooth_end=currenttime+1e-2
	
	duplicate/free/r=(Smooth_start,Smooth_end) ForceWave FtoSmth
	duplicate/free/r=(Smooth_start,Smooth_end) SepWave StoSmth
	make/o/n=0 ForceSm,SepSm

	DE_Filtering#FilterForceSep(FtoSmth,StoSmth,ForceSm,SepSm,"SVG",FilterPoints)

	
	
	variable FirstWLC_start=currenttime-FirstWLC_backtime
	variable FirstWLC_end=currenttime+FirstWLC_Fwdtime
	variable SecondWLC_start=currenttime-SecondWLC_backtime
	variable SecondWLC_end=currenttime+SecondWLC_Fwdtime



	duplicate/free/r=(FirstWLC_start,FirstWLC_end) ForceWave WLC1
	duplicate/free/r=(SecondWLC_start,SecondWLC_end) ForceWave WLC2

	duplicate/free/r=(FirstWLC_start,FirstWLC_end) ForceWave CutForceWave1
	duplicate/free/r=(SecondWLC_start,SecondWLC_end) ForceWave CutForceWave2
	duplicate/free/r=(FirstWLC_start,FirstWLC_end) SepWave CutSepWave1
	duplicate/free/r=(SecondWLC_start,SecondWLC_end) SepWave CutSepWave2
	
	duplicate/free/r=(FirstWLC_start,FirstWLC_end) ForceSm CutForceSm1
	duplicate/free/r=(SecondWLC_start,SecondWLC_end) ForceSm CutForceSm2
	duplicate/free/r=(FirstWLC_start,FirstWLC_end) SepSm CutsepSm1
	duplicate/free/r=(SecondWLC_start,SecondWLC_end) SepSm CutsepSm2

	duplicate/free/r=(FirstWLC_start-2e-2,FirstWLC_end+2e-2) ForceWave WideWLC1
	duplicate/free/r=(FirstWLC_start-2e-2,FirstWLC_end+2e-2) SepWave CutSepWaveWide

	if(TransType==1)//The first WLC is the unfolded state, the second is folded
		WLC1=WLC(CutSepWave1-WLCParms[3],.4e-9,WLCParms[1],298)-WLCParms[2]
		WLC2=WLC(CutSepWave2-WLCParms[3],.4e-9,WLCParms[0],298)-WLCParms[2]
		WideWLC1=WLC(CutSepWaveWide-WLCParms[3],.4e-9,WLCParms[1],298)-WLCParms[2]

	elseif(TransType==-1)//The first WLC is the folded state, the second is unfolded
		WLC1=WLC(CutSepWave1-WLCParms[3],.4e-9,WLCParms[0],298)-WLCParms[2]
		WLC2=WLC(CutSepWave2-WLCParms[3],.4e-9,WLCParms[1],298)-WLCParms[2]
		WideWLC1=WLC(CutSepWaveWide-WLCParms[3],.4e-9,WLCParms[0],298)-WLCParms[2]

	endif

	make/free/n=0 Garbage
	variable WLC1Slope=MakeLinearFitToWave(WLC1,Garbage,Garbage)
	variable WLC2Slope=MakeLinearFitToWave(WLC2,Garbage,Garbage)
	
//	DE_Filtering#FilterForceSep(CutForceWave2,CutSepWave2,CutForceSm2,CutsepSm2,"SVG",FilterPoints)

	duplicate/free CutForceSm1 ForcesShiftedByWLC1
	ForcesShiftedByWLC1-=wLC1
	duplicate/free CutForceSm2 ForcesShiftedByWLC2
	ForcesShiftedByWLC2-=wLC2

	make/free/n=0 WLC1Crossings,WLC2Crossings

	if(FindLevelsWithError(ForcesShiftedByWLC2,0,.00001,edge,WLC2Crossings)==-1)
		//print "WLC1 Error"
		FilterPoints=2*(round(FilterPoints/3))+1
		if(FilterPoints>5)
			Result=CalculateCrossingfromWLC(ForceWave,SepWave,WLCParms,FilterPoints,EstimatePoint,PrevPoint,NextPoint,TransType,Waveout,CopyWavesOUt=1)
			return Result
		else
		print "Filtering under 5"
		return -1
		endif

	endif
	duplicate/free WLC2Crossings,WLC2Forces
	WLC2Forces =CutForceSm2(WLC2Crossings)

	if(FindLevelsWithError(ForcesShiftedByWLC1,0,.00001,edge,WLC1Crossings)==-1)

		//print "WLC2 Error"
		FilterPoints=2*(round(FilterPoints/3))+1
		if(FilterPoints>5)
			Result=CalculateCrossingfromWLC(ForceWave,SepWave,WLCParms,FilterPoints,EstimatePoint,PrevPoint,NextPoint,TransType,Waveout,CopyWavesOUt=1)

			return Result
		else
		print "Filtering under 5"
		return -1
		endif
	endif
	duplicate/Free WLC1Crossings,WLC1Forces
	WLC1Forces =CutForceSm1(WLC1Crossings)

	EliminatePreCrossings(WLC1Crossings,WLC2Crossings)//Eliminate crossings in 2 that occur before
																		// any of the crossing in 1
																		
																		
	if(numpnts(WLC1Crossings)==0||numpnts(WLC2Crossings)==0)
		FilterPoints=2*(round(FilterPoints/3))+1
		if(FilterPoints>5)
			Result=CalculateCrossingfromWLC(ForceWave,SepWave,WLCParms,FilterPoints,EstimatePoint,PrevPoint,NextPoint,TransType,Waveout,CopyWavesOUt=1)

			return Result
		else
		print "Filtering under 5"
		return -1
		endif
	
	endif
	startingscan=max(LastPoint(WLC1Crossings,WLC2Crossings)-1e-4,currenttime) //start 0.1( ms before the last crossing in 1
																		//that occurs before any crossings in 2

	endingscan=NextPointinWave(WLC2Crossings,startingscan)+1e-4 //100 us after the first crossing of WLC2 
																				//after wherever we started

	duplicate/free/R=(startingscan,endingscan)  ForceWave RawForceShiftbyWLC1
	duplicate/free/R=(startingscan,endingscan)  WideWLC1 RawWLC

	RawForceShiftbyWLC1-=RawWLC
	make/free/n=0 RawWLCCrossings

	if(FindLevelsWithError(RawForceShiftbyWLC1,0,.00001,edge,RawWLCCrossings)==-1)
		print "Couldn't find crossing in raw wave"
		return 0
	endif
	
	RawForceShiftbyWLC1+=RawWLC
	
	variable Estimate 
	
	
	if(numpnts(RawWLCCrossings)==0)
	Estimate=Nan
	Result=cmplx(Estimate,NAN)
	else
	duplicate/free RawWLCCrossings RawWLCForce
	RawWLCForce=ForceWave(RawWLCCrossings)
	Estimate= RawWLCCrossings[numpnts(RawWLCCrossings)-1]

	Result=cmplx(Estimate,ForceWave(Estimate))
	endif
	make/D/free/n=4 ResultsOut
//	//	w1=(Estimate1-dimoffset(ForceWave,0))/dimdelta(ForceWave,0)
	ResultsOut={Estimate,ForceWave(Estimate),(Estimate-dimoffset(ForceWave,0))/dimdelta(ForceWave,0)}
	duplicate/o ResultsOut Waveout
	
	if(ParamisDefault(CopyWavesOUt)||CopyWavesOUt==0)
	
	else
		NewDataFolder/O root:CorrectRupture
		NewDataFolder/O root:CorrectRupture:WLC
		Duplicate/o WLC1 root:CorrectRupture:WLC:CR_WLC1
		Duplicate/o WLC2 root:CorrectRupture:WLC:CR_WLC2
		duplicate/o CutForceSm1 root:CorrectRupture:WLC:CR_CutForceSm1
		duplicate/o CutForceSm2 root:CorrectRupture:WLC:CR_CutForceSm2
		duplicate/o WLC2Crossings root:CorrectRupture:WLC:CR_WLC2Crossings
		duplicate/o WLC2Forces root:CorrectRupture:WLC:CR_WLC2Forces
		duplicate/o WLC1Crossings root:CorrectRupture:WLC:CR_WLC1Crossings
		duplicate/o WLC1Forces root:CorrectRupture:WLC:CR_WLC1Forces
		duplicate/o ForcesShiftedByWLC1 root:CorrectRupture:WLC:CR_ForcesShiftedByWLC1
		duplicate/o ForcesShiftedByWLC2 root:CorrectRupture:WLC:CR_ForcesShiftedByWLC2
		duplicate/o RawForceShiftbyWLC1 root:CorrectRupture:WLC:CR_RawShiftWLC
		duplicate/o CutForceWave1 root:CorrectRupture:WLC:CR_CutForceWave1
		duplicate/o CutForceWave2 root:CorrectRupture:WLC:CR_CutForceWave2
		make/o/n=(1,2) root:CorrectRupture:WLC:Final
		wave w1= root:CorrectRupture:WLC:Final
		w1={Estimate,ForceWave(Estimate),(Estimate-dimoffset(ForceWave,0))/dimdelta(ForceWave,0),currenttime,2*ForceWave(Estimate)}
	endif
	return Result

end


Static Function/C CalculateCrossingLinear(ForceWave,SepWave,WLCParms,FilterPoints,EstimatePoint,PrevPoint,NextPoint,TransType,Waveout,[CopyWavesOUt])
	wave ForceWave,SepWave,WLCParms,Waveout

	variable EstimatePoint,TransType,PrevPoint,NextPoint,FilterPoints

	variable CopyWavesOUt

	variable startingscan,endingscan
	variable/C Result

	variable edge //Depending on whether it's folding or unfolding we pick a direction to find levels
	if(TransType==1)
		edge=2
	elseif(transtype==-1)
		edge=1
	else
		return 0
	endif
	
	variable currenttime=pnt2x(ForceWave,EstimatePoint)
	variable nexttime=pnt2x(Forcewave,NextPoint)
	variable Prevtime=pnt2x(Forcewave,PrevPoint)

	variable SmoothStart=currenttime-.3
	variable SmoothEnd=currenttime+.3
	//		
	duplicate/free/r=(SmoothStart,SmoothEnd) ForceWave CutForceWave
	duplicate/free/r=(SmoothStart,SmoothEnd) SepWave CutSepWave
	//
	make/free/n=0 CutForceSm,CutsepSm,CutForceSuperSm,CutsepSuperSm
	DE_Filtering#FilterForceSep(CutForceWave,CutSepWave,CutForceSm,CutsepSm,"SVG",FilterPoints)
	DE_Filtering#FilterForceSep(CutForceWave,CutSepWave,CutForceSuperSm,CutsepSuperSm,"TVD",3e-9)

	variable startline1fit=currenttime-min(0.003,currenttime-prevtime)
	variable endline1fit=currenttime//-.00001
	variable startline1out=currenttime-2e-4//-.0005
	variable endline1out=currenttime+min(0.001,(nexttime-currenttime)*.3)

	duplicate/free/r=(startline1fit,endline1fit) CutForceSuperSm Line1FitRegion
	duplicate/free/r=(startline1out,endline1out) CutForceSm Line1

	variable startline2fit=currenttime+.0002
	variable endline2fit=currenttime+min(0.003,(nexttime-currenttime)*.8)
	variable startline2out=currenttime-2e-4//.0005
	variable endline2out=currenttime+min(1e-3,(nexttime-currenttime)*.6) //In the second WLC look ahead EITHER 2 ms, or 60% to the next transition
	duplicate/free/r=(startline2fit,endline2fit) CutForceSuperSm Line2FitRegion
	duplicate/free/r=(startline2out,endline2out) CutForceSm Line2
	
	variable WLC1Slope=ReturnWLCSlope(ForceWave,SepWave,WLCParms,EstimatePoint,NextPoint,TransType,0)
	variable WLC2Slope=ReturnWLCSlope(ForceWave,SepWave,WLCParms,EstimatePoint,NextPoint,TransType,1)


	variable slope1=MakeLinearFitToWaveLimited(WLC1Slope,Line1FitRegion,Line1,Line1)
	variable slope2=MakeLinearFitToWaveLimited(WLC2Slope,Line2FitRegion,Line2,Line2)

	duplicate/free/R=(startline2out,endline2out) CutForceSm ShiftForcebyLine2
	ShiftForcebyLine2-=Line2
	make/free/n=0 Line2Crossings
	if(FindLevelsWithError(ShiftForcebyLine2,0,.00001,edge,Line2Crossings)==-1)
		//	print "Line2 Error"
		FilterPoints=2*(round(FilterPoints/3))+1
		if(FilterPoints>5)
			Result=CalculateCrossingLinear(ForceWave,SepWave,WLCParms,FilterPoints,EstimatePoint,PrevPoint,NextPoint,TransType,Waveout,CopyWavesOUt=1)

			return Result
		else
			print "Filtering under 5"
			if(ParamisDefault(CopyWavesOUt)||CopyWavesOUt==0)
	
	else
		NewDataFolder/O root:CorrectRupture:Linear
		duplicate/o Line1 root:CorrectRupture:Linear:CR_Line1
		duplicate/o Line2 root:CorrectRupture:Linear:CR_Line2
		duplicate/o CutForceSm root:CorrectRupture:Linear:CR_CutForceSm
		duplicate/o Line1FitRegion root:CorrectRupture:Linear:CR_Line1FiTregion
		duplicate/o Line2FitRegion root:CorrectRupture:Linear:CR_Line2FitRegion
		duplicate/o CutForceWave root:CorrectRupture:Linear:CR_CutForceWave
		make/o/n=(1,2) root:CorrectRupture:Linear:Final
		wave w1= root:CorrectRupture:Linear:Final
		w1={NaN,NaN,NaN}
	endif
			return -1
		endif
	endif

	duplicate/free Line2Crossings,Line2Force

	if(numpnts(Line2Crossings)==0)
	else
		Line2Force =CutForceSm(Line2Crossings)

	endif

	ShiftForcebyLine2+=Line2

	duplicate/free/R=(startline1out,endline1out) CutForceSm ShiftForcebyLine1
	ShiftForcebyLine1-=Line1
	make/free/n=0 Line1Crossings
	if(FindLevelsWithError(ShiftForcebyLine1,0,.00001,edge,Line1Crossings)==-1)
		//print "Line1 Error"
		FilterPoints=2*(round(FilterPoints/3))+1
		if(FilterPoints>5)
			Result=CalculateCrossingLinear(ForceWave,SepWave,WLCParms,FilterPoints,EstimatePoint,PrevPoint,NextPoint,TransType,Waveout,CopyWavesOUt=1)

			return Result
		else
			print "Filtering under 5"
			if(ParamisDefault(CopyWavesOUt)||CopyWavesOUt==0)
	
	else
		NewDataFolder/O root:CorrectRupture:Linear
		duplicate/o Line1 root:CorrectRupture:Linear:CR_Line1
		duplicate/o Line2 root:CorrectRupture:Linear:CR_Line2
		Duplicate/o Line2Crossings root:CorrectRupture:Linear:CR_Line2Crossings
		Duplicate/o Line2Force root:CorrectRupture:Linear:CR_Line2Force
		duplicate/o CutForceSm root:CorrectRupture:Linear:CR_CutForceSm
		duplicate/o Line1FitRegion root:CorrectRupture:Linear:CR_Line1FiTregion
		duplicate/o Line2FitRegion root:CorrectRupture:Linear:CR_Line2FitRegion
		duplicate/o CutForceWave root:CorrectRupture:Linear:CR_CutForceWave
		make/o/n=(1,2) root:CorrectRupture:Linear:Final
		wave w1= root:CorrectRupture:Linear:Final
		w1={NaN,NaN,NaN}
	endif
			return -1
		endif
	endif

	duplicate/free Line1Crossings,Line1Force
	if(numpnts(Line1Crossings)==0)
	else
		Line1Force =CutForceSm(Line1Crossings)

	endif
	
	EliminatePreCrossings(Line1Crossings,Line2Crossings)//Eliminate crossings in 2 that occur before
																		// any of the crossing in 1
																		
																		
																		
	if(numpnts(Line1Crossings)==0||numpnts(Line2Crossings)==0)
		FilterPoints=2*(round(FilterPoints/3))+1
		if(FilterPoints>5)
			Result=CalculateCrossingLinear(ForceWave,SepWave,WLCParms,FilterPoints,EstimatePoint,PrevPoint,NextPoint,TransType,Waveout,CopyWavesOUt=1)

			return Result
		else
		print "Filtering under 5"
		if(ParamisDefault(CopyWavesOUt)||CopyWavesOUt==0)
	
	else
		NewDataFolder/O root:CorrectRupture:Linear
		duplicate/o Line1 root:CorrectRupture:Linear:CR_Line1
		duplicate/o Line2 root:CorrectRupture:Linear:CR_Line2
		Duplicate/o Line1Crossings root:CorrectRupture:Linear:CR_Line1Crossings
		Duplicate/o Line2Crossings root:CorrectRupture:Linear:CR_Line2Crossings
		Duplicate/o Line1Force root:CorrectRupture:Linear:CR_Line1Force
		Duplicate/o Line2Force root:CorrectRupture:Linear:CR_Line2Force
		duplicate/o CutForceSm root:CorrectRupture:Linear:CR_CutForceSm
		duplicate/o Line1FitRegion root:CorrectRupture:Linear:CR_Line1FiTregion
		duplicate/o Line2FitRegion root:CorrectRupture:Linear:CR_Line2FitRegion
		duplicate/o CutForceWave root:CorrectRupture:Linear:CR_CutForceWave
		make/o/n=(1,2) root:CorrectRupture:Linear:Final
		wave w1= root:CorrectRupture:Linear:Final
		w1={NaN,NaN,NaN}
	endif
		return -1
		endif
	
	endif
																		
	if(TransType==1)
	
		startingscan=LastPoint(Line1Crossings,Line2Crossings)-1e-4
		endingscan=Line2Crossings[0]+2e-6
	elseif(TransType==-1)
		startingscan=LastPoint(Line1Crossings,Line2Crossings)-1e-4
		endingscan=Line2Crossings[0]+2e-6
	endif

	if(endingscan>=endline1out)
		endingscan=endline1out
	endif
	make/free/n=0 FinalLine1
	duplicate/free/R=(startingscan,endingscan)  CutForceWave RawShiftLine
	//duplicate/o/R=(startingscan,endingscan)  CutForceSm RawShiftLine
	//
	MakeLinearFitToWaveLimited(WLC1Slope,Line1FitRegion,CutForceWave,FinalLine1)
	if(TransType==1)
		duplicate/free/R=(startingscan,endingscan)  FinalLine1 RawLine

	elseif(TransType==-1)
		duplicate/free/R=(startingscan,endingscan)  FinalLine1 RawLine

	endif

	RawShiftLine-=RawLine

	make/free/n=0 RawlineCrossings
	if(FindLevelsWithError(RawShiftLine,0,.00001,edge,RawlineCrossings)==-1)
		if(TransType==1)
	
			startingscan=LastPoint(Line1Crossings,Line2Crossings)-3e-4
			endingscan=Line2Crossings[0]+2e-6
		elseif(TransType==-1)
			startingscan=LastPoint(Line1Crossings,Line2Crossings)-3e-4
			endingscan=Line2Crossings[0]+2e-6
		endif
	
		if(endingscan>=endline1out)
			endingscan=endline1out
		endif
		make/free/n=0 FinalLine1
		duplicate/free/R=(startingscan,endingscan)  CutForceWave RawShiftLine
		//duplicate/o/R=(startingscan,endingscan)  CutForceSm RawShiftLine
		//
		MakeLinearFitToWaveLimited(WLC1Slope,Line1FitRegion,CutForceWave,FinalLine1)
		if(TransType==1)
			duplicate/free/R=(startingscan,endingscan)  FinalLine1 RawLine

		elseif(TransType==-1)
			duplicate/free/R=(startingscan,endingscan)  FinalLine1 RawLine

		endif

		RawShiftLine-=RawLine
		

		make/free/n=0 RawlineCrossings
		if(FindLevelsWithError(RawShiftLine,0,.00001,edge,RawlineCrossings)==-1)
			print "Final Line Failed"
			if(ParamisDefault(CopyWavesOUt)||CopyWavesOUt==0)
	
	else
		NewDataFolder/O root:CorrectRupture:Linear
		duplicate/o Line1 root:CorrectRupture:Linear:CR_Line1
		duplicate/o Line2 root:CorrectRupture:Linear:CR_Line2
		Duplicate/o Line1Crossings root:CorrectRupture:Linear:CR_Line1Crossings
		Duplicate/o Line2Crossings root:CorrectRupture:Linear:CR_Line2Crossings
		Duplicate/o Line1Force root:CorrectRupture:Linear:CR_Line1Force
		Duplicate/o Line2Force root:CorrectRupture:Linear:CR_Line2Force
		Duplicate/o RawShiftLine root:CorrectRupture:Linear:CR_RawShiftLine
		duplicate/o CutForceSm root:CorrectRupture:Linear:CR_CutForceSm
		duplicate/o Line1FitRegion root:CorrectRupture:Linear:CR_Line1FiTregion
		duplicate/o Line2FitRegion root:CorrectRupture:Linear:CR_Line2FitRegion
		duplicate/o CutForceWave root:CorrectRupture:Linear:CR_CutForceWave
		make/o/n=(1,2) root:CorrectRupture:Linear:Final
		wave w1= root:CorrectRupture:Linear:Final
		w1={NaN,NaN,NaN}
	endif
			return nAN
		endif
	endif
	RawShiftLine+=RawLine







	DUPLICATE/free RawlineCrossings RawlineForces
	variable Estimate
	if(numpnts(RawlineCrossings)==0)
		Estimate=Nan
		Result=cmplx(Estimate,NAN)
	else
		RawlineForces=ForceWave(RawlineCrossings)
		Estimate= RawLineCrossings[numpnts(RawLineCrossings)-1]
		Result=cmplx(Estimate,ForceWave(Estimate))
	endif
	make/D/free/n=4 ResultsOut
	ResultsOut={Estimate,ForceWave(Estimate),(Estimate-dimoffset(ForceWave,0))/dimdelta(ForceWave,0)}
	duplicate/o ResultsOut Waveout
	if(ParamisDefault(CopyWavesOUt)||CopyWavesOUt==0)
	
	else
			NewDataFolder/O root:CorrectRupture

		NewDataFolder/O root:CorrectRupture:Linear
		duplicate/o Line1 root:CorrectRupture:Linear:CR_Line1
		duplicate/o Line2 root:CorrectRupture:Linear:CR_Line2
		Duplicate/o Line1Crossings root:CorrectRupture:Linear:CR_Line1Crossings
		Duplicate/o Line2Crossings root:CorrectRupture:Linear:CR_Line2Crossings
		Duplicate/o Line1Force root:CorrectRupture:Linear:CR_Line1Force
		Duplicate/o Line2Force root:CorrectRupture:Linear:CR_Line2Force
		Duplicate/o RawShiftLine root:CorrectRupture:Linear:CR_RawShiftLine
		duplicate/o CutForceSm root:CorrectRupture:Linear:CR_CutForceSm
		duplicate/o Line1FitRegion root:CorrectRupture:Linear:CR_Line1FiTregion
		duplicate/o Line2FitRegion root:CorrectRupture:Linear:CR_Line2FitRegion
		duplicate/o CutForceWave root:CorrectRupture:Linear:CR_CutForceWave
		make/o/n=(1,2) root:CorrectRupture:Linear:Final
		wave w1= root:CorrectRupture:Linear:Final
		w1={Estimate,ForceWave(Estimate),(Estimate-dimoffset(ForceWave,0))/dimdelta(ForceWave,0),currenttime,2*ForceWave(Estimate)}
	endif

	return Result

end

Static Function ReturnWLCSlope(ForceWave,SepWave,WLCParms,EstimatePoint,NextPoint,TransType,transNumber)
	wave ForceWave,SepWave,WLCParms
	variable EstimatePoint,NextPoint,TransType,transNumber
	
	variable currenttime=pnt2x(ForceWave,EstimatePoint)

	variable nexttime=pnt2x(Forcewave,NextPoint)

	variable backtimeFirst=2e-4//min(.3e-3,(pnt2x(ForceWave,EstimatePoint)-Prevtime)*.9)
	variable forwardtimeFirst=min(2e-3,(nexttime-currenttime)*.9)
	variable backtimeSecond=-2e-4
	variable forwardtimeSecond=min(4e-3,(nexttime-currenttime)*.9)

	variable starttimeFirst=pnt2x(ForceWave,EstimatePoint)-backtimeFirst
	variable endtimeFirst=pnt2x(ForceWave,EstimatePoint)+forwardtimeFirst
	variable starttimeSecond=pnt2x(ForceWave,EstimatePoint)-backtimeSecond
	variable endtimeSecond=pnt2x(ForceWave,EstimatePoint)+forwardtimeSecond
	duplicate/free/r=(starttimeFirst,endtimeFirst) ForceWave WLC1
	duplicate/free/r=(starttimeSecond,endtimeSecond) ForceWave WLC2
	duplicate/free/r=(starttimeFirst,endtimeFirst) ForceWave CutForceWave1
	duplicate/free/r=(starttimeSecond,endtimeSecond) ForceWave CutForceWave2
	duplicate/free/r=(starttimeFirst,endtimeFirst) SepWave CutSepWave1
	duplicate/free/r=(starttimeSecond,endtimeSecond) SepWave CutSepWave2

	if(TransType==1)
		WLC1=WLC(CutSepWave1-WLCParms[3],.4e-9,WLCParms[1],298)-WLCParms[2]
		WLC2=WLC(CutSepWave2-WLCParms[3],.4e-9,WLCParms[0],298)-WLCParms[2]
	elseif(TransType==-1)
		WLC1=WLC(CutSepWave1-WLCParms[3],.4e-9,WLCParms[0],298)-WLCParms[2]

		WLC2=WLC(CutSepWave2-WLCParms[3],.4e-9,WLCParms[1],298)-WLCParms[2]
	endif

	make/free/n=0 Garbage
	variable WLC1Slope=MakeLinearFitToWave(WLC1,Garbage,Garbage)
	variable WLC2Slope=MakeLinearFitToWave(WLC2,Garbage,Garbage)
	variable result
	if(transNumber==0)
		result=WLC1Slope
	elseif(transNumber==1)
		result=WLC2Slope

	endif
	return result
end



Static Function FindLevelsWithError(WaveIn,Level,M,edge,Waveout)

	wave WaveIn,Waveout
	variable Level,M,edge
	FindLevels/Q/M=(M)/edge=(edge) WaveIn Level
	wave W_FindLevels
	if(numpnts(W_FindLevels)==0)
		killwaves W_FindLevels 

	return -1
	endif
	duplicate/o W_FindLevels WaveOut
	killwaves W_FindLevels 
	return 0

end
Static Function MakeLinearFitToWave(WaveIn,Overlap,Waveout)
	wave wavein,waveout,Overlap
	duplicate/free Overlap Resulting
	CurveFit/Q/W=2/NTHR=0 line Wavein
	
	wave w_coef,w_sigma,M_Jacobian

	Resulting=x*w_coef[1]+w_coef[0]
	variable result=w_coef[1]
	killwaves w_coef ,w_sigma,M_Jacobian
	duplicate/o Resulting WaveOut
	return result
end

Static Function MakeLinearFitToWaveFixed(SlopeIn,WaveIn,Overlap,Waveout)
	variable SlopeIn
	wave wavein,waveout,Overlap
	duplicate/free Overlap Resulting
	
	K1=SlopeIn
	CurveFit/H="01"/Q/W=2/NTHR=0 line Wavein
	
	
	wave w_coef,w_sigma,M_Jacobian

	Resulting=x*w_coef[1]+w_coef[0]
	variable result=w_coef[1]
	killwaves w_coef ,w_sigma,M_Jacobian
	duplicate/o Resulting WaveOut
	return result
end

Static Function MakeLinearFitToWaveLimited(SlopeIn,WaveIn,Overlap,Waveout)
	variable SlopeIn
	wave wavein,waveout,Overlap
	duplicate/free Overlap Resulting
	K1=SlopeIn
	CurveFit/H="00"/Q/W=2/NTHR=0 line Wavein
	wave w_coef
	Make/free/T/N=2 T_Constraints
	if(Slopein<0)
		T_Constraints[0] = {"K1>"+num2str(2*SlopeIn),"K1<"+num2str(-2*SlopeIn)}

	else
	T_Constraints[0] = {"K1<"+num2str(2*SlopeIn),"K1>"+num2str(-2*SlopeIn)}
	endif
	FuncFit/ODR=0/Q/W=2/NTHR=0 DE_linefit W_coef Wavein/C=T_Constraints
	if(abs(w_coef[1]-SlopeIn)<1e-12)
	FuncFit/ODR=2/Q/W=2/NTHR=0 DE_linefit W_coef Wavein/C=T_Constraints

	endif
	wave w_coef,w_sigma,M_Jacobian
	Resulting=x*w_coef[1]+w_coef[0]
	variable result=w_coef[1]
	killwaves w_coef ,w_sigma,M_Jacobian
	duplicate/o Resulting WaveOut
	return result
end

Function DE_linefit(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) =a+b*x
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 2
	//CurveFitDialog/ w[0] = a
	//CurveFitDialog/ w[1] = b

	return w[0]+w[1]*x
End

//Returns the LAST point in Crossings one that is before ANY of the crossings in 2
Static Function LastPoint(Crossings1,Crossings2)
	wave Crossings1, Crossings2
	
	variable n=0
	for(n=0;n<(numpnts(Crossings1)-1);n+=1)
		if(Crossings2[0]<Crossings1[n+1])
	
		return Crossings1[n]
		endif			
	endfor
	return Crossings1[numpnts(Crossings1)-1]
	


end

Static Function LastPointMod(Crossings1,Crossings2)
	wave Crossings1, Crossings2
	
	variable n=0
	for(n=0;n<(numpnts(Crossings1)-1);n+=1)
		if(Crossings2[numpnts(Crossings2)-1]<Crossings1[n])
	
		return Crossings1[n]
		endif			
	endfor
	return Crossings1[numpnts(Crossings1)-1]
	


end

Static Function NextPointinWave(Crossing,TimeIn)
	wave Crossing
	variable TimeIn
	variable n=0
	for(n=0;n<(numpnts(Crossing)-1);n+=1)
		if(Crossing[n]>TimeIn)
	
		return Crossing[n]
		endif			
	endfor
	return Crossing[numpnts(Crossing)-1]
	


end


//Eliminates any points in Crossings2 that occur BEFORE the first point in Crossings1
Static Function EliminatePreCrossings(Crossings1,Crossings2)
	wave Crossings1, Crossings2
	variable cut
	variable n=0
	for(n=0;n<(numpnts(Crossings2));n+=1)
		if(Crossings1[0]>Crossings2[n])
		cut=n+1
		endif			
	endfor
	deletepoints 0,cut, Crossings2
	


end

Static Function FixPicksSingle(PntWave,index,FilterPnts,Type)
	wave PntWave
	variable index //Down=-1, Up=1
	variable FilterPnts
	string type 
	variable prepnt,Direction,nextpnt
	//pH7M10006RupPntu_adj,1,0,pH7M10006Force_Align
	if(strsearch(NameofWave(PntWave),"PntD",0,2)!=-1)
		Direction=-1
	elseif(strsearch(NameofWave(PntWave),"PntU",0,2)!=-1)
		Direction=1
	else
		return 0
	endif
	
	
	if(Direction==-1)
		wave OtherPnt=$replacestring("PntD",nameofwave(pntWave),"PntU")
		wave ForceWave=$replacestring("RupPntD_adj",nameofwave(pntWave),"Force_Align")

		prepnt=OtherPnt[index]
		if(index==(numpnts(OtherPnt)-1))
			nextpnt=numpnts(ForceWave)
		else
			nextpnt=OtherPnt[index+1]

	endif

	elseif(Direction==1)
		wave OtherPnt=$replacestring("PntU",nameofwave(pntWave),"PntD")
		wave ForceWave=$replacestring("RupPntu_adj",nameofwave(pntWave),"Force_Align")

			if(index==0)
			prepnt=0
			nextpnt=OtherPnt[index]
		else
			prepnt=OtherPnt[index-1]
			nextpnt=OtherPnt[index]
		endif
		
	endif

	wave SepWave=$replacestring("Force_Align",nameofwave(ForceWave),"Sep_Align_Sm")
	wave WLC=$replacestring("Force_Align",nameofwave(ForceWave),"_WLCParms_Adj")
	duplicate/free ForceWave FWIN
	FWIN*=-1
	make/o/n=0 Test
	
	StrSwitch(Type)
	
	case "WLC":
		DE_CorrRup#CalculateCrossingfromWLC(FWIN,SepWave,WLC,FilterPnts,PntWave[index],prepnt,nextpnt,Direction,Test,CopyWavesOUt=1)
		MakeASinglePlot(ForceWave,"WLC")
	
	break

	case "Line":
		DE_CorrRup#CalculateCrossingLinear(FWIN,SepWave,WLC,FilterPnts,PntWave[index],prepnt,nextpnt,Direction,Test,CopyWavesOUt=1)
		MakeASinglePlot(ForceWave,"Line")

	break
	
	default:
	break
	endswitch

	//print PntWave[index];print prepnt;print nextpnt
	//DE_CorrRup#CalculateBestCrossingGuess(FWIN,SepWave,WLC,PntWave[index],prepnt,nextpnt,Direction,Test,CopyWavesOUt=1)

	//MakeASinglePlot(ForceWave)


	
end

Static Function MakeASinglePlot(ForceWave,Type)
	wave ForceWave
	string Type
		wave FSm=$(nameofwave(ForceWave)+"_Sm")
	DoWindow SinglePoint
		if(V_flag==1)
	killwindow SinglePoint
	endif
	display/N=SinglePoint ForceWave;appendtograph/W=SinglePoint FSM 
	ModifyGraph/W=SinglePoint rgb($nameofwave(ForceWave))=(65535,49151,49151)
	ModifyGraph/W=SinglePoint muloffset={0,-1}
	DFREF saveDFR = GetDataFolderDFR()


	StrSwitch(Type)
	
	case "WLC":
		SetDataFolder root:CorrectRupture:WLC
		AppendToGraph/W=SinglePoint $"CR_WLC1",$"CR_WLC2"
		ModifyGraph/W=SinglePoint lsize($"CR_WLC1")=2,rgb($"CR_WLC1")=(0,0,0),lsize($"CR_WLC2")=2,rgb($"CR_WLC2")=(0,0,0)
		AppendToGraph/W=SinglePoint $"CR_WLC1Forces" vs $"CR_WLC1Crossings"
		AppendToGraph/W=SinglePoint $"CR_WLC2Forces" vs $"CR_WLC2Crossings"
		ModifyGraph/W=SinglePoint mode($"CR_WLC1Forces")=3,marker($"CR_WLC1Forces")=19,rgb($"CR_WLC1Forces")=(19789,44975,19018),useMrkStrokeRGB($"CR_WLC1Forces")=1
		ModifyGraph/W=SinglePoint mode($"CR_WLC1Forces")=3,mode($"CR_WLC2Forces")=3,marker($"CR_WLC2Forces")=16,rgb($"CR_WLC2Forces")=(14906,32382,47288),useMrkStrokeRGB($"CR_WLC2Forces")=1
		AppendToGraph/W=SinglePoint $"CR_RawShiftWLC"
		ModifyGraph/W=SinglePoint lsize($"CR_RawShiftWLC")=1,rgb($"CR_RawShiftWLC")=(65535,30840,0)

	

	break
	
	case "Line":
		SetDataFolder root:CorrectRupture:Linear
		AppendToGraph/W=SinglePoint $"CR_Line1",$"CR_Line2"
		ModifyGraph/W=SinglePoint lsize($"CR_Line1")=2,rgb($"CR_Line1")=(23130,23130,23130),lsize($"CR_Line2")=2,rgb($"CR_Line2")=(23130,23130,23130)
		AppendToGraph/W=SinglePoint $"CR_Line1Force" vs $"CR_Line1Crossings"
		AppendToGraph/W=SinglePoint $"CR_Line2Force" vs $"CR_Line2Crossings"
		ModifyGraph/W=SinglePoint mode($"CR_Line1Force")=3,marker($"CR_Line1Force")=19,rgb($"CR_Line1Force")=(19789,44975,19018),useMrkStrokeRGB($"CR_Line1Force")=1
		ModifyGraph/W=SinglePoint mode($"CR_Line2Force")=3,mode($"CR_Line2Force")=3,marker($"CR_Line2Force")=16,rgb($"CR_Line2Force")=(14906,32382,47288),useMrkStrokeRGB($"CR_Line2Force")=1
		AppendToGraph/W=SinglePoint $"CR_RawShiftLine"
		ModifyGraph/W=SinglePoint lsize($"CR_RawShiftLine")=1,rgb($"CR_RawShiftLine")=(65535,30840,0)

	break
	
	default:
	break
	endswitch
		wave T1=$"Final"
		AppendToGraph/W=SinglePoint T1[1,1] vs T1[0,0]
		ModifyGraph/W=SinglePoint mode($nameofwave(T1))=3,marker($nameofwave(T1))=18,rgb($nameofwave(T1))=(36873,14755,58982),useMrkStrokeRGB($nameofwave(T1))=1
		//AppendToGraph/W=SinglePoint T1[3,3] vs T1[1,1]
		//ModifyGraph/W=SinglePoint mode($nameofwave(T1)+"#1")=3,marker($nameofwave(T1)+"#1")=26,rgb($nameofwave(T1)+"#1")=(65535,0,52428),useMrkStrokeRGB($nameofwave(T1)+"#1")=1
		SetAxis/W=SinglePoint bottom T1[0]-2e-3,T1[0]+3e-3
		SetAxis/W=SinglePoint/A=2 left
		AppendToGraph/W=SinglePoint T1[4,4] vs T1[3,3]
		ModifyGraph/W=SinglePoint mode($(nameofwave(T1)+"#1"))=1,rgb($(nameofwave(T1)+"#1"))=(26205,52428,1)


//	

//	


	SetDataFolder saveDFR

end


Static Function/C CalculateBestCrossingGuess(ForceWave,SepWave,WLCParms,EstimatePoint,PrevPoint,NextPoint,TransType,Waveout,[CopyWavesOUt])
	wave ForceWave,SepWave,WLCParms,Waveout

	variable EstimatePoint,TransType,PrevPoint,NextPoint

	variable CopyWavesOUt

	variable nexttime=pnt2x(Forcewave,NextPoint)
	variable Prevtime=pnt2x(Forcewave,PrevPoint)

	
	variable backtimeFirst=2e-4//min(.3e-3,(pnt2x(ForceWave,EstimatePoint)-Prevtime)*.9)
	variable forwardtimeFirst=min(2e-3,(nexttime-pnt2x(ForceWave,EstimatePoint))*.9)
	variable backtimeSecond=-2e-4
	variable forwardtimeSecond=min(4e-3,(nexttime-pnt2x(ForceWave,EstimatePoint))*.9)
	
	variable FilterPoints=31
		variable edge
	
	if(TransType==1)
		edge=2
	elseif(transtype==-1)
		edge=1
	endif

	
	variable starttimeFirst=pnt2x(ForceWave,EstimatePoint)-backtimeFirst
	variable endtimeFirst=pnt2x(ForceWave,EstimatePoint)+forwardtimeFirst
	variable starttimeSecond=pnt2x(ForceWave,EstimatePoint)-backtimeSecond
	variable endtimeSecond=pnt2x(ForceWave,EstimatePoint)+forwardtimeSecond

	duplicate/free/r=(starttimeFirst,endtimeFirst) ForceWave WLC1
	duplicate/free/r=(starttimeSecond,endtimeSecond) ForceWave WLC2
	duplicate/free/r=(starttimeFirst,endtimeFirst) ForceWave CutForceWave1
	duplicate/free/r=(starttimeSecond,endtimeSecond) ForceWave CutForceWave2

	duplicate/free/r=(starttimeFirst,endtimeFirst) SepWave CutSepWave1
	duplicate/free/r=(starttimeSecond,endtimeSecond) SepWave CutSepWave2
	duplicate/free/r=(starttimeFirst-1e-2,endtimeFirst+1e-2) ForceWave WideWLC1
	duplicate/free/r=(starttimeFirst-1e-2,endtimeFirst+1e-2) SepWave CutSepWaveWide


	if(TransType==1)
	WLC1=WLC(CutSepWave1-WLCParms[3],.4e-9,WLCParms[1],298)-WLCParms[2]
		WideWLC1=WLC(CutSepWaveWide-WLCParms[3],.4e-9,WLCParms[1],298)-WLCParms[2]

	WLC2=WLC(CutSepWave2-WLCParms[3],.4e-9,WLCParms[0],298)-WLCParms[2]
	
	elseif(TransType==-1)
	WLC1=WLC(CutSepWave1-WLCParms[3],.4e-9,WLCParms[0],298)-WLCParms[2]
		WideWLC1=WLC(CutSepWaveWide-WLCParms[3],.4e-9,WLCParms[0],298)-WLCParms[2]

	WLC2=WLC(CutSepWave2-WLCParms[3],.4e-9,WLCParms[1],298)-WLCParms[2]
	endif
	

	make/free/n=0 Garbage
	variable WLC1Slope=MakeLinearFitToWave(WLC1,Garbage,Garbage)
	variable WLC2Slope=MakeLinearFitToWave(WLC2,Garbage,Garbage)
	make/free/n=0 CutForceSm1,CutsepSm1,CutForceSm2,CutsepSm2
	DE_Filtering#FilterForceSep(CutForceWave1,CutSepWave1,CutForceSm1,CutsepSm1,"SVG",FilterPoints)
	DE_Filtering#FilterForceSep(CutForceWave2,CutSepWave2,CutForceSm2,CutsepSm2,"SVG",FilterPoints)

	//DE_Filtering#FilterForceSep(CutForceWave,CutSepWave,CutForceSm,CutsepSm,"TVD",1e-10)
	
	duplicate/free CutForceSm2 ShiftForcebyWLC2
	ShiftForcebyWLC2-=wLC2
	make/free/n=0 WLC2Crossings
	duplicate/free CutForceSm1 ShiftForcebyWLC1
	ShiftForcebyWLC1-=wLC1
	make/free/n=0 WLC1Crossings
	
	if(FindLevelsWithError(ShiftForcebyWLC2,0,.00001,edge,WLC2Crossings)==-1)

		//return nAN
	endif
	duplicate/free WLC2Crossings,WLC2Forces
	WLC2Forces =CutForceSm2(WLC2Crossings)

	if(FindLevelsWithError(ShiftForcebyWLC1,0,.00001,edge,WLC1Crossings)==-1)

		//return nAN
	endif

	duplicate/Free WLC1Crossings,WLC1Forces
	WLC1Forces =CutForceSm1(WLC1Crossings)
	
	variable startingscan
	variable endingscan
	
	
	if(TransType==1)
		EliminatePreCrossings(WLC1Crossings,WLC2Crossings)
		startingscan=LastPoint(WLC1Crossings,WLC2Crossings)-3e-4
		//endingscan=WLC2Crossings[0]-2e-5
				endingscan=NextPointinWave(WLC2Crossings,startingscan)+1e-4

		//endingscan=LastPointMod(WLC2Crossings,WLC1Crossings)+1e-4

	elseif(TransType==-1)
		EliminatePreCrossings(WLC1Crossings,WLC2Crossings)

		startingscan=LastPoint(WLC1Crossings,WLC2Crossings)-3e-4
		//endingscan=WLC2Crossings[0]-2e-5
		endingscan=NextPointinWave(WLC2Crossings,startingscan)+1e-4
		//endingscan=LastPointMod(WLC2Crossings,WLC1Crossings)+1e-4
	endif
	
	
	duplicate/free/R=(startingscan,endingscan)  ForceWave RawShiftWLC
	//duplicate/free/R=(startingscan,endingscan)  CutForceSm RawShiftWLC
	
//	if(endingscan>forwardtimeFirst)
//		duplicate/free/r=(starttimeFirst,endingscan) ForceWave WLC1
//
//
//		duplicate/o/r=(starttimeFirst,endingscan) SepWave CutSepWave1
//		duplicate/o/r=(starttimeSecond,endingscan) SepWave CutSepWave2
//
//		//duplicate/free/r=(Starttime,EndTime) SepWave CutSepWave
//		//duplicate/free/r=(Starttime-.05,EndTime+.05) SepWave CutSepWaveLong
//
//		if(TransType==1)
//			WLC1=WLC(CutSepWave1-WLCParms[3],.4e-9,WLCParms[1],298)-WLCParms[2]
//	
//		elseif(TransType==-1)
//			WLC1=WLC(CutSepWave1-WLCParms[3],.4e-9,WLCParms[0],298)-WLCParms[2]
//		endif
//	endif
	
	if(TransType==1)
		duplicate/free/R=(startingscan,endingscan)  WideWLC1 RawWLC


	elseif(TransType==-1)
		duplicate/free/R=(startingscan,endingscan)  WideWLC1 RawWLC

	endif
	RawShiftWLC-=RawWLC



	
	make/free/n=0 RawWLCCrossings
	if(FindLevelsWithError(RawShiftWLC,0,.00001,edge,RawWLCCrossings)==-1)
		//	return nAN
	endif
	//RawShiftWLC+=RawWLC
	DUPLICATE/free RawWLCCrossings RawWLCForces
	RawWLCForces=ForceWave(RawWLCCrossings)
	variable Estimate1= RawWLCCrossings[numpnts(RawWLCCrossings)-1]
	





	//Switch to the linear fit!

	variable SmoothStart=Estimate1-.3
	variable SmoothEnd=Estimate1+.3
		
	duplicate/free/r=(SmoothStart,SmoothEnd) ForceWave CutForceWave
	duplicate/free/r=(SmoothStart,SmoothEnd) SepWave CutSepWave

	make/free/n=0 CutForceSm,CutsepSm,CutForceSuperSm,CutsepSuperSm
	DE_Filtering#FilterForceSep(CutForceWave,CutSepWave,CutForceSm,CutsepSm,"SVG",5)
	DE_Filtering#FilterForceSep(CutForceWave,CutSepWave,CutForceSuperSm,CutsepSuperSm,"TVD",3e-9)



	variable startline1fit=Estimate1-min(0.003,Estimate1-prevtime)
	variable endline1fit=Estimate1-.00001
	variable startline1out=Estimate1-.0005
	variable endline1out=Estimate1+.0005
	
	
	duplicate/free/r=(startline1fit,endline1fit) CutForceSuperSm Line1FitRegion
	//duplicate/o/r=(startline1fit,endline1fit) CutForceSm Line1FitRegion
	duplicate/free/r=(startline1out,endline1out) CutForceSm Line1
	
	variable startline2fit=Estimate1+.0002
	variable endline2fit=Estimate1+min(0.003,(nexttime-Estimate1)*.8)
	variable startline2out=Estimate1-0//.0005
	variable endline2out=Estimate1+.0005

	//duplicate/o/r=(startline2fit,endline2fit) CutForceSm Line2FitRegion
	duplicate/free/r=(startline2fit,endline2fit) CutForceSuperSm Line2FitRegion
	duplicate/free/r=(startline2out,endline2out) CutForceSm Line2

	make/free/n=0 Disposal

	variable slope1=MakeLinearFitToWaveLimited(WLC1Slope,Line1FitRegion,Line1,Line1)
	variable slope2=MakeLinearFitToWaveLimited(WLC2Slope,Line2FitRegion,Line2,Line2)
	
	duplicate/free/R=(startline2out,endline2out) CutForceSm ShiftForcebyLine2
	ShiftForcebyLine2-=Line2
	make/free/n=0 Line2Crossings
	if(FindLevelsWithError(ShiftForcebyLine2,0,.00001,edge,Line2Crossings)==-1)
		//return nAN
	endif
	duplicate/free Line2Crossings,Line2Force
		if(numpnts(Line2Crossings)==0)
	else
		Line2Force =CutForceSm(Line2Crossings)

	endif
		ShiftForcebyLine2+=Line2

	//	
	duplicate/free/R=(startline1out,endline1out) CutForceSm ShiftForcebyLine1
	ShiftForcebyLine1-=Line1
	make/free/n=0 Line1Crossings
	if(FindLevelsWithError(ShiftForcebyLine1,0,.00001,edge,Line1Crossings)==-1)
		//return nAN
	endif
	duplicate/free Line1Crossings,Line1Force
	if(numpnts(Line1Crossings)==0)
	else
		Line1Force =CutForceSm(Line1Crossings)

	endif
	
		
	if(TransType==1)
	
		startingscan=LastPoint(Line1Crossings,Line2Crossings)-5e-4
		endingscan=Line2Crossings[0]-2e-6
	elseif(TransType==-1)
		startingscan=LastPoint(Line1Crossings,Line2Crossings)-5e-4
		endingscan=Line2Crossings[0]-2e-6
	endif

	if(endingscan>=endline1out)
	endingscan=endline1out
	endif

	duplicate/free/R=(startingscan,endingscan)  CutForceWave RawShiftLine
	//duplicate/o/R=(startingscan,endingscan)  CutForceSm RawShiftLine

	MakeLinearFitToWaveLimited(WLC1Slope,Line1FitRegion,CutForceWave,Line1)
	if(TransType==1)
		duplicate/free/R=(startingscan,endingscan)  line1 RawLine

	elseif(TransType==-1)
		duplicate/free/R=(startingscan,endingscan)  line1 RawLine

	endif

	RawShiftLine-=RawLine

	make/free/n=0 RawlineCrossings
	if(FindLevelsWithError(RawShiftLine,0,.00001,edge,RawlineCrossings)==-1)
		
	//	return nAN
	endif
	RawShiftLine+=RawLine

	DUPLICATE/free RawlineCrossings RawlineForces
	variable Estimate2
	variable/C Result
	if(numpnts(RawlineCrossings)==0)
	Estimate2=Nan
	Result=cmplx(Estimate2,NAN)
	else
	RawlineForces=ForceWave(RawlineCrossings)
	Estimate2= RawLineCrossings[numpnts(RawLineCrossings)-1]
	Result=cmplx(Estimate2,ForceWave(Estimate2))

	endif
	make/D/free/n=4 ResultsOut
	//	w1=(Estimate1-dimoffset(ForceWave,0))/dimdelta(ForceWave,0)
	ResultsOut={Estimate1,Estimate2,ForceWave(Estimate1),ForceWave(Estimate2),(Estimate1-dimoffset(ForceWave,0))/dimdelta(ForceWave,0),(Estimate2-dimoffset(ForceWave,0))/dimdelta(ForceWave,0)}
	duplicate/o ResultsOut Waveout
	if(ParamisDefault(CopyWavesOUt)||CopyWavesOUt==0)
	
	else
	NewDataFolder/O root:CorrectRupture
	Duplicate/o WLC1 root:CorrectRupture:CR_WLC1
	Duplicate/o WLC2 root:CorrectRupture:CR_WLC2
	duplicate/o Line1 root:CorrectRupture:CR_Line1
	duplicate/o Line2 root:CorrectRupture:CR_Line2
	Duplicate/o Line1Crossings root:CorrectRupture:CR_Line1Crossings
	Duplicate/o Line2Crossings root:CorrectRupture:CR_Line2Crossings
	Duplicate/o Line1Force root:CorrectRupture:CR_Line1Force
	Duplicate/o Line2Force root:CorrectRupture:CR_Line2Force
	Duplicate/o RawShiftLine root:CorrectRupture:CR_RawShiftLine
	duplicate/o CutForceSm root:CorrectRupture:CR_CutForceSm
	duplicate/o WLC2Crossings root:CorrectRupture:CR_WLC2Crossings
	duplicate/o WLC2Forces root:CorrectRupture:CR_WLC2Forces
	duplicate/o WLC1Crossings root:CorrectRupture:CR_WLC1Crossings
	duplicate/o WLC1Forces root:CorrectRupture:CR_WLC1Forces
	duplicate/o Line1FitRegion root:CorrectRupture:CR_Line1FiTregion
	duplicate/o Line2FitRegion root:CorrectRupture:CR_Line2FitRegion
	duplicate/o ShiftForcebyWLC1 root:CorrectRupture:CR_ShiftForcebyWLC1
	duplicate/o RawShiftWLC root:CorrectRupture:CR_RawShiftWLC
	duplicate/o CutForceWave root:CorrectRupture:CR_CutForceWave

	endif


	return Result
end
