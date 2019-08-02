#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma ModuleName=DE_CorrRup
Static Function FixPicks(UpWave,ForceWave,SmoothFW,WLC,SepWave)
	wave Upwave,ForceWave,SepWave,WLC,SmoothFW
	wave DownWave=$ReplaceString("PntU", nameofwave(Upwave), "PntD")
	variable n
	make/o/n=(numpnts(DownWave),4) DownCrap
	make/o/n=(numpnts(Upwave),4) UpCrap

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
		if(n==25)
	//	print DownWave[n]
	//	print nextpnt
	//	print prepnt
		endif
		DE_CorrRup#CalculateBestCrossingGuess(ForceWave,SepWave,WLC,DownWave[n],prepnt,nextpnt,-1,Test)
		
		DownCrap[n][]=Test[q]
	endfor
	
	for(n=0;n<numpnts(Upwave);n+=1)
			if(n==0)
				prepnt=0
				nextpnt=DownWave[n]
			else
				prepnt=DownWave[n-1]
				nextpnt=DownWave[n]
			endif
		DE_CorrRup#CalculateBestCrossingGuess(ForceWave,SepWave,WLC,Upwave[n],prepnt,nextpnt,1,Test)
		UpCrap[n][]=Test[q]

	endfor
	wavestats/Q Test
	if(V_numNans>0)
	print "We missed something"
	endif
	make/o/D/n=(numpnts(Upwave)) $(nameofwave(UpWave)+"_Mod")
	wave w1=$(nameofwave(UpWave)+"_Mod")
	w1=(UpCrap[p][1]-dimoffset(ForceWave,0))/dimdelta(ForceWave,0)
	//w1=x2pnt(ForceWave,UpCrap[p][0])
	make/o/D/n=(numpnts(DownWave)) $(nameofwave(DownWave)+"_Mod")
	wave w2=$(nameofwave(DownWave)+"_Mod")
	w2=(DownCrap[p][1]-dimoffset(ForceWave,0))/dimdelta(ForceWave,0)
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
	
	variable FilterPoints=5
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

	duplicate/o/r=(starttimeFirst,endtimeFirst) SepWave CutSepWave1
	duplicate/o/r=(starttimeSecond,endtimeSecond) SepWave CutSepWave2


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
		startingscan=LastPoint(WLC1Crossings,WLC2Crossings)-1e-4
		//endingscan=WLC2Crossings[0]-2e-5
				endingscan=NextPointinWave(WLC2Crossings,startingscan)

		//endingscan=LastPointMod(WLC2Crossings,WLC1Crossings)+1e-4

	elseif(TransType==-1)
		EliminatePreCrossings(WLC1Crossings,WLC2Crossings)

		startingscan=LastPoint(WLC1Crossings,WLC2Crossings)-1e-4
		//endingscan=WLC2Crossings[0]-2e-5
		endingscan=NextPointinWave(WLC2Crossings,startingscan)
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
		duplicate/free/R=(startingscan,endingscan)  WLC1 RawWLC


	elseif(TransType==-1)
		duplicate/free/R=(startingscan,endingscan)  WLC1 RawWLC

	endif
	RawShiftWLC-=RawWLC



	
	make/free/n=0 RawWLCCrossings
	if(FindLevelsWithError(RawShiftWLC,0,.00001,edge,RawWLCCrossings)==-1)
		//	return nAN
	endif
	RawShiftWLC+=RawWLC
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
	
	
	duplicate/o/r=(startline1fit,endline1fit) CutForceSuperSm Line1FitRegion
	//duplicate/o/r=(startline1fit,endline1fit) CutForceSm Line1FitRegion
	duplicate/o/r=(startline1out,endline1out) CutForceSm Line1
	
	variable startline2fit=Estimate1+.0002
	variable endline2fit=Estimate1+min(0.003,(nexttime-Estimate1)*.8)
	variable startline2out=Estimate1-0//.0005
	variable endline2out=Estimate1+.0005

	//duplicate/o/r=(startline2fit,endline2fit) CutForceSm Line2FitRegion
	duplicate/o/r=(startline2fit,endline2fit) CutForceSuperSm Line2FitRegion
	duplicate/o/r=(startline2out,endline2out) CutForceSm Line2

	make/free/n=0 Disposal

	variable slope1=MakeLinearFitToWaveLimited(WLC1Slope,Line1FitRegion,Line1,Line1)
	variable slope2=MakeLinearFitToWaveLimited(WLC2Slope,Line2FitRegion,Line2,Line2)
	
	duplicate/o/R=(startline2out,endline2out) CutForceSm ShiftForcebyLine2
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
	
		startingscan=LastPoint(Line1Crossings,Line2Crossings)-2e-4
		endingscan=Line2Crossings[0]-2e-6
	elseif(TransType==-1)
		startingscan=LastPoint(Line1Crossings,Line2Crossings)-2e-4
		endingscan=Line2Crossings[0]-2e-6
	endif

	if(endingscan>=endline1out)
	endingscan=endline1out
	endif

	duplicate/free/R=(startingscan,endingscan)  CutForceWave RawShiftLine
	//duplicate/o/R=(startingscan,endingscan)  CutForceSm RawShiftLine

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
	Duplicate/o WLC1 CR_WLC1
	Duplicate/o WLC2 CR_WLC2
	duplicate/o Line1 CR_Line1
	duplicate/o Line2 CR_Line2
	Duplicate/o Line1Crossings CR_Line1Crossings
	Duplicate/o Line2Crossings CR_Line2Crossings
	Duplicate/o Line1Force CR_Line1Force
	Duplicate/o Line2Force CR_Line2Force
	Duplicate/o RawShiftLine CR_RawShiftLine
	duplicate/o CutForceSm CR_CutForceSm
	duplicate/o WLC2Crossings CR_WLC2Crossings
	duplicate/o WLC2Forces CR_WLC2Forces
	duplicate/o WLC1Crossings CR_WLC1Crossings
	duplicate/o WLC1Forces CR_WLC1Forces
		duplicate/o Line1FitRegion CR_Line1FiTregion
	duplicate/o Line2FitRegion CR_Line2FitRegion

	duplicate/o ShiftForcebyWLC1 CR_ShiftForcebyWLC1
	duplicate/o RawShiftWLC CR_RawShiftWLC
endif


	return Result
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
	
	wave w_coef,w_sigma

	Resulting=x*w_coef[1]+w_coef[0]
	variable result=w_coef[1]
	killwaves w_coef ,w_sigma
	duplicate/o Resulting WaveOut
	return result
end

Static Function MakeLinearFitToWaveFixed(SlopeIn,WaveIn,Overlap,Waveout)
	variable SlopeIn
	wave wavein,waveout,Overlap
	duplicate/free Overlap Resulting
	
	K1=SlopeIn
	CurveFit/H="01"/Q/W=2/NTHR=0 line Wavein
	
	
	wave w_coef,w_sigma

	Resulting=x*w_coef[1]+w_coef[0]
	variable result=w_coef[1]
	killwaves w_coef ,w_sigma
	duplicate/o Resulting WaveOut
	return result
end

Static Function MakeLinearFitToWaveLimited(SlopeIn,WaveIn,Overlap,Waveout)
	variable SlopeIn
	wave wavein,waveout,Overlap
	duplicate/free Overlap Resulting
	K1=SlopeIn
	CurveFit/H="01"/Q/W=2/NTHR=0 line Wavein
	wave w_coef
	Make/free/T/N=2 T_Constraints
	if(Slopein<0)
		T_Constraints[0] = {"K1>"+num2str(2*SlopeIn),"K1<"+num2str(-2*SlopeIn)}

	else
	T_Constraints[0] = {"K1<"+num2str(2*SlopeIn),"K1>"+num2str(-2*SlopeIn)}
	endif

	FuncFit/ODR=1/Q/W=2/NTHR=0 DE_linefit W_coef Wavein/C=T_Constraints
	
	wave w_coef,w_sigma
	Resulting=x*w_coef[1]+w_coef[0]
	variable result=w_coef[1]
	killwaves w_coef ,w_sigma
	duplicate/o Resulting WaveOut
	return result
end



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

Static Function FixPicksSingle(PntWave,Direction,index,ForceWave)
	wave PntWave,ForceWave
	variable Direction //Down=-1, Up=1
	variable index
	variable prepnt
	variable nextpnt
	if(Direction==-1)
		wave OtherPnt=$replacestring("PntD",nameofwave(pntWave),"PntU")
	prepnt=OtherPnt[index]
		if(index==(numpnts(OtherPnt)-1))
			nextpnt=numpnts(ForceWave)
		else
						nextpnt=OtherPnt[index+1]

		endif

	elseif(Direction==1)
		wave OtherPnt=$replacestring("PntU",nameofwave(pntWave),"PntD")
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
	//print PntWave[index];print prepnt;print nextpnt
	DE_CorrRup#CalculateBestCrossingGuess(FWIN,SepWave,WLC,PntWave[index],prepnt,nextpnt,Direction,Test,CopyWavesOUt=1)
	MakeASinglePlot(ForceWave)


	
end

Static Function MakeASinglePlot(ForceWave)
	wave ForceWave

	DoWindow SinglePoint
	if(V_flag==1)
	killwindow SinglePoint
	endif
	wave FSm=$(nameofwave(ForceWave)+"_Sm")
	display/N=SinglePoint ForceWave;appendtograph/W=SinglePoint FSM 
	ModifyGraph/W=SinglePoint rgb($nameofwave(ForceWave))=(65535,49151,49151)
	ModifyGraph/W=SinglePoint muloffset={0,-1}
	AppendToGraph/W=SinglePoint $"CR_WLC1",$"CR_WLC2"
	ModifyGraph/W=SinglePoint lsize($"CR_WLC1")=2,rgb($"CR_WLC1")=(0,0,0),lsize($"CR_WLC2")=2,rgb($"CR_WLC2")=(0,0,0)
	
	AppendToGraph/W=SinglePoint $"CR_Line1",$"CR_Line2"
	ModifyGraph/W=SinglePoint lsize($"CR_Line1")=2,rgb($"CR_Line1")=(23130,23130,23130),lsize($"CR_Line2")=2,rgb($"CR_Line2")=(23130,23130,23130)
	
	AppendToGraph/W=SinglePoint $"CR_WLC1Forces" vs $"CR_WLC1Crossings"
	AppendToGraph/W=SinglePoint $"CR_WLC2Forces" vs $"CR_WLC2Crossings"
	ModifyGraph/W=SinglePoint mode($"CR_WLC1Forces")=3,marker($"CR_WLC1Forces")=19,rgb($"CR_WLC1Forces")=(19789,44975,19018),useMrkStrokeRGB($"CR_WLC1Forces")=1
	ModifyGraph/W=SinglePoint mode($"CR_WLC1Forces")=3,mode($"CR_WLC2Forces")=3,marker($"CR_WLC2Forces")=16,rgb($"CR_WLC2Forces")=(14906,32382,47288),useMrkStrokeRGB($"CR_WLC2Forces")=1
	wave T1=$"Test"
	AppendToGraph/W=SinglePoint T1[2,2] vs T1[0,0]
	ModifyGraph/W=SinglePoint mode($nameofwave(T1))=3,marker($nameofwave(T1))=18,rgb($nameofwave(T1))=(36873,14755,58982),useMrkStrokeRGB($nameofwave(T1))=1
	AppendToGraph/W=SinglePoint T1[3,3] vs T1[1,1]
	ModifyGraph/W=SinglePoint mode($nameofwave(T1)+"#1")=3,marker($nameofwave(T1)+"#1")=26,rgb($nameofwave(T1)+"#1")=(65535,0,52428),useMrkStrokeRGB($nameofwave(T1)+"#1")=1
	SetAxis/W=SinglePoint bottom T1[0]-5e-3,T1[0]+5e-3
	SetAxis/W=SinglePoint/A=2 left
end


