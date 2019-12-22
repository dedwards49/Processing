#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma ModuleName=DE_CountRates

Static Function/C FindCrossTimesForHistogram(HistogramIn,State,WLCPArms,WaveOut)

	wave HistogramIn,WLCPArms,Waveout
	String State
	variable index
		variable Folded
	StrSwitch(State)
		case "Folded":
			Folded=1
			break
		case "Unfolded":
			Folded=0
			break
		default:
			print "Bad String input: State"
			return -1
	endswitch
	
	variable startForce= pnt2x(HistogramIn,index)
	variable endforce=pnt2x(HistogramIn,index+1)
	
	variable FOldedLC=WLCPArms[0]
	variable UnfoldedLC=WLCPArms[1]
	Variable FOrceOff=WLCPArms[2]
	variable SepOff=WLCPArms[3]
	variable EXtension
	variable TimeOutgoing,TimeIngoing
	
	make/free/n=(numpnts(HistogramIn)+1) FreeHist
	variable n
	for(n=0;n<numpnts(HistogramIn)+1;n+=1)
		 startForce= pnt2x(HistogramIn,n)
		if(Folded==1)
		FreeHist[n][0]=	 DE_WLC#ReturnExtentionatForce(startForce+FOrceOff,.4e-9,FOldedLC,298)+SepOff
		elseif(Folded==0)
		FreeHist[n][0]=	DE_WLC#ReturnExtentionatForce(startForce+FOrceOff,.4e-9,UnFOldedLC,298)+SepOff

	endif
	
	endfor

	duplicate/o FreeHist WaveOut
end

Static Function FindDwellsforSingleRamp(FFilt,SFilt,StateWave,FoldedHistogram,UnfoldedHistogram,WLCParms,index,ReturnTimesFolded,ReturnTimesUnfolded)
	wave FoldedHistogram,UnfoldedHistogram,WLCParms,FFilt,SFilt,StateWave,ReturnTimesFolded,ReturnTimesUnfolded
	variable index
	make/n=0/free CrossExFolded,CrossExUNfolded
	
	make/free/n=(dimsize(StateWave,0)) Points,RupForce,Type,Trace
	Points=StateWave[p][0]
	RupForce=-StateWave[p][1]
	Type=StateWave[p][2]
	Trace=StateWave[p][3]
	Extract/INDX/Free Points, LocalIndex, Trace==index
	make/free/n=(numpnts(LocalIndex)+1) LocalPoints,LocalType,LocalTrace
	LocalPoints[]=Points[LocalIndex[0]+p][0]
	LocalType[]=Type[LocalIndex[0]+p][2]
	LocalTrace[]=Trace[LocalIndex[0]+p][2]
	FindValue/V=-2/T=.1 LocalType
	variable turnaround=v_value
	duplicate/free/R=[LocalPoints[0],LocalPoints[turnaround]] SFilt, OutgoingSepWave,OutgoingSepFolded,OutgoingSepUnFolded
	duplicate/free/R=[LocalPoints[turnaround],LocalPoints[numpnts(localpoints)-2]] SFilt, IncomingSepWave,IncomingSepFolded,IncomingSepUNfolded

	DE_NewDudko#GenerateSepLine(FFilt,SFilt,StateWave,index,0,1,OutgoingSepFolded)
	DE_NewDudko#GenerateSepLine(FFilt,SFilt,StateWave,index,0,0,OutgoingSepUnFolded)
	DE_NewDudko#GenerateSepLine(FFilt,SFilt,StateWave,index,1,1,IncomingSepFolded)
	DE_NewDudko#GenerateSepLine(FFilt,SFilt,StateWave,index,1,0,IncomingSepUNfolded)

	FindCrossTimesForHistogram(FoldedHistogram,"UnFolded",WLCPArms,CrossExFolded)
	FindCrossTimesForHistogram(UnFoldedHistogram,"folded",WLCPArms,CrossExUNfolded)
	make/free/n=(numpnts(CrossExFolded),2) CrossTimeFolded
	make/free/n=(numpnts(CrossExunFolded),2) CrossTimeUnFolded
	variable n
	variable/C Test
	for(n=0;n<dimsize(CrossExFolded,0);n+=1)
		
		Test=ExtensiontoTime(CrossExFolded[n],FFilt,SFilt,IncomingSepFolded,OutgoingSepFolded,LocalPoints,turnaround,"unFolded")
		CrossTimeFolded[n][0]=real(Test)
		CrossTimeFolded[n][1]=imag(Test)
	endfor
	
	for(n=0;n<dimsize(CrossExunFolded,0);n+=1)

		Test=ExtensiontoTime(CrossExunFolded[n],FFilt,SFilt,IncomingSepFolded,OutgoingSepFolded,LocalPoints,turnaround,"Folded")
		CrossTimeunFolded[n][0]=real(Test)
		CrossTimeunFolded[n][1]=imag(Test)
		

	endfor

	duplicate/o CrossTimeFolded ReturnTimesFolded
	duplicate/o  CrossTimeunFolded ReturnTimesUnfolded
end

Static Function/C ExtensiontoTime(Extension,ForceWave,SepWave,IncomingSepWave,OutgoingSepWave,LocalPoints,turnaround,State)
	variable Extension,turnaround
	wave ForceWave,SepWave,IncomingSepWave,OutgoingSepWave,LocalPoints
	string State
	variable Folded
	StrSwitch(State)
		case "Folded":
			Folded=1
			break
		case "Unfolded":
			Folded=0
			break
		default:
			print "Bad String input: State"
			return -1
	endswitch
	
	
	variable TimeOutgoing,TimeIngoing
	if(wavemin(OutgoingSepWave)>EXtension)
		TimeOutgoing=pnt2x(Sepwave,LocalPoints[0])
	elseif(wavemax(OutgoingSepWave)<EXtension)
		TimeOutgoing=pnt2x(Sepwave,LocalPoints[turnaround])
	else
		FindLevels/Q OutgoingSepWave EXtension
		wave W_FindLevels
		TimeOutgoing = mean(W_FindLevels)//-pnt2x(FirstSepWave,0)
	endif
	
	
	if(wavemax(IncomingSepWave)<EXtension)
		TimeIngoing=pnt2x(Sepwave,LocalPoints[turnaround])
	elseif(wavemin(IncomingSepWave)>EXtension)
		TimeIngoing=pnt2x(Sepwave,LocalPoints[numpnts(LocalPoints)-2])
	else
		FindLevels/Q IncomingSepWave EXtension
		TimeIngoing= mean(W_FindLevels)//-pnt2x(FirstSepWave,0)
	endif
	killwaves W_FindLevels
	return cmplx(TimeOutgoing,Timeingoing)
end

static Function AddUpTimes(Fsm,StateWave,index,FoldedTimes,UnfoldedTimes,FoldedHistogram,UnfoldedHistogram,FTimesOut,UTimesOut,FHistOut,UHistOut,ModStateWaveOUt)
	wave Fsm,StateWave, FoldedTimes,UnfoldedTimes,FoldedHistogram,UnfoldedHistogram,FHistOut,UHistOut,FTimesOut,UTimesOut,ModStateWaveOUt
	variable index
	make/free/n=(dimsize(StateWave,0)) Points,RupForce,Type,Trace
	Points=StateWave[p][0]
	RupForce=-StateWave[p][1]
	Type=StateWave[p][2]
	Trace=StateWave[p][3]
	
	duplicate/free FoldedHistogram FreeFHist
	duplicate/free UnFoldedHistogram FreeUFHist
	FreeFHist=0
	FreeUFHist=0
	Extract/INDX/Free Points, LocalIndex, Trace==index
	make/free/n=(numpnts(LocalIndex)+1) LocalPoints,LocalType,LocalTrace
	LocalPoints[]=Points[LocalIndex[0]+p][0]
	LocalType[]=Type[LocalIndex[0]+p][2]
	LocalTrace[]=Trace[LocalIndex[0]+p][2]
	FindValue/V=-2/T=.1 LocalType
	variable turnaround=v_value
	variable n=1,targettime
	make/free/n=(dimsize(FoldedTimes,0)) FoldedTimesOut,FoldedTimesIn
	FoldedTimesOut=FoldedTimes[p][0]
	FoldedTimesIn=FoldedTimes[p][1]

	make/free/n=(dimsize(unFoldedTimes,0)) unFoldedTimesOut,unFoldedTimesIn
	unFoldedTimesOut=unFoldedTimes[p][0]
		unFoldedTimesIn=unFoldedTimes[p][1]

	make/free/n=(dimsize(unFoldedTimes,0)-1) unFoldedTimeSpent
	make/free/n=(dimsize(FoldedTimes,0)-1) FoldedTimeSpent
	unFoldedTimeSpent=0
	FoldedTimeSpent=0
	
	variable currentlyFolded=1
	variable timeentered=max(unFoldedTimesOut[0],pnt2x(Fsm,LocalPoints[0]))
	variable pointentered,targetpoint
	make/free/n=(numpnts(LocalPoints)-4,20) FreeNewState



	for(n=1;n<turnaround;n+=1)
		targettime=pnt2x(Fsm,LocalPoints[n])
		
		FreeNewState[n-1][19]=targettime
		if(currentlyfolded==0)

			FindLevel/Q FoldedTimesOut, timeentered
			pointentered=v_levelx
			FindLevel/Q FoldedTimesOut, targettime
			targetpoint=v_levelx

			FreeFHist[floor(targetpoint)]+=1
			FreeNewState[n-1][0]=(-1)^currentlyFolded
			FreeNewState[n-1][1]=1
			FreeNewState[n-1][2]=pointentered
			FreeNewState[n-1][3]=targetpoint
			if(floor(pointentered)==floor(targetpoint))
				FoldedTimeSpent[floor(targetpoint)]+=targettime-timeentered
				FreeNewState[n-1][4]=targettime-timeentered

			elseif(floor(pointentered)+1==floor(targetpoint))
				FoldedTimeSpent[floor(pointentered)]+=FoldedTimesOut[floor(pointentered)+1]-timeentered
				FoldedTimeSpent[floor(targetpoint)]+=targettime-FoldedTimesOut[floor(targetpoint)]
				FreeNewState[n-1][4]=FoldedTimesOut[floor(pointentered)+1]-timeentered
				FreeNewState[n-1][5]=targettime-FoldedTimesOut[floor(targetpoint)]

			else
				FoldedTimeSpent[floor(pointentered)]+=FoldedTimesOut[floor(pointentered)+1]-timeentered
				FoldedTimeSpent[floor(pointentered)+1,floor(targetpoint)-1]+=FoldedTimesOut[p+1]-FoldedTimesOut[p]
				FoldedTimeSpent[floor(targetpoint)]+=targettime-FoldedTimesOut[floor(targetpoint)]
				FreeNewState[n-1][4]=FoldedTimesOut[floor(pointentered)+1]-timeentered
				FreeNewState[n-1][5]=1
				FreeNewState[n-1][6]=targettime-FoldedTimesOut[floor(targetpoint)]

			endif

			currentlyFolded=	1
			timeentered=targettime	

		elseif(currentlyfolded==1)
			if(timeentered<=wavemin(unFoldedTimesOut))
				pointentered=0
			else
				FindLevel/Q UnFoldedTimesOut, timeentered
				pointentered=v_levelx
			endif
			FindLevel/Q UnFoldedTimesOut, targettime
			targetpoint=v_levelx


			FreeUFHist[floor(targetpoint)]+=1
			FreeNewState[n-1][0]=(-1)^currentlyFolded
			FreeNewState[n-1][1]=1
			FreeNewState[n-1][2]=pointentered
			FreeNewState[n-1][3]=targetpoint
			if(floor(pointentered)==floor(targetpoint))
				unFoldedTimeSpent[floor(targetpoint)]+=targettime-timeentered   
				FreeNewState[n-1][4]=targettime-timeentered


			elseif(floor(pointentered)+1==floor(targetpoint))
				unFoldedTimeSpent[floor(pointentered)]+=unFoldedTimesOut[floor(pointentered)+1]-timeentered
				unFoldedTimeSpent[floor(targetpoint)]+=targettime-unFoldedTimesOut[floor(targetpoint)]

				FreeNewState[n-1][4]=unFoldedTimesOut[floor(pointentered)+1]-timeentered
				FreeNewState[n-1][5]=targettime-unFoldedTimesOut[floor(targetpoint)]
			else
		
				unFoldedTimeSpent[floor(pointentered)]+=unFoldedTimesOut[floor(pointentered)+1]-timeentered
				unFoldedTimeSpent[floor(pointentered)+1,floor(targetpoint)-1]+=unFoldedTimesOut[p+1]-unFoldedTimesOut[p]
				unFoldedTimeSpent[floor(targetpoint)]+=targettime-unFoldedTimesOut[floor(targetpoint)]
				FreeNewState[n-1][4]=unFoldedTimesOut[floor(pointentered)+1]-timeentered
				FreeNewState[n-1][5]=1
				FreeNewState[n-1][6]=targettime-unFoldedTimesOut[floor(targetpoint)]

			endif

			currentlyFolded=	0
			timeentered=targettime	

		endif
		

	endfor
	
	//	//////////////////////
	FindLevel/Q FoldedTimesOut, timeentered
	pointentered=v_levelx
	targettime=pnt2x(FSm,LocalPoints[turnaround])
	
	if(targettime>=wavemax(FoldedTimesOut))
		targettime=wavemax(FoldedTimesOut)
		targetpoint=numpnts(FoldedTimesOut)-2
	else
		FindLevel/Q FoldedTimesOut, targettime
		targetpoint=v_levelx
	endif
		FreeNewState[turnaround-1][7]=pointentered
		FreeNewState[turnaround-1][8]=targetpoint
		FreeNewState[turnaround-1][12]=-1

	if(floor(pointentered)==floor(targetpoint))
		FoldedTimeSpent[floor(targetpoint)]+=targettime-timeentered
		FreeNewState[turnaround-1][9]=targettime-timeentered

	elseif(floor(pointentered)+1==floor(targetpoint))
		FoldedTimeSpent[floor(pointentered)]+=FoldedTimesOut[floor(pointentered)+1]-timeentered
		FoldedTimeSpent[floor(targetpoint)]+=targettime-FoldedTimesOut[floor(targetpoint)]
		FreeNewState[turnaround-1][9]=FoldedTimesOut[floor(pointentered)+1]-timeentered
		FreeNewState[turnaround-1][10]=targettime-FoldedTimesOut[floor(targetpoint)]
	else
				
		FoldedTimeSpent[floor(pointentered)]+=FoldedTimesOut[floor(pointentered)+1]-timeentered
		FoldedTimeSpent[floor(pointentered)+1,floor(targetpoint)-1]+=FoldedTimesOut[p+1]-FoldedTimesOut[p]
		FoldedTimeSpent[floor(targetpoint)]+=targettime-FoldedTimesOut[floor(targetpoint)]
		FreeNewState[turnaround-1][9]=FoldedTimesOut[floor(pointentered)+1]-timeentered
		FreeNewState[turnaround-1][10]=1
		FreeNewState[turnaround-1][11]=targettime-FoldedTimesOut[floor(targetpoint)]
	endif
	
	
	////////////////////////////////////////////////////////////////////
	
			
			
	currentlyFolded=0
	timeentered=max(wavemin(FoldedTimesIn),pnt2x(Fsm,LocalPoints[turnaround]))
	
	for(n=turnaround+1;n<numpnts(LocalPoints)-2;n+=1)
		targettime=pnt2x(Fsm,LocalPoints[n])
				FreeNewState[n-2][19]=targettime

		if(currentlyfolded==0)
			if(timeentered<=wavemin(FoldedTimesIn))
				pointentered=numpnts(FoldedTimesIn)-2
			else
				FindLevel/Q FoldedTimesIn, timeentered
				pointentered=v_levelx
			endif
			FindLevel/Q FoldedTimesIn, targettime
			targetpoint=v_levelx

			FreeFHist[floor(targetpoint)]+=1
			FreeNewState[n-2][0]=(-1)^currentlyFolded
			FreeNewState[n-2][1]=-1
			FreeNewState[n-2][2]=pointentered
			FreeNewState[n-2][3]=targetpoint
			if(floor(pointentered)==floor(targetpoint))
				FoldedTimeSpent[floor(targetpoint)]+=targettime-timeentered
				FreeNewState[n-2][4]=targettime-timeentered

			elseif(floor(pointentered)-1==floor(targetpoint))
				FoldedTimeSpent[floor(pointentered)]+=FoldedTimesIn[floor(pointentered)]-timeentered
				FoldedTimeSpent[floor(targetpoint)]+=targettime-FoldedTimesIn[floor(targetpoint)+1]
				FreeNewState[n-2][4]=FoldedTimesIn[floor(pointentered)]-timeentered
				FreeNewState[n-2][5]=targettime-FoldedTimesIn[floor(targetpoint)+1]
			else
				FoldedTimeSpent[floor(pointentered)]+=FoldedTimesIn[floor(pointentered)]-timeentered
				FoldedTimeSpent[floor(targetpoint)+1,floor(pointentered)-1]+=FoldedTimesin[p]-FoldedTimesin[p+1]
				FoldedTimeSpent[floor(targetpoint)]+=targettime-FoldedTimesIn[floor(targetpoint)+1]
				FreeNewState[n-2][4]=FoldedTimesIn[floor(pointentered)]-timeentered
				FreeNewState[n-2][5]=1
				FreeNewState[n-2][6]=targettime-FoldedTimesIn[floor(targetpoint)+1]
			endif

			currentlyFolded=	1
			timeentered=targettime		
			
		elseif(currentlyfolded==1)
			FindLevel/Q UnFoldedTimesIn, timeentered
			pointentered=v_levelx

			FindLevel/Q UnFoldedTimesIn, targettime

			targetpoint=v_levelx

			FreeUFHist[floor(targetpoint)]+=1
			FreeNewState[n-2][0]=(-1)^currentlyFolded
			FreeNewState[n-2][1]=-1
			FreeNewState[n-2][2]=pointentered
			FreeNewState[n-2][3]=targetpoint
			if(floor(pointentered)==floor(targetpoint))
				unFoldedTimeSpent[floor(targetpoint)]+=targettime-timeentered
				FreeNewState[n-2][4]=targettime-timeentered

			elseif(floor(pointentered)-1==floor(targetpoint))
				unFoldedTimeSpent[floor(pointentered)]+=unFoldedTimesIn[floor(pointentered)]-timeentered
				unFoldedTimeSpent[floor(targetpoint)]+=targettime-unFoldedTimesIn[floor(targetpoint)+1]

				FreeNewState[n-2][4]=unFoldedTimesIn[floor(pointentered)]-timeentered
				FreeNewState[n-2][5]=targettime-unFoldedTimesIn[floor(targetpoint)+1]
			else
				unFoldedTimeSpent[floor(pointentered)]+=unFoldedTimesIn[floor(pointentered)]-timeentered
				unFoldedTimeSpent[floor(targetpoint)+1,floor(pointentered)-1]+=unFoldedTimesin[p]-unFoldedTimesin[p+1]
				unFoldedTimeSpent[floor(targetpoint)]+=targettime-unFoldedTimesIn[floor(targetpoint)+1]
				FreeNewState[n-2][4]=unFoldedTimesIn[floor(pointentered)]-timeentered

				FreeNewState[n-2][5]=1
				FreeNewState[n-2][6]=targettime-unFoldedTimesIn[floor(targetpoint)+1]
			endif

			currentlyFolded=	0
			timeentered=targettime	
		endif
	endfor
	
	FindLevel/Q unFoldedTimesin, timeentered
	pointentered=v_levelx
	targettime=pnt2x(FSm,LocalPoints[numpnts(LocalPoints)-2])

	if(targettime>=wavemax(unFoldedTimesin))
		targettime=wavemax(unFoldedTimesin)
		targetpoint=0
	else
		FindLevel/Q unFoldedTimesin, targettime
		targetpoint=v_levelx
	endif
		FreeNewState[numpnts(LocalPoints)-5][13]=pointentered
		FreeNewState[numpnts(LocalPoints)-5][14]=targetpoint
		FreeNewState[numpnts(LocalPoints)-5][18]=1
	
	if(floor(pointentered)==floor(targetpoint))
		unFoldedTimeSpent[floor(targetpoint)]+=targettime-timeentered
		FreeNewState[numpnts(LocalPoints)-5][15]=targettime-timeentered
	elseif(floor(pointentered)-1==floor(targetpoint))
		unFoldedTimeSpent[floor(pointentered)]+=unFoldedTimesIn[floor(pointentered)]-timeentered
		unFoldedTimeSpent[floor(targetpoint)]+=targettime-unFoldedTimesIn[floor(targetpoint)+1]
		FreeNewState[numpnts(LocalPoints)-5][15]=unFoldedTimesIn[floor(pointentered)]-timeentered
		FreeNewState[numpnts(LocalPoints)-5][16]=targettime-unFoldedTimesIn[floor(targetpoint)+1]
	else
		unFoldedTimeSpent[floor(pointentered)]+=unFoldedTimesIn[floor(pointentered)]-timeentered
		unFoldedTimeSpent[floor(targetpoint)+1,floor(pointentered)-1]+=unFoldedTimesin[p]-unFoldedTimesin[p+1]
		unFoldedTimeSpent[floor(targetpoint)]+=targettime-unFoldedTimesIn[floor(targetpoint)+1]
		FreeNewState[numpnts(LocalPoints)-5][15]=unFoldedTimesIn[floor(pointentered)]-timeentered
		FreeNewState[numpnts(LocalPoints)-5][10]=1
		FreeNewState[numpnts(LocalPoints)-5][17]=targettime-unFoldedTimesIn[floor(targetpoint)+1]
	endif

	duplicate/o FoldedTimeSpent FTimesOut
	duplicate/o unFoldedTimeSpent UTimesOut
	duplicate/o FreeFHist FHistOut
	duplicate/o FreeUFHist UHistOut
	duplicate/o FreeNewState ModStateWaveOUt
end

Static Function PushThrough(ForceWave,SepWave,Smoothed,StateWave,InFHistograms,InUHistograms,WLCParms,Shifting)
	wave ForceWave,SepWave,StateWave,InFHistograms,InUHistograms,WLCParms
	string Shifting
	variable smoothed
	if(smoothed==1)
		duplicate/free ForceWave FSM
		duplicate/free SepWave SSm
	else
	make/o/n=0 FSm,SSM;DE_Filtering#FilterForceSep(ForceWave,SepWave,FSm,SSm,"svg",51)
	endif
	variable/c Shifts=CorrectShift(ForceWave,Shifting)
	SSm+=imag(Shifts)
	Fsm+=real(Shifts)
	make/free/n=0 ReturnTimesFolded,ReturnTimesUnfolded
	make/free/n=(dimsize(StateWave,0)) Points,RupForce,Type,Trace
	Points=StateWave[p][0]
	RupForce=-StateWave[p][1]
	Type=StateWave[p][2]
	Trace=StateWave[p][3]
	variable top=wavemax(trace)
	variable n
	duplicate/o InFHistograms, FAddHist,FAddtime
	FAddHist=0
	FAddtime=0
	duplicate/o InUHistograms, UAddHist,UAddTime
	UAddHist=0
	UAddTime=0
	variable bottom=0
	//top=11
	for(n=bottom;n<top;n+=1)
		make/free/n=0 Ftimes,UTimes,FHistOut,UHistOut,ModStateWaveOut
		FindDwellsforSingleRamp(FSm,SSM,StateWave,InFHistograms,InUHistograms,WLCParms,n,ReturnTimesFolded,ReturnTimesUnfolded)
		AddUpTimes(Fsm,StateWave,n,ReturnTimesFolded,ReturnTimesUnfolded,InFHistograms,InUHistograms,Ftimes,UTimes,FHistOut,UHistOut,ModStateWaveOut)
		make/free/n=(dimsize(ModStateWaveOut,0)) Nthing
		NThing=n
		concatenate/o/NP=1 {ModStateWaveOut,NThing}, FreeAdded
	

		if(n==bottom)
			duplicate/o FreeAdded ResultstoHold
			duplicate/o ReturnTimesFolded ReturnTimesFoldedMassive
			duplicate/o ReturnTimesUnfolded ReturnTimesUnfoldedMassive
		else 
			FreeAdded[0][7,12]=ResultstoHold[dimsize(ResultstoHold,0)-1][q+6]
			FreeAdded[0][12]=2
			ResultstoHold[dimsize(ResultstoHold,0)-1][13,18]=0
			concatenate/o/NP=0 {ResultstoHold,FreeAdded}, FreedHold
			concatenate/o/NP=1 {ReturnTimesFoldedMassive,ReturnTimesFolded}, ReturnTimesFoldedFree
			concatenate/o/NP=1 {ReturnTimesUnfoldedMassive,ReturnTimesUnfolded}, ReturnTimesUnfoldedFree

			duplicate/o FreedHold ResultstoHold
			duplicate/o ReturnTimesFoldedFree ReturnTimesFoldedMassive
			duplicate/o ReturnTimesUnfoldedFree ReturnTimesUnfoldedMassive
		endif
		FAddHist+=FHistOut
		FAddtime+=FTimes
		UAddHist+=UHistOut
		UAddTime+=UTimes
	endfor
	duplicate/o UAddHist URate
	duplicate/o FAddHist FRate


	killwaves FreeAdded,FreedHold,ReturnTimesFoldedFree,ReturnTimesUnfoldedFree
	SetScale/P x pnt2x(InFHistograms,0),dimdelta(InFHistograms,0),"", ReturnTimesFolded;SetScale/P x pnt2x(InUHistograms,0),dimdelta(InUHistograms,0),"", ReturnTimesUnfolded
	URate=UAddHist[p]/UAddTime[p]
	FRate=FAddHist[p]/FAddtime[p]
	DE_Filtering#FilterForceSep(ForceWave,SepWave,FSm,SSm,"svg",51)
	SSm+=imag(Shifts)
	Fsm+=real(Shifts)
end


Static Function/C CorrectShift(ForceWave,Type)
	wave ForceWave
	String Type
	String UsedAllowed
	variable ForceAdj,SepAdj,ShiftedFOff,ShiftedSOff,UnShiftedFOff,UnShiftedSOff,UsedFoff,UsedSoff,AltFOff,AltSOff
	wave w1=ForceWave
	if(cmpstr(Stringbykey("UsedAlignmentFShift",note(w1),":","\r"),"")==0)
	print "Shifting Problem"
		return cmplx(0,0)
	else
		UsedFoff=str2num(Stringbykey("UsedAlignmentFShift",note(w1),":","\r"))
		UsedSOff=str2num(Stringbykey("UsedAlignmentSShift",note(w1),":","\r"))
		AltFOff=str2num(Stringbykey("AltAlignmentFShift",note(w1),":","\r"))
		AltSOff=str2num(Stringbykey("AltAlignmentSShift",note(w1),":","\r"))
		if(str2num(Stringbykey("UsedAlignmentSShift",note(w1),":","\r"))==0)
			UsedAllowed="No"
			ShiftedFOff=str2num(Stringbykey("AltAlignmentFShift",note(w1),":","\r"))
			ShiftedSOff=str2num(Stringbykey("AltAlignmentSShift",note(w1),":","\r"))
			UnShiftedFOff=str2num(Stringbykey("UsedAlignmentFShift",note(w1),":","\r"))
			UnShiftedSOff=str2num(Stringbykey("UsedAlignmentSShift",note(w1),":","\r"))
				
		else
			UsedAllowed="Yes"
			unShiftedFOff=str2num(Stringbykey("AltAlignmentFShift",note(w1),":","\r"))
			unShiftedSOff=str2num(Stringbykey("AltAlignmentSShift",note(w1),":","\r"))
			ShiftedFOff=str2num(Stringbykey("UsedAlignmentFShift",note(w1),":","\r"))
			ShiftedSOff=str2num(Stringbykey("UsedAlignmentSShift",note(w1),":","\r"))
		endif
			
		Strswitch(Type)
			case "Used":
				ForceAdj=0
				SepAdj=0
				break
			case "Shifted":
				ForceAdj=UsedFoff-ShiftedFOff
				SepAdj=UsedSOff-ShiftedSOff
				break	
			case "UnShifted":
				ForceAdj=UsedFoff-unShiftedFOff
				SepAdj=UsedSOff-UnShiftedSOff
				break	
			case "Alt":
				ForceAdj=UsedFoff-AltFOff
				SepAdj=UsedSOff-AltSOff
				break	
			default:
				ForceAdj=0
				SepAdj=0
		endswitch
	endif
		
	
	
	return cmplx(ForceAdj,SepAdj)
end

Function RunJustFromFolder(FolderString,Shifting)
	String FolderString,Shifting
	wave ForceWave=$(FolderString+":"+stringfromlist(0,DE_CountRates#ListWaves(FolderString,"*Force_Align")))
	wave SepWave=$(FolderString+":"+stringfromlist(0,DE_CountRates#ListWaves(FolderString,"*Sep_Align")))
	wave WLC=$(FolderString+":"+stringfromlist(0,DE_CountRates#ListWaves(FolderString,"*WLCParms_Adj")))
	string UFHIstString=(FolderString+":UFHIst")
	string FHIstString=(FolderString+":FHIst")
	make/o/n=12 $UFHIstString/Wave=UFHIst
	make/o/n=16 $FHIstString/Wave=FHIst

	
	SetScale/P x 6e-12,.5e-12,"", FHist;SetScale/P x 8e-12,1e-12,"", UFhist
	
	wave State=$(FolderString+":"+stringfromlist(0,DE_CountRates#ListWaves(FolderString,"*2States*")))
	make/free/n=0 WLCOut;DE_CountRates#CorrectWLC(WLC,WLCOut,ForceWave,SepWave,Shifting)
	PushThrough(ForceWave,SepWave,0,State,FHist,UFHist,WLCOut,Shifting)
	string NameString=stringfromlist(1,FolderString,":")
	wave FAddHist,FaddTime,Frate,URate,UAddHist,UAddTime,FSm,SSm,ResultstoHold,ReturnTimesFoldedMassive,ReturnTimesUnfoldedMassive
	duplicate/o FAddHist $(NameString+"_FoldHist")
	duplicate/o FAddTime $(NameString+"_FoldTimes")
	duplicate/o Frate $(NameString+"_FoldRate")
	duplicate/o UAddHist $(NameString+"_UNfoldHist")
	duplicate/o UAddTime $(NameString+"_UnfoldTimes")
	duplicate/o Urate $(NameString+"_UnfoldRate")
	duplicate/o UAddTime $(NameString+"_UnfoldTimes")
	duplicate/o ResultstoHold $(NameString+"_AllTrans")
	duplicate/o ReturnTimesFoldedMassive $(NameString+"_FoldTimes")
	duplicate/o ReturnTimesUnfoldedMassive $(NameString+"_UnfoldTimes")

	killwaves FAddHist,FaddTime,Frate,URate,UAddHist,UAddTime,FSm,SSm,ResultstoHold,ReturnTimesFoldedMassive,ReturnTimesUnfoldedMassive
	killwaves UFHIst,FHISt
	//PlotForceandFit(FolderString,Shifting,10e-9)
	//PlotThingsNiceWithaNameString(NameString)
end

Function RunOnCombinedData(BaseString,ForceShifting,WLCShifting,[AltWLC])
	String BaseString,ForceShifting,WLCShifting
	wave AltWLC
	if(ParamisDefault(AltWLC))
		wave WLC=$(BaseString+"_WLCParms")
	else
		wave WLC=AltWLC
	
	endif
	variable WLCTYPE
	wave ForceWave=$(BaseString+"_Force")
	wave SepWave=$(BaseString+"_Sep")
	wave StateWave=$(BaseString+"_CompState")
	//stupidly I wrote the code so that it takes the WLC of the form {UnfoldedLC,FoldedLC,ForceShift,SepShift}
	//this just fixes that a bit if the WLC fits are from the older Dudko analysis (of the form {FoldedLC,UnfoldedLC,Forcehift,SepShift}
	if(strsearch(nameofwave(WLC),"_WLCParms",0)!=0)
		duplicate/free WLC WLCOut

		WLCOut[2]*=-1
	elseif(strsearch(nameofwave(WLC),"_WLCParms_Adj",0)!=0)
		make/free/n=0 WLCOut;DE_CountRates#CorrectWLC(WLC,WLCOut,ForceWave,SepWave,WLCShifting)

	endif
	string UFHIstString=(BaseString+"_Count_UFHIst")
	string FHIstString=(BaseString+"_Count_FHIst")
	make/o/n=16 $UFHIstString/Wave=UFHIst
	make/o/n=22 $FHIstString/Wave=FHIst
	SetScale/P x 6e-12,.5e-12,"", FHist;SetScale/P x 6e-12,1e-12,"", UFhist

	PushThrough(ForceWave,SepWave,0,StateWave,FHist,UFhist,WLCOut,ForceShifting)
	string NameString=stringfromlist(0,nameofwave(ForceWave),"_")+"_"+stringfromlist(1,nameofwave(ForceWave),"_")
	wave FAddHist,FaddTime,Frate,URate,UAddHist,UAddTime,FSm,SSm,ResultstoHold,ReturnTimesFoldedMassive,ReturnTimesUnfoldedMassive
	duplicate/o FAddHist $(NameString+"_FoldHist")
	duplicate/o FAddTime $(NameString+"_FoldTimes")
	duplicate/o Frate $(NameString+"_FoldRate")
	duplicate/o UAddHist $(NameString+"_UNfoldHist")
	duplicate/o UAddTime $(NameString+"_UnfoldTimes")
	duplicate/o Urate $(NameString+"_UnfoldRate")
	duplicate/o UAddTime $(NameString+"_UnfoldTimes")
	duplicate/o ResultstoHold $(NameString+"_AllTrans")
	duplicate/o ReturnTimesFoldedMassive $(NameString+"_FoldTimes")
	duplicate/o ReturnTimesUnfoldedMassive $(NameString+"_UnfoldTimes")
		killwaves FAddHist,FaddTime,Frate,URate,UAddHist,UAddTime,ResultstoHold,ReturnTimesFoldedMassive,ReturnTimesUnfoldedMassive
		killwaves FSm,SSM
	killwaves UFHIst,FHISt
end





Static Function PlotThingsNiceWithaNameString(BaseNameString)
	
	String BaseNameString
	
	wave UnfoldedHist=$(BaseNameString+"_UNfoldHist")
	wave FoldedHist=$(BaseNameString+"_FoldHist")
	wave UnfoldedNum= $(BaseNameString+"_UnfoldTimes")
	wave FoldedNum=$(BaseNameString+"_FoldTimes")
	wave UnfoldedRate=$(BaseNameString+"_UnfoldRate")
	wave FoldedRate=$(BaseNameString+"_FoldRate")

	DoWindow $BaseNameString
	if(V_Flag==1)
	killwindow $BaseNameString
	endif
	Display/N=$BaseNameString UnfoldedRate,FoldedRate
	AppendToGraph/L=L1 FoldedNum,UnfoldedNum
	AppendToGraph/L=L2 FoldedHist,UnfoldedHist
	ModifyGraph freePos(L1)=0
	ModifyGraph freePos(L2)=0
	ModifyGraph axisEnab(left)={0,0.30},axisEnab(L1)={0.35,0.65}
	ModifyGraph axisEnab(L2)={0.75,1}
	ModifyGraph log(left)=1
	ModifyGraph mode($nameofwave(UnfoldedRate))=3,marker($nameofwave(UnfoldedRate))=19, rgb($nameofwave(UnfoldedRate))=(58368,6656,7168);DelayUpdate
	ModifyGraph useMrkStrokeRGB($nameofwave(UnfoldedRate))=1;
	ModifyGraph mode($nameofwave(FoldedRate))=3,marker($nameofwave(FoldedRate))=19;DelayUpdate
	ModifyGraph rgb($nameofwave(FoldedRate))=(14848,32256,47104), useMrkStrokeRGB($nameofwave(FoldedRate))=1

	ModifyGraph mode($nameofwave(FoldedNum))=4,marker($nameofwave(FoldedNum))=29;DelayUpdate
	ModifyGraph rgb($nameofwave(FoldedNum))=(14848,32256,47104),useMrkStrokeRGB($nameofwave(FoldedNum))=1;
	ModifyGraph mode($nameofwave(UnfoldedNum))=4,marker($nameofwave(UnfoldedNum))=29;DelayUpdate
	ModifyGraph useMrkStrokeRGB($nameofwave(UnfoldedNum))=1
	//ModifyGraph lsize($nameofwave(UnfoldedSlope))=1.5, lsize($nameofwave(foldedSlope))=1.5;DelayUpdate
	//ModifyGraph rgb($nameofwave(foldedSlope))=(14848,32256,47104), rgb($nameofwave(UnfoldedSlope))=(58368,6656,7168)
	ModifyGraph mode($nameofwave(FoldedHist))=5,rgb($nameofwave(FoldedHist))=(14848,32256,47104)
	ModifyGraph useBarStrokeRGB($nameofwave(FoldedHist))=1,hbFill($nameofwave(FoldedHist))=2
	ModifyGraph mode($nameofwave(unFoldedHist))=5,rgb($nameofwave(unFoldedHist))=(58368,6656,7168);
	ModifyGraph useBarStrokeRGB($nameofwave(unFoldedHist))=1,hbFill($nameofwave(unFoldedHist))=2
	
	//ErrorBars $nameofwave(UnfoldedRate) Y,wave=(UnfoldedErrorP,UnfoldedErrorM);DelayUpdate
	//ErrorBars $nameofwave(FoldedRate) Y,wave=(FoldedErrorP,FoldedErrorM)


	Label bottom "\\f01Force(pN)"
	ModifyGraph fSize=9,font="Arial"
	ModifyGraph prescaleExp(bottom)=12, muloffset={0,0},prescaleExp(L2)=12,lblPosMode(L1)=1,lblPosMode(L2)=1
	Label left "\\f01Rate\r (1/s)"
	Label L1 "\\f01Number"
	Label L2 "\\f01Hist"
	•ModifyGraph margin(bottom)=29,margin(top)=14,margin(right)=14,margin(left)=43



end



Function CalculateErrors(FolderString,BaseNameString)
	String FolderString,BaseNameString


	CompileTransitionWaves(FolderString,BaseNameString)
	IterativeResampling(FolderString,BaseNameString,500)
end


Static Function CompileTransitionWaves(FolderString,BaseNameString)
	String FolderString,BaseNameString
	Wave InputWave=$(folderstring+":"+BaseNameString+"_AllTrans")
	wave OrigFoldHist=$(folderstring+":"+BaseNameString+"_FoldHist")
	wave OrigUnFoldHist=$(folderstring+":"+BaseNameString+"_UnFoldHist")
	wave FoldedTimes=$(folderstring+":"+BaseNameString+"_FoldTimes")
	wave unFoldedTimes=$(folderstring+":"+BaseNameString+"_unFoldTimes")

	variable n
	variable top=dimsize(InputWave,0)
	duplicate/o OrigFoldHist NewFoldingHist,NewFoldingTimes,NewFoldingRate
	NewFoldingHist=0;NewFoldingTimes=0;NewFoldingRate=0
	duplicate/o OrigUnFoldHist NewUnFoldingHist,NewUnFoldingTimes,NewUnFoldingRate 
	duplicate/o OrigUnFoldHist CurrentUHist,CurrentUTime
	duplicate/o OrigFoldHist CurrentFHist,CurrentFTime

	NewUnFoldingHist=0;NewUnFoldingTimes=0;NewUnFoldingRate=0
	variable targetpoint,pointentered
	variable trace
	for(n=0;n<top;n+=1)
		trace=InputWave[n][20]
		CurrentFHist=0
		CurrentFTime=0
		currentUHist=0
		CurrentUTime=0
		
		make/free/n=(dimsize(FoldedTimes,0)) FoldedTimesOut,FoldedTimesIn
		FoldedTimesOut=FoldedTimes[p][2*trace]
		FoldedTimesIn=FoldedTimes[p][2*trace+1]

		make/free/n=(dimsize(unFoldedTimes,0)) unFoldedTimesOut,unFoldedTimesIn
		unFoldedTimesOut=unFoldedTimes[p][2*trace]
		unFoldedTimesIn=unFoldedTimes[p][2*trace+1]
		
	
	
		targetpoint=InputWave[n][3]
		pointentered=InputWave[n][2]
		if(InputWave[n][0]==-1&&InputWave[n][1]==1)//unfolding OutBound
			NewUnFoldingHist[floor(targetpoint)]+=1
			CurrentUHist[floor(targetpoint)]+=1

			if(floor(pointentered)==floor(targetpoint))
				NewUnFoldingTimes[floor(targetpoint)]+=InputWave[n][4]                                                                                                                                                               
				CurrentUTime[floor(targetpoint)]+=InputWave[n][4]                                                                                                                                                               


			elseif(floor(pointentered)+1==floor(targetpoint))
				NewUnFoldingTimes[floor(pointentered)]+=InputWave[n][4]
				NewUnFoldingTimes[floor(targetpoint)]+=InputWave[n][5]
				CurrentUTime[floor(pointentered)]+=InputWave[n][4]
				CurrentUTime[floor(targetpoint)]+=InputWave[n][5]
			else
				NewUnFoldingTimes[floor(pointentered)]+=InputWave[n][4]
				NewUnFoldingTimes[floor(pointentered)+1,floor(targetpoint)-1]+=unFoldedTimesout[p+1]-unFoldedTimesout[p]
				NewUnFoldingTimes[floor(targetpoint)]+=InputWave[n][6]
				CurrentUTime[floor(pointentered)]+=InputWave[n][4]
				CurrentUTime[floor(pointentered)+1,floor(targetpoint)-1]+=unFoldedTimesout[p+1]-unFoldedTimesout[p]
				CurrentUTime[floor(targetpoint)]+=InputWave[n][6]

			endif
			if(InputWave[n][12]==2)
				

				make/free/n=(dimsize(unFoldedTimes,0)) PrevunFoldedTimesOut,PrevunFoldedTimesIn
				PrevunFoldedTimesOut=unFoldedTimes[p][2*(trace-1)]
				PrevunFoldedTimesIn=unFoldedTimes[p][2*(trace-1)+1]
				
				targetpoint=InputWave[n][8]
				pointentered=InputWave[n][7]
				if(floor(pointentered)==floor(targetpoint))
					NewUnFoldingTimes[floor(targetpoint)]+=InputWave[n][9]
					CurrentUTime[floor(targetpoint)]+=InputWave[n][9]

				elseif(floor(pointentered)-1==floor(targetpoint))
					NewUnFoldingTimes[floor(pointentered)]+=InputWave[n][9]
					NewUnFoldingTimes[floor(targetpoint)]+=InputWave[n][10]
					CurrentUTime[floor(pointentered)]+=InputWave[n][9]
					CurrentUTime[floor(targetpoint)]+=InputWave[n][10]

				else
					NewUnFoldingTimes[floor(pointentered)]+=InputWave[n][9]
					NewUnFoldingTimes[floor(targetpoint)+1,floor(pointentered)-1]+=PrevunFoldedTimesIn[p]-PrevunFoldedTimesIn[p+1]
					NewUnFoldingTimes[floor(targetpoint)]+=InputWave[n][11]
					CurrentUTime[floor(pointentered)]+=InputWave[n][9]
					CurrentUTime[floor(targetpoint)+1,floor(pointentered)-1]+=PrevunFoldedTimesIn[p]-PrevunFoldedTimesIn[p+1]
					CurrentUTime[floor(targetpoint)]+=InputWave[n][11]
				endif
		
			endif
			wave FreeUTimeMaster
			if(waveexists(FreeUTimeMaster)==0)
				duplicate/o CurrentUTime FreeUTimeMaster
				duplicate/o CurrentUHist FreeUHistMaster
			else
				Concatenate/o {FreeUTimeMaster,CurrentUTime}, HoldUtime
				duplicate/o HoldUtime FreeUTimeMaster
				Concatenate/o {FreeUHistMaster,CurrentUHist}, HoldUHist
				duplicate/o HoldUHist FreeUHistMaster
			endif
		elseif(InputWave[n][0]==1&&InputWave[n][1]==1)//folding OutBound
			NewFoldingHist[floor(targetpoint)]+=1
			CurrentFHist[floor(targetpoint)]+=1

			if(floor(pointentered)==floor(targetpoint))
				NewFoldingTimes[floor(targetpoint)]+=InputWave[n][4]
				CurrentFTime[floor(targetpoint)]+=InputWave[n][4]

			elseif(floor(pointentered)+1==floor(targetpoint))
				NewFoldingTimes[floor(pointentered)]+=InputWave[n][4]
				NewFoldingTimes[floor(targetpoint)]+=InputWave[n][5]
				CurrentFTime[floor(pointentered)]+=InputWave[n][4]
				CurrentFTime[floor(targetpoint)]+=InputWave[n][5]

			else
				NewFoldingTimes[floor(pointentered)]+=InputWave[n][4]
				NewFoldingTimes[floor(pointentered)+1,floor(targetpoint)-1]+=FoldedTimesOut[p+1]-FoldedTimesOut[p]
				NewFoldingTimes[floor(targetpoint)]+=InputWave[n][6]
				CurrentFTime[floor(pointentered)]+=InputWave[n][4]
				CurrentFTime[floor(pointentered)+1,floor(targetpoint)-1]+=FoldedTimesOut[p+1]-FoldedTimesOut[p]
				CurrentFTime[floor(targetpoint)]+=InputWave[n][6]
			endif
			wave FreeFTimeMaster
			if(waveexists(FreeFTimeMaster)==0)
				duplicate/o CurrentFTime FreeFTimeMaster
				duplicate/o CurrentFHist FreeFHistMaster

			else
				Concatenate/o {FreeFTimeMaster,CurrentFTime}, HoldFtime
				duplicate/o HoldFtime FreeFTimeMaster
			
				Concatenate/o {FreeFHistMaster,CurrentFHist}, HoldFHist
				duplicate/o HoldFHist FreeFHistMaster
			endif
		
		elseif(InputWave[n][0]==-1&&InputWave[n][1]==-1)//unfolding InBound
			NewunFoldingHist[floor(targetpoint)]+=1
			CurrentUHist[floor(targetpoint)]+=1

			if(floor(pointentered)==floor(targetpoint))
				NewUnFoldingTimes[floor(targetpoint)]+=InputWave[n][4]
				CurrentUTime[floor(targetpoint)]+=InputWave[n][4]

			elseif(floor(pointentered)-1==floor(targetpoint))
				NewUnFoldingTimes[floor(pointentered)]+=InputWave[n][4]
				NewUnFoldingTimes[floor(targetpoint)]+=InputWave[n][5]
				CurrentUTime[floor(pointentered)]+=InputWave[n][4]
				CurrentUTime[floor(targetpoint)]+=InputWave[n][5]
			else
				NewUnFoldingTimes[floor(pointentered)]+=InputWave[n][4]
				NewUnFoldingTimes[floor(targetpoint)+1,floor(pointentered)-1]+=unFoldedTimesIn[p]-unFoldedTimesIn[p+1]
				NewUnFoldingTimes[floor(targetpoint)]+=InputWave[n][6]
				CurrentUTime[floor(pointentered)]+=InputWave[n][4]
				CurrentUTime[floor(targetpoint)+1,floor(pointentered)-1]+=unFoldedTimesIn[p]-unFoldedTimesIn[p+1]
				CurrentUTime[floor(targetpoint)]+=InputWave[n][6]
			endif
			wave FreeUTimeMaster
			if(waveexists(FreeUTimeMaster)==0)
				duplicate/o CurrentUTime FreeUTimeMaster
				duplicate/o CurrentUHist FreeUHistMaster
			else
				Concatenate/o {FreeUTimeMaster,CurrentUTime}, HoldUtime
				duplicate/o HoldUtime FreeUTimeMaster
				Concatenate/o {FreeUHistMaster,CurrentUHist}, HoldUHist
				duplicate/o HoldUHist FreeUHistMaster
			endif
		
		elseif(InputWave[n][0]==1&&InputWave[n][1]==-1)//folding InBound
			NewFoldingHist[floor(targetpoint)]+=1
			CurrentFHist[floor(targetpoint)]+=1

			if(floor(pointentered)==floor(targetpoint))
				NewFoldingTimes[floor(targetpoint)]+=InputWave[n][4]
				CurrentFTime[floor(targetpoint)]+=InputWave[n][4]

			elseif(floor(pointentered)-1==floor(targetpoint))
				NewFoldingTimes[floor(pointentered)]+=InputWave[n][4]
				NewFoldingTimes[floor(targetpoint)]+=InputWave[n][5]
				CurrentFTime[floor(pointentered)]+=InputWave[n][4]
				CurrentFTime[floor(targetpoint)]+=InputWave[n][5]

			else
				NewFoldingTimes[floor(pointentered)]+=InputWave[n][4]
				NewFoldingTimes[floor(targetpoint)+1,floor(pointentered)-1]+=FoldedTimesIn[p]-FoldedTimesIn[p+1]
				NewFoldingTimes[floor(targetpoint)]+=InputWave[n][6]
				CurrentFTime[floor(pointentered)]+=InputWave[n][4]
				CurrentFTime[floor(targetpoint)+1,floor(pointentered)-1]+=FoldedTimesIn[p]-FoldedTimesIn[p+1]
				CurrentFTime[floor(targetpoint)]+=InputWave[n][6]

			endif
			
			if(InputWave[n][12]==-1)
				targetpoint=InputWave[n][8]
				pointentered=InputWave[n][7]
				if(floor(pointentered)==floor(targetpoint))
					NewFoldingTimes[floor(targetpoint)]+=InputWave[n][9]
					CurrentFTime[floor(targetpoint)]+=InputWave[n][9]
				elseif(floor(pointentered)+1==floor(targetpoint))
					NewFoldingTimes[floor(pointentered)]+=InputWave[n][9]
					NewFoldingTimes[floor(targetpoint)]+=InputWave[n][10]
					CurrentFTime[floor(pointentered)]+=InputWave[n][9]
					CurrentFTime[floor(targetpoint)]+=InputWave[n][10]
				else
				
					NewFoldingTimes[floor(pointentered)]+=InputWave[n][9]
					NewFoldingTimes[floor(pointentered)+1,floor(targetpoint)-1]+=FoldedTimesOut[p+1]-FoldedTimesOut[p]
					NewFoldingTimes[floor(targetpoint)]+=InputWave[n][11]
					CurrentFTime[floor(pointentered)]+=InputWave[n][9]
					CurrentFTime[floor(pointentered)+1,floor(targetpoint)-1]+=FoldedTimesOut[p+1]-FoldedTimesOut[p]
					CurrentFTime[floor(targetpoint)]+=InputWave[n][11]
				endif
			endif
			wave FreeFTimeMaster
			if(waveexists(FreeFTimeMaster)==0)
				duplicate/o CurrentFTime FreeFTimeMaster
				duplicate/o CurrentFHist FreeFHistMaster

			else
				Concatenate/o {FreeFTimeMaster,CurrentFTime}, HoldFtime
				duplicate/o HoldFtime FreeFTimeMaster
			
				Concatenate/o {FreeFHistMaster,CurrentFHist}, HoldFHist
				duplicate/o HoldFHist FreeFHistMaster
			endif
		else
		endif
		
	endfor
	variable timeentered=InputWave[top-1][19]
	FindLevel/Q unFoldedTimesIn timeentered
	pointentered=v_levelx
	variable targettime=unFoldedTimesIn[0]
	targetpoint=0
	
	if(floor(pointentered)==floor(targetpoint))
		NewUnFoldingTimes[floor(targetpoint)]+=targettime-timeentered
		FreeFTimeMaster[floor(targetpoint)][dimsize(FreeFTimeMaster,1)-1]+=targettime-timeentered

	elseif(floor(pointentered)-1==floor(targetpoint))
		NewUnFoldingTimes[floor(pointentered)]+=unFoldedTimesIn[floor(pointentered)]-timeentered
		NewUnFoldingTimes[floor(targetpoint)]+=targettime-unFoldedTimesIn[floor(targetpoint)+1]
		
		FreeFTimeMaster[floor(pointentered)][dimsize(FreeFTimeMaster,1)-1]+=unFoldedTimesIn[floor(pointentered)]-timeentered
		FreeFTimeMaster[floor(targetpoint)][dimsize(FreeFTimeMaster,1)-1]+=targettime-unFoldedTimesIn[floor(targetpoint)+1]
	else
		NewUnFoldingTimes[floor(pointentered)]+=unFoldedTimesIn[floor(pointentered)]-timeentered
		NewUnFoldingTimes[floor(targetpoint)+1,floor(pointentered)-1]+=unFoldedTimesIn[p]-unFoldedTimesIn[p+1]
		NewUnFoldingTimes[floor(targetpoint)]+=targettime-unFoldedTimesIn[floor(targetpoint)+1]
		FreeFTimeMaster[floor(pointentered)][dimsize(FreeFTimeMaster,1)-1]+=unFoldedTimesIn[floor(pointentered)]-timeentered
		FreeFTimeMaster[floor(targetpoint)+1,floor(pointentered)-1][dimsize(FreeFTimeMaster,1)-1]+=unFoldedTimesIn[p]-unFoldedTimesIn[p+1]
		FreeFTimeMaster[floor(targetpoint)][dimsize(FreeFTimeMaster,1)-1]+=targettime-unFoldedTimesIn[floor(targetpoint)+1]
	endif
	NewUnFoldingRate=NewUnFoldingHist/NewUnfoldingTimes
	NewFoldingRate=NewFoldingHist/NewfoldingTimes

	SetScale/P x pnt2x(OrigFoldHist,0),dimdelta(OrigFoldHist,0),"", FreeFTimeMaster,FreeFHistMaster
	SetScale/P x pnt2x(OrigunFoldHist,0),dimdelta(OrigunFoldHist,0),"", FreeUTimeMaster,FreeUHISTMaster
	Duplicate/o FreeFTimeMaster $(FolderString+":"+BaseNameString+"_MasterFDwell")
	Duplicate/o FreeFHistMaster $(FolderString+":"+BaseNameString+"_MasterFHist")
	Duplicate/o FreeUTimeMaster $(FolderString+":"+BaseNameString+"_MasterUDwell")
	Duplicate/o FreeUHISTMaster $(FolderString+":"+BaseNameString+"_MasterUHist")
	killwaves FreeFTimeMaster,FreeUTimeMaster,FreeUHistMaster,FreeFHistMaster,CurrentFHist,CurrentFTime,CurrentUTime,CurrentUHist,HoldFHist,HoldFTime,Holduhist,HoldUtime
	killwaves NewFoldingHist,NewFoldingRate,NewFoldingTimes,NewunFoldingRate,NewunFoldingTimes,NewunFoldingHist
end

Static Function ResampleTheTransWaves(FolderString,BaseNameString,[OutPutString])
	string FolderString,BaseNameString,OutPutString
//	wave UTimeMaster,FTimeMaster,UHistMaster,FHistMaster,NewUTimeMaster,NewFTimeMaster,NewUHistMaster,NewFHistMaster
	wave UTimeMaster=$(FolderString+":"+BaseNameString+"_MasterUDwell")
	wave FTimeMaster=$(FolderString+":"+BaseNameString+"_MasterFDwell")
	wave UHistMaster=$(FolderString+":"+BaseNameString+"_MasterUHist")
	wave FHistMaster=$(FolderString+":"+BaseNameString+"_MasterFHist")
		If(paramisdefault(OutputString))
		OutputString=BaseNameString+"_Resample"
	
	endif
	make/free/n=(dimsize(UTimemaster,1)) UNums
	UNums=p
	make/free/n=(dimsize(UTimemaster,1)) FNums
	FNums=p
	
	StatsResample/N=(dimsize(UTimemaster,1)) UNums
	wave 	W_Resampled
	duplicate/free W_Resampled UNumsResampled
	duplicate/free UTimeMaster ResUTimeMaster
	duplicate/free UHistMaster ResUHistMaster
	ResUTimeMaster[][]=UTimeMaster[p][UNumsResampled[q]]
	ResUHistMaster[][]=UHistMaster[p][UNumsResampled[q]]
	
	StatsResample/N=(dimsize(UTimemaster,1)) FNums
	wave 	W_Resampled
	duplicate/free W_Resampled FNumsResampled
	duplicate/free FTimeMaster ResFTimeMaster
	duplicate/free FHistMaster ResFHistMaster
	ResFTimeMaster[][]=FTimeMaster[p][FNumsResampled[q]]
	ResFHistMaster[][]=FHistMaster[p][FNumsResampled[q]]
	
	duplicate/o ResFTimeMaster $(FolderString+":"+OutputString+"_MasterFDwell")
	duplicate/o ResFHistMaster $(FolderString+":"+OutputString+"_MasterFHist")

	duplicate/o ResUTimeMaster $(FolderString+":"+OutputString+"_MasterUDwell")
	duplicate/o ResUHistMaster $(FolderString+":"+OutputString+"_MasterUHist")
	killwaves W_Resampled

end

Static Function IterativeResampling(FolderString,BaseNameString,iterations)
variable iterations
	String FolderString,BaseNameString

	variable n=0
	for(n=0;n<iterations;n+=1)
		ResampleTheTransWaves(FolderString,BaseNameString,outPutString="Resample")
	ProcessingTransitionWave(FolderString,"Resample",OutPutString="Result")
	wave FRate=$(FolderString+":Result_FRate")
	wave URate=$(FolderString+":ReSult_URate")
	RemoveNan(FRate)
		RemoveNan(URate)

	if(n==0)
		duplicate/o FRate TotFRate
		duplicate/o URate TotURate

	else
		Concatenate/o {TotFRate,Frate}, HoldFRate
		Concatenate/o {TotURate,Urate}, HoldURate
		duplicate/o HoldFRate TotFRate
		duplicate/o HoldURate TotURate

	endif
	endfor
	SetScale/P x pnt2x(FRate,0),dimdelta(FRate,0), TotFRate
	SetScale/P x pnt2x(URate,0),dimdelta(URate,0), TotURate
	duplicate/o FRate, FMean,Fstd
	duplicate/o URate, UMean,Ustd
	make/free/n=(dimsize(TotFRate,1)) Fortho
	make/free/n=(dimsize(TotURate,1)) Uortho
	for(n=0;n<(dimsize(FRate,0));n+=1)
	
		Fortho=TotFRate[n][p]
		wavestats/Q Fortho
		FMean[n]=v_avg
		Fstd[n]=v_sdev
	endfor
	for(n=0;n<(dimsize(URate,0));n+=1)
	
		Uortho=TotURate[n][p]
		wavestats/Q Uortho
		UMean[n]=v_avg
		Ustd[n]=v_sdev
	endfor
	RemoveZeros(Ustd)
	RemoveZeros(Fstd)
		RemoveZeros(UMean)
	RemoveZeros(FMean)
	duplicate/o Ustd $(FolderString+":"+BaseNameString+"_UErr")
		duplicate/o Fstd $(FolderString+":"+BaseNameString+"_FErr")

	killwaves Umean,FMEan,Fstd,Ustd
	killwaves $(FolderString+":ReSult_UHist"),$(FolderString+":Result_UDwell"),$(FolderString+":Result_FDwell"),$(FolderString+":Result_FHist")
	killwaves $(FolderString+":Resample_MasterUDwell"),$(FolderString+":Resample_MasterFDwell"),$(FolderString+":Resample_MasterUHist"),$(FolderString+":Resample_MasterFHist")
	killwaves FRate,URate,HoldFRate,HoldURate,TotURate,TotFRate
end
	
Static Function RemoveNan(WaveIn)
	wave WaveIn
	
	
	variable n=0
	for(n=0;n<numpnts(WaveIn);n+=1)
		If(numtype(WaveIn[n])==2)
			WaveIn[n]=0
		endif
	
	
	endfor
end	
Static Function RemoveZeros(WaveIn)
	wave WaveIn
	
	
	variable n=0
	for(n=0;n<numpnts(WaveIn);n+=1)
		If((WaveIn[n])==0)
			WaveIn[n]=Nan
		endif
	
	
	endfor
end	
	
Static Function ProcessingTransitionWave(FolderString,BaseNameString,[OutputString])
	String FolderString,BaseNameString,OutputString
	wave UTimeMaster=$(FolderString+":"+BaseNameString+"_MasterUDwell")
	wave FTimeMaster=$(FolderString+":"+BaseNameString+"_MasterFDwell")
	wave UHistMaster=$(FolderString+":"+BaseNameString+"_MasterUHist")
	wave FHistMaster=$(FolderString+":"+BaseNameString+"_MasterFHist")
	If(paramisdefault(OutputString))
		OutputString=BaseNameString+"_Redone"
	
	endif
	make/free/n=(dimsize(FTimeMaster,0)) FinalFrate,FinalFTime,FinalFHist
	make/free/n=(dimsize(UTimeMaster,0)) FinalUrate,FinalUTime,FinalUHist
	FinalFrate=0
	FinalFTime=0
	FinalFHist=0
	FinalUrate=0
	FinalUTime=0
	FinalUHist=0
	variable FoldingNum=dimsize(FTimeMaster,1)
	variable UnFoldingNum=dimsize(FTimeMaster,1)

	variable n
	for(n=0;n<FoldingNum;n+=1)
	FinalFTime+=(FTimeMaster[p][n])
		FinalFHist+=(FHistMaster[p][n])

	endfor
	
		for(n=0;n<UnFoldingNum;n+=1)
	FinalUTime+=(UTimeMaster[p][n])
		FinalUHist+=(UHistMaster[p][n])

	endfor
	FinalFrate=FinalFHist/FinalFTime
	FinalUrate=FinalUHist/FinalUTime
	SetScale/P x pnt2x(FTimeMaster,0),dimdelta(FTimeMaster,0),"", FinalFrate,FinalFTime,FinalFHist
	SetScale/P x pnt2x(UTimeMaster,0),dimdelta(UTimeMaster,0),"", FinalUrate,FinalUHist,FinalUTime
	duplicate/o FinalFrate $(FolderString+":"+OutputString+"_FRate")
	duplicate/o FinalFHist $(FolderString+":"+OutputString+"_FHist")
	duplicate/o FinalFTime $(FolderString+":"+OutputString+"_FDwell")
	duplicate/o FinalUrate $(FolderString+":"+OutputString+"_URate")
	duplicate/o FinalUHist $(FolderString+":"+OutputString+"_UHist")
	duplicate/o FinalUTime $(FolderString+":"+OutputString+"_UDwell")

end

Static Function/S ListWaves(FolderString,SearchString)
	string SearchString,FolderString
	String saveDF=GetDataFolder(1)
	
	SetDataFolder FolderString
	String list = WaveList(SearchString, ";", "")
	SetDataFolder saveDF
	return list
end

Static Function CorrectWLC(WLCIn,WLCOut,ForceWave,SepWave,Shifting)
	wave WLCIn,WLCOut,ForceWave,SepWave
	String Shifting
	duplicate/free WLCIn WLCFree
	WLCFree[0]= WLCIn[1]
		WLCFree[1]= WLCIn[0]
		WLCFree[2]= WLCIn[2]
		WLCFree[3]= WLCIn[3]
		variable/c Shift= CorrectShift(ForceWave,Shifting)


	WLCFree[3]	+=imag(shift)
	WLCFree[2]	+=real(shift)
//	make/o/n=1000 TS
//	TS=wavemin(SepWave)+(Wavemax(SepWave)-Wavemin(SepWave))/999*x
//	duplicate/o TS ForceOut1,ForceOut2
//	ForceOut1=WLC(TS-WLCFree[3],.4e-9,WLCFree[0],298)-WLCFree[2]
//	ForceOut2=WLC(TS-WLCFree[3],.4e-9,WLCFree[1],298)-WLCFree[2]
	duplicate/o WLCFree WLCOut

end

Static Function PLotEverythingwithARate()

	string FoldRatesList=wavelist("*_FoldRate",";","")
	string UnFoldRatesList=wavelist("*UnFoldRate",";","")
	string FoldedString,BaseName,unFoldedString
	variable n,loc
	DoWindow RatesComparison
	if(V_Flag==1)
	
	Killwindow RatesComparison
	endif
	Display/N=RatesComparison
	for(n=0;n<itemsinlist(FoldRatesList);n+=1)
		FoldedString=stringfromlist(n,FoldRatesList)

		BaseName=Stringfromlist(0,FoldedString,"_")
		loc= WhichListItem(BaseName+"_UnfoldRate",UnFoldRatesList)
		unFoldedString=stringfromlist(loc,UnFoldRatesList)
		appendtograph/W=RatesComparison $FoldedString
		appendtograph/W=RatesComparison $unFoldedString
		ModifyGraph/W=RatesComparison mode($FoldedString)=3,marker($FoldedString)=19,useMrkStrokeRGB($FoldedString)=1
		ModifyGraph/W=RatesComparison mode($unFoldedString)=3,marker($unFoldedString)=23,useMrkStrokeRGB($unFoldedString)=1
		ModifyGraph/W=RatesComparison rgb($FoldedString)=(ColorList(n+1,0),ColorList(n+1,1),ColorList(n+1,2))
		ModifyGraph/W=RatesComparison rgb($unFoldedString)=(ColorList(n+1,0),ColorList(n+1,1),ColorList(n+1,2))

		
	endfor
	ModifyGraph/W=RatesComparison log(left)=1
end

Static Function ColorList(index,column)

	variable index,column
	make/free/n=(3,7) ColorWave
	ColorWave[][1]={58596,6682,7196}
	ColorWave[][2]={14906,32382,47288}
	ColorWave[][3]={39064,20046,41891}
	ColorWave[][4]={19789,44975,19018}
	ColorWave[][5]={65535,32639,0}
	ColorWave[][6]={26214,26212,0}
	ColorWave[][7]={65535,16385,55749}
	if(index>7)
	return 0
	endif
	return ColorWave[column][index]



end

Static Function PlotForceandFit(FolderString,Type,Filtering)

	variable Filtering
	string FolderString,Type
	wave ForceWave=$(FolderString+":"+stringfromlist(0,DE_CountRates#ListWaves(FolderString,"*Force_Align")))
	wave SepWave=$(FolderString+":"+stringfromlist(0,DE_CountRates#ListWaves(FolderString,"*Sep_Align")))
	wave ForceWaveSm=$(FolderString+":"+stringfromlist(0,DE_CountRates#ListWaves(FolderString,"*Force_Align"))+"_Sm")
	wave SepWaveSm=$(FolderString+":"+stringfromlist(0,DE_CountRates#ListWaves(FolderString,"*Sep_Align"))+"_Sm")
	wave WLCWave=$(FolderString+":"+stringfromlist(0,DE_CountRates#ListWaves(FolderString,"*WLCParms_Adj")))
	
	make/o/n=1000 $(FolderString+":WLCSeparation")/Wave=WLCSep
	make/o/n=1000 $(FolderString+":WLCForce1")/Wave=WLCForce1
	make/o/n=1000 $(FolderString+":WLCForce2")/Wave=WLCForce2
	duplicate/o ForceWave $(FolderString+":FSm")
	Wave FSm=$(FolderString+":FSm")
	duplicate/o SepWave $(FolderString+":SSm")
	Wave SSm=$(FolderString+":SSm")	
	variable/c Shifts=CorrectShift(ForceWave,Type)
	if(Filtering>5)
	DE_Filtering#FilterForceSep(ForceWave,SepWave,FSm,SSm,"SVG",filtering)
	elseif(filtering<1)
	DE_Filtering#FilterForceSep(ForceWave,SepWave,FSm,SSm,"Tvd",filtering)

	else
	print "BadFilering"
	return 0
	endif
	
	
	
	SSm+=imag(Shifts)
	Fsm+=real(Shifts)
	make/free/n=0 WLCFree
	CorrectWLC(WLCWave,WLCFree,ForceWave,SepWave,Type)


	WLCSep=wavemin(SepWave)+(Wavemax(SepWave)-Wavemin(SepWave))/999*x
	WLCForce1=WLC(WLCSep-WLCFree[3],.4e-9,WLCFree[0],298)-WLCFree[2]
	WLCForce2=WLC(WLCSep-WLCFree[3],.4e-9,WLCFree[1],298)-WLCFree[2]
	dowindow WLCQuality
	
	if(V_flag==1)
	killwindow WLCQuality
	endif
	display/N=WLCQuality
	Appendtograph/W=WLCQuality FSm vs SSm

	Appendtograph/W=WLCQuality WLCForce1 vs WLCSep
	Appendtograph/W=WLCQuality WLCForce2 vs WLCSep
	ModifyGraph/W=WLCQuality muloffset($nameofwave($nameofwave(FSM)))={0,-1}
	ModifyGraph/W=WLCQuality rgb($nameofwave(WLCForce1))=(0,0,0),rgb($nameofwave(WLCForce2))=(0,0,0)
	print WLCFree[0]-WLCFree[1]
	//print WLCFree
end

Function QuickFit(SepWave,WLCParms)
	wave SepWave,WLCParms

	make/o/n=1000 $("WLCSeparation")/Wave=WLCSep
	make/o/n=1000 $(":WLCForce1")/Wave=WLCForce1
	make/o/n=1000 $(":WLCForce2")/Wave=WLCForce2
	
	
//	duplicate/o ForceWave $("MyFSm")
//	Wave FSm=$("MyFSm")
//	duplicate/o SepWave $("MySSm")
//	Wave SSm=$("MySSm")	
//	variable/c Shifts=CorrectShift(ForceWave,Type)
//	if(Filtering>5)
//	DE_Filtering#FilterForceSep(ForceWave,SepWave,FSm,SSm,"SVG",filtering)
//	elseif(filtering<1)
//	DE_Filtering#FilterForceSep(ForceWave,SepWave,FSm,SSm,"Tvd",filtering)
//
//	else
//	print "BadFilering"
//	return 0
//	endif
	
	
	
	//SSm-=imag(Shifts)
	//Fsm-=real(Shifts)
	//make/free/n=0 WLCFree
	//CorrectWLC(WLCParms,WLCFree,ForceWave,SepWave,Type)
	duplicate/free WLCParms WLCFRee
	WLCSep=wavemin(SepWave)+(Wavemax(SepWave)-Wavemin(SepWave))/999*x
	WLCForce1=WLC(WLCSep-WLCFree[3],.4e-9,WLCFree[0],298)-WLCFree[2]
	WLCForce2=WLC(WLCSep-WLCFree[3],.4e-9,WLCFree[1],298)-WLCFree[2]

end


Static Function CombineAllRates()


	string FoldHistList=wavelist("*_FoldHist",";","")

	string FoldedHistStr
	variable n,loc
	//	DoWindow RatesComparison
	//	if(V_Flag==1)
	
	//	Killwindow RatesComparison
	//	endif
	//	Display/N=RatesComparisond
	make/o/n=(numpnts($stringfromlist(0,FoldHistList))) Combined_FoldHist,Combined_FoldTimes,Combined_FoldRates
	SetScale/P x pnt2x($stringfromlist(0,FoldHistList),0),dimdelta($stringfromlist(0,FoldHistList),0),"", Combined_FoldHist,Combined_FoldTimes,Combined_FoldRates
	
	make/o/n=(numpnts($replacestring("FoldHist",stringfromlist(0,FoldHistList),"UnfoldHist"))) Combined_UnFoldHist,Combined_UnFoldTimes,Combined_unFoldRates
	SetScale/P x pnt2x($replacestring("FoldHist",stringfromlist(0,FoldHistList),"UnfoldHist"),0),dimdelta($replacestring("FoldHist",stringfromlist(0,FoldHistList),"UnfoldHist"),0),"", Combined_unFoldHist,Combined_unFoldTimes,Combined_unFoldRates

			Combined_FoldHist=0
		Combined_FoldTimes=0
		Combined_FoldRates=0
		Combined_unFoldHist=0
		Combined_unFoldTimes=0
		Combined_unFoldRates=0
	
	for(n=0;n<itemsinlist(FoldHistList);n+=1)
		FoldedHistStr=stringfromlist(n,FoldHistList)
		wave FoldedHist=$FoldedHistStr
		wave UnfoldedHist=$replacestring("FoldHist",nameofwave(FoldedHist),"UnfoldHist")
		wave foldedTimes=$replacestring("FoldHist",nameofwave(FoldedHist),"FoldTimes")
		wave UnfoldedTimes=$replacestring("FoldHist",nameofwave(FoldedHist),"UnFoldTimes")
		Combined_FoldHist+=FoldedHist
		Combined_FoldTimes+=foldedTimes
		Combined_FoldRates=Combined_FoldHist/Combined_FoldTimes
		Combined_unFoldHist+=UnfoldedHist
		Combined_unFoldTimes+=UnfoldedTimes
		Combined_unFoldRates=Combined_unFoldHist/	Combined_unFoldTimes
		
	endfor



end