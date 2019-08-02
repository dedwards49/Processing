#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma modulename=DE_OverlapRamps

Static Function ExtractForceCurveWithGaps(ForceIn,SepIn,StateWave,num,CutWave,ForceOut,SepOut)
	Wave Forcein,ForceOut,StateWave,SepIn,SepOut,CutWave
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

	duplicate/free/r=[startRet,EndExt] ForceIn FO
	duplicate/free/r=[startRet,EndExt] SepIn SO
	SetScale/P x 0,dimdelta(ForceIn,0),"s", FO
	SetScale/P x 0,dimdelta(ForceIn,0),"s", SO

	variable start1= CutWave[0]
	variable end1= CutWave[1]
	variable start2= CutWave[2]
	variable end2= CutWave[3]

	deletepoints start2,(end2-start2), FO
	deletepoints start1,(end1-start1), FO
	deletepoints start2,(end2-start2), SO
	deletepoints start1,(end1-start1), SO
	duplicate/o FO ForceOut
	duplicate/o SO SepOut

end

Static Function ReturnGaps(ForceIn,SepIn,StateWave,num,Results)
	Wave Forcein,SepIn,StateWave,Results
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
	
	duplicate/free/r=[startRet,EndExt] ForceIn FO
	duplicate/free/r=[startRet,EndExt] SepIn SO
	make/free/n=0 ReturnPoints
	ExtractPointsFromThisState(StateWave,num,ReturnPoints)
	FindLevel/Q ReturnPoints, Mid
	variable start1= ReturnPoints[1]-250
	start1-=startRet
	variable end1= ReturnPoints[v_levelx-1]+250
	end1-=startRet
	
variable start2
	if(V_levelx-1==numpnts(ReturnPoints))
		start2= ReturnPoints[v_levelx+1]	-250

	else
	start2= ReturnPoints[v_levelx+1]	-250
	endif
	start2-=startRet

	variable end2= ReturnPoints[numpnts(returnpoints)-1]+250
	end2-=startRet

	make/free/n=4 Garbage
	Garbage={start1,end1,start2,end2}
	duplicate/o Garbage Results

end

Static Function ReturnForcesWidestGaps(ForceWave1,SepWave1,StateWave1,num1,num2,OutForce1,OutSep1,OutForce2,OutSep2)
	wave ForceWave1,SepWave1,StateWave1,OutForce1,OutSep1,OutForce2,OutSep2
	variable num1,num2
	make/free/N=0 Points1,Points2,RF1,RF2,RS1,RS2
	make/free/n=4 FinalPoints
	ReturnGaps(ForceWave1,SepWave1,StateWave1,num1,Points1)
	ReturnGaps(ForceWave1,SepWave1,StateWave1,num2,Points2)
	FinalPoints[0]=min(Points1[0],Points2[0])
	FinalPoints[1]=max(Points1[1],Points2[1])
	FinalPoints[2]=min(Points1[2],Points2[2])
	FinalPoints[3]=max(Points1[3],Points2[3])
	ExtractForceCurveWithGaps(ForceWave1,SepWave1,StateWave1,num1,FinalPoints,RF1,RS1)
	ExtractForceCurveWithGaps(ForceWave1,SepWave1,StateWave1,num2,FinalPoints,RF2,RS2)
	duplicate/o RF1 OutForce1
	duplicate/o RS1 OutSep1
	duplicate/o RF2 OutForce2
	duplicate/o RS2 OutSep2

end



Static Function CalculateMaxIndexforForceWave(ForceWave)
	wave ForceWave
	
	
	string FullInd=Stringbykey("DE_Ind",note(ForceWave),":","\r")
	variable numberofpulls=itemsinlist(FullINd,";")/2
	return numberofpulls
end



Static Function CalculateForceOffsetIndex(ForceWave,SepWave,Statewave,index1,index2)
	wave ForceWave,SepWave,Statewave
	variable index1,index2


	
	make/free/n=0 F1,F2,S1,S2
	ReturnForcesWidestGaps(ForceWave,SepWave,Statewave,index1,index2,F1,S1,F2,S2)
	make/o/n=2 W_Coef
	W_Coef={0,1}
	variable lengthdif=numpnts(F1)-numpnts(F2)
	print lengthdif
	if(lengthdif==0)
	elseif(lengthdif<=2)
		deletepoints numpnts(F1)-lengthdif,1e3,F1
		else
	endif
	CurveFit/Q/H="01"/NTHR=0 line  kwCWave=W_Coef F1 /X=F2 
	wave W_Coef,w_Sigma

	variable returns=W_Coef[0]
	killwaves W_Coef,w_Sigma
	return returns
end


Static Function ExtractPointsFromThisState(StateWave,Index,ReturnPoints)
	wave StateWave,ReturnPoints
	variable Index
	
	make/free/n=(dimsize(StateWave,0)) Points,Type,States
	Points=StateWave[p][0]
	Type=StateWave[p][2]
	States=StateWave[p][3]
	Extract/Free/Indx States, HiDevin, States==index
	make/free/n=(numpnts(HiDevin)) Results
	Results=Points[HiDevin]
	duplicate/o Results ReturnPoints

end

Static Function OffsetEachStepinForceWave(ForceIn,ForceOut,States)
	wave ForceIn,ForceOut,States
	
variable n,prevcut=0,currcut,FOFF
	duplicate/o ForceIn ForceOut
	for(n=1;n<dimsize(States,0);n+=1)
		if(States[n][2]==2)
			currcut=States[n][0]
			FOFF=ReturnOffForceForState(ForceIn,States[n][3]-1)
			ForceOut[prevcut,currcut]+=FOFF
			prevcut=States[n][0]
		endif
	
	endfor
end

Static Function/S AddForceOffsetstoForceWave(ForceWave,SepWave,StateWave)
	wave ForceWave,SepWave,StateWave
	variable numsteps=CalculateMaxIndexforForceWave(ForceWave)
	print numsteps
	variable n=0
	string Offsets="0;"
	variable RefOff=CalculateForceOffsetFromPause(ForceWave,0)
	for(n=1;n<numsteps;n+=1)
		//Offsets+=num2str(CalculateForceOffsetIndex(ForceWave,SepWave,StateWave,0,n))+";"
		Offsets+=num2str(CalculateForceOffsetFromPause(ForceWave,n)-RefOff)+";"
	endfor
	note/K ForceWave,ReplaceStringByKey("DE_FOff", note(ForceWave), Offsets,":","\r" )
	note/K SepWave,ReplaceStringByKey("DE_FOff", note(ForceWave), Offsets,":","\r" )
	return (Offsets)
end

Static Function ReturnOffForceForState(ForceWave,index)

		wave ForceWave
		variable index
		return str2num(stringfromlist(index,stringbykey("DE_FOff",note(ForceWave),":","\r")))
end

Static Function ReturnForcewithOffset(ForceWave,StateWave,Point)
	wave ForceWave,StateWave	
	variable Point
	
	make/free/n=(dimsize(StateWave,0)) Points
	Points=StateWave[p][0]
	FindLevel/Q Points Point
	variable index=StateWave[ceil(V_levelx)][3]
	return (ForceWave[Point]+ReturnOffForceForState(ForceWave,index))
end

Static Function AddShiftToStates(ForceWave,StateWave)
	wave ForceWave,StateWave
	make/free/n=(dimsize(StateWave,0),6) AboutToFinish
	AboutToFinish[][0]=StateWave[p][0]
	AboutToFinish[][1]=StateWave[p][1]
	AboutToFinish[][2]=StateWave[p][2]
	AboutToFinish[][3]=StateWave[p][3]
	AboutToFinish[][4]=StateWave[p][4]
	AboutToFinish[][5]=StateWave[p][1]+ReturnOffForceForState(ForceWave,StateWave[p][3])

	duplicate/o AbouttoFinish StateWave
end


Static Function CalculateForceOffsetFromPause(ForceWave,RampNumber)
	wave ForceWave
	variable RampNumber
	String PauseLocs=stringbykey("DE_PauseLoc",note(ForceWave),":","\r")
	variable startpnt=str2num(stringfromlist(2*RampNumber,PauseLocs))
	variable endpnt=str2num(stringfromlist(2*RampNumber+1,PauseLocs))
	variable diff=endpnt-startpnt
	duplicate/o/R=[startpnt+.1*diff,endpnt-.1*diff] ForceWave ForceWaveCut

	wavestats/Q ForceWaveCut
//	
//	make/free/n=0 F1,F2,S1,S2
//	ReturnForcesWidestGaps(ForceWave,SepWave,Statewave,index1,index2,F1,S1,F2,S2)
//	make/o/n=2 W_Coef
//	W_Coef={0,1}
//	variable lengthdif=numpnts(F1)-numpnts(F2)
//	print lengthdif
//	if(lengthdif==0)
//	elseif(lengthdif<=2)
//		deletepoints numpnts(F1)-lengthdif,1e3,F1
//		else
//	endif
//	CurveFit/Q/H="01"/NTHR=0 line  kwCWave=W_Coef F1 /X=F2 
//	wave W_Coef,w_Sigma
//
//	variable returns=W_Coef[0]
//	killwaves W_Coef,w_Sigma
	return v_avg
end