#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma modulename=DE_Transitionpath
#Include "DE_FiLtering"
#Include "Panel Progs"
Static Function OneClickTransitionPathTimes(InputWave,InitialSmoothing)

	wave InputWave
	variable InitialSmoothing
	make/o/n=0 CrudeTransitions,Lifetime1,Lifetime2
	make/o/n=0 SmoothedWave
	make/free/n=0 IndexMatrix
	DE_Filtering#TVD1D_denoise(InputWave,InitialSmoothing,SmoothedWave)
	Make/N=50/free Hist
	Histogram/C/B=1 SmoothedWave,Hist
	DE_Transitionpath#FindTransitionsSmoothed(SmoothedWave,pnt2x(Hist,24),CrudeTransitions)
	
	if(numpnts(CrudeTransitions)==0)
	return -1
	endif
	DE_TransitionPath#CalculateTwoStateLifeTimes(CrudeTransitions,Lifetime1,Lifetime2)
	DE_transitionpath#LogisticFitAllTransitions(InputWave,CrudeTransitions,IndexMatrix)
	make/free/n=0 WaveOut
	FindAllFinalTPs(InputWave,IndexMatrix,1,0,WaveOut)
	make/free/n=(dimsize(WaveOut,0)) TransitionTimes
	TransitionTimes=WaveOut[p][3]
	SaveOut(InputWave,IndexMatrix,WaveOut,TransitionTimes,Lifetime1,Lifetime2)
end

static function FindTransitions(INputWave,Smoothing,Level,OutWave)
	wave InputWave,Outwave
	variable smoothing,Level
	make/free/n=0 SmoothedWave
	DE_Filtering#TVD1D_denoise(InputWave,smoothing,SmoothedWave)
	FindLevels/Q  SmoothedWave Level
	wave W_FindLevels
	duplicate/o W_FindLevels OutWave
end

static function FindTransitionsSmoothed(INputWave,Level,OutWave)
	wave InputWave,Outwave
	variable Level

	FindLevels/Q  INputWave Level
	wave W_FindLevels
	duplicate/o W_FindLevels OutWave
end

Static Function FindAllFinalTPs(InputWave,IndexMatrix,smth,Select,WaveOut)

	wave InputWave,IndexMatrix,WaveOut
	variable smth,Select
	make/free/n=(dimsize(Indexmatrix,0)) IndexWave,LWaveA,LWave
	make/free/n=0 StartTime,EndTime,TransitionTime,StartLevel,EndLevel,CenterLevelWave,FilterSet
	IndexWave=IndexMatrix[p][0]
	LWave=(IndexMatrix[p][1])
	LWaveA=abs(LWave)
	variable n,TimeGuess,CenterGuess,Index,centerlevel,edge,Result,v1,v2,Keep,MedSm,loops,CenterValue
	variable Lavg=mean(LWaveA)	
	String ResString
	for(n=0;n<numpnts(IndexWave);n+=1)
		Index=n
		TimeGuess=Indexmatrix[n][2]
		CenterGuess=Indexmatrix[n][3]
		make/n=0/o Segment
		SliceTransitionbyIndex(INputWave,IndexWave,Index,max(10*TimeGuess,.2e-3),Segment)
		Duplicate/o Segment SegmentSM
		MedSm= round(smth*TimeGuess/2/dimdelta(Segment,0))
		//print "Filtering: "+num2str(round(smth*TimeGuess/2/dimdelta(Segment,0)))
		Smooth/M=0 round(smth*TimeGuess/2/dimdelta(Segment,0)), SegmentSM 
				centerlevel=Indexmatrix[n][7]

		//centerlevel=SegmentSM(CenterGuess)
		if(LWave[index]>0)//Up
			FindLevels/Q/edge=(1) /D= Lower SegmentSM, (centerlevel-LAvg/4)
			FindLevels/Q/edge=(1) /D= Upper SegmentSM, (centerlevel+LAvg/4)
			if(numpnts(Lower)==0||numpnts(Upper)==0)
				result=NaN
			else
				ResString=DE_TransitionPath#FindProductiveTransitions(Lower,Upper)
				sscanf (Stringfromlist(0,ResString,",")),"%10f",v1
				sscanf (Stringfromlist(1,ResString,",")),"%10f",v2
				Result=v2-v1
				duplicate/o/r=(v1,v2) SegmentSM CutSeg
			endif
		else //Down
			FindLevels/Q/edge=(2) /D= Lower SegmentSM, (centerlevel-LAvg/4)
			FindLevels/Q/edge=(2) /D= Upper SegmentSM, (centerlevel+LAvg/4)
			if(numpnts(Lower)==0||numpnts(Upper)==0)
				result=NaN
	

			else
				ResString=DE_TransitionPath#FindProductiveTransitions(Upper,Lower)
				sscanf (Stringfromlist(0,ResString,",")),"%10f",v1
				sscanf (Stringfromlist(1,ResString,",")),"%10f",v2
				Result=v2-v1
				duplicate/o/r=(v1,v2) SegmentSM CutSeg
			endif
		endif
				loops=0
		do
		//if(Result<MedSm)
			if(loops==0)
			else
			Duplicate/o Segment SegmentSM
			TimeGuess=Result
					MedSm= round(smth*TimeGuess/2/dimdelta(Segment,0))

			print "Revise Filtering: "+num2str(round(smth*TimeGuess/2/dimdelta(Segment,0)))
			Smooth/M=0 round(smth*Result/2/dimdelta(Segment,0)), SegmentSM 
			
			centerlevel=SegmentSM(CenterGuess)
			if(LWave[index]>0)//Up
				FindLevels/Q/edge=(1) /D= Lower SegmentSM, (centerlevel-LAvg/4)
				FindLevels/Q/edge=(1) /D= Upper SegmentSM, (centerlevel+LAvg/4)
				if(numpnts(Lower)==0||numpnts(Upper)==0)
					result=NaN
				else
					ResString=DE_TransitionPath#FindProductiveTransitions(Lower,Upper)
					sscanf (Stringfromlist(0,ResString,",")),"%10f",v1
					sscanf (Stringfromlist(1,ResString,",")),"%10f",v2
					Result=v2-v1
					duplicate/o/r=(v1,v2) SegmentSM CutSeg
				endif
			else //Down
				FindLevels/Q/edge=(2) /D= Lower SegmentSM, (centerlevel-LAvg/4)
				FindLevels/Q/edge=(2) /D= Upper SegmentSM, (centerlevel+LAvg/4)
				if(numpnts(Lower)==0||numpnts(Upper)==0)
					result=NaN
	

				else
					ResString=DE_TransitionPath#FindProductiveTransitions(Upper,Lower)
					sscanf (Stringfromlist(0,ResString,",")),"%10f",v1
					sscanf (Stringfromlist(1,ResString,",")),"%10f",v2
					Result=v2-v1
					duplicate/o/r=(v1,v2) SegmentSM CutSeg
				endif
			endif
		endif
		loops+=1
		while((Result<TimeGuess/2||Result>TimeGuess*2)&&loops<=5)
		
		if(numpnts(Lower)==0)
		else
			Duplicate/o LOWer LowerLV;
			LowerLv=SegmentSm(Lower)
		endif
		if(numpnts(Upper)==0)
		else
			Duplicate/o Upper UpperLV;
			UpperLv=SegmentSm(Upper)
		endif
		if((Result>100e-6||Result<10e-6||loops>5)&&Select==1)
			print Result
			print Index
			DoWindow DTE
			if( V_Flag==0 )
				display/N=DTE 	InputWave
				appendtograph Segment
				appendtograph SegmentSm
				appendtograph CutSeg
				appendtograph UpperLV vs Upper 
				appendtograph LowerLV vs Lower
				ModifyGraph rgb($nameofwave(InputWave))=(65280,48896,48896),rgb(SegmentSM)=(0,0,0);DelayUpdate
				ModifyGraph mode(UpperLV)=3,marker(UpperLV)=19,rgb(UpperLV)=(0,52224,26368);DelayUpdate
				ModifyGraph mode(LowerLV)=3,marker(LowerLV)=19,rgb(LowerLV)=(14848,32256,47104)
				ModifyGraph lsize(CutSeg)=2,rgb(CutSeg)=(29440,0,58880)
				SetAxis bottom (CEnterGuess-max(10*TimeGuess,.2e-3)/2),(CEnterGuess+max(10*TimeGuess,.2e-3)/2)
			else
				killwindow DTE
				display/N=DTE 	InputWave
				appendtograph Segment
				appendtograph SegmentSm
				appendtograph CutSeg
				appendtograph UpperLV vs Upper 
				appendtograph LowerLV vs Lower
				ModifyGraph rgb($nameofwave(InputWave))=(65280,48896,48896),rgb(SegmentSM)=(0,0,0);DelayUpdate
				ModifyGraph mode(UpperLV)=3,marker(UpperLV)=19,rgb(UpperLV)=(0,52224,26368);DelayUpdate
				ModifyGraph mode(LowerLV)=3,marker(LowerLV)=19,rgb(LowerLV)=(14848,32256,47104)
				ModifyGraph lsize(CutSeg)=2,rgb(CutSeg)=(29440,0,58880)
				SetAxis bottom (CEnterGuess-max(10*TimeGuess,.2e-3)/2),(CEnterGuess+max(10*TimeGuess,.2e-3)/2)
			endif
			Keep=UserCursorAdjust("DTE",0)
			
		else
			if(loops>5)
			Keep=1
			else
			keep=0
			endif
		endif
		if(Keep==0)
			insertpoints n,1, StartTime,EndTime,TransitionTime,StartLevel,EndLevel,CenterLevelWave,FilterSet
			StartTime[n]=v1
			EndTime[n]=v2
			StartLevel[n]=(centerlevel-LAvg/4)
			EndLevel[n]=(centerlevel+LAvg/4)
			CenterLevelWave[n]=centerlevel
			TransitionTime[n]=Result
			FilterSet[n]=MedSm
		else
			insertpoints n,1, StartTime,EndTime,TransitionTime,StartLevel,EndLevel,CenterLevelWave,FilterSet
			StartTime[n]=NaN
			EndTime[n]=NaN
			StartLevel[n]=NaN
			EndLevel[n]=NaN
			CenterLevelWave[n]=NaN
			TransitionTime[n]=NaN
			FilterSet[n]=NaN
		endif
	endfor
		

	make/free/n=(numpnts(StartTime),8) Combined
	Combined[][0]=IndexWave

	Combined[][1]=StartTime[p]
	Combined[][2]=EndTime[p]

	Combined[][3]=TransitionTime[p]
	Combined[][4]=StartLevel[p]
	Combined[][5]=EndLevel[p]
	Combined[][6]=CenterLevelWave[p]
		Combined[][7]=FilterSet[p]

	duplicate/o Combined WaveOut
	wave W_FINDLEVELS
	killwaves Segment,SegmentSm,CutSeg,Lower,Upper,LowerLV,UpperLV,W_FindLevels

end

Static Function SliceTransitionbyIndex(INputWave,IndexWave,Index,width,OUtputWave)
	wave INputWave,IndexWave,OUtputWave
	variable index,width

	if(index>numpnts(indexwave))
		print "NOOOOOO"
		return -1
	endif

	duplicate/o/r=(IndexWave[index]-width/2,IndexWave[index]+width/2) Inputwave Outputwave

end

Function Logistic(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) =y0+L/(1+exp(-k*(x-x0)))
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 4
	//CurveFitDialog/ w[0] = L
	//CurveFitDialog/ w[1] = k
	//CurveFitDialog/ w[2] = x0
	//CurveFitDialog/ w[3] = y0

	return w[3]+w[0]/(1+exp(-w[1]*(x-w[2])))
End

Static Function LogisticFitAllTransitions(InputWave,IndexWave,OutputLs)
	wave InputWave,IndexWave,OutputLs
	make/free/n=0 Ls,TransitionWidth,TransitionCenter,W_coef1,W_coef2,W_coef3,CenterValue
	variable n
	variable start,ends
	
	for(n=0;n<numpnts(IndexWave);n+=1)
		make/o/n=0 Segment
		SliceTransitionbyIndex(InputWave,IndexWave,n,1e-3,Segment)
		make/d/o/n=4 w_coef
		wavestats/Q/r=[0,2*numpnts(Segment)/10] Segment
		start=v_avg
		wavestats/Q/r=[8/10*numpnts(Segment),numpnts(Segment)-1] Segment
		ends=v_avg
		w_coef={ends-start,1e5,pnt2x(Segment,numpnts(Segment)/2),start}

		FuncFit/Q/q/w=2/NTHR=0 Logistic W_coef  Segment
		insertpoints n,1, Ls,TransitionWidth,TransitionCenter,CenterValue,W_coef1,W_coef2,W_coef3
		TransitionWidth[n]=(ln(3)-ln(1/3))/w_coef[1]
		TransitionCenter[n]=w_coef[2]-ln(1)/w_coef[1]
		Ls[n]=w_coef[0]
		W_coef1[n]=w_coef[1]
		W_coef2[n]=w_coef[2]
		W_coef3[n]=w_coef[3]
		CenterValue[n]=w_coef[3]+w_coef[0]/2
	endfor

	wave w_coef,W_sigma
	make/free/n=(numpnts(Ls),8) Combined
	Combined[][0]=IndexWave

	Combined[][1]=(Ls[p])
	Combined[][2]=TransitionWidth[p]

	Combined[][3]=TransitionCenter[p]
	Combined[][4]=W_coef1[p]
	Combined[][5]=W_coef2[p]
	Combined[][6]=W_coef3[p]
Combined[][7]=CenterValue[p]
	duplicate/o Combined OutputLs
	killwaves w_coef,W_sigma



	
end


Static Function/C RefitTransitionByIndex(InputWave,IndexMatrix,index,smth)

	wave InputWave,IndexMatrix
	variable index,smth
	make/free/n=(dimsize(Indexmatrix,0)) IndexWave,LWaveA,LWave
	IndexWave=IndexMatrix[p][0]
	LWave=(IndexMatrix[p][1])
	LWaveA=abs(LWave)
	variable TimeGuess=Indexmatrix[index][2]
	variable CenterGuess=Indexmatrix[index][3]
	variable Lavg=mean(LWaveA)	
	make/n=0/o Segment
	SliceTransitionbyIndex(INputWave,IndexWave,Index,max(10*TimeGuess,.2e-3),Segment)
	Duplicate/o Segment SegmentSM

	Smooth/M=0 round(smth*TimeGuess/2/dimdelta(Segment,0)), SegmentSM 
	variable centerlevel=SegmentSM(CenterGuess)

	variable edge
	variable Result,n	
	make/o/n=0 CutSeg
	string ResString
	variable v1,v2
	if(LWave[index]>0)//Up
		FindLevels/Q/edge=(1) /D= Lower SegmentSM, (centerlevel-LAvg/4)
		FindLevels/Q/edge=(1) /D= Upper SegmentSM, (centerlevel+LAvg/4)
		if(numpnts(Lower)==0||numpnts(Upper)==0)
			result=NaN
//		elseif(Upper[0]>Lower[numpnts(Lower)-1])
//			Result= Upper[0]-Lower[numpnts(Lower)-1]
//			duplicate/o/r=(Lower[numpnts(Lower)-1],Upper[0]) Segment CutSeg
//		elseif(Lower[0]>Upper[numpnts(Upper)-1])
//			Result= NaN
//		elseif(Upper[0]<Lower[0])
//		
//		
		else
			ResString=DE_TransitionPath#FindProductiveTransitions(Lower,Upper)
			sscanf (Stringfromlist(0,ResString,",")),"%10f",v1
			sscanf (Stringfromlist(1,ResString,",")),"%10f",v2
			Result=v2-v1

//			FindLevel/p Lower, Upper[0]
//			result=Upper[0]-Lower[floor(v_levelx)]
			duplicate/o/r=(v1,v2) Segment CutSeg
//
		endif
	else //Down
		FindLevels/Q/edge=(2) /D= Lower SegmentSM, (centerlevel-LAvg/4)
		FindLevels/Q/edge=(2) /D= Upper SegmentSM, (centerlevel+LAvg/4)
		if(numpnts(Lower)==0||numpnts(Upper)==0)
			result=NaN
	
//		elseif(Lower[0]>Upper[numpnts(Upper)-1])
//			Result= Lower[0]-Upper[numpnts(Upper)-1]
//			duplicate/o/r=(Upper[numpnts(Upper)-1],Lower[0]) Segment CutSeg
//
//		elseif(Upper[0]<Lower[numpnts(Lower)-1])
//			Result= NaN
		else
			ResString=DE_TransitionPath#FindProductiveTransitions(Upper,Lower)
			sscanf (Stringfromlist(0,ResString,",")),"%10f",v1
			sscanf (Stringfromlist(1,ResString,",")),"%10f",v2
			Result=v2-v1
//			FindLevel/Q Upper, Lower[0]
//			result=Lower[0]-Upper[floor(v_levelx)]
			duplicate/o/r=(v1,v2) Segment CutSeg
//
		endif
	endif

	if(numpnts(Lower)==0)
	else
		Duplicate/o LOWer LowerLV;
		LowerLv=SegmentSm(Lower)
	endif
	if(numpnts(Upper)==0)
	else
		Duplicate/o Upper UpperLV;
		UpperLv=SegmentSm(Upper)
	endif
	if(Result>1000e-6||Result<.1e-6)
	print Result
	print Index
		DoWindow DTE
	if( V_Flag==0 )
		display/N=DTE 	InputWave
		appendtograph Segment
		appendtograph SegmentSm
		appendtograph CutSeg
		appendtograph UpperLV vs Upper 
		appendtograph LowerLV vs Lower
		ModifyGraph rgb($nameofwave(InputWave))=(65280,48896,48896),rgb(SegmentSM)=(0,0,0);DelayUpdate
		ModifyGraph mode(UpperLV)=3,marker(UpperLV)=19,rgb(UpperLV)=(0,52224,26368);DelayUpdate
		ModifyGraph mode(LowerLV)=3,marker(LowerLV)=19,rgb(LowerLV)=(14848,32256,47104)
		ModifyGraph lsize(CutSeg)=2,rgb(CutSeg)=(29440,0,58880)
		SetAxis bottom (CEnterGuess-max(10*TimeGuess,.2e-3)/2),(CEnterGuess+max(10*TimeGuess,.2e-3)/2)
	else
		killwindow DTE
		display/N=DTE 	InputWave
		appendtograph Segment
		appendtograph SegmentSm
		appendtograph CutSeg
		appendtograph UpperLV vs Upper 
		appendtograph LowerLV vs Lower
		ModifyGraph rgb($nameofwave(InputWave))=(65280,48896,48896),rgb(SegmentSM)=(0,0,0);DelayUpdate
		ModifyGraph mode(UpperLV)=3,marker(UpperLV)=19,rgb(UpperLV)=(0,52224,26368);DelayUpdate
		ModifyGraph mode(LowerLV)=3,marker(LowerLV)=19,rgb(LowerLV)=(14848,32256,47104)
		ModifyGraph lsize(CutSeg)=2,rgb(CutSeg)=(29440,0,58880)
		SetAxis bottom (CEnterGuess-max(10*TimeGuess,.2e-3)/2),(CEnterGuess+max(10*TimeGuess,.2e-3)/2)
	endif
	if (UserCursorAdjust("DTE",5) != 0)
		//print "aborted"
	endif
			killwindow DTE

	endif
	wave W_FINDLEVELS
	killwaves Segment,SegmentSm,CutSeg,Lower,Upper,LowerLV,UpperLV,W_FindLevels
	return Result 
end

Static Function/S FindProductiveTransitions(Locations1,Locations2)
	wave Locations1, Locations2
	variable n,m
	String IDs=""
	string Insert
	for(n=0;n<(numpnts(Locations1)-1);n+=1)
	
		for(m=0;m<numpnts(Locations2);m+=1)
			if(Locations2[m]>Locations1[n+1])
				break
			elseif(Locations2[m]>Locations1[n])
			sprintf  Insert "%f15,%f15", Locations1[n],Locations2[m]
				IDs+=Insert+";"

				Return IDs
				endif
		endfor
	
	endfor
		for(m=0;m<numpnts(Locations2);m+=1)
			if(Locations2[m]>Locations1[numpnts(Locations1)-1])

			sprintf  Insert "%f15,%f15", Locations1[numpnts(Locations1)-1],Locations2[m]
				IDs+=Insert+";"
				return IDS

				break
			elseif(Locations2[m])
				endif
		endfor

end



Static Function CalculateTwoStateLifeTimes(IndexWave,Lifetime1,Lifetime2)
	wave IndexWave,Lifetime1,Lifetime2
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
duplicate/o State1 Lifetime1
duplicate/o State2 Lifetime2

end
Function UserCursorAdjust(graphName,autoAbortSecs)
	String graphName
	Variable autoAbortSecs
	DoWindow/F $graphName // Bring graph to front
	if (V_Flag == 0) // Verify that graph exists
		Abort "UserCursorAdjust: No such graph."
		return -1
	endif
	
	
	NewDataFolder/O root:tmp_PauseforCursorDF
	Variable/G root:tmp_PauseforCursorDF:canceled= 0



	
	NewPanel /K=2 /W=(187,368,437,531) as "Pause for Cursor"
	DoWindow/C tmp_PauseforCursor // Set to an unlikely name
	AutoPositionWindow/E/M=1/R=$graphName // Put panel near the graph
	DrawText 21,20,"Adjust the cursors and then"
	DrawText 21,40,"Click Continue."
	Button button0,pos={80,58},size={92,20},title="Continue"
	Button button0,proc=UserCursorAdjust_ContButtonProc
	Button button1,pos={80,88},size={92,20},title="Kill"
	Button button1,proc=UserCursorAdjust_KillButtonProc
	Variable didAbort= 0
	if( autoAbortSecs == 0 )
		PauseForUser tmp_PauseforCursor,$graphName
	else
		SetDrawEnv textyjust= 1
		DrawText 162,103,"sec"
		SetVariable sv0,pos={48,97},size={107,15},title="Aborting in "
		SetVariable sv0,limits={-inf,inf,0},value= _NUM:10
		Variable td= 10,newTd
		Variable t0= ticks
		Do
			newTd= autoAbortSecs - round((ticks-t0)/60)
			if( td != newTd )
				td= newTd
				SetVariable sv0,value= _NUM:newTd,win=tmp_PauseforCursor
				if( td <= 10 )
					SetVariable sv0,valueColor= (65535,0,0),win=tmp_PauseforCursor
				endif
			endif
			if( td <= 0 )
				DoWindow/K tmp_PauseforCursor
				didAbort= 1
				break
			endif
			PauseForUser/C tmp_PauseforCursor,$graphName
		while(V_flag)
	endif
	NVAR gCaneled= root:tmp_PauseforCursorDF:canceled
	Variable canceled= gCaneled			// Copy from global to local before global is killed
	KillDataFolder root:tmp_PauseforCursorDF

	return canceled
End
Function UserCursorAdjust_ContButtonProc(ctrlName) : ButtonControl
	String ctrlName
	DoWindow/K tmp_PauseforCursor // Kill panel
End

Function UserCursorAdjust_KillButtonProc(ctrlName) : ButtonControl
	String ctrlName
	Variable/G root:tmp_PauseforCursorDF:canceled= 1

	DoWindow/K tmp_PauseforCursor // Kill panel
End

Static Function SaveOut(BaseWave,IndexMatric,WaveMatrix,TransitionPaths,Lifetime1,Lifetime2)
	wave BaseWave,IndexMatric,TransitionPaths,Lifetime1,Lifetime2,WaveMatrix
	String Base=nameofwave(BaseWave)
	String x=Base+"_Matrix"
	controlinfo de_Viewer_popup3
	String y=Base+"_TransitionPath"
	String z=Base+"_Lifetime1"
	String i=Base+"_Lifetime2"
	String j=Base+"_Segments"

	Prompt x, "Transition Matrix Name: "		// Set prompt for x param
	Prompt y, "Enter Transition path Wave Name : "		// Set prompt for y param
		Prompt z, "Enter First Lifetime Wave Name : "		// Set prompt for y param
	Prompt i, "Enter Second Lifetime path Wave Name : "		// Set prompt for y param
	Prompt j, "Enter Segments path Wave Name : "		// Set prompt for y param

	DoPrompt "Enter Names", x, y,z,i,j
	if (V_Flag)
		return -1								// User canceled
	endif
	duplicate/o IndexMatric $X
	duplicate/o TransitionPaths $Y
	duplicate/o Lifetime1 $Z
	duplicate/o Lifetime2 $i
	duplicate/o WaveMatrix $j

end

Static Function PlotASegment(ForceWave,LogisticWave,TransitionWave,index)
	wave ForceWave,LogisticWave,TransitionWave
	variable index
	variable TimeGuess,CenterGuess,centerlevel,v1,v2,C1,C2,C3,C4
	TimeGuess=LogisticWave[index][2]
	CenterGuess=LogisticWave[index][3]
	make/n=0/o Segment
	SliceTransitionbyIndex(ForceWave,LogisticWave,Index,max(50*TimeGuess,2e-3),Segment)
	Duplicate/o Segment SegmentSM
	print "Filtering: "+num2str(round(TimeGuess/2/dimdelta(Segment,0)))
	Smooth/M=0 round(TimeGuess/2/dimdelta(Segment,0)), SegmentSM 
	centerlevel=SegmentSM(CenterGuess)
	v1=TransitionWave[index][1]
	v2=TransitionWave[index][2]
	duplicate/o/r=(v1,v2) SegmentSM CutSeg
	Make/o/n=1 Lower,LowerLv,Upper,UpperLv
	make/o/n=0 LogFit
	C1=LogisticWave[index][1]
	C2=LogisticWave[index][4]
	C3=LogisticWave[index][5]
	C4=LogisticWave[index][6]
	
	LogFit=MakeLogisticFit(Segment,C1,C2,C3,C4,LogFit)
	Lower=v1
	LowerLv=SegmentSm(Lower)
	Upper=v2
	UpperLv=SegmentSm(Upper)
	DoWindow DTE
	if( V_Flag==0 )
		display/N=DTE 	ForceWave
		appendtograph Segment
		appendtograph SegmentSm
		appendtograph CutSeg
		appendtograph UpperLV vs Upper 
		appendtograph LowerLV vs Lower
		appendtograph LogFit
		ModifyGraph rgb($nameofwave(ForceWave))=(65280,48896,48896)
		ModifyGraph rgb(Segment)=(65280,48896,48896)
		ModifyGraph rgb(SegmentSM)=(58368,6656,7168)
		ModifyGraph mode(UpperLV)=3,marker(UpperLV)=19,rgb(UpperLV)=(0,52224,26368);DelayUpdate
		ModifyGraph mode(LowerLV)=3,marker(LowerLV)=19,rgb(LowerLV)=(14848,32256,47104)
		ModifyGraph lsize(CutSeg)=2,rgb(CutSeg)=(29440,0,58880)
		ModifyGraph lstyle(LogFit)=2,rgb(LogFit)=(0,0,0),lsize(LogFit)=1.5
		SetAxis bottom (CEnterGuess-max(10*TimeGuess,.2e-3)/2),(CEnterGuess+max(10*TimeGuess,.2e-3)/2)
		SetAxis/A=2 left
		TextBox/C/N=text0/F=0/A=LT "Time="+num2str(TransitionWave[index][3]*1e6)+" µs"
	else
	killwindow DTE
		display/N=DTE 	ForceWave
		appendtograph Segment
		appendtograph SegmentSm
		appendtograph CutSeg
		appendtograph UpperLV vs Upper 
		appendtograph LowerLV vs Lower
		appendtograph LogFit
		ModifyGraph rgb($nameofwave(ForceWave))=(65280,48896,48896)
		ModifyGraph rgb(Segment)=(65280,48896,48896)
		ModifyGraph rgb(SegmentSM)=(58368,6656,7168)
		ModifyGraph mode(UpperLV)=3,marker(UpperLV)=19,rgb(UpperLV)=(0,52224,26368);DelayUpdate
		ModifyGraph mode(LowerLV)=3,marker(LowerLV)=19,rgb(LowerLV)=(14848,32256,47104)
		ModifyGraph lsize(CutSeg)=2,rgb(CutSeg)=(29440,0,58880)
		ModifyGraph lstyle(LogFit)=7,rgb(LogFit)=(0,0,0),lsize(LogFit)=1.5
		SetAxis bottom (CEnterGuess-max(10*TimeGuess,.2e-3)/2),(CEnterGuess+max(10*TimeGuess,.2e-3)/2)
		SetAxis/A=2 left
			TextBox/C/N=text0/F=0/A=LT "Time="+num2str(TransitionWave[index][3]*1e6)+" µs"

	endif

end

Static Function LBP(ctrlName,row,col,event) : ListBoxControl

	String ctrlName
	Variable row
	Variable col
	Variable event	//1=mouse down, 2=up, 3=dbl click, 4=cell select with mouse or keys
	//5=cell select with shift key, 6=begin edit, 7=end
	//
	switch (event)
		case -1:
			break
		case 1:
			break
		case 3:
			DE_TransitionPath#SelectRamp()
			//UpdateLocalPoints()
			//DE_MultiRampViewer#MakeSlopeFits()

			break			
		case 4:	
			DE_TransitionPath#SelectRamp()
			//UpdateLocalPoints()
		//	DE_MultiRampViewer#MakeSlopeFits()
			break
	endswitch
	
	return 0
End

Static Function SelectRamp()

	string saveDF = GetDataFolder(1)
	controlinfo de_TP_popup0
	SetDataFolder s_value

	

	wave ForceWave=root:DE_TP:FullForceTrace
	wave LogisticWave=root:DE_TP:Matrix
	wave TransitionWave=root:DE_TP:Segments

	controlinfo de_TP_list0
	
	variable index=V_Value

	variable TimeGuess,CenterGuess,centerlevel,v1,v2,C1,C2,C3,C4
	TimeGuess=LogisticWave[index][2]
	CenterGuess=LogisticWave[index][3]
	make/n=0/free Segment
	SliceTransitionbyIndex(ForceWave,LogisticWave,index,max(50*TimeGuess,2e-3),Segment)
	Duplicate/free Segment SegmentSM
	variable Filtering=TransitionWave[index][7]
	print "Filtering: "+num2str(Filtering)
	Smooth/M=0 Filtering, SegmentSM 
	centerlevel=SegmentSM(CenterGuess)
	v1=TransitionWave[index][1]
	v2=TransitionWave[index][2]
	duplicate/free/r=(v1,v2) SegmentSM CutSeg
	Make/free/n=1 Lower,LowerLv,Upper,UpperLv
	make/free/n=0 LogFit
	C1=LogisticWave[index][1]
	C2=LogisticWave[index][4]
	C3=LogisticWave[index][5]
	C4=LogisticWave[index][6]
	LogFit=MakeLogisticFit(Segment,C1,C2,C3,C4,LogFit)
	Lower=v1
	LowerLv=SegmentSm(Lower)
	Upper=v2
	UpperLv=SegmentSm(Upper)
	duplicate/o Segment root:DE_TP:Unfilt
	duplicate/o SegmentSM root:DE_TP:smoothed
	duplicate/o CutSeg root:DE_TP:TransitionPath
	duplicate/o LogFit root:DE_TP:LogisticFit
duplicate/o Lower root:DE_TP:Lower
duplicate/o LowerLV root:DE_TP:LowerLV
duplicate/o Lower root:DE_TP:Upper
duplicate/o Upper	LV root:DE_TP:UpperLV

		SetAxis/W=ScanTP#Data bottom (CEnterGuess-max(20*TimeGuess,.2e-3)/2),(CEnterGuess+max(20*TimeGuess,.2e-3)/2)
		SetAxis/W=ScanTP#Data/A=2 left
		TextBox/W=ScanTP#Data/C/N=text0/F=0/A=LT "Time="+num2str(TransitionWave[index][3]*1e6)+" µs"
	SetDataFolder saveDF
	
end


Static Function PopMenuProc(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa
	Variable popNum = pa.popNum
	String popStr = pa.popStr
	string saveDF
	strswitch(pa.ctrlName)
	
	case "de_TP_popup1":
	switch( pa.eventCode )
		case 2: // mouse up
			DE_TransitionPath#LoadForce()
			break
		case -1: // control being killed
			break
	endswitch
	break 
	
	endswitch
	return 0
End

Static Function LoadForce()

	string saveDF = GetDataFolder(1)
	controlinfo de_TP_popup0
	SetDataFolder s_value
	
	
	controlinfo de_TP_popup1
	wave ForceStart=$S_value
	wave Matrix=$(S_Value+"_Matrix")
	wave Segments=$(S_Value+"_Segments")

	duplicate/o ForceStart root:DE_TP:FullForceTrace
	duplicate/o Matrix root:DE_TP:Matrix
	duplicate/o Segments root:DE_TP:Segments

	
	DE_TransitionPath#UpdateSegmentList()
end

Static Function UpdateSegmentList()	
	string saveDF = GetDataFolder(1)
	controlinfo de_TP_popup0
	SetDataFolder s_value

	wave ForceWave=root:DE_TP:FullForceTrace
	wave Matrix=root:DE_TP:Matrix
	wave Segments=root:DE_TP:Segments

//	SetDataFolder saveDF
	variable Num=dimsize(Segments,0)
	make/T/o/n=(num) root:DE_TP:ListWave1
	wave/T/z LW=root:DE_TP:ListWave1
	LW=num2str(p)
	make/o/n=(num) root:DE_TP:SelWave1
	wave SW=root:DE_TP:SelWave1
	SW=0
	SW[0]=1

end

Macro ScanTransitions() : Panel
	//
	PauseUpdate; Silent 1		// building window...
	NewPanel/N=ScanTP /W=(50,50,1300,750)
	NewDataFolder/o root:DE_TP
	make/T/o/n=0 root:DE_TP:ListWave1
	make/o/n=0 root:DE_TP:SelWave1
	PopupMenu de_TP_popup0,pos={250,2},size={129,21}
	PopupMenu de_TP_popup0,mode=1,popvalue="X",value= #"DE_PanelProgs#ListFolders()"
	PopupMenu de_TP_popup1,pos={250,40},size={129,21},proc=DE_TransitionPath#PopMenuProc
	PopupMenu de_TP_popup1,mode=1,popvalue="X",value= #"DE_TransitionPath#ListWaves(\"de_TP_popup0\")"
	
	make/o/n=0 root:DE_TP:Unfilt, root:DE_TP:Smoothed,root:DE_TP:LogisticFit, root:DE_TP:TransitionPath
//
	display/host=ScanTP/N=Data/W=(250,100,700,400) root:DE_TP:Unfilt
	appendtograph/W=ScanTP#Data root:DE_TP:Smoothed
	appendtograph/W=ScanTP#Data root:DE_TP:LogisticFit
	appendtograph/W=ScanTP#Data root:DE_TP:TransitionPath
	
	•ModifyGraph/W=ScanTP#Data rgb(Unfilt)=(63232,45824,45824),rgb(Smoothed)=(58368,6656,7168);DelayUpdate
•ModifyGraph/W=ScanTP#Data lsize(LogisticFit)=1.5,rgb(LogisticFit)=(0,0,0);DelayUpdate
•ModifyGraph/W=ScanTP#Data lsize(TransitionPath)=1.5,rgb(TransitionPath)=(29440,0,58880)
//	ModifyGraph/W=MRViewer#Data rgb(YDisp)=(63232,45824,45824),rgb(YDisp#1)=(49152,54784,59648);DelayUpdate
//	ModifyGraph/W=MRViewer#Data rgb(YDispSm)=(58368,6656,7168),rgb(YDispSm#1)=(14848,32256,47104)
//	
	ListBox de_Tp_list0,pos={10,40},size={100,100},proc=DE_TransitionPath#LBP,listWave=root:DE_TP:ListWave1
	ListBox de_TP_list0,row= 0,mode=2,selRow= 0//,selWave=root:DE_Viewer:SelWave1
//	
//	SetVariable de_Viewer_setvar0,pos={310,10},size={100,20},proc=DE_MultiRampViewer#SVP,value= _NUM:0,title="Smoothing"
//	SetVariable de_Viewer_setvar0 value= _NUM:5e-9
//	SetVariable de_Viewer_setvar1,pos={710,10},size={180,20},proc=DE_MultiRampViewer#SVP,value= _NUM:0,title="Starting BackGround"
//	SetVariable de_Viewer_setvar1 value= _NUM:1000,limits={1,inf,1}
//	
//	
//	make/o/n=0 root:DE_Viewer:AllUpPoints,root:DE_Viewer:AllDownPoints,root:DE_Viewer:CurrUpF,root:DE_Viewer:CurrUpX,root:DE_Viewer:CurrDownF,root:DE_Viewer:CurrDownX
//	make/o/n=0 root:DE_Viewer:SelUpX,root:DE_Viewer:SelUpY,root:DE_Viewer:SelDownX,root:DE_Viewer:SelDownY
//	appendtograph/W=MRViewer#TimeData  root:DE_Viewer:CurrUpF vs root:DE_Viewer:CurrUpX
//	appendtograph/W=MRViewer#TimeData root:DE_Viewer:CurrDownF vs root:DE_Viewer:CurrDownX
//	appendtograph/W=MRViewer#TimeData root:DE_Viewer:SelUpY vs root:DE_Viewer:SelUpX
//	appendtograph/W=MRViewer#TimeData root:DE_Viewer:SelDownY vs root:DE_Viewer:SelDownX
//	ModifyGraph mode(SelUpY)=3,marker(SelUpY)=17,mode(SelDownY)=3,marker(SelDownY)=17
//	•ModifyGraph msize(SelUpY)=4,rgb(SelUpY)=(14848,32256,47104),msize(SelDownY)=4;DelayUpdate
//	ModifyGraph rgb(SelDownY)=(19712,44800,18944)
//	ModifyGraph/W=MRViewer#TimeData mode(CurrUpF)=3,mode(CurrDownF)=3
//	ModifyGraph/W=MRViewer#TimeData marker(CurrUpF)=19,rgb(CurrUpF)=(2816,5632,8192),marker(CurrDownF)=19;DelayUpdate
//	ModifyGraph/W=MRViewer#TimeData rgb(CurrDownF)=(29440,0,58880)
//	Cursor/W=MRViewer#TimeData A, YDispSm, 0	// cursor A on first point of myWave
//
//	make/t/o/n=(0) root:DE_Viewer:TCurrUpF,root:DE_Viewer:TCurrDownF,root:DE_Viewer:TCurrUpL,root:DE_Viewer:TCurrDownL
//	make/o/n=0 root:DE_Viewer:TCurrUpL_Sel,root:DE_Viewer:TCurrDownL_Sel
//	ListBox de_Viewer_list1,pos={10,200},size={100,100},proc=DE_MultiRampViewer#LBP1,listWave=root:DE_Viewer:TCurrUpF
//	ListBox de_Viewer_list1,row= 0,mode=2,selRow= 0//selWave=root:DE_Viewer:SelWave1
//	ListBox de_Viewer_list2,pos={10,400},size={100,100},proc=DE_MultiRampViewer#LBP1,listWave=root:DE_Viewer:TCurrDownF
//	ListBox de_Viewer_list2,row= 0,mode=2,selRow= 0//selWave=root:DE_Viewer:SelWave1
//	ListBox de_Viewer_list3,pos={110,200},size={100,100},proc=DE_MultiRampViewer#LBP2,listWave=root:DE_Viewer:TCurrUpL
//	ListBox de_Viewer_list3,row= 0,mode=2,selRow= 0,selWave=root:DE_Viewer:TCurrUpL_Sel
//	ListBox de_Viewer_list4,pos={110,400},size={100,100},proc=DE_MultiRampViewer#LBP2,listWave=root:DE_Viewer:TCurrDownL
//	ListBox de_Viewer_list4,row= 0,mode=2,selRow= 0,selWave=root:DE_Viewer:TCurrDownL_Sel
//
//	make/o/n=0 root:DE_Viewer:AllUpForceHist,root:DE_Viewer:AllDownForceHist,root:DE_Viewer:AllUpForce,root:DE_Viewer:AllDownForce
//	display/host=MRViewer/N=Hist/W=(750,100,950,400) root:DE_Viewer:AllUpForceHist
//	appendtograph/W=MRViewer#Hist root:DE_Viewer:AllDownForceHist
//	ModifyGraph/W=MRViewer#Hist rgb(AllUpForceHist)=(19712,44800,18944)
//	ModifyGraph/W=MRViewer#Hist rgb(AllDownForceHist)=(0,0,65280)
//	ModifyGraph/W=MRViewer#Hist mode=5,hbFill=5
//	make/o/n=0  root:DE_viewer:UpLoad,root:DE_viewer:DownLoad,root:DE_viewer:CurrUpLoad,root:DE_viewer:CurrDownLoad
//	display/host=MRViewer/N=DFS/W=(1000,100,1200,400) root:DE_Viewer:AllUpForce vs root:DE_viewer:UpLoad
//	appendtograph/W=MRViewer#DFS root:DE_Viewer:AllDownForce vs root:DE_viewer:DownLoad
//	appendtograph/W=MRViewer#DFS root:DE_Viewer:CurrUpF vs root:DE_viewer:CurrUpLoad
//	appendtograph/W=MRViewer#DFS root:DE_Viewer:CurrDownF vs root:DE_viewer:CurrDownLoad
//	ModifyGraph/W=MRViewer#DFS mode=3,marker=8
//	ModifyGraph/W=MRViewer#DFS msize=4,mrkThick=1.5
//	ModifyGraph/W=MRViewer#DFS rgb(AllUpForce)=(2816,5632,8192),rgb(AllDownForce)=(29440,0,58880)
//
//	
//	Button de_viewer_but0,pos={750,500},size={150,25},proc=DE_MultiRampViewer#ButtonProc,title="Delete This Trace"
//	Button de_viewer_but1,pos={910,540},size={150,25},proc=DE_MultiRampViewer#ButtonProc,title="Modify Selected Unfolding"
//	Button de_viewer_but6,pos={1070,540},size={150,25},proc=DE_MultiRampViewer#ButtonProc,title="Delete Selected Unfolding"
//	Button de_viewer_but8,pos={750,540},size={150,25},proc=DE_MultiRampViewer#ButtonProc,title="Add Unfolding"
//
//	Button de_viewer_but2,pos={910,580},size={150,25},proc=DE_MultiRampViewer#ButtonProc,title="Modify Selected Folding"
//	Button de_viewer_but7,pos={1070,580},size={150,25},proc=DE_MultiRampViewer#ButtonProc,title="Delete Selected Folding"
//	Button de_viewer_but9,pos={750,580},size={150,25},proc=DE_MultiRampViewer#ButtonProc,title="Add Folding"
//
//	Button de_viewer_but3,pos={750,620},size={150,25},proc=DE_MultiRampViewer#ButtonProc,title="Reset"
//	Button de_viewer_but4,pos={750,660},size={150,25},proc=DE_MultiRampViewer#ButtonProc,title="Save Working"
//	//Button de_viewer_but5,pos={100,750},size={150,25},proc=DE_MultiRampViewer#ButtonProc,title="Update Lengths"
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


Static Function MakeLogisticFit(Segment,C1,C2,C3,C4,OutWave)
	wave	OutWave,Segment
	variable C1,C2,C3,C4
	make/free/n=4 Coefs
	Coefs={C1,C2,C3,C4}
	duplicate/o Segment OutWave
	OutWave=Logistic(Coefs,x)

end