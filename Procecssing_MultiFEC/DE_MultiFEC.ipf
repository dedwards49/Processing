#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma modulename=DE_MultiFEC
#include "SimpleWLCPrograms"
#include "DE_Filtering"
#include "DE_NewFeather"
Function FindAllForces(filtering,[bottom,top])
	variable filtering,bottom,top
	string AllForceRet= wavelist("*Force_Ret",";","")
	String ForceWaveList="",SepWaveList=""
	if(ParamisDefault(Bottom))
	Bottom=0
	endif
	if(ParamisDefault(top))
	top=itemsinlist(AllFOrceRet)
	endif
	variable n

	for(n=bottom;n<top;n+=1)
	print "N:"+num2str(n)
		//for(n=0;n<itemsinlist(AllFOrceRet);n+=1)
		wave ForceRetWave=$stringfromlist(n,AllForceRet)
		wave ForceExtWave=$replacestring("Ret",nameofwave(ForceRetWave),"Ext")
		wave SepRetWave=$replacestring("Force",nameofwave(ForceRetWave),"Sep")
		wave SepExtWave=$replacestring("Force",nameofwave(ForceExtWave),"Sep")
		make/free/n=0 ForceAll,SepAll
		Concatenate/o {ForceExtWave,ForceRetWave},ForceAll
		Concatenate/o {SepExtWave,SepRetWave},SepAll

		make/free/n=0 FSm,SSm
		make/o/n=0 FRetSm,SRetSm
		if(filtering==1)
			duplicate/free ForceAll FSM
			duplicate/free SepAll SSM
			duplicate/o ForceRetWave FRetSm
			duplicate/o SepRetWave SRetSm
		elseif(filtering>5)
			DE_Filtering#FilterForceSep(ForceAll,SepAll,FSm,SSm,"SVG",filtering)
			DE_Filtering#FilterForceSep(ForceRetWave,SepRetWave,FRetSm,SRetSm,"SVG",filtering)

		elseif(filtering<1)
			DE_Filtering#FilterForceSep(ForceAll,SepAll,FSm,SSm,"TVD",filtering)
			DE_Filtering#FilterForceSep(ForceRetWave,SepRetWave,FRetSm,SRetSm,"TVD",filtering)

		
		endif

		duplicate/o FSm $(Replacestring("Force_Ret",nameofwave(ForceRetWave),"Force"))
		wave ForceWave=$(Replacestring("Force_Ret",nameofwave(ForceRetWave),"Force"))
		duplicate/o SSm $(Replacestring("force_Ret",nameofwave(ForceRetWave),"Sep"))
		wave SepWave=$(Replacestring("Force_Ret",nameofwave(ForceRetWave),"Sep"))
		make/free/n=5 OptionsWave
		OptionsWave={4e-7,2e-3,str2num(StringbyKey("TriggerTime",note(ForceWave),":","\r")),0.0,str2num(StringbyKey("SpringConstant",note(ForceWave),":","\r"))}
		DE_NewFeather#OutportForce(ForceWave,SepWave)
		DE_NewFeather#RunFeatheronOutput(OptionsWave)
		wave event_starts

		duplicate/o event_starts $replacestring("Force_ext",nameofwave(ForceExtWave),"Starts")
		duplicate/o FRetSm $replacestring("Force_ext",nameofwave(ForceExtWave),"FSm")
		duplicate/o SRetSm $replacestring("Force_ext",nameofwave(ForceExtWave),"SSm")
		killwaves event_starts,ForceAll,SepAll,FretSm,SRetSm
		
	endfor
	print 	top

//	
//	for(n=bottom;n<top;n+=1)
//		wave ForceRetWave=$stringfromlist(n,AllForceRet)
//		wave ForceExtWave=$replacestring("Ret",nameofwave(ForceRetWave),"Ext")
//		wave SepRetWave=$replacestring("Force",nameofwave(ForceRetWave),"Sep")
//		wave SepExtWave=$replacestring("Force",nameofwave(ForceExtWave),"Sep")
//		wave FRetSm=$replacestring("Force_ext",nameofwave(ForceExtWave),"FSm")
//		wave SRetSm=$replacestring("Force_ext",nameofwave(ForceExtWave),"SSm")
//		wave ThisEvent=$replacestring("Force_ext",nameofwave(ForceExtWave),"Starts")
//		duplicate/o ThisEvent FreeRupPnts
//		FreeRupPnts-=NUMPNTS(SepExtWave)
//		MakeNicePlot(ForceRetWave,sEPRetWave,FRetSm,SRetSm,FreeRupPnts)
//		FreeRupPnts+=NUMPNTS(SepExtWave)
//		duplicate/o FreeRupPnts ThisEvent
//		killwaves FreeRupPnts
//		//killwaves FRetSm,SRetSm
//
//	endfor
//	
	//ExportWaveLists(ForceWaveList,SepWaveList)
end

Static Function TestZero(ForceRetWave,[distance])
	wave ForceRetWave
	variable distance
	if(ParamisDefault(distance))
		distance=40e-9
	endif
	wave ForceExtWave=$replacestring("Ret",nameofwave(ForceRetWave),"Ext")
	wave SepRetWave=$replacestring("Force",nameofwave(ForceRetWave),"Sep")
	wave SepExtWave=$replacestring("Force",nameofwave(ForceExtWave),"Sep")
	wave FRetSm=$replacestring("Force_ext",nameofwave(ForceExtWave),"FSm")
	wave SRetSm=$replacestring("Force_ext",nameofwave(ForceExtWave),"SSm")
	wave ThisEvent=$replacestring("Force_ext",nameofwave(ForceExtWave),"Starts")
	if(numpnts(ThisEvent)<4)
		//print nameofwave(ForceExtWave)+":Booted"

	return 0
	endif
//	variable FinalEvent=ThisEVent[numpnts(ThisEVent)-1]-numpnts(ForceExtWave)
//	variable offsetpoints=distance/str2num(stringbykey("Velocity",note(ForceRetWave),":","\r"))/dimdelta(ForceRetWave,0)
//	wavestats/Q/r=[FinalEvent+offsetpoints,FinalEvent+2*offsetpoints] ForceRetWave
	variable FinalEvent=ThisEVent[numpnts(ThisEVent)-1]-numpnts(ForceExtWave)
	variable offsetpoints=(numpnts(ForceRetWave)-FinalEvent)/2
	wavestats/Q/r=[FinalEvent+offsetpoints,] ForceRetWave

	variable postrup= v_avg
	wavestats/q/r=[numpnts(ForceRetWave)-1-offsetpoints,numpnts(ForceRetWave)-1] ForceRetWave
	variable traceend= v_avg
	return postrup
end


Function TestCondition(Number,filtering,tol,temporal)
	variable Number,filtering,tol,temporal
	string AllForceRet= wavelist("*Force_Ret",";","")
	String ForceWaveList="",SepWaveList=""
	variable n
	variable bottom=Number
	variable top=Number+1//itemsinlist(AllFOrceRet)
	for(n=bottom;n<top;n+=1)
		//for(n=0;n<itemsinlist(AllFOrceRet);n+=1)
		wave ForceRetWave=$stringfromlist(n,AllForceRet)
		wave ForceExtWave=$replacestring("Ret",nameofwave(ForceRetWave),"Ext")
		wave SepRetWave=$replacestring("Force",nameofwave(ForceRetWave),"Sep")
		wave SepExtWave=$replacestring("Force",nameofwave(ForceExtWave),"Sep")
		make/free/n=0 ForceAll,SepAll
		Concatenate/o {ForceExtWave,ForceRetWave},ForceAll
		Concatenate/o {SepExtWave,SepRetWave},SepAll

		make/free/n=0 FSm,SSm
		make/o/n=0 FRetSm,SRetSm
		if(filtering==1)
			duplicate/free ForceAll FSM
			duplicate/free SepAll SSM
			duplicate/o ForceRetWave FRetSm
			duplicate/o SepRetWave SRetSm
		elseif(filtering>5)
			DE_Filtering#FilterForceSep(ForceAll,SepAll,FSm,SSm,"SVG",filtering)
			DE_Filtering#FilterForceSep(ForceRetWave,SepRetWave,FRetSm,SRetSm,"SVG",filtering)

		elseif(filtering<1)
			DE_Filtering#FilterForceSep(ForceAll,SepAll,FSm,SSm,"TVD",filtering)
			DE_Filtering#FilterForceSep(ForceRetWave,SepRetWave,FRetSm,SRetSm,"TVD",filtering)

		
		endif

		duplicate/o FSm $(Replacestring("Force_Ret",nameofwave(ForceRetWave),"Force"))
		wave ForceWave=$(Replacestring("Force_Ret",nameofwave(ForceRetWave),"Force"))
		duplicate/o SSm $(Replacestring("force_Ret",nameofwave(ForceRetWave),"Sep"))
		wave SepWave=$(Replacestring("Force_Ret",nameofwave(ForceRetWave),"Sep"))
		make/free/n=5 OptionsWave
		OptionsWave={tol,temporal,str2num(StringbyKey("TriggerTime",note(ForceWave),":","\r")),0.0,str2num(StringbyKey("SpringConstant",note(ForceWave),":","\r"))}
		DE_NewFeather#OutportForce(ForceWave,SepWave)
		DE_NewFeather#RunFeatheronOutput(OptionsWave)
		wave event_starts

		duplicate/o event_starts $replacestring("Force_ext",nameofwave(ForceExtWave),"Starts")
		duplicate/o FRetSm $replacestring("Force_ext",nameofwave(ForceExtWave),"FSm")
		duplicate/o SRetSm $replacestring("Force_ext",nameofwave(ForceExtWave),"SSm")
		killwaves event_starts,ForceAll,SepAll,FretSm,SRetSm
		
	endfor
	print 	top

//	
	for(n=bottom;n<top;n+=1)
		wave ForceRetWave=$stringfromlist(n,AllForceRet)
		wave ForceExtWave=$replacestring("Ret",nameofwave(ForceRetWave),"Ext")
		wave SepRetWave=$replacestring("Force",nameofwave(ForceRetWave),"Sep")
		wave SepExtWave=$replacestring("Force",nameofwave(ForceExtWave),"Sep")
		wave FRetSm=$replacestring("Force_ext",nameofwave(ForceExtWave),"FSm")
		wave SRetSm=$replacestring("Force_ext",nameofwave(ForceExtWave),"SSm")
		wave ThisEvent=$replacestring("Force_ext",nameofwave(ForceExtWave),"Starts")
		duplicate/o ThisEvent FreeRupPnts
		FreeRupPnts-=NUMPNTS(SepExtWave)
		MakeNicePlot(ForceRetWave,sEPRetWave,FRetSm,SRetSm,FreeRupPnts)
		FreeRupPnts+=NUMPNTS(SepExtWave)
		duplicate/o FreeRupPnts ThisEvent
		killwaves FreeRupPnts
		//killwaves FRetSm,SRetSm

	endfor
//	
	//ExportWaveLists(ForceWaveList,SepWaveList)
end
Function PlotaBlock(Start,EndNum)
	variable start,endnum
	variable n
	for(n=start;n<endnum;n+=1)
	
	TweakRupturesandAssignValue(n)
	print n
	endfor

end
Function TweakRupturesandAssignValue(Number)
	variable Number
	string AllForceRet= wavelist("*Force_Ret",";","")
	String ForceWaveList="",SepWaveList=""
	variable n
	variable top,start
	if(number==-1)
		start=0
		top=itemsinlist(AllFOrceRet)
	else
		start=number 
		top=number+1
	endif
	for(n=start;n<top;n+=1)
		wave ForceRetWave=$stringfromlist(n,AllForceRet)
		wave ForceExtWave=$replacestring("Ret",nameofwave(ForceRetWave),"Ext")
		wave SepRetWave=$replacestring("Force",nameofwave(ForceRetWave),"Sep")
		wave SepExtWave=$replacestring("Force",nameofwave(ForceExtWave),"Sep")
		wave FRetSm=$replacestring("Force_Ret",nameofwave(ForceRetWave),"FSm")
		wave SRetSm=$replacestring("Force_Ret",nameofwave(ForceRetWave),"SSm")
		wave ThisEvent=$replacestring("Force_Ret",nameofwave(ForceRetWave),"Starts")
		duplicate/o ThisEvent FreeRupPnts
		FreeRupPnts-=NUMPNTS(SepExtWave)
		MakeNicePlot(ForceRetWave,sEPRetWave,FRetSm,SRetSm,FreeRupPnts)
		
		make/free/n=0 LCSBack,SLopesBack
		duplicate/free FreeRupPnts ForcesBack
		//ForcesBack=FRetSm[FreeRupPnts[p]]

	FreeRupPnts+=numpnts(ForceExtWave)
				CalculateForcesFromPoints(ForceRetWave,SepRetWave,FreeRupPnts,ForcesBack,500)

		FreeRupPnts-=numpnts(ForceExtWave)
		

		variable offset=-1*str2num(stringbykey("DE_SChollOffset",note(ForceRetWave),":","\r"))-5e-9
		CalcAllLCs(ForceRetWave,SepRetWave,FreeRupPnts,offset,LCsBack)
		CalculateSlopes(ForceRetWave,SepRetWave,FreeRupPnts,offset,SlopesBack)
		AddNotes(ForceRetWave,SepRetWave,FreeRupPnts,ForcesBack,LCSBack,SLopesBack)
		FreeRupPnts+=NUMPNTS(SepExtWave)
		duplicate/o FreeRupPnts ThisEvent
		killwaves FreeRupPnts
		//killwaves FRetSm,SRetSm

	endfor

end




Function ReassignValue(Number)
	variable Number
	string AllForceRet= wavelist("*Force_Ret",";","")
	String ForceWaveList="",SepWaveList=""
	variable n
	variable top,start
	if(number==-1)
		start=0
		top=itemsinlist(AllFOrceRet)
	else
		start=number 
		top=number+1
	endif
	for(n=start;n<top;n+=1)
		string A=stringfromlist(n,AllForceRet)
		if(StrSearch(A,"fit",0)==-1)
		wave ForceRetWave=$stringfromlist(n,AllForceRet)
		wave ForceExtWave=$replacestring("Ret",nameofwave(ForceRetWave),"Ext")
		wave SepRetWave=$replacestring("Force",nameofwave(ForceRetWave),"Sep")
		wave SepExtWave=$replacestring("Force",nameofwave(ForceExtWave),"Sep")
		wave FRetSm=$replacestring("Force_Ret",nameofwave(ForceRetWave),"FSm")
		wave SRetSm=$replacestring("Force_Ret",nameofwave(ForceRetWave),"SSm")
		wave ThisEvent=$replacestring("Force_Ret",nameofwave(ForceRetWave),"Starts")
		Sort ThisEvent ThisEvent

		duplicate/o ThisEvent FreeRupPnts
	
		make/free/n=0 LCSBack,SLopesBack
		duplicate/free FreeRupPnts ForcesBack
		//ForcesBack=FRetSm[FreeRupPnts[p]]
		CalculateForcesFromPoints(ForceRetWave,SepRetWave,FreeRupPnts,ForcesBack,500)
		FreeRupPnts-=NUMPNTS(SepExtWave)

		variable offset=-1*str2num(stringbykey("DE_SChollOffset",note(ForceRetWave),":","\r"))-5e-9
		CalcAllLCs(ForceRetWave,SepRetWave,FreeRupPnts,offset,LCsBack)
		CalculateSlopes(ForceRetWave,SepRetWave,FreeRupPnts,offset,SlopesBack)
		AddNotes(ForceRetWave,SepRetWave,FreeRupPnts,ForcesBack,LCSBack,SLopesBack)
		FreeRupPnts+=NUMPNTS(SepExtWave)
		duplicate/o FreeRupPnts ThisEvent
		killwaves FreeRupPnts
		//killwaves FRetSm,SRetSm
		else
		endif

	endfor

end

Function RunAll()
	ReassignValue(-1)
	AddZeroOffsetsToNotes()
	ProcessForceCurves()
	ReturnRupForces()
	ReturnRupForcesZero()

	ReturnSlopes()
	ReturnRelevantLCs()
end

Static Function/D DE_Median(w) // Returns median value of wave w
Wave w
Variable result
Duplicate/o w, tempMedianWave // Make a clone of wave
Sort tempMedianWave, tempMedianWave // Sort clone
SetScale/P x 0,1,tempMedianWave
result = tempMedianWave((numpnts(tempMedianWave)-1)/2)
KillWaves tempMedianWave
return Result
end

Static Function MakeSingleContoursAndDisplay(ForceRetwave,SepRetWave,PointWave,index)
	
	wave ForceRetwave,SepRetWave,PointWave
	variable index
	variable tot=numpnts(PointWave)
	variable prevpnt
	variable result
	variable offset=-1*str2num(stringbykey("DE_SChollOffset",note(ForceRetwave),":","\r"))-5e-9
	FindLevels/P/Q SepRetWave, -1*offset
	wave w_FindLevels
	variable SurfacePnt=DE_MultiFEC#DE_Median(W_FindLevels)
	make/o/n=0 LCHistOut
	if(index==0)
		prevpnt=SurfacePnt
	else
		prevpnt=PointWave[index-1]+10//+max(100,PointWave[1]-PointWave[n-1]

	endif

	result= CalcCurrLC(ForceRetwave,SepRetWave,prevpnt,PointWave[index],offset,LCHistOut)
	DoWindow TempHistPlot
	if(V_flag==1)
		killwindow TempHistPlot
	else
	endif
	//Display/W=(500,500,700,700)/N=TempHistPlot LCHistOut
	//AutoPositionWindow/E/M=0/R=Test TempHistPlot // Put panel near the graph
	//TextBox/W=TempHistPlot/C/N=text0/F=0/A=RT "LC = "+num2str(result)
	return result

end


Static Function MakeSlopesAndAddtoPlot(ForceRetwave,SepRetWave,ForceRetSm,SepRetWaveSm,FreeRupPnts,offset,backcalc,SlopesBack)
	wave ForceRetwave,SepRetWave,ForceRetSm,SepRetWaveSm,FreeRupPnts,SlopesBack
	variable offset,backcalc
		FindLevels/P/Q SepRetWave, -1*offset
	wave w_FindLevels
	variable SurfacePnt=DE_MultiFEC#DE_Median(W_FindLevels)
	make/free/n=0 TempSLopesBack
	variable tot=numpnts(FreeRupPnts)
	variable n=0
	variable prevpnt
	duplicate/free FreeRupPnts TempSlopes
	for(n=0;n<tot;n+=1)
		if(n==0)
				prevpnt=Max(FreeRupPnts[n]-backcalc,SurfacePnt+10)

		else
				prevpnt=Max(FreeRupPnts[n]-backcalc,FreeRupPnts[n-1]+10)

		endif
		
		duplicate/free/R=[prevpnt,FreeRupPnts[n]] ForceRetwave, FFit 
		duplicate/o FFit $("FSlopefit_"+num2str(n))
		CurveFit/Q/W=2 line FFit/D=$("FSlopefit_"+num2str(n))
		wave w_coef,w_sigma
		TempSlopes[n]= w_coef[1]
		duplicate/free/R=[prevpnt,FreeRupPnts[n]] SepRetWave, SFit
	//	duplicate/o FFit $("FSlopefit_"+num2str(n))
	//	duplicate/o SFit $("SSlopefit_"+num2str(n))
//
	//	CurveFit/Q/W=2 line FFit /X=SFit /D=$("FSlopefit_"+num2str(n))
//
//		
//	
	endfor
	duplicate/o TempSlopes SlopesBack
	killwaves w_coef,w_sigma,w_FindLevels
		

end


Static Function CalcAllLCs(ForceWave,SepWave,PointWave,SepOff,LCsBack)
	wave ForceWave,SepWave,PointWave,LCsBack
	variable Sepoff
	variable tot=numpnts(PointWave)
	variable n=0
	variable prevpnt
	variable backcalc=5000
	FindLevels/P/Q SepWave, -1*Sepoff
	wave w_FindLevels
	variable SurfacePnt=DE_MultiFEC#DE_Median(W_FindLevels)
	duplicate/free PointWave TempLcs
	for(n=0;n<tot;n+=1)
		make/o/n=0 HistOUt
		if(n==0)
		prevpnt=SurfacePnt
		else
				prevpnt=PointWave[n-1]+10//+max(100,PointWave[1]-PointWave[n-1]

		endif
		TempLcs[n]= CalcCurrLC(ForceWave,SepWave,prevpnt,PointWave[n],SepOff,HistOut)
	
	endfor
	duplicate/o TempLCs LCsBack
	Killwaves HistOUt,W_FindLevels

end

Static Function CalculateSlopes(ForceWave,SepWave,PointWave,SepOff,SlopesBack)
	wave ForceWave,SepWave,PointWave,SlopesBack
	variable Sepoff
	variable tot=numpnts(PointWave)
	variable n=0
	variable prevpnt
	variable backdist=5e-9
	variable backcalc=backdist/str2num(stringbykey("retractvelocity",note(ForceWave),":","\r"))/dimdelta(ForceWave,0)
	FindLevels/P/Q SepWave, -1*Sepoff
	wave w_FindLevels
	variable SurfacePnt=DE_MultiFEC#DE_Median(W_FindLevels)
	
	duplicate/free PointWave TempSlopes
	for(n=0;n<tot;n+=1)
		if(n==0)
				prevpnt=max(PointWave[n]-backcalc,SurfacePnt)

		else
			prevpnt=max(PointWave[n]-backcalc,PointWave[n-1]+10)
		endif
		make/o/n=0 HistOUt
		duplicate/free/R=[prevpnt,PointWave[n]] ForceWave, FFit
		CurveFit/Q/W=2 line FFit
		wave w_coef,w_sigma
		TempSlopes[n]= w_coef[1]
		
	
	endfor
	TempSlopes*=-1
	duplicate/o TempSlopes SlopesBack
	killwaves w_coef,w_sigma,W_FindLevels
	Killwaves HistOUt

end
	
	
Static Function CalcCurrLC(ForceWave,SepWave,StartPnt,EndPnt,SepOff,HistOut)
	wave ForceWave,SepWave,HistOut
	variable StartPnt,EndPnt,SepOff
	
	duplicate/o/r=[startpnt,endpnt] FOrceWave TempForce,TempLC
	duplicate/o/r=[startpnt,endpnt] SepWave TempSep
	endpnt-=startpnt

	startpnt=0
	TempForce*=-1
	TempSep+=SepOff
	TempLC=DE_WLC#ContourTransform(TempForce,TempSep,.4e-9,298)
	variable lastLC=DE_WLC#ContourTransform(TempForce[endpnt],TempSep[endpnt],.4e-9,298)
	variable startLC=TempSep[endpnt]
	variable stepLC=2*(lastlc-startlc)/49
	make/o/n=50 TempHist 
	variable Q=numpnts(TempLC)
	Histogram/C/B={startLC,stepLC,50} TempLC,TempHist;
	CurveFit/Q/W=2 gauss TempHist 
	wave w_coef,w_sigma
	variable result=w_coef[2]
	killwaves w_coef,w_sigma
	duplicate/o TempHist HistOUt
	return result
end


Static Function MakeNicePlot(ForceWave,SepWave,ForceRetSm,SepRetSm,FreeRupPnts)
	wave ForceWave,SepWave,FreeRupPnts,ForceRetSm,SepRetSm
	Dowindow Test
	if(V_Flag==1)
		killwindow Test
	endif
	make/o/n=(numpnts(FreeRupPnts),2) RupTimes
	RupTimes[][0]=SepWave[FreeRupPnts]
	RupTimes[][1]=ForceWave[FreeRupPnts[p]]
		make/free/n=0 LCSBack,SLopesBack
		wave ForceExtWave=$replacestring("Ret",nameofwave(ForceWave),"Ext")
	variable offset=-1*str2num(stringbykey("DE_SChollOffset",note(ForceWave),":","\r"))-5e-9
	variable PointsToFit=1000
	//MakeSlopesAndAddtoPlot(ForceWave,SepWave,ForceRetSm,SepRetSm,FreeRupPnts,offset,PointsToFit,SlopesBack)
	variable zeroforce=TestZero(ForceWave)
	if(abs(zeroforce)>2e-12)
	print "Beware of the 0"
	print	zeroforce
	endif
	display/N=Test ForceWave vs SepWave
	Appendtograph/W=Test RupTimes[][1] vs RupTimes[][0]
	Appendtograph/W=Test ForceRetSm vs SepRetSm
	ModifyGraph/W=Test rgb($nameofwave(ForceWave))=(65535,49151,49151)
	ModifyGraph/W=Test mode($nameofwave(RupTimes))=3,marker($nameofwave(RupTimes))=16,rgb($nameofwave(RupTimes))=(0,0,0)
	ModifyGraph/W=Test muloffset={0,-1}
	if (DE_MultiFEC#UserCursorAdjust("Test",0) != 0)
		return -1
	endif
	killwindow Test
	DoWindow TempHistPlot
	if(V_Flag==1)
	killwindow TempHistPlot
	endif
	
	duplicate/free FreeRupPnts ForcesBack

	FreeRupPnts+=numpnts(ForceExtWave)
	CalculateForcesFromPoints(ForceWave,SepWave,FreeRupPnts,ForcesBack,500)
		FreeRupPnts-=numpnts(ForceExtWave)

	//ForcesBack=ForceRetSm[FreeRupPnts[p]]
	CalcAllLCs(ForceWave,SepWave,FreeRupPnts,offset,LCsBack)
	CalculateSlopes(ForceWave,SepWave,FreeRupPnts,offset,SlopesBack)
	AddNotes(ForceWave,SepWave,FreeRupPnts,ForcesBack,LCSBack,SLopesBack,FOffset=zeroforce)
	killwaves RupTimes
end

Static function CalculateForcesFromPoints(ForceWave,SepWave,PointWave,ForcesBack,backcalc)


	wave ForceWave,SepWave,PointWave,ForcesBack
	variable backcalc
	wave ExtForce=$ReplaceString("Ret",nameofwave(ForceWave),"Ext")
	variable correctPoints=numpnts(ExtForce)
	variable offset=str2num(stringbykey("DE_SchollOffset",note(ForceWave),":","\r"))
	FindLevels/P/Q SepWave, offset
	wave w_FindLevels
	variable SurfacePnt=DE_MultiFEC#DE_Median(W_FindLevels)
	variable n,prevpnt,CriticalTime
	duplicate/free PointWave TempForces
	variable tot=dimsize(PointWave,0)
	for(n=0;n<tot;n+=1)
		if(n==0)
			prevpnt=max(PointWave[n]-backcalc-correctPoints,SurfacePnt)

		else
			prevpnt=max(PointWave[n]-backcalc-correctPoints,PointWave[n-1]-correctPoints+10)
		endif
		make/o/n=0 HistOUt
	 
		duplicate/free/R=[prevpnt,PointWave[n]-correctPoints] ForceWave, FFit
		if(numpnts(FFit)>50)
		
		CurveFit/Q/W=2 line FFit
		CriticalTime=pnt2x(ForceWave,PointWave[n]-correctPoints)
		wave w_coef,w_sigma
		TempForces[n]= w_coef[1]*CriticalTime+w_coef[0]
		else
		TempForces[n]=ForceWave[PointWave[n]-correctPoints]
		endif
	
	endfor
	TempForces*=-1
	duplicate/o TempForces ForcesBack
	killwaves w_coef,w_sigma,W_FindLevels
	Killwaves HistOUt


end

Static Function AddNotes(ForceWave,SepWave,RupPnts,ForcesBack,LCSBack,SLopesBack,[Foffset])
	wave ForceWave,SepWave,ForcesBack,LCSBack,SLopesBack,RupPnts
	variable Foffset
	String RupForces=""
	String Lcs=""
	String Slopes=""
	String RupPntString=""
	variable tot=numpnts(ForcesBack)
	variable n
	for(n=0;n<tot;n+=1)
		RupPntString+=num2str(RupPnts[n])+";"
		RupForces+=num2str(ForcesBack[n])+";"
		Lcs+=num2str(LCSBack[n])+";"
		Slopes+=num2str(SLopesBack[n])+";"
	
	endfor
	
	String StartingNote=note(ForceWave)
	StartingNote=ReplaceStringbyKey("RupPnts",StartingNote,RupPntString,":","\r")

	StartingNote=ReplaceStringbyKey("RupForce",StartingNote,RupForces,":","\r")
	StartingNote=ReplaceStringbyKey("ContourLengths",StartingNote,Lcs,":","\r")
	StartingNote=ReplaceStringbyKey("Slopes",StartingNote,Slopes,":","\r")
	if(ParamisDefault(Foffset))
	else
		StartingNote=ReplaceStringbyKey("DE_FOff",StartingNote,num2str(Foffset),":","\r")

	endif
	note/K ForceWave, StartingNote
	note/K SepWave, StartingNote

end

Static Function AddZeroOffsetsToNotes()
	string AllForceRet= wavelist("*Force_Ret",";","")
	String ForceWaveList="",SepWaveList=""
	variable n
	variable bottom=0
	variable top=itemsinlist(AllFOrceRet)
	variable ZeroForce
	String StartingNote=""
	for(n=bottom;n<top;n+=1)
		wave ForceRetWave=$stringfromlist(n,AllForceRet)
		ZeroForce=TestZero(ForceRetWave)
		print nameofwave(ForceRetWave)+":"+num2str(ZeroForce)
		StartingNote=note(ForceRetWave)
		StartingNote=ReplaceStringbyKey("DE_FOff",StartingNote,num2str(ZeroForce),":","\r")
		note/K ForceRetWave, StartingNote

	endfor

end

Static Function CheckAllZeros()

	string AllForceRet= wavelist("*Force_Ret",";","")
	String ForceWaveList="",SepWaveList=""
	variable n
	variable bottom=0
	variable top=itemsinlist(AllFOrceRet)
		variable ZeroForce

	for(n=bottom;n<top;n+=1)
//		//for(n=0;n<itemsinlist(AllFOrceRet);n+=1)
		wave ForceRetWave=$stringfromlist(n,AllForceRet)
		ZeroForce=TestZero(ForceRetWave,distance=70e-9)
		if(abs(ZeroForce)>3e-12)
		print "Huge:"+nameofwave(ForceRetWave)+"_"+num2str(ZeroForce)
		elseif(abs(ZeroForce)>2e-12)
		print "Biggish:"+nameofwave(ForceRetWave)
				elseif(abs(ZeroForce)>1e-12)
		print "NonZero:"+nameofwave(ForceRetWave)
		endif

		
	endfor

end

Static Function ReFilterandCalculateAll()
	string AllForceRet= wavelist("*Force_Ret",";","")
	String ForceWaveList="",SepWaveList=""
	variable n
	variable bottom=0
	variable top=itemsinlist(AllFOrceRet)
	variable filtering
	wave FRetSm=$replacestring("Force_Ret", stringfromlist(0,wavelist("*Force_Ret",";","")),"FSm")

	
		if(cmpstr(stringbykey("DE_Filtering",note(FRetSm),":","\r"),"")==0)
			Prompt filtering, "Enter Filtering"
			DoPrompt "Enter Filtering" filtering
		else
			filtering=str2num(stringbykey("DE_Filtering",note(FRetSm),":","\r"))
		endif
	
	
	for(n=bottom;n<top;n+=1)
		wave ForceRetWave=$stringfromlist(n,AllForceRet)
		
		UpdatedRawWaves(ForceRetWave,filtering=filtering)
	endfor

end

Static Function UpdatedRawWaves(ForceRetWave,[filtering])
	wave ForceRetWave
	variable filtering


	wave ForceExtWave=$replacestring("Ret",nameofwave(ForceRetWave),"Ext")
	wave SepRetWave=$replacestring("Force",nameofwave(ForceRetWave),"Sep")
	wave SepExtWave=$replacestring("Force",nameofwave(ForceExtWave),"Sep")
	wave FRetSm=$replacestring("Force_ext",nameofwave(ForceExtWave),"FSm")
	wave SRetSm=$replacestring("Force_ext",nameofwave(ForceExtWave),"SSm")
	wave ThisEvent=$replacestring("Force_ext",nameofwave(ForceExtWave),"Starts")
	
	if(ParamisDefault(filtering))
	
		if(cmpstr(stringbykey("DE_Filtering",note(FRetSm),":","\r"),"")==0)
			Prompt filtering, "Enter Filtering"
			DoPrompt "Enter Filtering" filtering
		else
			filtering=str2num(stringbykey("DE_Filtering",note(FRetSm),":","\r"))
		endif
	endif
	
	if(filtering==1)

		duplicate/o ForceRetWave FRetSm
		duplicate/o SepRetWave SRetSm
	elseif(filtering>5)
		DE_Filtering#FilterForceSep(ForceRetWave,SepRetWave,FRetSm,SRetSm,"SVG",filtering)

	elseif(filtering<1)
		DE_Filtering#FilterForceSep(ForceRetWave,SepRetWave,FRetSm,SRetSm,"TVD",filtering)

	endif
	
	duplicate/o ThisEvent FreeRupPnts
	FreeRupPnts-=NUMPNTS(SepExtWave)
	string Yuck=nameofwave(ThisEvent)
	make/free/n=0 LCSBack,SLopesBack
	duplicate/free FreeRupPnts ForcesBack
//	ForcesBack=FRetSm[FreeRupPnts[p]]
		CalculateForcesFromPoints(ForceRetWave,SepRetWave,FreeRupPnts,ForcesBack,500)

	variable offset=-1*str2num(stringbykey("DE_SChollOffset",note(ForceRetWave),":","\r"))-5e-9
	CalcAllLCs(ForceRetWave,SepRetWave,FreeRupPnts,offset,LCsBack)
	CalculateSlopes(ForceRetWave,SepRetWave,FreeRupPnts,offset,SlopesBack)
	AddNotes(ForceRetWave,SepRetWave,FreeRupPnts,ForcesBack,LCSBack,SLopesBack)
	FreeRupPnts+=NUMPNTS(SepExtWave)
	//duplicate/o FreeRupPnts ThisEvent
	killwaves FreeRupPnts
		
		
end

Static Function UserCursorAdjust(graphName,autoAbortSecs)
	String graphName
	Variable autoAbortSecs
	DoWindow/F $graphName // Bring graph to front
	if (V_Flag == 0) // Verify that graph exists
		Abort "UserCursorAdjust: No such graph."
		return -1
	endif
	ShowInfo/W=Test
	NewPanel /K=2 /W=(187,368,637,531) as "Pause for Cursor"
	DoWindow/C tmp_PauseforCursor // Set to an unlikely name
	AutoPositionWindow/E/M=1/R=$graphName // Put panel near the graph
	DrawText 21,20,"Adjust the cursors and then"
	DrawText 21,40,"Click Continue."
	Button button0,pos={80,58},size={92,20},title="Continue"
	Button button0,proc=DE_MultiFEC#UserCursorAdjust_ContButtonProc
	
	PopupMenu pop0,pos={250,58},size={92,20},title="Garbage"
	PopupMenu pop0,proc=DE_MultiFEC#PopMenuProc,value= DE_MultiFEC#MakeStringList()
	
	Button button1,pos={250,88},size={92,20},title="Fix That"
	Button button1,proc=DE_MultiFEC#UpdateAPoint
	Button button2,pos={250,110},size={92,20},title="Delete That"
	Button button2,proc=DE_MultiFEC#DeleteButton
	Button button3,pos={250,135},size={92,20},title="Add Here"
	Button button3,proc=DE_MultiFEC#AddButton
	
	Button button4,pos={80,135},size={92,20},title="Resize"
	Button button4,proc=DE_MultiFEC#ResizeButton
	
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
	return didAbort
End

Static Function ResizeButton(ctrlName) : ButtonControl
	String ctrlName

	SetAxis/W=Test/A

end

Static Function UserCursorAdjust_ContButtonProc(ctrlName) : ButtonControl
	String ctrlName
	DoWindow/K tmp_PauseforCursor // Kill panel
End

Static Function UpdateAPoint(ctrlName) : ButtonControl
	String ctrlName
		wave FreeRupPnts

	string CursorString= CsrInfo(A,"Test")
	if(cmpstr(CursorString,"")==0)
		print "No Cursor Dummy"
	else
		Controlinfo/W=tmp_PauseforCursor  pop0
		variable index=V_Value-1
		variable newlocation=pcsr(A,"Test")

		FreeRupPnts[index]=newlocation
	ReCalcDependWaves()


	endif
End

Static Function DeleteButton(ctrlName) : ButtonControl
	String ctrlName

	wave FreeRupPnts
	controlinfo/W=tmp_PauseforCursor pop0
	deletepoints (V_Value-1),1, FreeRupPnts
	ReCalcDependWaves()
end

Static Function AddButton(ctrlName) : ButtonControl
	String ctrlName

	wave FreeRupPnts
	String StringCsr=CsrInfo(A,"Test")
	if(cmpstr(StringCsr,"")==0)

		print "No Cursor A Dummy"
		return -1
	endif

	variable spot=pcsr(A,"Test")
	variable wheretoinsert
	FindLevel/P/Q FreeRupPnts,spot

	if(V_Flag==1)
		if(spot<wavemin(FreeRupPnts))
			wheretoinsert=0
		elseif(spot>wavemax(FreeRupPnts))
			wheretoinsert=numpnts(FreeRupPnts)
	
		endif
	else
		wheretoinsert=ceil(v_levelx)
	endif
	InsertPoints wheretoinsert,1,FreeRupPnts
	FreeRupPnts[wheretoinsert]=spot
	ReCalcDependWaves()
end

Static Function/S MakeStringList()
	wave FreeRupPnts
	variable	Num=dimsize(FreeRupPnts,0)
	variable n
	string Result=""
	for(n=0;n<Num;n+=1)
		Result+=num2str(n)+";"
		
	endfor
	return Result


end

Static Function PopMenuProc(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa

	switch( pa.eventCode )
		case 2: // mouse up
			Variable popNum = pa.popNum-1
			String popStr = pa.popStr
			wave FreeRupPnts
			wave ForceWave=$stringfromlist(0,TraceNameList("Test",";",1))
			Wave SepWave=$replacestring("Force",nameofwave(ForceWave),"Sep")
			variable ThisHerePoint=FreeRupPnts[popNum]
			variable ThisHereSep=SepWave[ThisHerePoint]
			variable ForceMax=-1*ForceWave[ThisHerePoint]+25e-12
			variable ForceMin=-1*ForceWave[ThisHerePoint]-25e-12
			SetAxis/A=2/W=Test bottom ThisHereSep-10e-9,ThisHereSep+10e-9
			SetAxis/A=2/W=Test left ForceMin,ForceMax

			MakeSingleContoursAndDisplay(Forcewave,SepWave,FreeRupPnts,popNum)
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End


Static Function ReCalcDependWaves()
	wave FreeRupPnts
	wave FRetSm=$StringFromList(2,TraceNameList("Test",";",1))
	wave SRetSm=$ReplaceString("FSm",nameofwave(FRetSm),"SSm")
	wave FRet=$StringFromList(0,TraceNameList("Test",";",1))
	wave SRet=$ReplaceString("Force",nameofwave(FRet),"Sep")
	PopupMenu pop0,proc=DE_MultiFEC#PopMenuProc,value= DE_MultiFEC#MakeStringList()

	make/o/n=(numpnts(FreeRupPnts),2) RupTimes
	RupTimes[][0]=SRet[FreeRupPnts]
	RupTimes[][1]=FRet[FreeRupPnts[p]]
end

Static Function CorrectAllWavesForOffsetError()
string AllForceRet= wavelist("*Force_Ret",";","")
	String ForceWaveList="",SepWaveList=""
	variable n
	variable top=itemsinlist(AllFOrceRet)
	variable Entries
	String RupForces,Contours,Slopes
	variable ZeroForce
	make/o/n=(0,5) First,Second,Third,Fourth,Fifth
	for(n=1;n<top;n+=1)
		Wave ForceRetWave=$StringFromList(n,wavelist("*Force_Ret",";",""))
		wave ForceExtWave=$replacestring("Ret",nameofwave(ForceRetWave),"Ext")
		wave ThisEvent=$replacestring("Force_ext",nameofwave(ForceExtWave),"Starts")
	//	if(cmpstr(nameofwave(ForceRetWave),"Best500021Force_Ret")==0)
	//	return 0
	//	endif
		print nameofwave(ForceRetWave)
		ThisEvent+=numpnts(ForceExtWave)
		//wave SepExtWave
	endfor

end

Static Function ProcessForceCurves()
	string AllForceRet= wavelist("*Force_Ret",";","")
	String ForceWaveList="",SepWaveList=""
	variable n
	variable top=itemsinlist(AllFOrceRet)
	variable Entries
	String RupForces,Contours,Slopes
	variable ZeroForce
	make/o/n=(0,5) First,Second,Third,Fourth,Fifth
	for(n=0;n<top;n+=1)
		Wave ForceWave=$StringFromList(n,wavelist("*Force_Ret",";",""))
		Entries=CountEntries(ForceWave)
		RupForces=Stringbykey("RupForce",note(ForceWave),":","\r")
		Contours=Stringbykey("ContourLengths",note(ForceWave),":","\r")
		ZeroForce=str2num(Stringbykey("DE_FOff",note(ForceWave),":","\r"))

		Slopes=Stringbykey("Slopes",note(ForceWave),":","\r")
		if(Entries==4)
			InsertPoints/M=0 0,1, First,Third,Fourth,Fifth

			First[0][0]=n
			First[0][1]=str2num(Stringfromlist(0,RupForces))
			First[0][2]=First[0][1]-ZeroForce
			First[0][3]=str2num(Stringfromlist(0,Contours))
			First[0][4]=str2num(Stringfromlist(0,Slopes))
			

			
			Third[0][0]=n
			Third[0][1]=str2num(Stringfromlist(1,RupForces))
			Third[0][2]=Third[0][1]-ZeroForce
			Third[0][3]=str2num(Stringfromlist(1,Contours))
			Third[0][4]=str2num(Stringfromlist(1,Slopes))
			
			Fourth[0][0]=n
			Fourth[0][1]=str2num(Stringfromlist(2,RupForces))
			Fourth[0][2]=Fourth[0][1]-ZeroForce
			Fourth[0][3]=str2num(Stringfromlist(2,Contours))
			Fourth[0][4]=str2num(Stringfromlist(2,Slopes))
			
			Fifth[0][0]=n
			Fifth[0][1]=str2num(Stringfromlist(3,RupForces))
			Fifth[0][2]=Fifth[0][1]-ZeroForce
			Fifth[0][3]=str2num(Stringfromlist(3,Contours))
			Fifth[0][4]=str2num(Stringfromlist(3,Slopes))
		
		elseif(Entries==5)
		InsertPoints/M=0 0,1, First,Second,Third,Fourth,Fifth
			First[0][0]=n
			First[0][1]=str2num(Stringfromlist(0,RupForces))
			First[0][2]=First[0][1]-ZeroForce
			First[0][3]=str2num(Stringfromlist(0,Contours))
			First[0][4]=str2num(Stringfromlist(0,Slopes))
			
			Second[0][0]=n
			Second[0][1]=str2num(Stringfromlist(1,RupForces))
			Second[0][2]=Second[0][1]-ZeroForce
			Second[0][3]=str2num(Stringfromlist(1,Contours))
			Second[0][4]=str2num(Stringfromlist(1,Slopes))
			
			Third[0][0]=n
			Third[0][1]=str2num(Stringfromlist(2,RupForces))
			Third[0][2]=Third[0][1]-ZeroForce
			Third[0][3]=str2num(Stringfromlist(2,Contours))
			Third[0][4]=str2num(Stringfromlist(2,Slopes))
			
			Fourth[0][0]=n
			Fourth[0][1]=str2num(Stringfromlist(3,RupForces))
			Fourth[0][2]=Fourth[0][1]-ZeroForce
			Fourth[0][3]=str2num(Stringfromlist(3,Contours))
			Fourth[0][4]=str2num(Stringfromlist(3,Slopes))
		
			Fifth[0][0]=n
			Fifth[0][1]=str2num(Stringfromlist(4,RupForces))
			Fifth[0][2]=Fifth[0][1]-ZeroForce
			Fifth[0][3]=str2num(Stringfromlist(4,Contours))
			Fifth[0][4]=str2num(Stringfromlist(4,Slopes))

		
		else
		

		endif
	endfor
	
	
	//ExportWaveLists(ForceWaveList,SepWaveList)
end

Static Function CountEntries(ForceWave)
	wave ForceWave
	String ForceRups=Stringbykey("RupForce",note(ForceWave),":","\r")
	return itemsinlist(ForceRups)
end

Static Function ReturnRelevantLCs()
	wave First,Second,Third,Fourth,Fifth

	make/free/n=(dimsize(First,0)) Firstn,Thirdn
	FirstN=First[p][0]
	thirdN=Third[p][0]

	variable IntRupNum=dimsize(Second,0)
	variable TotalTraces=dimsize(First,0)

	make/o/n=(TotalTraces) LC1,LC3,LC4
	make/o/n=(IntRupNum) LC2first,Lc2Second

	variable n,m
	for(n=0;n<IntRupNum;n+=1)
		m=Second[n][0]
		FindValue/V=(m) FirstN
		LC2First[n]=Second[n][3]-First[v_value][3]
		FindValue/V=(m) ThirdN
		LC2Second[n]=Third[v_value][3]-Second[n][3]
	endfor

	for(n=0;n<TotalTraces;n+=1)
		LC1[n]=Third[n][2]-First[n][3]
		LC3[n]=Fourth[n][2]-Third[n][3]
		LC4[n]=Fifth[n][2]-Fourth[n][3]


	endfor

end

Static Function ReturnRupForces()
	wave First,Second,Third,Fourth,Fifth

	make/o/n=(dimsize(First,0)) Rup1
	Rup1=First[p][1]
		make/o/n=(dimsize(Second,0)) Rup2
	Rup2=Second[p][1]
		make/o/n=(dimsize(Third,0)) Rup3
	Rup3=Third[p][1]
		make/o/n=(dimsize(Fourth,0)) Rup4
	Rup4=Fourth[p][1]
		make/o/n=(dimsize(Fifth,0)) Rup5
	Rup5=Fifth[p][1]
	Rup5*=-1;Rup4*=-1;Rup3*=-1;Rup2*=-1;Rup1*=-1

end
Static Function ReturnRupForcesZero()
	wave First,Second,Third,Fourth,Fifth

	make/o/n=(dimsize(First,0)) ZRup1
	ZRup1=First[p][2]
		make/o/n=(dimsize(Second,0)) ZRup2
	ZRup2=Second[p][2]
		make/o/n=(dimsize(Third,0)) ZRup3
	ZRup3=Third[p][2]
		make/o/n=(dimsize(Fourth,0)) ZRup4
	ZRup4=Fourth[p][2]
		make/o/n=(dimsize(Fifth,0)) ZRup5
	ZRup5=Fifth[p][2]
	ZRup5*=-1;ZRup4*=-1;ZRup3*=-1;ZRup2*=-1;ZRup1*=-1

end

Static Function ReturnSlopes()
	wave First,Second,Third,Fourth,Fifth

	make/o/n=(dimsize(First,0)) Slope1
	Slope1=First[p][4]
		make/o/n=(dimsize(Second,0)) Slope2
	Slope2=Second[p][4]
		make/o/n=(dimsize(Third,0)) Slope3
	Slope3=Third[p][4]
		make/o/n=(dimsize(Fourth,0)) Slope4
	Slope4=Fourth[p][4]
		make/o/n=(dimsize(Fifth,0)) Slope5
	Slope5=Fifth[p][4]
	//Slope5*=-1;Slope4*=-1;Slope3*=-1;Slope2*=-1;Slope1*=-1


end

Static Function PlotOne(ForceWaveNumber)
	variable ForceWaveNumber
		string AllForceRet= wavelist("*Force_Ret",";","")
	Wave ForceRetWave=$stringfromlist(ForceWaveNumber,AllForceRet)
	wave ForceExtWave=$replacestring("Ret",nameofwave(ForceRetWave),"Ext")
	wave SepRetWave=$replacestring("Force",nameofwave(ForceRetWave),"Sep")
	wave SepExtWave=$replacestring("Force",nameofwave(ForceExtWave),"Sep")
	wave FRetSm=$replacestring("Force_ext",nameofwave(ForceExtWave),"FSm")
	wave SRetSm=$replacestring("Force_ext",nameofwave(ForceExtWave),"SSm")
	wave ThisEvent=$replacestring("Force_ext",nameofwave(ForceExtWave),"Starts")
	duplicate/free ThisEvent AdjPoints
	AdjPoints-=numpnts(SepExtWave)
	
	DoWindow PlotOne
	if(V_Flag==1)
	KillWindow PlotOne
	endif
	display/N=PlotOne ForceRetWave //vs SepRetWave
	Appendtograph/W=PlotOne FRetSm// vs SRetSm
	ModifyGraph/W=PlotOne rgb($Nameofwave(ForceRetWave))=(65535,49151,49151)
	
	
//	wave ForceWave,SepWave,PointWave,SlopesBack
//	variable Sepoff
//	variable tot=numpnts(PointWave)
//	variable n=0
//	variable prevpnt
	variable offset=-1*str2num(stringbykey("DE_SChollOffset",note(ForceRetWave),":","\r"))-5e-9

	variable backdist=3e-9
	variable backcalc=backdist/str2num(stringbykey("retractvelocity",note(ForceRetWave),":","\r"))/dimdelta(ForceRetWave,0)
	FindLevels/P/Q SepRetWave, -1*offset
	wave w_FindLevels
	variable SurfacePnt=DE_MultiFEC#DE_Median(W_FindLevels)
	variable n,prevpnt
//	duplicate/free PointWave TempSlopes
	for(n=0;n<numpnts(AdjPoints);n+=1)
		if(n==0)
				prevpnt=max(AdjPoints[n]-backcalc,SurfacePnt)

		else
			prevpnt=max(AdjPoints[n]-backcalc,AdjPoints[n-1]+10)
		endif
		make/o/n=0 HistOUt
		duplicate/free/R=[prevpnt,AdjPoints[n]] ForceRetWave, FFit
		duplicate/o FFit $("Fit"+num2str(n))
		CurveFit/Q/W=2 line FFit /D=$("Fit"+num2str(n))
		wave FitWave=$("Fit"+num2str(n))
		appendtograph/W=PlotOne FitWave
		ModifyGraph/W=PlotOne rgb($nameofwave(FitWave) )=(0,0,0)
		endfor
//		wave w_coef,w_sigma
//		TempSlopes[n]= w_coef[1]
//		
//	
//	endfor
//	TempSlopes*=-1
//	duplicate/o TempSlopes SlopesBack
//	killwaves w_coef,w_sigma,W_FindLevels
//	Killwaves HistOUt
	
	
end

//Static Function 	ExportWaveLists(ForceWaveList,SepWaveList)
//
//	string ForceWaveList,SepWaveList
//	variable n
//	display/n=TMP_D
//	for(n=0;n<itemsinlist(ForceWaveList);n+=1)
//		wave Forcewave=$(stringfromlist(n,FOrceWaveList))
//		wave Sepwave=$(stringfromlist(n,SepWaveList))
//		duplicate/o ForceWave $(replaceString("Force",nameofwave(ForceWave),"Time"))
//		wave TimeWave=$(replaceString("Force",nameofwave(ForceWave),"Time"))
//		TimeWave=pnt2x(ForceWave,p)
//		appendtograph/w=TMP_D Forcewave vs SepWave 
//		Appendtograph/w=TMP_D  TimeWave
//		
//
//	endfor
//	String Path="D:\Data\Feather\Hold.pxp"
//	SaveGraphCopy/o as Path
//	KillWindow TMP_D
//
//
//end