#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma modulename=DE_MultiFEC
#include "SimpleWLCPrograms"
#include "DE_Filtering"
#include "DE_NewFeather"

Function BatchProcess(Threshold,tau,[bottom,top])

	variable bottom,top,Threshold,tau
	string AllForceRet= wavelist("*Force_Ret",";","")
	String ForceWaveList="",SepWaveList=""
	if(ParamisDefault(Bottom))
		Bottom=0
	endif
	if(ParamisDefault(top))
		top=itemsinlist(AllFOrceRet)
	endif
	variable n
	string WaveNote
	//
	for(n=bottom;n<top;n+=1)
		print "N:"+num2str(n)
		wave ForceRetWave=$stringfromlist(n,AllForceRet)
		wave ForceExtWave=$replacestring("Ret",nameofwave(ForceRetWave),"Ext")
		wave SepRetWave=$replacestring("Force",nameofwave(ForceRetWave),"Sep")
		wave SepExtWave=$replacestring("Force",nameofwave(ForceExtWave),"Sep")
		make/free/n=0 ForceAll,SepAll
		Concatenate/NP/o {ForceExtWave,ForceRetWave},ForceAll
		Concatenate/NP/o {SepExtWave,SepRetWave},SepAll

		WaveNote=ReplaceStringbyKey("DwellTime",note(ForceAll),"0",":","\r")
		note/K ForceAll, WaveNote
		note/K SepAll, WaveNote
		//
		duplicate/o ForceAll $(Replacestring("Force_Ret",nameofwave(ForceRetWave),"Force"))
		wave ForceWave=$(Replacestring("Force_Ret",nameofwave(ForceRetWave),"Force"))
		duplicate/o SepAll $(Replacestring("force_Ret",nameofwave(ForceRetWave),"Sep"))
		wave SepWave=$(Replacestring("Force_Ret",nameofwave(ForceRetWave),"Sep"))
	endfor
	DE_NEwFeather#SaveOutAllWaves("*Force")
	make/free/n=2 OptionsWave
	OptionsWave={Threshold,tau}
	DE_NewFeather#RunFeatheronOutputFolder(OptionsWave)
	DE_NewFeather#LoadTheWaves("*Force")
	killwaves ForceaLl,SepAll
		
	for(n=bottom;n<top;n+=1)
		wave ForceRetWave=$stringfromlist(n,AllForceRet)
		wave ForceExtWave=$replacestring("Ret",nameofwave(ForceRetWave),"Ext")
		wave SepRetWave=$replacestring("Force",nameofwave(ForceRetWave),"Sep")
		wave SepExtWave=$replacestring("Force",nameofwave(ForceExtWave),"Sep")
		wave ForceWave=$(Replacestring("Force_Ret",nameofwave(ForceRetWave),"Force"))
		wave SepWave=$(Replacestring("Force_Ret",nameofwave(ForceRetWave),"Sep"))
		killwaves ForceWave,SepWave
	endfor
		
end


Static Function CleartheAllWave()

	string AllForceRet= wavelist("*Force_ret",";","")
	
	
	variable n
	variable	top=itemsinlist(AllFOrceRet)

//
	for(n=0;n<top;n+=1)
//		
		wave ForceRetWave=$stringfromlist(n,AllForceRet)
		wave ForceExtWave=$replacestring("Ret",nameofwave(ForceRetWave),"Ext")
		wave SepRetWave=$replacestring("Force",nameofwave(ForceRetWave),"Sep")
		wave SepExtWave=$replacestring("Force",nameofwave(ForceExtWave),"Sep")
		wave ForceWave=$(Replacestring("Force_Ret",nameofwave(ForceRetWave),"Force"))
		wave SepWave=$(Replacestring("Force_Ret",nameofwave(ForceRetWave),"Sep"))
		killwaves ForceWave,SepWave
	endfor


end

Static Function CleartheStarts()

	string AlLStarts= wavelist("*STarts",";","")
	
	
	variable n
	variable	top=itemsinlist(AlLStarts)

//
	for(n=0;n<top;n+=1)
	print "N:"+num2str(n)
//		
		wave StartWave=$stringfromlist(n,AlLStarts)
		
		killwaves StartWave
	endfor


end


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

		wave ThisEvent=$replacestring("Force_Ret",nameofwave(ForceRetWave),"Starts")
		make/o/n=0 $replacestring("Force_Ret",nameofwave(ForceRetWave),"FSm"),$replacestring("Force_Ret",nameofwave(ForceRetWave),"SSm")
		wave FRetSm=$replacestring("Force_Ret",nameofwave(ForceRetWave),"FSm")
		wave SRetSm=$replacestring("Force_Ret",nameofwave(ForceRetWave),"SSm")
		DE_Filtering#FilterForceSep(ForceRetWave,SepRetWave,FRetSm,SRetSm,"SVG",51)
		duplicate/o ThisEvent FreeRupPnts
		FreeRupPnts-=NUMPNTS(SepExtWave)
		MakeNicePlot(ForceRetWave,sEPRetWave,FRetSm,SRetSm,FreeRupPnts)
		
		make/free/n=0 LCSBack,SLopesBack
		duplicate/free FreeRupPnts ForcesBack
		//ForcesBack=FRetSm[FreeRupPnts[p]]

	FreeRupPnts+=numpnts(ForceExtWave)
				CalculateForcesFromPoints(ForceRetWave,SepRetWave,FreeRupPnts,ForcesBack)

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
		CalculateForcesFromPoints(ForceRetWave,SepRetWave,FreeRupPnts,ForcesBack)
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
	ProcessFCs_RLCCohDock()
	//ProcessForceCurves()
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

Static Function MakeSingleContoursAndDisplay(ForceRetwave,SepRetWave,PointWave,index,[Plot])
	
	wave ForceRetwave,SepRetWave,PointWave
	variable index,Plot
	if(ParamisDefault(Plot))
	
	Plot=0
	endif
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
	if(Plot==1)
	Display/W=(500,500,700,700)/N=TempHistPlot LCHistOut
	//AutoPositionWindow/E/M=0/R=Test TempHistPlot // Put panel near the graph
	TextBox/W=TempHistPlot/C/N=text0/F=0/A=RT "LC = "+num2str(result)
	endif
	return result

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
				if(prevpnt>PointWave[n]-10)
			prevpnt=PointWave[n]-10
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
	variable backdist=4e-9
		variable velocitysynch=str2num(stringbykey("VelocitySynch",note(ForceWave),":","\r"))
	variable velocity
	if(velocitysynch==1)
	
	 velocity=str2num(stringbykey("velocity",note(ForceWave),":","\r"))
	else
		 velocity=str2num(stringbykey("retractvelocity",note(ForceWave),":","\r"))

	
	endif
	variable backcalc=backdist/velocity/dimdelta(ForceWave,0)
		//backcalc=.05/dimdelta(ForceWave,0)

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
		
		if(prevpnt>PointWave[n]-10)
			prevpnt=PointWave[n]-10
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
	
Static Function ReturnWaveNamesForaNum(number,Trans)

	variable number,Trans
	variable prevpnt
	string AllForceRet= wavelist("*Force_Ret",";","")
	print "N:"+num2str(number)
	print stringfromlist(number,AllForceRet)
	wave ForceRet=$stringfromlist(number,AllForceRet)
	print replacestring("Ret",stringfromlist(number,AllForceRet),"Ext")
	wave ForceExt=$replacestring("Ret",stringfromlist(number,AllForceRet),"Ext")

	print replacestring("Force",stringfromlist(number,AllForceRet),"Sep")
	wave SepRet=$replacestring("Force",stringfromlist(number,AllForceRet),"Sep")

	print replacestring("Force",stringfromlist(number,AllForceRet),"Sep")
	print replacestring("Force",stringfromlist(number,AllForceRet),"Starts")
	wave PointWave=$replacestring("Force_ret",stringfromlist(number,AllForceRet),"Starts")
	variable offset=-1*str2num(stringbykey("DE_SChollOffset",note(ForceRet),":","\r"))-5e-9
	print "Offset: "+num2str(offset)
	FindLevels/P/Q SepRet, -1*offset
	wave w_FindLevels
	variable SurfacePnt=DE_MultiFEC#DE_Median(W_FindLevels)
	if(Trans==0)
		prevpnt=SurfacePnt
	else
		prevpnt=PointWave[Trans-1]+10-numpnts(ForceExt)//+max(100,PointWave[1]-PointWave[n-1]

	endif
	print "Startpnt: " +num2str(prevpnt)
	print "Endpnt: " +num2str(PointWave[Trans]-numpnts(ForceExt))


end
	
Static Function CalcCurrLC(ForceWave,SepWave,StartPnt,EndPnt,SepOff,HistOut)
	wave ForceWave,SepWave,HistOut
	variable StartPnt,EndPnt,SepOff

	duplicate/free/r=[startpnt,endpnt] FOrceWave TempForce,TempLC
	duplicate/free/r=[startpnt,endpnt] SepWave TempSep
	endpnt-=startpnt

	startpnt=0
	TempForce*=-1
	TempSep+=SepOff
	TempLC=DE_WLC#ContourTransform_Per_QMC(TempForce,TempSep,.4e-9,298)
	//variable lastLC=DE_WLC#ContourTransform_Per_QMC(TempForce[endpnt],TempSep[endpnt],.4e-9,298)
	//variable lastLC=DE_WLC#FindLC_Woodside(TempForce[endpnt],TempSep[endpnt],.4e-9,25e-9,298)
	//variable startLC=wavemin(TempLC)
	//variable startLC=TempSep[endpnt]
		variable startLC=max(TempSep[numpnts(TempSep)-1],wavemin(TempLC))-3e-9

	variable lastLC=min(wavemax(TempLC),2*TempSep[numpnts(TempSep)-1])+3e-9
		variable stepLC,numsteps

	if((lastlc-startlc)/49<1e-9)
	stepLC=1e-9
	numsteps=ceil((lastlc-startlc)/1e-9)
	else
	stepLC=1*(lastlc-startlc)/49
	numsteps=50
	endif
	
	make/o/n=(numsteps) TempHist 
	variable Q=numpnts(TempLC)
	Histogram/C/B={startLC,stepLC,numsteps} TempLC,TempHist;
	CurveFit/Q/W=2 gauss TempHist /d
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
	display/W=(50,50,1000,600 )/N=Test ForceWave vs SepWave
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
	CalculateForcesFromPoints(ForceWave,SepWave,FreeRupPnts,ForcesBack)
		FreeRupPnts-=numpnts(ForceExtWave)

	//ForcesBack=ForceRetSm[FreeRupPnts[p]]
	CalcAllLCs(ForceWave,SepWave,FreeRupPnts,offset,LCsBack)
	CalculateSlopes(ForceWave,SepWave,FreeRupPnts,offset,SlopesBack)
	AddNotes(ForceWave,SepWave,FreeRupPnts,ForcesBack,LCSBack,SLopesBack,FOffset=zeroforce)
	killwaves RupTimes
end

Static function CalculateForcesFromPoints(ForceWave,SepWave,PointWave,ForcesBack)


	wave ForceWave,SepWave,PointWave,ForcesBack
	wave ExtForce=$ReplaceString("Ret",nameofwave(ForceWave),"Ext")
	variable correctPoints=numpnts(ExtForce)
	variable offset=str2num(stringbykey("DE_SchollOffset",note(ForceWave),":","\r"))
	FindLevels/P/Q SepWave, offset
	wave w_FindLevels
	variable SurfacePnt=DE_MultiFEC#DE_Median(W_FindLevels)
	variable n,prevpnt,CriticalTime
	
	variable backdist=4e-9
	variable backcalc=backdist/str2num(stringbykey("retractvelocity",note(ForceWave),":","\r"))/dimdelta(ForceWave,0)
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
		
		CurveFit/Q/W=2 line FFit /D
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




Static Function ProcessFCs_RLCCohDock()
	string AllForceRet= wavelist("*Force_Ret",";","")
	String ForceWaveList="",SepWaveList=""
	variable n
	variable top=itemsinlist(AllFOrceRet)
	variable Entries
	String RupForces,Contours,Slopes
	variable ZeroForce
	make/o/n=(0,5) FirstRLC,SecondRLC,FirstDD,SecondDD,SoloCoh,FirstCoh,SecondCoh
	for(n=0;n<top;n+=1)
		Wave ForceWave=$StringFromList(n,wavelist("*Force_Ret",";",""))
		Entries=CountEntries(ForceWave)
		RupForces=Stringbykey("RupForce",note(ForceWave),":","\r")
		Contours=Stringbykey("ContourLengths",note(ForceWave),":","\r")
		ZeroForce=str2num(Stringbykey("DE_FOff",note(ForceWave),":","\r"))
		Slopes=Stringbykey("Slopes",note(ForceWave),":","\r")
		Switch(Entries)
			case 6:
			
				AddanEntrytoWave(FirstRLC,n,0,ZeroForce,RupForces,Contours,Slopes)
				AddanEntrytoWave(SecondRLC,n,1,ZeroForce,RupForces,Contours,Slopes)
				AddanEntrytoWave(FirstDD,n,2,ZeroForce,RupForces,Contours,Slopes)
				AddanEntrytoWave(SecondDD,n,3,ZeroForce,RupForces,Contours,Slopes)
				AddanEntrytoWave(FirstCoh,n,4,ZeroForce,RupForces,Contours,Slopes)
				AddanEntrytoWave(SecondCoh,n,5,ZeroForce,RupForces,Contours,Slopes)

			
				break
			
			case 5:
			
				if(str2num(Stringfromlist(4,RupForces))>str2num(Stringfromlist(3,RupForces)))
					AddanEntrytoWave(FirstRLC,n,0,ZeroForce,RupForces,Contours,Slopes)
					AddanEntrytoWave(SecondRLC,n,1,ZeroForce,RupForces,Contours,Slopes)
					AddanEntrytoWave(FirstDD,n,2,ZeroForce,RupForces,Contours,Slopes)
					AddanEntrytoWave(SecondDD,n,3,ZeroForce,RupForces,Contours,Slopes)
					AddanEntrytoWave(SoloCoh,n,4,ZeroForce,RupForces,Contours,Slopes)

				
				
				
				else
					AddanEntrytoWave(FirstRLC,n,0,ZeroForce,RupForces,Contours,Slopes)
					AddanEntrytoWave(FirstDD,n,1,ZeroForce,RupForces,Contours,Slopes)
					AddanEntrytoWave(SecondDD,n,2,ZeroForce,RupForces,Contours,Slopes)
					AddanEntrytoWave(FirstCoh,n,3,ZeroForce,RupForces,Contours,Slopes)
					AddanEntrytoWave(SecondCoh,n,4,ZeroForce,RupForces,Contours,Slopes)
				
				endif
				break
			
			case 4:
				AddanEntrytoWave(FirstRLC,n,0,ZeroForce,RupForces,Contours,Slopes)
				AddanEntrytoWave(FirstDD,n,1,ZeroForce,RupForces,Contours,Slopes)
				AddanEntrytoWave(SecondDD,n,2,ZeroForce,RupForces,Contours,Slopes)
				AddanEntrytoWave(SoloCoh,n,3,ZeroForce,RupForces,Contours,Slopes)
				
				break
			
			default :
				break
	
		endswitch 

	endfor
	
	
end

Static Function AddanEntrytoWave(WaveIn,index,location,ZeroForce,RupFStr,ContourStr,SlopesStr)

	wave WaveIn
	variable index,location,ZeroForce
	String RupFStr,ContourStr,SlopesStr
	InsertPoints/M=0 0,1, WaveIn
	
	WaveIn[0][0]=index
	WaveIn[0][1]=str2num(Stringfromlist(location,RupFStr))
	WaveIn[0][2]=WaveIn[0][1]+ZeroForce
	WaveIn[0][3]=str2num(Stringfromlist(location,ContourStr))
	WaveIn[0][4]=str2num(Stringfromlist(location,SlopesStr))
	
end

Static Function CountEntries(ForceWave)
	wave ForceWave
	String ForceRups=Stringbykey("RupForce",note(ForceWave),":","\r")
	return itemsinlist(ForceRups)
end

Static Function ReturnRelevantLCs()
	
	
	variable n,CurrentNum
	wave FirstRLC,SecondRLC,FirstDD,SecondDD,SoloCoh,FirstCoh,SecondCoh
	make/o/n=(dimsize(FirstRLC,0)) LC_RLCTotal,LC_DDOne,LC_DDTwo
	make/o/n=(dimsize(FirstCoh,0)) LC_CoH
	make/o/n=(dimsize(SecondRLC,0)) LC_RLC1,LC_RLC2


	LC_RLCTotal=FirstDD[p][3]-FirstRLC[p][3]
	LC_DDOne=SecondDD[p][3]-FirstDD[p][3]
	LC_CoH=SecondCoh[p][3]-SoloCoh[p][3]
	make/free/n=(dimsize(SoloCoh,0)) SoloCohNums
	make/free/n=(dimsize(FirstCoh,0)) FirstCohNums

	FirstCohNums=FirstCoh[p][0]
	SoloCohNums=SoloCoh[p][0]

	for(n=0;n<dimsize(SecondDD,0);n+=1)
		CurrentNum=SecondDD[n][0]
		FindValue/V=(CurrentNum) SoloCohNums
		if(v_value!=-1)
			LC_DDTwo[n]=SoloCoh[v_value][3]-SecondDD[n][3]
		
		else
			FindValue/V=(CurrentNum) FirstCohNums
			LC_DDTwo[n]=FirstCoh[v_value][3]-SecondDD[n][3]

		
		endif

	endfor
	
	make/free/n=(dimsize(FirstRLC,0)) FirstRLCNums
	make/free/n=(dimsize(FirstDD,0)) FirstDDNums

	FirstRLCNums=FirstRLC[p][0]
	FirstDDNums=FirstDD[p][0]

	for(n=0;n<dimsize(SecondRLC,0);n+=1)
		CurrentNum=SecondRLC[n][0]
	
		FindValue/V=(CurrentNum) FirstRLCNums
		LC_RLC1[n]=SecondRLC[n][3]-FirstRLC[v_value][3]
		FindValue/V=(CurrentNum) FirstDDNums
		LC_RLC2[n]=FirstDD[n][3]-SecondRLC[n][3]
		
		
	endfor

end

Static Function ReturnRupForces()
	
	wave FirstRLC,SecondRLC,FirstDD,SecondDD,SoloCoh,FirstCoh,SecondCoh

	make/o/n=(dimsize(FirstRLC,0)) Rup_FirstRLC
	Rup_FirstRLC=FirstRLC[p][1]
	make/o/n=(dimsize(SecondRLC,0)) Rup_SecondRLC
	Rup_SecondRLC=SecondRLC[p][1]
	make/o/n=(dimsize(FirstDD,0)) Rup_FirstDD
	Rup_FirstDD=FirstDD[p][1]
	make/o/n=(dimsize(SecondDD,0)) Rup_SecondDD
	Rup_SecondDD=SecondDD[p][1]
	make/o/n=(dimsize(SoloCoh,0)) Rup_SoloCoh
	Rup_SoloCoh=SoloCoh[p][1]
	make/o/n=(dimsize(FirstCoh,0)) Rup_FirstCoh
	Rup_FirstCoh=FirstCoh[p][1]	
	make/o/n=(dimsize(SecondCoh,0)) Rup_SecondCoh
	Rup_SecondCoh=SecondCoh[p][1]
	
	Rup_FirstRLC*=-1;Rup_SecondRLC*=-1;Rup_FirstDD*=-1;
	Rup_SecondDD*=-1;Rup_SoloCoh*=-1;Rup_FirstCoh*=-1;Rup_SecondCoh*=-1;
	
	
//	wave First,Second,Third,Fourth,Fifth
//
//	make/o/n=(dimsize(First,0)) Rup1
//	Rup1=First[p][1]
//	make/o/n=(dimsize(Second,0)) Rup2
//	Rup2=Second[p][1]
//	make/o/n=(dimsize(Third,0)) Rup3
//	Rup3=Third[p][1]
//	make/o/n=(dimsize(Fourth,0)) Rup4
//	Rup4=Fourth[p][1]
//	make/o/n=(dimsize(Fifth,0)) Rup5
//	Rup5=Fifth[p][1]
//	Rup5*=-1;Rup4*=-1;Rup3*=-1;Rup2*=-1;Rup1*=-1

end

Static Function ReturnRupForcesZero()


	wave FirstRLC,SecondRLC,FirstDD,SecondDD,SoloCoh,FirstCoh,SecondCoh

	make/o/n=(dimsize(FirstRLC,0)) zRup_FirstRLC
	zRup_FirstRLC=FirstRLC[p][2]
	make/o/n=(dimsize(SecondRLC,0)) ZRup_SecondRLC
	ZRup_SecondRLC=SecondRLC[p][2]
	make/o/n=(dimsize(FirstDD,0)) ZRup_FirstDD
	ZRup_FirstDD=FirstDD[p][2]
	make/o/n=(dimsize(SecondDD,0)) ZRup_SecondDD
	ZRup_SecondDD=SecondDD[p][2]
	make/o/n=(dimsize(SoloCoh,0)) ZRup_SoloCoh
	ZRup_SoloCoh=SoloCoh[p][2]
	make/o/n=(dimsize(FirstCoh,0)) ZRup_FirstCoh
	zRup_FirstCoh=FirstCoh[p][2]	
	make/o/n=(dimsize(SecondCoh,0)) zRup_SecondCoh
	zRup_SecondCoh=SecondCoh[p][2]
	
//	ZRup_FirstRLC*=-1;ZRup_SecondRLC*=-1;ZRup_FirstDD*=-1;
//	ZRup_SecondDD*=-1;ZRup_SoloCoh*=-1;ZRup_FirstCoh*=-1;ZRup_SecondCoh*=-1;	
//	
//	wave First,Second,Third,Fourth,Fifth
//
//	make/o/n=(dimsize(First,0)) ZRup1
//	ZRup1=First[p][2]
//		make/o/n=(dimsize(Second,0)) ZRup2
//	ZRup2=Second[p][2]
//		make/o/n=(dimsize(Third,0)) ZRup3
//	ZRup3=Third[p][2]
//		make/o/n=(dimsize(Fourth,0)) ZRup4
//	ZRup4=Fourth[p][2]
//		make/o/n=(dimsize(Fifth,0)) ZRup5
//	ZRup5=Fifth[p][2]
//	ZRup5*=-1;ZRup4*=-1;ZRup3*=-1;ZRup2*=-1;ZRup1*=-1

end

Static Function ReturnSlopes()
	
	wave FirstRLC,SecondRLC,FirstDD,SecondDD,SoloCoh,FirstCoh,SecondCoh

	make/o/n=(dimsize(FirstRLC,0)) Slope_FirstRLC
	Slope_FirstRLC=FirstRLC[p][4]
	make/o/n=(dimsize(SecondRLC,0)) Slope_SecondRLC
	Slope_SecondRLC=SecondRLC[p][4]
	make/o/n=(dimsize(FirstDD,0)) Slope_FirstDD
	Slope_FirstDD=FirstDD[p][4]
	make/o/n=(dimsize(SecondDD,0)) Slope_SecondDD
	Slope_SecondDD=SecondDD[p][4]
	make/o/n=(dimsize(SoloCoh,0)) Slope_SoloCoh
	Slope_SoloCoh=SoloCoh[p][4]
	make/o/n=(dimsize(FirstCoh,0)) Slope_FirstCoh
	Slope_FirstCoh=FirstCoh[p][4]	
	make/o/n=(dimsize(SecondCoh,0)) Slope_SecondCoh
	Slope_SecondCoh=SecondCoh[p][4]

	
	
//	wave First,Second,Third,Fourth,Fifth
//
//	make/o/n=(dimsize(First,0)) Slope1
//	Slope1=First[p][4]
//		make/o/n=(dimsize(Second,0)) Slope2
//	Slope2=Second[p][4]
//		make/o/n=(dimsize(Third,0)) Slope3
//	Slope3=Third[p][4]
//		make/o/n=(dimsize(Fourth,0)) Slope4
//	Slope4=Fourth[p][4]
//		make/o/n=(dimsize(Fifth,0)) Slope5
//	Slope5=Fifth[p][4]
//	//Slope5*=-1;Slope4*=-1;Slope3*=-1;Slope2*=-1;Slope1*=-1


end

Static Function PlotOne(ForceWaveNumber,[Trans])
	variable ForceWaveNumber,Trans
	string AllForceRet= wavelist("*Force_Ret",";","")
	variable n,prevpnt
	Wave ForceRetWave=$stringfromlist(ForceWaveNumber,AllForceRet)
	wave ForceExtWave=$replacestring("Ret",nameofwave(ForceRetWave),"Ext")
	wave SepRetWave=$replacestring("Force",nameofwave(ForceRetWave),"Sep")
	wave SepExtWave=$replacestring("Force",nameofwave(ForceExtWave),"Sep")
	wave FRetSm=$replacestring("Force_ext",nameofwave(ForceExtWave),"FSm")
	wave SRetSm=$replacestring("Force_ext",nameofwave(ForceExtWave),"SSm")
	wave ThisEvent=$replacestring("Force_ext",nameofwave(ForceExtWave),"Starts")
	duplicate/free ThisEvent AdjPoints
	AdjPoints-=numpnts(SepExtWave)
	String WindowName
	if(ParamisDefault(Trans))
	WindowName="PlotWide"
	else
	WindowName="PlotOne"
	endif
	
	DoWindow $WindowName
	if(V_Flag==1)
	KillWindow $WindowName
	endif
	
	if(numpnts(ThisEvent)<3||Trans>=numpnts(ThisEvent))

	return 0
	endif
	
	Variable ZeroForce=str2num(Stringbykey("DE_FOff",note(ForceRetWave),":","\r"))
	Duplicate/o ForceRetWave ForceRetZero
	ForceRetZero=ZeroForce
	
	display/N=$WindowName ForceRetWave
	Appendtograph/W=$WindowName FRetSm
	Appendtograph/W=$WindowName ForceRetZero
	ModifyGraph/W=$WindowName rgb($Nameofwave(ForceRetWave))=(65535,49151,49151)
	ModifyGraph/W=$WindowName rgb(ForceRetZero)=(0,0,0)
	variable offset=-1*str2num(stringbykey("DE_SChollOffset",note(ForceRetWave),":","\r"))-5e-9

	variable backdist=4e-9
	variable velocitysynch=str2num(stringbykey("VelocitySynch",note(ForceRetWave),":","\r"))
	variable velocity
	if(velocitysynch==1)
	
	 velocity=str2num(stringbykey("velocity",note(ForceRetWave),":","\r"))
	else
		 velocity=str2num(stringbykey("retractvelocity",note(ForceRetWave),":","\r"))

	
	endif
	variable backcalc=backdist/velocity/dimdelta(ForceRetWave,0)
	variable backtime=backdist/velocity
	
	FindLevels/P/Q SepRetWave, -1*offset
	wave w_FindLevels
	variable SurfacePnt=DE_MultiFEC#DE_Median(W_FindLevels)
	for(n=0;n<numpnts(AdjPoints);n+=1)
		if(n==0)
			prevpnt=max(AdjPoints[n]-backcalc,SurfacePnt)

		else
			prevpnt=max(AdjPoints[n]-backcalc,AdjPoints[n-1]+25)
		endif
		make/o/n=0 HistOUt
		duplicate/free/R=[prevpnt,AdjPoints[n]] ForceRetWave, FFit
		duplicate/o FFit $("Fit"+num2str(n))
		CurveFit/Q/W=2 line FFit /D=$("Fit"+num2str(n))
		wave W_coef
		wave FitWave=$("Fit"+num2str(n))
		appendtograph/W=$WindowName FitWave
		ModifyGraph/W=$WindowName rgb($nameofwave(FitWave) )=(0,0,0)
		Execute "LinLinPlot_1()"
		ModifyGraph/W=$WindowName fSize=9
		endfor

	if(ParamisDefault(Trans))
		ModifyGraph/W=$WindowName width=300,height=200
		MoveWindow/W=$WindowName 700,50,800,300
	else
	for(n=0;n<numpnts(AdjPoints);n+=1)
		if(n==Trans)

		else
		ModifyGraph hideTrace($("Fit"+num2str(n)))=1
		
		endif
		
		endfor
		ThisEvent-=numpnts(ForceExtWave)
		MakeSingleContoursAndDisplay(ForceRetWave,SepRetWave,ThisEvent,Trans,Plot=1)
				ThisEvent+=numpnts(ForceExtWave)

		wave ThisFit=$("Fit"+num2str(Trans))
		print "Slope: "+ num2str((ThisFit[numpnts(Thisfit)-1]-ThisFit[0])/(pnt2x(Thisfit,numpnts(Thisfit)-1)-pnt2x(Thisfit,0)))
		SetAxis/W=$WindowName bottom pnt2x(ForceRetWave,AdjPoints[trans])-3*backtime,pnt2x(ForceRetWave,AdjPoints[trans])+backtime
		SetAxis/W=$WindowName/A=2 left
		print "Force: "+num2str(str2num((stringfromlist(Trans,Stringbykey("RupForce",note(ForceRetWave),":","\r"))))+str2num(Stringbykey("DE_FOff",note(ForceRetWave),":","\r")))
		make/o/n=1 RupPnt
		RupPnt=(-str2num(stringfromlist(Trans,Stringbykey("RupForce",note(ForceRetWave),":","\r"))))//-str2num(Stringbykey("DE_FOff",note(ForceRetWave),":","\r"))
		SetScale/P x pnt2x(ForceRetWave,str2num(stringfromlist(Trans,Stringbykey("RupPnts",note(ForceRetWave),":","\r")))),0.0001,"s", RupPnt

		Appendtograph/W=$WindowName RupPnt
		ModifyGraph/W=$WindowName mode(RupPnt)=3,marker(RupPnt)=16,rgb(RupPnt)=(29524,1,58982),useMrkStrokeRGB(RupPnt)=1
		ModifyGraph/W=$WindowName width=300,height=200
		MoveWindow/W=$WindowName 350,50,800,300

	endif
	
end
Window ViewFECs() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel/N=ViewFECs /W=(150,77,450,277)
	SetVariable DE_ViewFECs0,pos={30.00,40.00},size={189.00,18.00},title="TraceNumber"
	SetVariable DE_ViewFECs0,proc=DE_MultiFEC#TraceSVA,value=_NUM:0

	SetVariable DE_ViewFECs1,pos={33.00,66.00},size={183.00,18.00},title="RupNumber"
	SetVariable DE_ViewFECs1,proc=DE_MultiFEC#TraceSVA,value=_NUM:0

EndMacro

Static Function TraceSVA(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
	
	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
//AutoPositionWindow/E/M=1/R=$graphName
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
				StrSwitch(sva.ctrlName)
				case "DE_ViewFECs0":
					Controlinfo/W=ViewFEcs DE_ViewFECs1
					•DE_MultiFEC#PlotOne(dval)
					•DE_MultiFEC#PlotOne(dval,Trans=V_Value)
				break
				
				case "DE_ViewFECs1":
				
					Controlinfo/W=ViewFEcs DE_ViewFECs0
					PlotOne(V_Value,Trans=dval)
				break
				default:
				break
				endswitch
				
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Static Function FakeApproach(InputForceWave,InputSepWave,Ratio,OutputForceWave,OutputSepWave)

	wave InputForceWave,InputSepWave,OutputForceWave,OutputSepWave
	variable Ratio
	variable Total=Ratio*numpnts(InputForceWave)
	make/free/n=(Total) FreeForce,FreeSep
	Interpolate2/T=1/N=(Total)/Y=FreeForce InputForceWave
	Interpolate2/T=1/N=(Total)/Y=FreeSep InputSepWave
	duplicate/o FreeForce OutputForceWave
	duplicate/o FreeSep OutputSepWave
	note/K OutputForceWave,note (InputForceWave)
		note/K OutputSepWave,note (InputSepWave)

end
//
//
//Function iterateUnit(Number)
//variable Number
//	string AllForceRet= wavelist("*Force_Ret",";","")
//	String ForceWaveList="",SepWaveList=""
//	variable n,i,j
//	variable tolnum=20
//variable Timenum=20
//	variable bottom=Number
//	variable top=Number+1//itemsinlist(AllFOrceRet)
//	make/o/n=(tolnum) TolValues
//	make/o/n=(Timenum) TimeValues
//	make/o/n=(tolnum,Timenum) Grid
//	TolValues=1e-7*(3^p)
//	TimeValues=1e-6*(3^p)
//	for(n=bottom;n<top;n+=1)
//		
//		Grid[][]=UnitTest(Number,1,TolValues[p],TimeValues[q],0)
//	endfor
//
//
//
//end
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


//Static Function MakeSlopesAndAddtoPlot(ForceRetwave,SepRetWave,ForceRetSm,SepRetWaveSm,FreeRupPnts,offset,backcalc,SlopesBack)
//	wave ForceRetwave,SepRetWave,ForceRetSm,SepRetWaveSm,FreeRupPnts,SlopesBack
//	variable offset,backcalc
//		FindLevels/P/Q SepRetWave, -1*offset
//	wave w_FindLevels
//	variable SurfacePnt=DE_MultiFEC#DE_Median(W_FindLevels)
//	make/free/n=0 TempSLopesBack
//	variable tot=numpnts(FreeRupPnts)
//	variable n=0
//	variable prevpnt
//	duplicate/free FreeRupPnts TempSlopes
//	for(n=0;n<tot;n+=1)
//		if(n==0)
//				prevpnt=Max(FreeRupPnts[n]-backcalc,SurfacePnt+10)
//
//		else
//				prevpnt=Max(FreeRupPnts[n]-backcalc,FreeRupPnts[n-1]+10)
//
//		endif
//		
//		duplicate/free/R=[prevpnt,FreeRupPnts[n]] ForceRetwave, FFit 
//		duplicate/o FFit $("FSlopefit_"+num2str(n))
//		CurveFit/Q/W=2 line FFit/D=$("FSlopefit_"+num2str(n))
//		wave w_coef,w_sigma
//		TempSlopes[n]= w_coef[1]
//		duplicate/free/R=[prevpnt,FreeRupPnts[n]] SepRetWave, SFit
//	//	duplicate/o FFit $("FSlopefit_"+num2str(n))
//	//	duplicate/o SFit $("SSlopefit_"+num2str(n))
////
//	//	CurveFit/Q/W=2 line FFit /X=SFit /D=$("FSlopefit_"+num2str(n))
////
////		
////	
//	endfor
//	duplicate/o TempSlopes SlopesBack
//	killwaves w_coef,w_sigma,w_FindLevels
//		
//
//end

//Static Function CorrectAllWavesForOffsetError()
//string AllForceRet= wavelist("*Force_Ret",";","")
//	String ForceWaveList="",SepWaveList=""
//	variable n
//	variable top=itemsinlist(AllFOrceRet)
//	variable Entries
//	String RupForces,Contours,Slopes
//	variable ZeroForce
//	make/o/n=(0,5) First,Second,Third,Fourth,Fifth
//	for(n=1;n<top;n+=1)
//		Wave ForceRetWave=$StringFromList(n,wavelist("*Force_Ret",";",""))
//		wave ForceExtWave=$replacestring("Ret",nameofwave(ForceRetWave),"Ext")
//		wave ThisEvent=$replacestring("Force_ext",nameofwave(ForceExtWave),"Starts")
//	//	if(cmpstr(nameofwave(ForceRetWave),"Best500021Force_Ret")==0)
//	//	return 0
//	//	endif
//		print nameofwave(ForceRetWave)
//		ThisEvent+=numpnts(ForceExtWave)
//		//wave SepExtWave
//	endfor
//
//end

//Function FindAllForces(filtering,[bottom,top])
//	variable filtering,bottom,top
//	string AllForceRet= wavelist("*Force_Ret",";","")
//	String ForceWaveList="",SepWaveList=""
//	if(ParamisDefault(Bottom))
//	Bottom=0
//	endif
//	if(ParamisDefault(top))
//	top=itemsinlist(AllFOrceRet)
//	endif
//	variable n
//
//	for(n=bottom;n<top;n+=1)
//	print "N:"+num2str(n)
//		//for(n=0;n<itemsinlist(AllFOrceRet);n+=1)
//		wave ForceRetWave=$stringfromlist(n,AllForceRet)
//		wave ForceExtWave=$replacestring("Ret",nameofwave(ForceRetWave),"Ext")
//		wave SepRetWave=$replacestring("Force",nameofwave(ForceRetWave),"Sep")
//		wave SepExtWave=$replacestring("Force",nameofwave(ForceExtWave),"Sep")
//		make/free/n=0 ForceAll,SepAll
//		Concatenate/o {ForceExtWave,ForceRetWave},ForceAll
//		Concatenate/o {SepExtWave,SepRetWave},SepAll
//
//		make/free/n=0 FSm,SSm
//		make/o/n=0 FRetSm,SRetSm
//		if(filtering==1)
//			duplicate/free ForceAll FSM
//			duplicate/free SepAll SSM
//			duplicate/o ForceRetWave FRetSm
//			duplicate/o SepRetWave SRetSm
//		elseif(filtering>5)
//			DE_Filtering#FilterForceSep(ForceAll,SepAll,FSm,SSm,"SVG",filtering)
//			DE_Filtering#FilterForceSep(ForceRetWave,SepRetWave,FRetSm,SRetSm,"SVG",filtering)
//
//		elseif(filtering<1)
//			DE_Filtering#FilterForceSep(ForceAll,SepAll,FSm,SSm,"TVD",filtering)
//			DE_Filtering#FilterForceSep(ForceRetWave,SepRetWave,FRetSm,SRetSm,"TVD",filtering)
//
//		
//		endif
//
//		duplicate/o FSm $(Replacestring("Force_Ret",nameofwave(ForceRetWave),"Force"))
//		wave ForceWave=$(Replacestring("Force_Ret",nameofwave(ForceRetWave),"Force"))
//		duplicate/o SSm $(Replacestring("force_Ret",nameofwave(ForceRetWave),"Sep"))
//		wave SepWave=$(Replacestring("Force_Ret",nameofwave(ForceRetWave),"Sep"))
//		make/free/n=5 OptionsWave
//		OptionsWave={1e-5,1e-3,str2num(StringbyKey("TriggerTime",note(ForceWave),":","\r")),0.0,str2num(StringbyKey("SpringConstant",note(ForceWave),":","\r"))}
//		DE_NewFeather#OutportForce(ForceWave,SepWave)
//		DE_NewFeather#RunFeatheronOutput(OptionsWave)
//		wave event_starts
//
//		duplicate/o event_starts $replacestring("Force_ext",nameofwave(ForceExtWave),"Starts")
//		duplicate/o FRetSm $replacestring("Force_ext",nameofwave(ForceExtWave),"FSm")
//		duplicate/o SRetSm $replacestring("Force_ext",nameofwave(ForceExtWave),"SSm")
//		killwaves event_starts,ForceAll,SepAll,FretSm,SRetSm
//		
//	endfor
//
//
//end
//
//Function UnitTest(Number,filtering,tol,temporal,compensate)
//	variable Number,filtering,tol,temporal,compensate
//	string AllForceRet= wavelist("*Force_Ret",";","")
//	String ForceWaveList="",SepWaveList=""
//	variable n,approachvelocity,retractvelocity,ratio
//	variable bottom=Number
//	variable top=Number+1//itemsinlist(AllFOrceRet)
//	for(n=bottom;n<top;n+=1)
//		//for(n=0;n<itemsinlist(AllFOrceRet);n+=1)
//		wave ForceRetWave=$stringfromlist(n,AllForceRet)
//		wave ForceExtWave=$replacestring("Ret",nameofwave(ForceRetWave),"Ext")
//		wave SepRetWave=$replacestring("Force",nameofwave(ForceRetWave),"Sep")
//		wave SepExtWave=$replacestring("Force",nameofwave(ForceExtWave),"Sep")
//		approachvelocity=str2num(stringbykey("ApproachVelocity",note(ForceRetWave),":","\r"))
//		retractvelocity=str2num(stringbykey("RetractVelocity",note(ForceRetWave),":","\r"))
//		ratio=approachvelocity/retractvelocity
//		if(compensate==1)
//		
//			print ratio
//			make/o/n=0  OutputForceWave,OutputSepWave
//			FakeApproach(ForceExtWave,SepExtWave,Ratio,OutputForceWave,OutputSepWave)
//			SetScale/P x 0,(dimdelta(ForceRetWave,0)),"s", OutputForceWave,OutputSepWave
//			make/free/n=0 ForceAll,SepAll
//			Concatenate/o {OutputForceWave,ForceRetWave},ForceAll
//			Concatenate/o {OutputSepWave,SepRetWave},SepAll
//		else
//			Concatenate/o {ForceExtWave,ForceRetWave},ForceAll
//			Concatenate/o {SepExtWave,SepRetWave},SepAll
//		endif
//	
//		make/free/n=0 FSm,SSm
//		make/o/n=0 FRetSm,SRetSm
//		if(filtering==1)
//			duplicate/free ForceAll FSM
//			duplicate/free SepAll SSM
//			duplicate/o ForceRetWave FRetSm
//			duplicate/o SepRetWave SRetSm
//		elseif(filtering>5)
//			DE_Filtering#FilterForceSep(ForceAll,SepAll,FSm,SSm,"SVG",filtering)
//			DE_Filtering#FilterForceSep(ForceRetWave,SepRetWave,FRetSm,SRetSm,"SVG",filtering)
//
//		elseif(filtering<1)
//			DE_Filtering#FilterForceSep(ForceAll,SepAll,FSm,SSm,"TVD",filtering)
//			DE_Filtering#FilterForceSep(ForceRetWave,SepRetWave,FRetSm,SRetSm,"TVD",filtering)
//
//		
//		endif
//
//		duplicate/o FSm $(Replacestring("Force_Ret",nameofwave(ForceRetWave),"Force"))
//		wave ForceWave=$(Replacestring("Force_Ret",nameofwave(ForceRetWave),"Force"))
//		duplicate/o SSm $(Replacestring("force_Ret",nameofwave(ForceRetWave),"Sep"))
//		wave SepWave=$(Replacestring("Force_Ret",nameofwave(ForceRetWave),"Sep"))
//		make/free/n=5 OptionsWave
//		OptionsWave={tol,temporal,str2num(StringbyKey("TriggerTime",note(ForceWave),":","\r")),0.0,str2num(StringbyKey("SpringConstant",note(ForceWave),":","\r"))}
//		DE_NewFeather#OutportForce(ForceWave,SepWave)
//		DE_NewFeather#RunFeatheronOutput(OptionsWave)
//		wave event_starts
//		variable result=numpnts(event_starts)
//		if(compensate==1)
//			event_starts-=numpnts(OutputSepWave)-numpnts(SepExtWave)
//		endif
//
//	
//		killwaves event_starts,ForceAll,SepAll,FretSm,SRetSm
//		
//	endfor
//	return result
//
//	//	
//	
//	//ExportWaveLists(ForceWaveList,SepWaveList)
//end

//
//Function TestCondition(Number,filtering,tol,temporal)
//	variable Number,filtering,tol,temporal
//	string AllForceRet= wavelist("*Force_Ret",";","")
//	String ForceWaveList="",SepWaveList=""
//	variable n,approachvelocity,retractvelocity,ratio
//	variable bottom=Number
//	variable top=Number+1//itemsinlist(AllFOrceRet)
//	for(n=bottom;n<top;n+=1)
//		//for(n=0;n<itemsinlist(AllFOrceRet);n+=1)
//		wave ForceRetWave=$stringfromlist(n,AllForceRet)
//		wave ForceExtWave=$replacestring("Ret",nameofwave(ForceRetWave),"Ext")
//		wave SepRetWave=$replacestring("Force",nameofwave(ForceRetWave),"Sep")
//		wave SepExtWave=$replacestring("Force",nameofwave(ForceExtWave),"Sep")
//		approachvelocity=str2num(stringbykey("ApproachVelocity",note(ForceRetWave),":","\r"))
//		retractvelocity=str2num(stringbykey("RetractVelocity",note(ForceRetWave),":","\r"))
//		ratio=approachvelocity/retractvelocity
//		print ratio
//		make/o/n=0  OutputForceWave,OutputSepWave
//		FakeApproach(ForceExtWave,SepExtWave,Ratio,OutputForceWave,OutputSepWave)
//		SetScale/P x 0,(dimdelta(ForceRetWave,0)),"s", OutputForceWave,OutputSepWave
//		make/free/n=0 ForceAll,SepAll
//	//	Concatenate/o {ForceExtWave,ForceRetWave},ForceAll
//	//	Concatenate/o {SepExtWave,SepRetWave},SepAll
//		Concatenate/o {OutputForceWave,ForceRetWave},ForceAll
//		Concatenate/o {OutputSepWave,SepRetWave},SepAll
//		make/free/n=0 FSm,SSm
//		make/o/n=0 FRetSm,SRetSm
//		if(filtering==1)
//			duplicate/free ForceAll FSM
//			duplicate/free SepAll SSM
//			duplicate/o ForceRetWave FRetSm
//			duplicate/o SepRetWave SRetSm
//		elseif(filtering>5)
//			DE_Filtering#FilterForceSep(ForceAll,SepAll,FSm,SSm,"SVG",filtering)
//			DE_Filtering#FilterForceSep(ForceRetWave,SepRetWave,FRetSm,SRetSm,"SVG",filtering)
//
//		elseif(filtering<1)
//			DE_Filtering#FilterForceSep(ForceAll,SepAll,FSm,SSm,"TVD",filtering)
//			DE_Filtering#FilterForceSep(ForceRetWave,SepRetWave,FRetSm,SRetSm,"TVD",filtering)
//
//		
//		endif
//
//		duplicate/o FSm $(Replacestring("Force_Ret",nameofwave(ForceRetWave),"Force"))
//		wave ForceWave=$(Replacestring("Force_Ret",nameofwave(ForceRetWave),"Force"))
//		duplicate/o SSm $(Replacestring("force_Ret",nameofwave(ForceRetWave),"Sep"))
//		wave SepWave=$(Replacestring("Force_Ret",nameofwave(ForceRetWave),"Sep"))
//		make/free/n=5 OptionsWave
//		OptionsWave={tol,temporal,str2num(StringbyKey("TriggerTime",note(ForceWave),":","\r")),0.0,str2num(StringbyKey("SpringConstant",note(ForceWave),":","\r"))}
//		DE_NewFeather#OutportForce(ForceWave,SepWave)
//		DE_NewFeather#RunFeatheronOutput(OptionsWave)
//		wave event_starts
//	event_starts-=numpnts(OutputSepWave)-numpnts(SepExtWave)
//		duplicate/o event_starts $replacestring("Force_ext",nameofwave(ForceExtWave),"Starts")
//		
//		duplicate/o FRetSm $replacestring("Force_ext",nameofwave(ForceExtWave),"FSm")
//		duplicate/o SRetSm $replacestring("Force_ext",nameofwave(ForceExtWave),"SSm")
//		killwaves event_starts,ForceAll,SepAll,FretSm,SRetSm
//		
//	endfor
//	print 	top
//
////	
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
////	
//	//ExportWaveLists(ForceWaveList,SepWaveList)
//end

//Static Function ReFilterandCalculateAll()
//	string AllForceRet= wavelist("*Force_Ret",";","")
//	String ForceWaveList="",SepWaveList=""
//	variable n
//	variable bottom=0
//	variable top=itemsinlist(AllFOrceRet)
//	variable filtering
//	wave FRetSm=$replacestring("Force_Ret", stringfromlist(0,wavelist("*Force_Ret",";","")),"FSm")
//
//	
//		if(cmpstr(stringbykey("DE_Filtering",note(FRetSm),":","\r"),"")==0)
//			Prompt filtering, "Enter Filtering"
//			DoPrompt "Enter Filtering" filtering
//		else
//			filtering=str2num(stringbykey("DE_Filtering",note(FRetSm),":","\r"))
//		endif
//	
//	
//	for(n=bottom;n<top;n+=1)
//		wave ForceRetWave=$stringfromlist(n,AllForceRet)
//		
//		RefilterRawWaves(ForceRetWave,filtering=filtering)
//	endfor
//
//end


//Static Function RefilterRawWaves(ForceRetWave,[filtering])
//	wave ForceRetWave
//	variable filtering
//
//
//	wave ForceExtWave=$replacestring("Ret",nameofwave(ForceRetWave),"Ext")
//	wave SepRetWave=$replacestring("Force",nameofwave(ForceRetWave),"Sep")
//	wave SepExtWave=$replacestring("Force",nameofwave(ForceExtWave),"Sep")
//	wave FRetSm=$replacestring("Force_ext",nameofwave(ForceExtWave),"FSm")
//	wave SRetSm=$replacestring("Force_ext",nameofwave(ForceExtWave),"SSm")
//	wave ThisEvent=$replacestring("Force_ext",nameofwave(ForceExtWave),"Starts")
//	
//	if(ParamisDefault(filtering))
//	
//		if(cmpstr(stringbykey("DE_Filtering",note(FRetSm),":","\r"),"")==0)
//			Prompt filtering, "Enter Filtering"
//			DoPrompt "Enter Filtering" filtering
//		else
//			filtering=str2num(stringbykey("DE_Filtering",note(FRetSm),":","\r"))
//		endif
//	endif
//	
//	if(filtering==1)
//
//		duplicate/o ForceRetWave FRetSm
//		duplicate/o SepRetWave SRetSm
//	elseif(filtering>5)
//		DE_Filtering#FilterForceSep(ForceRetWave,SepRetWave,FRetSm,SRetSm,"SVG",filtering)
//
//	elseif(filtering<1)
//		DE_Filtering#FilterForceSep(ForceRetWave,SepRetWave,FRetSm,SRetSm,"TVD",filtering)
//
//	endif
//	
//	duplicate/o ThisEvent FreeRupPnts
//	FreeRupPnts-=NUMPNTS(SepExtWave)
//	string Yuck=nameofwave(ThisEvent)
//	make/free/n=0 LCSBack,SLopesBack
//	duplicate/free FreeRupPnts ForcesBack
////	ForcesBack=FRetSm[FreeRupPnts[p]]
//		CalculateForcesFromPoints(ForceRetWave,SepRetWave,FreeRupPnts,ForcesBack)
//
//	variable offset=-1*str2num(stringbykey("DE_SChollOffset",note(ForceRetWave),":","\r"))-5e-9
//	CalcAllLCs(ForceRetWave,SepRetWave,FreeRupPnts,offset,LCsBack)
//	CalculateSlopes(ForceRetWave,SepRetWave,FreeRupPnts,offset,SlopesBack)
//	AddNotes(ForceRetWave,SepRetWave,FreeRupPnts,ForcesBack,LCSBack,SLopesBack)
//	FreeRupPnts+=NUMPNTS(SepExtWave)
//	//duplicate/o FreeRupPnts ThisEvent
//	killwaves FreeRupPnts
//		
//		
//end

//Static Function ProcessForceCurves()
//	string AllForceRet= wavelist("*Force_Ret",";","")
//	String ForceWaveList="",SepWaveList=""
//	variable n
//	variable top=itemsinlist(AllFOrceRet)
//	variable Entries
//	String RupForces,Contours,Slopes
//	variable ZeroForce
//	make/o/n=(0,5) First,Second,Third,Fourth,Fifth
//	for(n=0;n<top;n+=1)
//		Wave ForceWave=$StringFromList(n,wavelist("*Force_Ret",";",""))
//		Entries=CountEntries(ForceWave)
//		RupForces=Stringbykey("RupForce",note(ForceWave),":","\r")
//		Contours=Stringbykey("ContourLengths",note(ForceWave),":","\r")
//		ZeroForce=str2num(Stringbykey("DE_FOff",note(ForceWave),":","\r"))
//
//		Slopes=Stringbykey("Slopes",note(ForceWave),":","\r")
//		if(Entries==4)
//			InsertPoints/M=0 0,1, First,Third,Fourth,Fifth
//
//			First[0][0]=n
//			First[0][1]=str2num(Stringfromlist(0,RupForces))
//			First[0][2]=First[0][1]+ZeroForce
//			First[0][3]=str2num(Stringfromlist(0,Contours))
//			First[0][4]=str2num(Stringfromlist(0,Slopes))
//			
//
//			
//			Third[0][0]=n
//			Third[0][1]=str2num(Stringfromlist(1,RupForces))
//			Third[0][2]=Third[0][1]+ZeroForce
//			Third[0][3]=str2num(Stringfromlist(1,Contours))
//			Third[0][4]=str2num(Stringfromlist(1,Slopes))
//			
//			Fourth[0][0]=n
//			Fourth[0][1]=str2num(Stringfromlist(2,RupForces))
//			Fourth[0][2]=Fourth[0][1]+ZeroForce
//			Fourth[0][3]=str2num(Stringfromlist(2,Contours))
//			Fourth[0][4]=str2num(Stringfromlist(2,Slopes))
//			
//			Fifth[0][0]=n
//			Fifth[0][1]=str2num(Stringfromlist(3,RupForces))
//			Fifth[0][2]=Fifth[0][1]+ZeroForce
//			Fifth[0][3]=str2num(Stringfromlist(3,Contours))
//			Fifth[0][4]=str2num(Stringfromlist(3,Slopes))
//		
//		elseif(Entries==5)
//		InsertPoints/M=0 0,1, First,Second,Third,Fourth,Fifth
//			First[0][0]=n
//			First[0][1]=str2num(Stringfromlist(0,RupForces))
//			First[0][2]=First[0][1]+ZeroForce
//			First[0][3]=str2num(Stringfromlist(0,Contours))
//			First[0][4]=str2num(Stringfromlist(0,Slopes))
//			
//			Second[0][0]=n
//			Second[0][1]=str2num(Stringfromlist(1,RupForces))
//			Second[0][2]=Second[0][1]+ZeroForce
//			Second[0][3]=str2num(Stringfromlist(1,Contours))
//			Second[0][4]=str2num(Stringfromlist(1,Slopes))
//			
//			Third[0][0]=n
//			Third[0][1]=str2num(Stringfromlist(2,RupForces))
//			Third[0][2]=Third[0][1]+ZeroForce
//			Third[0][3]=str2num(Stringfromlist(2,Contours))
//			Third[0][4]=str2num(Stringfromlist(2,Slopes))
//			
//			Fourth[0][0]=n
//			Fourth[0][1]=str2num(Stringfromlist(3,RupForces))
//			Fourth[0][2]=Fourth[0][1]+ZeroForce
//			Fourth[0][3]=str2num(Stringfromlist(3,Contours))
//			Fourth[0][4]=str2num(Stringfromlist(3,Slopes))
//		
//			Fifth[0][0]=n
//			Fifth[0][1]=str2num(Stringfromlist(4,RupForces))
//			Fifth[0][2]=Fifth[0][1]+ZeroForce
//			Fifth[0][3]=str2num(Stringfromlist(4,Contours))
//			Fifth[0][4]=str2num(Stringfromlist(4,Slopes))
//
//		
//		else
//		
//
//		endif
//	endfor
//	
//	
//	//ExportWaveLists(ForceWaveList,SepWaveList)
//end
