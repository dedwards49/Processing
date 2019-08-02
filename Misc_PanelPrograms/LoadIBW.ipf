#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma modulename=DE_IBWLoad
#include ":\Misc_PanelPrograms\AsylumNaming"
#include "SimpleWLCPrograms"
#include "FEC_Fast_Corr"


Static function/S LoadAsylumIBWbyIndex(Index,PathString)
	variable Index
	String PathString
	NewPath/o/q Test PathString
	string Filename=IndexedFile(Test,Index,".ibw")
	string Waven=IndexedFile(Test,Index,".ibw")[0,strlen(Filename)-5]
	LoadWave/O/P=Test/Q/A  IndexedFile(Test,Index,".ibw")
	wave w1=$waveN	
	return nameofwave(w1)
end

Static function GlideLoad(index,Pathstring)
	variable Index
	String PathString
	wave w1=$LoadAsylumIBWbyIndex(Index,PathString)
	string origname=nameofwave(w1)
	variable npnt=DimSize(w1,0)
	make/free/n=(npnt,2)	Test
	Test[][0]=w1[p][%ZSnsr]
	Test[][1]=w1[p][%Defl]
	Note/K Test note(w1)
	setscale/P x,dimoffset(w1,0),dimdelta(w1,0),Test
	SetDimLabel 1,0,ZSnsr,Test 
	SetDimLabel 1,1,Defl,Test 
	duplicate/o Test $origname

end

Static Function LoadAllinFolder(type,PathString)
	String type,PathString
	NewPath/o/q Test PathString
	variable index=0
	do
		GlideLoad(index,Pathstring)
		Index+=1
	while(cmpstr(IndexedFile(Test,Index,".ibw"),"")!=0)
	
end

Static Function RunonAll(FunctionString)
	String FunctionString
	variable n
	do
		wave w1=WaveRefIndexedDFR(root:,n )
		strswitch(FunctionString)
		
		case "BG":
			CorrectBG(w1)
			break
		case "Shift Sep":
			SepShift(w1)
			break
		endswitch
		n+=1
	while(cmpstr(nameofwave(WaveRefIndexedDFR(root:,n )),"")!=0)
end


Static Function CorrectBG(input)
	wave input
	variable offset
	string Speed=stringbykey("Pull Speed (Pair)",note(input),":","\r")

	strswitch(Speed)
	
		case "Fast":
			CorrectDeflection(input,input)
			input[][%Force]*=-1
			break
		
		case "Slow":
			struct ForceWave Name
			DE_Naming#WavetoStruc(nameofwave(input),Name)
			string BaseWave=DE_Naming#StringCreate(Name.Name,Name.VNum-1,Name.SMeas,Name.SDirec)
			wave w1=$BaseWave
			CorrectDeflection(w1,input)
			input[][%Force]*=-1


			break
		
		
		default:
			return -1
	
	endswitch
end

Static Function SepShift(input)
	wave input
	variable offset
	string Speed=stringbykey("Pull Speed (Pair)",note(input),":","\r")
	strswitch(Speed)
	
		case "Fast":
			offset=SepAtTrigger(input)
			input[][%Sep]-=offset
			break
		
		case "Slow":
			struct ForceWave Name
			DE_Naming#WavetoStruc(nameofwave(input),Name)
			string BaseWave=DE_Naming#StringCreate(Name.Name,Name.VNum-1,Name.SMeas,Name.SDirec)
			wave w1=$BaseWave
			offset=SepAtTrigger(w1)
			input[][%Sep]-=offset
			break
		
		
		default:
			return -1
	
	endswitch


end

Static Function ContourLength(IBWWave)
	Wave IBWWave
	
	if(dimsize(IBWWave,1)==5)
	else
		print "Incorrect wave"
		return -1
	endif
	
	make/free/n=(dimsize(IBWWave,0)) CL
	
	CL=DE_WLC#ContourTransform(IBWWave[p][%Force],IBWWave[p][%Sep],0.4e-9,298)
	Make/N=100/o CL_Hist
	Histogram/B={1e-009,2e-009,100} CL,CL_Hist
	
	if(FindDimLabel(IBWWave, 1, "CL" )==-2)
		Insertpoints/M=1 dimsize(IBWWave,1),1, IBWWave
		SetDimLabel 1,(dimsize(IBWWave,1)-1),Sep,IBWWave 

	else
	endif
	IBWWave[][ FindDimLabel(IBWWave, 1, "CL" )]=0
	IBWWave[][ FindDimLabel(IBWWave, 1, "CL" )]=0
 end


Static Function CorrectDeflection(BaseWave,CorrWave)
	wave BaseWave,CorrWave
	//	DE_Naming#WavetoStruc(nameofwave(IBWWave),Name)
	//	string FinalName
	string Direction= stringbykey("Direction",note(BaseWave),":","\r")
	string Indexes= stringbykey("Indexes",note(BaseWave),":","\r")
	make/free/n=(ItemsInList(Direction,",")) DirecWave
	variable n,startpnt,endpnt
	for(n=0;n<ItemsInList(Direction,",");n+=1)
		DirecWave[n]=str2num(StringFromList(n, Direction,","))
	endfor
	FindLevels/Q DirecWave, 1
	wave w1=W_FindLevels
	deletepoints 0,1, w1
	//
	startpnt=str2num(StringFromList(w1[0]-1, Indexes,","))
	endpnt=str2num(StringFromList(w1[0], Indexes,","))
	duplicate/free/r=[startpnt,endpnt] BaseWave Approach
	variable Cutoff=FindCutoff(Approach)

	make/free/n=(dimsize(BaseWave,0)) WDef_fit,WZsen_fit
	WDef_fit=BaseWave[p][%Defl]
	WZsen_fit=BaseWave[p][%ZSnsr]
	make/free/n=(dimsize(CorrWave,0)) WDef_corr,WZsen_corr
	WDef_corr=CorrWave[p][%Defl]
	WZsen_corr=CorrWave[p][%ZSnsr]
	note/K WDef_fit,note(BaseWave)
	note/K WZsen_fit,note(BaseWave)

	note/K WDef_corr,note(BaseWave)
	note/K WZsen_corr,note(CorrWave)

	CorrFEC(WDef_fit,WZsen_fit,WDef_corr,WZsen_corr,0,Cutoff,NewName="TestDeflCorr",ForceName="TestForce",SepName="TestSep",FitType="LinSin",ResName="Res")
	wave TestDeflCorr,TestForce,TestSep,Res
	
	if(FindDimLabel(CorrWave, 1, "DeflCorr" )==-2)
	Insertpoints/M=1 2,1, CorrWave
	SetDimLabel 1,2,DeflCorr,Corrwave 

	else
	endif
	
	if(FindDimLabel(CorrWave, 1, "Force" )==-2)
	Insertpoints/M=1 3,1, CorrWave
	SetDimLabel 1,3,Force,Corrwave 

	else
	endif
	
	if(FindDimLabel(CorrWave, 1, "Sep" )==-2)
	Insertpoints/M=1 4,1, CorrWave
	SetDimLabel 1,4,Sep,Corrwave 

	else
	endif
	
	Corrwave[][%DeflCorr]=TestDeflCorr[p]
	Corrwave[][%Force]=TestForce[p]
	Corrwave[][%Sep]=TestSep[p]

	
	killwaves TestDeflCorr,TestForce,TestSep,Res,w1
end

Static Function FindCutoff(Approach)
	wave Approach
	
	variable triggerDistance=str2num(stringbykey("TriggerPoint",note(approach),":","\r"))/str2num(stringbykey("SpringConstant",note(approach),":","\r"))
	variable timetotrigger=1.25* triggerDistance/str2num(stringbykey("ApproachVelocity",note(approach),":","\r"))*1e6
	
	variable finaltime=DimOffset(Approach, 0) + (dimsize(Approach,0)-1) *DimDelta(Approach,0)
	variable crossingpnt=((finaltime-timetotrigger) - DimOffset(Approach, 0))/DimDelta(Approach,0)
	return floor(crossingpnt)
end

Static Function SepAtTrigger(BaseWave)
	wave BaseWave
	
	variable trigtime= str2num(stringbykey("TriggerTime1",note(BaseWave),":","\r"))
	variable pnt=(trigtime - DimOffset(BaseWave, 0))/DimDelta(BaseWave,0)
	variable res
	if(dimsize(BaseWave,1)==5)
		res=-BaseWave[pnt][%ZSnsr]+BaseWave[pnt][%DeflCorr]
	else
		res=0 
	endif
	return res
end

Static function  GenerateCurves(IBWWave,Meas,Direc,[number])
	wave IBWWave
	string Meas,Direc
	variable number
	variable npnt=DimSize(IBWWave,0)
	struct ForceWave Name
	DE_Naming#WavetoStruc(nameofwave(IBWWave),Name)
	string FinalName
	string Direction= stringbykey("Direction",note(IBWWave),":","\r")
	string Indexes= stringbykey("Indexes",note(IBWWave),":","\r")
	variable n
	make/free/n=(ItemsInList(Direction,",")) DirecWave
	for(n=0;n<ItemsInList(Direction,",");n+=1)
		DirecWave[n]=str2num(StringFromList(n, Direction,","))
	endfor
	
	If(ParamisDefault(number))
		number=0
	endif
	
	strswitch(Meas)
		case "Defl":
			make/free/n=(npnt)	Test
			Test=	IBWWave[p][%Defl]
			Note/K Test note(IBWWave)
			setscale/P x,dimoffset(IBWWAVE,0),dimdelta(IBWWAVE,0),Test
			break
		default:
		
		case "Force":
			make/free/n=(npnt) Test
			Test=	IBWWave[p][%Defl]
			Note/K Test note(IBWWave)
			setscale/P x,dimoffset(IBWWAVE,0),dimdelta(IBWWAVE,0),Test
			variable k=str2num(Stringbykey("SpringConstant",note(IBWWave),":","\r"))
			Test*=k
			break
		case "ZSnsr":
			make/free/n=(npnt) Test
			Test=	IBWWave[p][%ZSnsr]
			Note/K Test note(IBWWave)
			setscale/P x,dimoffset(IBWWAVE,0),dimdelta(IBWWAVE,0),Test

			break
		
		default:
			print "Unknown Measure"
			return -1
	endswitch
	
	variable startpnt,endpnt
	
	strswitch(Direc)
		case "All":
			duplicate/free/o Test Final

			break
			
		case "Ret":
			FindLevels/Q DirecWave, -1
			wave w1=W_FindLevels
			startpnt=str2num(StringFromList(w1[number]-1, Indexes,","))
			endpnt=str2num(StringFromList(w1[number], Indexes,","))
			duplicate/free/r=[startpnt,endpnt] Test Final
			 
			break
		case "Ext":
			FindLevels/Q DirecWave, 1
			wave w1=W_FindLevels
			deletepoints 0,1, w1

			startpnt=str2num(StringFromList(w1[number]-1, Indexes,","))
			endpnt=str2num(StringFromList(w1[number], Indexes,","))
			duplicate/free/r=[startpnt,endpnt] Test Final
			break
			
		case "Towd":
			FindLevels/Q DirecWave, 0
			wave w1=W_FindLevels
			startpnt=str2num(StringFromList(w1[2*number]-1, Indexes,","))
			endpnt=str2num(StringFromList(w1[number], Indexes,","))
			duplicate/free/r=[startpnt,endpnt] Test Final
			break
			
		case "Away":
			FindLevels/Q DirecWave, 0
			wave w1=W_FindLevels
			startpnt=str2num(StringFromList(w1[2*number+1]-1, Indexes,","))
			endpnt=str2num(StringFromList(w1[number], Indexes,","))
			duplicate/free/r=[startpnt,endpnt] Test Final
			break
		default:
			print "Unknown Direction"
			return -1
	endswitch
	FinalName=DE_Naming#StringCreate(Name.Name,Name.VNum,Meas,Direc)
	killwaves/z W_FindLevels
	duplicate/o Final $FinalName
end
