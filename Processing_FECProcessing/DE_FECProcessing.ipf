#pragma rtGlobals=3		// Use modern global access method and strict wave access.
//The goal of this software suite is to complement calculations 

#pragma ModuleName = DE_FECProcessing
#Include "SimpleWLCPrograms"


Static Function ContourLengthChanges(Force,Sep,OutWave,[PersistanceLength,Temp])
	Wave Force,Sep,OutWave
	variable PersistanceLength,Temp
	
	if(ParamIsDefault(PersistanceLength))	//Default
		PersistanceLength=0.4e-9
			
	endif
	if(ParamIsDefault(Temp))	//Default
		Temp=298
			
	endif
	
	if(itemsinlist(stringbykey("PRHEvents",note(Force),":","\r"),":")!=0)
		variable n,Starts,Ends
		String EventList=StringbyKey("PRHEvents",note(Force),":","\r")
		make/free/n=(itemsinlist(EventList,":")) ContourLengthPeak
		for(n=0;n<(itemsinlist(EventList,":"));n+=1)
			Ends=str2num(StringFromList(n, EventList ,":"))+1
			if(n==0)
				Starts=0
			else
				Starts=str2num(StringFromList((n-1), EventList ,":"))
			endif
			duplicate/free/r=[Starts,Ends] Force CutForce
			duplicate/free/r=[Starts,Ends] Sep CutSep
			variable shift=565e-9
			CutForce*=-1
			CutSep-=shift
			ContourLengthPeak[n]=ContourLengthSingle(CutForce,CutSep,PersistanceLength=PersistanceLength,Temp=Temp)

		endfor
		make/free/n=(numpnts(ContourLengthPeak)-1) DeltaLC
		DeltaLC=ContourLengthPeak[p+1]-ContourLengthPeak[p]
	endif
	duplicate/o DeltaLC OutWave
end

Static Function ContourLengthSingle(Force,Sep,[PersistanceLength,Temp,Method])
	Wave Force,Sep
	variable PersistanceLength,Temp
	string Method
	if(ParamIsDefault(PersistanceLength))	//Default
		PersistanceLength=0.4e-9
			
	endif
	if(ParamIsDefault(Temp))	//Default
		Temp=298
			
	endif
	
	if(ParamIsDefault(Method))	//Default
		Method="Median"
			
	endif


	duplicate/free Force ContourLength
	ContourLength=ContourTransform(Force[p],Sep[p],PersistanceLength,Temp)
	
	strswitch(Method)
		case "Mean":
			wavestats/q ContourLength
			return v_avg
			break	
		
		case "Median":
			StatsQuantiles/q ContourLength
			return v_median
			break
		case "Gauss":
			variable number=100
			StatsQuantiles/Q ContourLength
			variable Lower=max(1e-9,v_median-4*V_IQR)
			variable Steps=3*(v_median-Lower)/number
			make /free/n=(number) Histo
			Histogram/B={Lower,Steps,number} ContourLength,Histo
			CurveFit/NTHR=0 gauss  Histo 
			wave w_coef
			return w_coef[2]
			break
	endswitch

end

Static Function HistogramSomeShit(InputWave,OutputWave)

	Wave InputWave,OutputWave

	variable number=100
	StatsQuantiles/Q InputWave
	
	variable Lower=max(1e-9,v_median-4*V_IQR)
	variable Steps=3*(v_median-Lower)/number
	make /free/n=(number) Histo
	Histogram/B={Lower,Steps,number} InputWave,Histo
	CurveFit/NTHR=0 gauss  Histo 
		duplicate/o Histo OutputWave

	wave w_coef
	return w_coef[2]
end

Static Function EventContourTransformAll(Force,Sep,PersistanceLength,Temp,[Histograms])
	Wave Force,Sep
	variable PersistanceLength,Temp,Histograms
	
	if(ParamIsDefault(Histograms ))	//I assume you just want the histograms. 
		Histograms=1
			
	endif
	
	if (Histograms==1)
		variable Lower,Steps			
		variable number=40

	else

		
	endif
		
	variable shift=565e-9
	if(itemsinlist(stringbykey("PRHEvents",note(Force),":","\r"),":")!=0)
		variable n,Starts,Ends
		String EventList=StringbyKey("PRHEvents",note(Force),":","\r")
		//n=0
		for(n=0;n<(itemsinlist(EventList,":"));n+=1)
			Ends=str2num(StringFromList(n, EventList ,":"))+1
			if(n==0)
				Starts=0
			else
				Starts=str2num(StringFromList((n-1), EventList ,":"))
			endif
			duplicate/free/r=[Starts,Ends] Force CutForce
			duplicate/free/r=[Starts,Ends] Sep CutSep
			CutForce*=-1
			CutSep-=shift
			duplicate/free/r=[Starts,Ends] Force ContourLength

			ContourLength=ContourTransform(CutForce[p],CutSep[p],PersistanceLength,Temp)
		
			if (Histograms==1)
				wavestats/Q ContourLength
				StatsQuantiles/Q ContourLength
				Lower=max(1e-9,v_median-4*V_IQR)
				Steps=2*(v_median-Lower)/number
				make /free/n=(number) Histo
				Histogram/B={Lower,Steps,number} ContourLength,Histo

				duplicate/o Histo $(nameofwave(Force)+"_CLH_"+num2str(n))

			else
				duplicate/o ContourLength $(nameofwave(Force)+"_CL_"+num2str(n))

		
			endif
		endfor
			
	else
			
	endif
	
	

end