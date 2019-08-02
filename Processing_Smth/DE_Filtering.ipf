#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma ModuleName = DE_Filtering
#include ":\Misc_PanelPrograms\AsylumNaming"


Static Function FilterByDefAuto(Number,Direction,FilterType,FilterPoints)
	variable Number,FilterPoints
	String Direction,FilterType

	String NameString="Image"

	String DefWaveName=DE_Naming#StringCreate(NameString,Number,"DeflCor",Direction)
	String ZSnsrWaveName=DE_Naming#StringCreate(NameString,Number,"ZSnsr",Direction)
	wave DefWave=$DefWaveName
	if( waveexists(DefWave)==0)
		Print "Error, something invalid in the Deflection wave naming"
	
	endif
	wave ZWave=$ZSnsrWaveName
	if(waveexists(ZWave)==0)
		Print "Error, something invalid in the ZSnsr wave naming"
	
	endif

	if(waveexists(ZWave)==0||waveexists(DefWave)==0)
		return -1
	endif

	FilterFECs(DefWave,ZWave,FilterType,FilterPoints)
	
	wave Force_Smth,Sep_Smth,Defl_Smth
	
	String ForceWaveName=DE_Naming#StringCreate(NameString,Number,"Force",Direction)
	String SepWaveName=DE_Naming#StringCreate(NameString,Number,"Sep",Direction)

	duplicate/o Defl_Smth $(DefWaveName+"_sm")
	duplicate/o Sep_Smth $(SepWaveName+"_sm")
	duplicate/o Force_Smth $(ForceWaveName+"_sm")
	killwaves Defl_Smth,Sep_Smth,Force_Smth
	//display $ForceWaveName vs $SepWaveName
	//appendtograph $(ForceWaveName+"_sm")
	//ModifyGraph rgb($ForceWaveName)=(65280,43520,43520)
end


Static Function FilterForceByDefAuto(Number,Direction,FilterType,FilterPoints)
	variable Number,FilterPoints
	String Direction,FilterType

	String NameString="Image"

	String ForceWaveName=DE_Naming#StringCreate(NameString,Number,"Force",Direction)
	String SepWaveName=DE_Naming#StringCreate(NameString,Number,"Sep",Direction)
	String DefWaveName=DE_Naming#StringCreate(NameString,Number,"DeflCor",Direction)
	String ZSnsrWaveName=DE_Naming#StringCreate(NameString,Number,"ZSnsr",Direction)
	
	wave ForWave=$ForceWaveName
	if( waveexists(ForWave)==0)
		Print "Error, something invalid in the Force wave naming"
	
	endif
		
	wave SepWave=$SepWaveName
	if(waveexists(SepWave)==0)
		Print "Error, something invalid in the Sep wave naming"
	
	endif
	
	
	duplicate/o ForWave $(DefWaveName)
	duplicate/o SepWave $(ZSnsrWaveName)
	//duplicate/o Force_Smth $(ForceWaveName+"_sm")
	
	wave DefWave= $(DefWaveName)
	wave ZWave=$(ZSnsrWaveName)

	variable spring=str2num(stringbykey("SpringConstant", note(DefWave),":","\r"))

	if(numtype(spring)==2)
		spring=str2num(stringbykey("Spring Constant", note(DefWave),":","\r"))
	endif
	if(numtype(spring)==2)
		print "Error: Couldn't find a reasonable spring constant"
	endif
	
	DefWave/=spring
	ZWave-=DefWave
	ZWave*=-1
	if(waveexists(ZWave)==0||waveexists(DefWave)==0)
		return -1
	endif

	FilterFECs(DefWave,ZWave,FilterType,FilterPoints)
	
	wave Force_Smth,Sep_Smth,Defl_Smth
	

	duplicate/o Defl_Smth $(DefWaveName+"_sm")
	duplicate/o Sep_Smth $(SepWaveName+"_sm")
	duplicate/o Force_Smth $(ForceWaveName+"_sm")
	//killwaves Defl_Smth,Sep_Smth,Force_Smth
	//display $ForceWaveName vs $SepWaveName
	//appendtograph $(ForceWaveName+"_sm")
	//ModifyGraph rgb($ForceWaveName)=(65280,43520,43520)
end


Static Function FilterFECs(DeflWave,ZSnsrWave,FilterType,FilterPoints)

	wave DeflWave,ZSnsrWave
	string FilterType
	variable FilterPoints
	variable spring=str2num(stringbykey("SpringConstant", note(DeflWave),":","\r"))

	if(numtype(spring)==2)
		spring=str2num(stringbykey("Spring Constant", note(DeflWave),":","\r"))
	endif
	if(numtype(spring)==2)
		spring=str2num(stringbykey("K", note(DeflWave),"=","\r"))
	endif
	if(numtype(spring)==2)
		print "Error: Couldn't find a reasonable spring constant"
	endif
	if(numpnts(DEFlWave)<=FilterPoints)
		duplicate/o DeflWave Sep_Smth, Force_Smth
		duplicate/o DeflWave Defl_Smth
		Sep_Smth=-ZSnsrWave+Defl_Smth
		Force_Smth=Defl_Smth*spring
		return -1
	endif
	duplicate/o DeflWave Defl_Smth
	if(cmpstr(FilterType,"SVG")==0)
		
		Smooth/S=2 FilterPoints, Defl_Smth
	elseif(cmpstr(FilterType,"TVD")==0)
		TVD1D_denoise(DeflWave,FilterPoints,Defl_Smth)
	//	TVD1D_denoise(ZSnsrWave,FilterPoints,ZSnsrWave)

	endif	
	duplicate/o DeflWave Sep_Smth, Force_Smth
	
	
	Sep_Smth=-ZSnsrWave+Defl_Smth
	Force_Smth=Defl_Smth*spring
end

Static Function FilterForceSep(ForceIn,SepIn,ForceOut,SepOut,FilterType,FilterPoints)
	Wave ForceIn,SepIn,ForceOut,SepOut
	string FilterType
	Variable FilterPoints

	
	variable spring=str2num(stringbykey("SpringConstant", note(ForceIn),":","\r"))

	if(numtype(spring)==2)
		spring=str2num(stringbykey("Spring Constant", note(ForceIn),":","\r"))
	endif
	
	if(numtype(spring)==2)
		spring=str2num(stringbykey("K", note(ForceIn),"=","\r"))
	endif
	
	if(numtype(spring)==2)
		print "Error: Couldn't find a reasonable spring constant"
	endif
	duplicate/free ForceIn DeflIn
	DeflIn/=spring
	
	duplicate/free SepIn ZSnsrIn
	ZSnsrIn-=DeflIn
	ZSnsrIn*=-1
	
	FilterFECs(DeflIn,ZSnsrIn,FilterType,FilterPoints)
	wave Force_Smth,Sep_Smth,Defl_Smth
	string FNote=note(Force_Smth)
	FNote=replacestringbykey("DE_Filtering",Fnote,num2str(FilterPoints),":","\r")
	note/K Force_Smth FNote
		string SNote=note(Sep_Smth)
	SNote=replacestringbykey("DE_Filtering",SNote,num2str(FilterPoints),":","\r")
	note/K Sep_Smth SNote
	duplicate/o Sep_Smth SepOut
	duplicate/o Force_Smth ForceOut
	killwaves Defl_Smth,Sep_Smth,Force_Smth

end

Static Function TVD1D_denoise(input1,lambda2,anotherwave)
	wave input1,anotherwave
	variable lambda2
	variable width=numpnts(input1)
	make/free/n=(numpnts(input1)) output
	variable k=0,k0=0
	variable umin=lambda2,umax=-lambda2
	variable vmin=input1[0]-lambda2, vmax=input1[0]+lambda2
	variable  kplus=0, kminus=0
	variable twolambda=2.0*lambda2	
	variable minlambda=-lambda2

	do
		if((k==(width-1)))
			do
				if(umin<0.0)
					do
						output[k0] = vmin
						k0 += 1
					while(k0<=kminus)
					k=k0
					kminus=k
					vmin=input1[k0]
					umin=lambda2
					umax=(vmin)+(lambda2)-vmax
				elseif(umax > 0.0)
					do
						output[k0] = vmax
						k0 += 1
					while(k0 <= kplus)					
					k=k0
					kplus=k
					vmax=input1[k0]
					umax=minlambda
					umin=(vmax)+(minlambda)-vmin
				else
					vmin += umin/(k-k0+1)
					do
						output[k0]=vmin
						k0 += 1
					while(k0 <= k)
					SetScale/P x,pnt2x(input1,0),dimdelta(input1,0),"s", output
					Note/K output, note(input1)
					duplicate/o output anotherwave	
				
					return 1 
				endif
			while (k==(width-1))
		endif
		umin += input1[k+1]-vmin
		umax += input1[k+1]-vmax
		if(umin < minlambda)
			do
				output[k0]=vmin
				k0 += 1
			while(k0 <= kminus)
			k=k0
			kminus=k
			kplus=kminus
			vmin=input1[k0]
			vmax=(vmin)+twolambda
			umin=lambda2
			umax=minlambda
		elseif(umax > lambda2)
			do
				output[k0]=vmax
				k0 += 1
			while(k0 <= kplus)
			k=k0
			kminus=k
			kplus=kminus
			vmax=input1[k0]
			vmin=(vmax)-twolambda
			umin=lambda2
			umax=minlambda
		else
			k += 1
			if(umin >= lambda2)
				kminus=k
				vmin+=(umin-lambda2)/((k)-k0+1)
				umin=lambda2
			endif			
			if (umax <= minlambda)
				kplus = k
				vmax+=(umax+lambda2)/((k)-k0+1)
				umax=minlambda
			endif
		endif


	
	while(1>-1)


end

Menu "Filtering"
	"Smooth data with TVD",TVDTop()
//	"Add lines for white noise", AllanWhite()
End
//Menu "Filtering"
//	//SubMenu "Processing"
//	"SmoothSomething", TVD()
//
//
//	//end
//	
//end

Function TVDTop()
//	Variable timerRefNum
//	Variable microSeconds
//	Variable n
//	timerRefNum = startMSTimer
	String wName
	Variable nSmooth = 1
	Prompt wName, "Smooth:",popup, WaveList("*", ";","WIN:")
	Prompt nSmooth, "Lambda for TVD"
	DoPrompt "TVD Smoothing",wName, nSmooth
	if(V_flag)	//user canceled
		Abort "User canceled TVD."
	endif
	WAVE w = $wName
	make/o/n=0 SmoothedWave
	TVD1D_denoise(w,nSmooth,SmoothedWave)

End


ThreadSafe Static Function AcausalLPFilter(Data,FilterRate,[StartIndex,StopIndex])
	Wave Data
	Variable FilterRate, StartIndex, StopIndex
	Variable UseSubRange = 2
	
	//overwrites the wave data
	//uses a 2 pass 2 pole butter LP filter on data at FilterRate
	//gets the SampleRate from the X scaling of the wave data
	//can provide startIndex and / or StopIndex to do a part of the wave.
	
	
	if (ParamIsDefault(StartIndex))
		StartIndex = 0
		UseSubRange -= 1
	endif
	if (ParamIsDefault(StopIndex))
		StopIndex = DimSize(Data,0)-1
		UseSubRange -= 1
	endif

	Variable SampleRate = 1/DimDelta(Data,0)
	Variable IgorFilterRate = FilterRate/SampleRate*(sqrt(2)-1)^.25
	if (IgorFilterRate <= 0)
		return(0)
	endif
	Variable NewFilterRate
	if (IgorFilterRate >= .5)
		IgorFilterRate = .49
		NewFilterRate = IgorFilterRate/((sqrt(2)-1)^.25)*SampleRate
		Print "LimitFilterRate in AcausalLPFilter from: "+num2str(FilterRate)+" to: "+Num2str(NewFilterRate)		//debug level, but can't do that from threads...
		//I am thinking that all thread print statements from threads are by definition, debug
	endif
	
	
	
	
	if (UseSubRange)
		Duplicate/O/R=[StartIndex,StopIndex] Data,AcausalFilterData
		Wave OrgData = Data
		Wave Data = AcausalFilterData
	endif

	Variable RampPoints, RampFraction = 1
	RampPoints = 1/FilterRate/RampFraction*SampleRate		//this is how many points needs to be added to each end.
	RampPoints = floor(RampPoints)
	//print RampPoints
	Variable Offset, Offset2, TempOffset

	Reverse/P Data
	offset = Data[0]*2-Data[RampPoints]
	FastOp Data = Data+(-Offset)
	InsertPoints/M=0 0,RampPoints,Data

	TempOffset = Data[RampPoints*2]-Data[RampPoints]
	Data[0,RampPoints-1] = Data[P+RampPoints]-TempOffset

//DoUpdate





	FilterIIR/CASC/LO=(IgorFilterRate)/ORD=2 Data
//DoUpdate
	DeletePoints/M=0 0,RampPoints,Data
	Reverse/P Data
//DoUpdate
	Offset2 = Data[0]*2-Data[RampPoints]
	Offset += Offset2
//DoUpdate
	FastOp Data = Data+(-Offset2)
	InsertPoints/M=0 0,RampPoints,Data
	TempOffset = Data[RampPoints*2]-Data[RampPoints]
	Data[0,RampPoints-1] = Data[P+RampPoints]-TempOffset
//DoUpdate
	FilterIIR/CASC/LO=(IgorFilterRate)/ORD=2 Data
	DeletePoints/M=0 0,RampPoints,Data
//DoUpdate
//DoUpdate
	FastOp Data = Data+(Offset)
//DoUpdate


	if (UseSubRange)
		OrgData[StartIndex,StopIndex] = Data[P-StartIndex]
		KillWaves Data
	endif



End //AcausalLPFilter

ThreadSafe Static Function ZPosRateFunc(SampleRate,ForceDist,MaxTime,nop)//,DirInfo)
	Variable SampleRate, ForceDist, MaxTime, nop
	//String DirInfo



	
	Variable MaxFilterRate = .49/((sqrt(2)-1)^.25)*SampleRate
	Variable MagicRatio = 185184*20
	Variable Divisor
//	if (NumType(ForceDist) == 1)
//		Divisor = 2300
//	else
		Divisor = 5555
//	endif
//	if (SampleRate > 15000)
		MagicRatio *= (SampleRate/Divisor-1)/2+1
//	endif
	
	if (NumType(ForceDist) == 1)		//triggered?
		if (MaxTime > 5)
			MagicRatio *= MaxTime/5
		endif
	endif

	Variable Rate
	Rate = MagicRatio/nop
	Rate = limit(Rate,20,MaxFilterRate/2)
	//Drop the Rate for short force plots.
	ForceDist = Log(Abs(ForceDist))
	if ((ForceDist < -6.5) && (SampleRate < 1e4))
//	if (ForceDist < -6.5)
		//-7 = 2.5
		//-8 = 4.5
		Rate /= (ForceDist*-2-11.5)
	elseif ((ForceDist > -5.2)	&& !numtype(ForceDist))	//6.3 um
		ForceDist = 10^(ForceDist)
		Rate *= (ForceDist)/4e-6
//		Rate *= (-5.2-ForceDist+1)
//		Rate *= 2
//		Rate
		
	endif
//	if (Strlen(DirInfo) && (Numtype(Str2num(StringFromList(0,DirInfo,","))) == 1))
	//	Rate = max(Rate,200)
	//endif
		
	Rate = limit(Rate,20,MaxFilterRate/2)
	

	return(Rate)
End //ZPosRateFunc


Static Function MakeZPositionFinal(Zposition,ForceDist,[Rate])
	Wave Zposition
	Variable ForceDist,Rate
	if(ParamisDefault(Rate))
		Variable Nop = DimSize(ZPosition,0)
	Variable SampleRate = 1/DimDelta(ZPosition,0)
	 Rate = ZPosRateFunc(SampleRate,ForceDist,RightX(ZPosition),nop)
	else
	
	endif


	



	MakeZPositionNew(Zposition,Rate)
	
	return(Rate)	
End //MakeZPositionFinal

Static Function MakeZPositionNew(Zposition,Rate)
	Wave Zposition
	Variable Rate

	AcausalLPFilter(ZPosition,Rate)
	return(0)


	
	
End //MakeZPositionNew

Static Function PointAtWaveandSmooth(InputForceWave,points)

	wave InputForceWave
	variable points
	Struct ForceWave Original 
	DE_Naming#TrueWavetoStruc(InputForceWave,Original)
	if(WaveExists($DE_Naming#StringCreate(Original.Name,Original.VNum,"Sep",Original.SDirec))==1)
		wave InputSepWave=$replacestring("Force",nameofwave(InputForceWave),"Sep")
	else
		print "Invalid Sep wave Associated"
		return -1
	endif
	make/free/n=0 FSm,SSm
	if(points>5)
		FilterForceSep(InputForcewave,InputSepWave,FSm,SSm,"SVG",points)	
	else 
		FilterForceSep(InputForcewave,InputSepWave,FSm,SSm,"TVD",points)	
	endif
	string ForceString
	DFREF Folder=Original.Folder
	duplicate/o FSm Folder:$DE_Naming#StringCreate(Original.Name,Original.VNum,"Force",Original.SDirec,modifier="_Sm")
	duplicate/o SSm Folder:$DE_Naming#StringCreate(Original.Name,Original.VNum,"Sep",Original.SDirec,modifier="_Sm")

end

Static Function/S AutofilterAForceSep(ForceWave,Type,Pnts,[WindowString])
	Wave ForceWave
	String Type,WindowString
	Variable Pnts
	
	wave SepWave=$replacestring("Force",nameofwave(ForceWave),"Sep")
	make/o/n=0 $replacestring("Force",nameofwave(ForceWave),"FSm"),$replacestring("Force",nameofwave(ForceWave),"Ssm")
	wave SmForceWave=$replacestring("Force",nameofwave(ForceWave),"FSm")
	wave SmSepWave=$replacestring("Force",nameofwave(ForceWave),"Ssm")
	FilterForceSep(ForceWave,SepWave,SmForceWave,SmSepWave,Type,Pnts)
	if(ParamisDefault(WindowString))
	else
	appendtograph/W=$WindowString SmForceWave vs SmSepWave
	endif
end