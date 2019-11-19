#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma Modulename=DE_AdjustInvols


Static Function ExtractNoise(waveout1,waveout2,waveout3,waveout4)

	wave waveout1,waveout2,waveout3,waveout4
	string AllForceRet= wavelist("*Force_Ret",";","")
	String ForceWaveList="",SepWaveList=""
	variable	Bottom=0
	variable	top=itemsinlist(AllFOrceRet)
	make/free/n=(top) FreeNoise,FreePntsperSec,FreeBW,FreeInvols
	variable n,offsetpoints
	variable desiredfreq=1e3
	string WaveNote
	//
	for(n=bottom;n<top;n+=1)
		wave ForceRetWave=$stringfromlist(n,AllForceRet)
		wave ForceExtWave=$replacestring("Force_ret",nameofwave(ForceRetWave),"Force_ext")
		wave Points=$replacestring("Force_ret",nameofwave(ForceRetWave),"Starts")
		offsetpoints=numpnts(ForceExtWave)
		FreeNoise[n]=ExtractNoisefromWave(ForceRetWave,Points,offsetpoints,desiredfreq)
		FreePntsperSec[n]=str2num(Stringbykey("NumPtsPerSec",note(ForceRetWave),":","\r"))
		FreeBW[n]=str2num(Stringbykey("ForceFilterBW",note(ForceRetWave),":","\r"))
		FreeInvols[n]=str2num(Stringbykey("InvOLS",note(ForceRetWave),":","\r"))

	endfor
	duplicate/o FreeNoise waveout1
	duplicate/o FreePntsperSec waveout2
	duplicate/o FreeBW waveout3
	duplicate/o FreeInvols waveout4


end

Static Function ExtractNoisefromWave(ForceWave,PntWave,offsetpoints,desiredfreq)
	wave ForceWave,PntWave
	variable offsetpoints,desiredfreq
	variable currentfreq=1/dimdelta(ForceWave,0)
	variable cutdown=desiredfreq/currentfreq
	variable finalpoint=PntWave[numpnts(PntWave)-1]-offsetpoints
	variable endofwave=numpnts(ForceWave)-1
	duplicate/o/r=[endofwave-(endofwave-finalpoint)+300,endofwave] ForceWave Garbage
	Make/O/D/N=0 coefs; DelayUpdate
	Duplicate/free Garbage, filtered; DelayUpdate

	FilterFIR/DIM=0/LO={cutdown,1.2*cutdown,101}/COEF coefs, filtered
	wavestats/q filtered
	return v_sdev


end

Static Function RecalcWaveswithNewInvols(Force,Sep,NewInvols,ForceOut,SepOut)

	Wave Force,Sep,ForceOut,SepOut
	variable NewInvols
	
	duplicate/free Force FreeOrigForce,FreeCalcDefl,FreeCalcDefV
	duplicate/free Sep FreeOrigSep,FreeCalcZSnsr
	variable spring=str2num(stringbykey("SpringConstant",note(Force),":","\r"))
	variable OldInvols=str2num(stringbykey("Invols",note(Force),":","\r"))
	
	FreeCalcDefl=FreeOrigForce/spring
	FreeCalcDefV=FreeCalcDefl/OldInvols
	FreeCalcZSnsr=FreeOrigSep-FreeCalcDefl
	FreeCalcDefl=FreeCalcDefV*NewInvols
	FreeOrigForce=FreeCalcDefl*spring
	FreeOrigSep=FreeCalcZSnsr+FreeCalcDefl
	duplicate/o FreeOrigForce ForceOut
	duplicate/o FreeOrigSep SepOut

end

Static Function DetermineInvolsAdj(ForceWave,PntWave,offsetpoints,TargetNoise,Frequency)
	wave ForceWave,PntWave
	variable TargetNoise,Frequency,offsetpoints
	
	variable noise=ExtractNoisefromWave(ForceWave,PntWave,offsetpoints,Frequency)
	variable InvAdj=TargetNoise/noise
	return InvAdj
end