#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
Function PareandAverage(Rup,Slope,mAXfORCE)
	wave Rup,Slope
	VARIABLE mAXfORCE	
	
	duplicate/free Rup FreeRup
	duplicate/free Slope FreeSlope
	variable n
	for(n=0;n<numpnts(FreeSlope);n+=1)
		if(FreeSlope[n]<0)
		
			deletepoints n,1, FreeSlope,FreeRup
			n-=1
		endif
		
	endfor
		for(n=0;n<numpnts(FreeSlope);n+=1)
		if(FreerUP[n]<0||FreerUP[n]>mAXfORCE)
		
			deletepoints n,1, FreeSlope,FreeRup
			n-=1
		endif
		
	endfor
	variable AvgSlope,AvgRup,StdSlope,StdRup,SEMSlope,SEMRup
		wavestats/Q FreeRup
AvgRup=v_avg
StdRup=V_sdev
SEMRup=V_Sem
	wavestats/Q FreeSlope
	AvgSlope=v_avg
StdSlope=V_sdev
SEMSlope=v_sem
print AvgSlope
print StdSlope
print SEMSlope
print AvgRup
print StdRup
print SEMRup

end

Function BellEvans(RuptureForce,Slope)
	wave RuptureForce,Slope

	duplicate/free Slope FreeSlope
		duplicate/free RuptureForce FreeForce
FreeForce[2]=Nan
FreeSlope[2]=Nan

	FreeSlope=ln(Slope)
	CurveFit/Q/W=2 Line FreeForce/X=FreeSlope 
	wave w_coef
		variable KT=4.11e-21

	variable DX=kT/w_coef[1]
	variable K0=DX*(exp(w_coef[0]*DX/KT)*KT)^-1
	print DX
	print K0
	make/o/n=0 $(nameofwave(RuptureForce)+"_fit")
	wave w1=$(nameofwave(RuptureForce)+"_fit")
	makeTheFit(w_coef,Slope,w1)

end

Function makeTheFit(CoefIn,SlopeIn,WaveOut)
	wave CoefIn,SlopeIn,Waveout
	variable slope=Coefin[1]
	variable Inter=Coefin[0]
	
	make/free/n=100 FreeResult
	SetScale/P x SlopeIn[0],SlopeIn[numpnts(SlopeIn)-1],"", FreeResult
	FreeResult=slope*ln(x)+Inter
		duplicate/o FreeResult WaveOut
end