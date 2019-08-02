#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
Function CleanUpRates(RateIn,RateOut,MinNum)
	wave RateIn,RateOut
	variable MinNum
	duplicate/free Ratein RateHold
	variable n=0
	wave Hist=$replacestring("Rate",nameofwave(Ratein),"Hist")
	for(n=0;n<numpnts(Ratein);n+=1)
		if(Hist[n]<MinNum)
		RateHold[n]=NaN
		
		endif	
	endfor
	duplicate/o RateHold RateOut

end

Function AutoClean(WaveIn,Divisor)
	wave Wavein
	variable Divisor
	make/o/n=0 $Replacestring("Rate",nameofwave(Wavein),"CutRate")
	wave OutWave=$Replacestring("Rate",nameofwave(Wavein),"CutRate")
		wave Hist=$replacestring("Rate",nameofwave(Wavein),"Hist")


	CleanUpRates(Wavein,OutWave,sum(Hist)/divisor)
	print sum(Hist)/divisor
	print nameofwave(OutWave)
end