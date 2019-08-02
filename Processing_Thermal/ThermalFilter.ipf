#pragma rtGlobals=3		// Use modern global access method and strict wave access.

Function/C DE_ARC2EWMAFilterResp(freq, f3dB)
	// Response of a digital EWMA filter
	
	Variable freq, f3dB
	 variable cMasterSampleRate=5e4

	f3dB = abs(f3dB/cMasterSampleRate)

	Variable alpha = (1 - sin(2*pi*f3dB))/cos(2*pi*f3dB)
	
	Variable/C fresp, cz
	cz = exp(cmplx(0,1)*2*pi*freq/cMasterSampleRate)	

	fresp = ((1 - alpha)/2)*(cz + 1)/(cz - alpha)
	
	if (f3dB < 0)
		fresp = 1 - fresp	// highpass
	endif
	
	return fresp

End

Function/C BackpackFilterResp(freq, f3dB, selectivity)
	// Response of Backpack filters
	// Inputs: freq - frequency in Hz; f3db - filter 3dB cutoff frequency in Hz, selectivity - 1 or >> 1
	// Output: Complex response of the filter @ freq
	
	Variable freq, f3dB, selectivity
	
	Variable/C fresp, ca
	
	ca = cmplx(0,1)*(freq/f3dB)	// s domain becuase the backpack sampling rate (5 MHz) is quite high
	
	if (selectivity == 1)	
		fresp = 1/(ca + 1)	// order 1
	else
		fresp = 1/((ca + 1)*(ca^2 - 2*ca*cos((2 + 3 - 1)/(2*3)*pi) + 1))	// order 3
	endif
		
	return (fresp)

End	//BackpackFilterResp