#pragma rtGlobals=1	// Use modern global access method.


Function Cypher_Thermal(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = a*x
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 4
//CurveFitDialog/ w[0] = Amplitude at DC
	//CurveFitDialog/ w[1] = Q
	//CurveFitDialog/ w[2] = omega0, in either radians or Hz
	//CurveFitDialog/ w[3] = white noise. Units are m/rtHz
	
	
	variable DC = w[0]							//using local variables is faster than using wave points
	variable QTerm = w[1]
	variable Omega0 =w[2]
	variable WNTerm = w[3]^2




	return sqrt(WNTerm+DC^2*Omega0^4/QTerm^2*(((x)^2-Omega0^2 )^2+Omega0 ^2*x^2/QTerm^2)^(-1))
End

Function Cypher_ThermalSq(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = a*x
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 4
//CurveFitDialog/ w[0] = Amplitude at DC
	//CurveFitDialog/ w[1] = Q
	//CurveFitDialog/ w[2] = omega0, in either radians or Hz
	//CurveFitDialog/ w[3] = white noise. Units are m/rtHz
	
	
	variable DC = w[0]							//using local variables is faster than using wave points
	variable QTerm = w[1]
	variable Omega0 =w[2]
	variable WNTerm = w[3]^2




	return abs(WNTerm+DC^2*Omega0^4/QTerm^2*(((x)^2-Omega0^2 )^2+Omega0 ^2*x^2/QTerm^2)^(-1))
End

Function Pirzer_Thermal(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = a*x
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 4
//CurveFitDialog/ w[0] = Amplitude at DC
	//CurveFitDialog/ w[1] = Q
	//CurveFitDialog/ w[2] = omega0, in either radians or Hz
	//CurveFitDialog/ w[3] = white noise. Units are m/rtHz
	
	
	variable DC = w[0]							//using local variables is faster than using wave points
	variable QTerm = 1/w[1]^2
	variable Omega0 =w[2]
	variable WNTerm = w[3]^2




	return (WNTerm+DC^2*Omega0^2*Qterm/4*(((x)-Omega0 )^2+Omega0^2*QTerm/4)^(-1))
End




Function SHOAmpWhiteTest(w,xwave)	: FitFunc		//this is an all at once fit function
	Wave w
	Variable xwave							//the XWave does nothing at the moment
	
	
	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = a*x
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 4
	//CurveFitDialog/ w[0] = Amplitude at DC
	//CurveFitDialog/ w[1] = Q
	//CurveFitDialog/ w[2] = omega0, in either radians or Hz
	//CurveFitDialog/ w[3] = white noise. Units are m/rtHz
	
	
	variable DC = w[0]							//using local variables is faster than using wave points
	variable QTerm = 1/w[1]^2
	variable Omega0 =w[2]
	variable WNTerm = w[3]^2
	
	return Sqrt((DC*Omega0/XWave)^2/ ((Omega0/XWave - XWave/Omega0)^2 + QTerm) + WNTerm)	//calculate the wave
	

end //SHOAmpWhite

function Thermalk_Cypher(w)
	wave w
	//returns the thermal spring constant calculated from one of the SHO functions
	//This assumes you have fit a linear Power Spectral Density of the cantilever motion (units of m/rtHz)
	//using a function such as SHOAmp or SHOAmpPinkWhite
	//w[0] = Amplitude at DC
	//w[1] = Q
	//w[2] = omega0, in Hz
	//x is the variable omega
	
	
	//ref:
	//Short cantilevers for atomic force microscopy.  Walters DA, Cleveland
	//JP, Thomson NH, Hansma PK, Wendman MA, Gurley G, Elings V, Rev. Sci.
	//Instrum. 67 (10) 1996, p. 3583, endnote 39
	
	//variable energy = 1.38e-23 * GV("ThermalTemperature")
	variable energy = 1.38e-23 * 300

	
	return(2*energy/(Pi*w[2]*w[0]^2))*abs(w[1])

end //Thermalk

function Thermalk_Pirzer(w)
	wave w
	//returns the thermal spring constant calculated from one of the SHO functions
	//This assumes you have fit a linear Power Spectral Density of the cantilever motion (units of m/rtHz)
	//using a function such as SHOAmp or SHOAmpPinkWhite
	//w[0] = Amplitude at DC
	//w[1] = Q
	//w[2] = omega0, in Hz
	//x is the variable omega
	
	
	//ref:
	//Short cantilevers for atomic force microscopy.  Walters DA, Cleveland
	//JP, Thomson NH, Hansma PK, Wendman MA, Gurley G, Elings V, Rev. Sci.
	//Instrum. 67 (10) 1996, p. 3583, endnote 39
	
	//variable energy = 1.38e-23 * GV("ThermalTemperature")
	variable energy = 1.38e-23 * 300

	
	return(energy/(w[2]*w[0]^2)*2*w[1]/(pi/2+atan(2*w[1])))

end //Thermalk