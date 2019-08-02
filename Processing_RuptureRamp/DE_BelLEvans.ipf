#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma modulename=DE_BellEvans

Function/C Intercepts(m1,m2,b1,b2)
	variable m1,m2,b1,b2
	variable x= (b2-b1)/(m1-m2)
	variable y=m2*x+b2
	return cmplx(x,exp(y))

end

Function Bell(w,F) : FitFunc
	Wave w
	Variable F

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(F) = ln(k0)+Folding*F*Dx/4.11e-21
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ F
	//CurveFitDialog/ Coefficients 3
	//CurveFitDialog/ w[0] = k0
	//CurveFitDialog/ w[1] = Dx
	//CurveFitDialog/ w[2] = Folding

	return ln(w[0])+w[2]*F*w[1]/4.11e-21
End

Static Function/C CalcImportantParms(FoldedFit,UnfoldedFit)
	wave FoldedFit,UnfoldedFit
	variable FoldedIntercept=str2num(StringbyKey("Intercept",note(FoldedFit),":","\r"))
	variable FoldedSlope=str2num(StringbyKey("Slope",note(FoldedFit),":","\r"))
	variable UnFoldedIntercept=str2num(StringbyKey("Intercept",note(UnFoldedFit),":","\r"))
	variable UnFoldedSlope=str2num(StringbyKey("Slope",note(UnFoldedFit),":","\r"))
	
	return Intercepts(FoldedSlope,UnFoldedSlope,FoldedIntercept,UnFoldedIntercept)
end

Static Function/C BellEvansFit(WaveIn,FitOut)
	wave WaveIn,FitOut
	duplicate/free wavein Logged,TestOut
	Logged=ln(wavein)
	CurveFit/Q/W=2  line Logged
	wave w_coef,w_sigma
	TestOut=w_coef[0]+w_coef[1]*x
	TestOut=exp(testOut)
	duplicate/o TestOut FitOut
	variable/C result= cmplx(w_Coef[0],w_coef[1])
	killwaves w_Coef,w_sigma
	return result	
end

Static Function/C BellEvansDistance(WaveIn)
	wave WaveIn
	variable Slope=str2num(StringbyKey("Slope",note(WaveIn),":","\r"))
	variable kbt=4.11433e-21
	return slope*kbt
end