#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma ModuleName= DE_Convol

Static Function ConvolveSeries(ConvolutedHist,OutWave,[OptStarting])
	wave ConvolutedHist,OutWave,OptStarting
	variable NumELPoints=numpnts(ConvolutedHist)
	variable PointCount=0
	variable	MaxIterations=  300
	make/o/n=(MaxIterations) ChiOut
	variable	mu=2
	Variable DeconvCount,ChiLast=1,ChiNow=1

	Duplicate/o ConvolutedHist, StartSoln, CurrentSoln, Conv, DeconvEL,Relaxn,constructprob,PreviousSoln,DiffK,StartDist
	Variable AStartSoln=Sum(ConvolutedHist)
	StartDist=ConvolutedHist/AStartSoln
	if(ParamIsDefault(OptStarting))
		StartSoln=StartDist
	else
		variable AOptStart=sum(OptStarting)
		duplicate/o OptStarting StartSoln
		StartSoln/=AOptStart
	endif
	PreviousSoln=StartSoln
	CurrentSoln=PreviousSoln
	For(DeconvCount=0;DeconvCount<MaxIterations;DeconvCount+=1)
		ChiLast=ChiNow;

		For(PointCount=0;PointCount<NumELPoints;PointCount+=1)
			// Get PSF wave
			Wave PSF=$("psf2_"+num2str(PointCount))
			CurrentSoln=PreviousSoln
			Duplicate/O CurrentSoln,Conv
					
			// Convolve current guess with the PSF
			Convolve PSF, Conv				
			// Remove the extra points from deconvolution
		
			// Trying a new method to get rid of extra points.  Worked well in my conv test.
			Variable OriginalLength=DimSize(StartSoln,0)
					
			Variable ConvLength=DimSize(Conv,0)
			Variable NumPntsRemove=((ConvLength-OriginalLength)/2)
			DeletePoints floor(ConvLength-NumPntsRemove), NumPntsRemove, Conv
			DeletePoints 0, ceil(NumPntsRemove), Conv				//Return wave to original length, normalise
			// Normalize convoluted guess
			Variable AConv=Sum(Conv)
			Conv/=AConv
			constructprob[PointCount]=Conv[PointCount]

		endfor
		duplicate/o constructprob conv,NewSoln,convtest

		Relaxn = mu*(1 - 2*abs(CurrentSoln - 0.5))
		//Relaxn = mu*CurrentSoln*(1-CurrentSoln )

		DiffK = (StartDist-Conv)
		Duplicate/O DiffK,DiffK_smth
		•Smooth 3, DiffK_smth
		NewSoln=CurrentSoln+Relaxn*DiffK_smth
		Variable ANewSoln=Sum(NewSoln)
		NewSoln/=ANewSoln
		PreviousSoln=NewSoln

		Duplicate/O DiffK ChiSqr				
		ChiSqr = (DiffK)^2
		ChiNow = (AStartSoln-AConv)*Sum(ChiSqr)
		ChiOut[DeconvCount]= ChiNow
		If (ChiNow > ChiLast)
			//chinow=chilast
		//	Print "error =",ChiNow,ChiLast, PointCount, DeconvCount+1

		endif

	endfor
	duplicate/o constructprob OutWave

	//killwaves StartSoln, CurrentSoln, Conv, DeconvEL,Relaxn,constructprob,PreviousSoln,DiffK
end

Static Function MakingMatchingHistogram(RawTrace,Smtrace)


	wave RawTrace,Smtrace

	Variable ForceBinSize=50e-3
	String RHistogramName=nameofwave(RawTrace)+"_Hist"
	String SHistogramName=nameofwave(SMTrace)+"_Hist"

	WaveStats/Q RawTrace
	Variable ForceMin=V_min
	Variable NumForceBins=Ceil((V_max-ForceMin)/ForceBinSize)
	Make/O/N=1 $RHistogramName
	Wave RHistogram=$RHistogramName
	// Make probability distribution function wave using Histogram.
	// Flags being used are /C for moving x location half a bin width over
	// /P to give probability instead of raw counts
	// /B sets the min, bin size and number of bins.
	Histogram/B={ForceMin,ForceBinSize,NumForceBins}/C RawTrace, $RHistogramName
	Variable TotalCount=Sum(RHistogram)
	RHistogram/=TotalCount

	Make/O/N=1 $SHistogramName
	Wave SHistogram=$SHistogramName
	// Make probability distribution function wave using Histogram.
	// Flags being used are /C for moving x location half a bin width over
	// /P to give probability instead of raw counts
	// /B sets the min, bin size and number of bins.
	Histogram/B={ForceMin,ForceBinSize,NumForceBins}/C SmTrace, $SHistogramName
	TotalCount=Sum(SHistogram)
	SHistogram/=TotalCount



end

Static Function MakePSF(ConvolutedHist)
	wave ConvolutedHist
	Variable NumELPoints=DimSize(ConvolutedHist,0),PSFCount=0
	Variable ELMidPoint=pnt2x(ConvolutedHist,NumELPoints/2)
	variable PSFWidth_Slope=17e-3
	variable PSFWidth_Offset=.074
	//variable PSFWidth_Slope=0
	//variable PSFWidth_Offset=.500/sqrt(2)
	For(PSFCount=0;PSFCount<NumELPoints;PSFCount+=1)
		String PSFWaveName="psf2_"+num2str(PSFCount)
		Duplicate/O ConvolutedHist, $PSFWaveName
		Wave PSFforPoint=$PSFWaveName
		//Variable PSFWidth=Abs(PSFWidth_Slope*ConvolutedEL[PSFCount]+PSFWidth_Offset)
		Variable PSFWidth=Abs(PSFWidth_Slope*pnt2x(ConvolutedHist,PSFCount)+PSFWidth_Offset)
		PSFforPoint=Gauss(x,ELMidPoint,PSFWidth)
		Note PSFforPoint, num2str(PSFWidth)
		//PSFforPoint=Gauss(x,pnt2x(ConvolutedEL,PSFCount),PSFWidth)
		Variable APSFforPoint=Sum(PSFForPoint)
		PSFforPoint/=APSFforPoint
	EndFor


End