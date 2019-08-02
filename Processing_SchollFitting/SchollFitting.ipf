#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma ModuleName = DE_SchollFitting
#include <Peak AutoFind> version>=5.07

Static function Scholl(w1,w2,fftdist,stepdist,smoothing,percent,[TotalNumber])
	wave w1,w2
	variable fftdist,stepdist ,smoothing,percent,TotalNumber
	variable rate= 1/str2num(stringbykey("RetractVelocity", note(w1),":","\r"))*str2num(stringbykey("NumPtsPerSec", note(w1),":","\r"))/1000
	variable fftLength
	variable stepWindow=stepdist
	//variable stepWindow=  round(rate*stepdist)
	//variable windowlength=round(rate*fftdist);
	variable windowlength=2^fftdist
	if(mod(windowlength,2)==0)
	else
	windowlength-=1
	endif
	variable TimeBack=.005

	wavestats/q w1
	variable datalength=v_npnts
	variable windowedDataLength = ceil(dataLength/stepWindow)
	if (windowLength/2 == round(windowLength/2))
		fftLength = (windowLength/2) + 1;
	else
		fftLength = ceil(windowLength/2);
	endif
	make/o/n=(windowedDataLength+1,fftLength) wfftsignal
	make/o/n=(windowlength+datalength+1) PadData
	PadData=0
	wfftsignal=0
	PadData[windowlength/2+1,windowlength/2+datalength]=w1[p-windowlength/2-1]
	make/o/n=(windowlength) dataframe
	make/o/n=(fftLength) fftValues
	fftvalues=0
	dataframe=0
	variable i,ii
	for (i = 0;i<=datalength;i+=stepWindow)

		duplicate/O/R=[i,i+windowlength-1] PadData,dataFrame
		FFT /Z /out=3  /DEST=fftValues  dataFrame
		wfftSignal[ii][] = fftValues[q]
		ii = ii + 1;

	endfor

	make/o/n=(windowedDataLength+1) coefficientsum
	coefficientsum=0
	i=1
	for (i=0;i<=(fftLength-1);i+=2)

		coefficientsum[]=coefficientsum[p]+wfftsignal[p][i]
	endfor
	coefficientSum = coefficientSum/windowLength;
	variable/c Smooths=EstPeakNoiseAndSmfact(coefficientsum, 0, (windowedDataLength))
	killwaves PadData,dataFrame,fftValues,wfftsignal
	//Variable peaksFound= AutoFindPeaks(coefficientsum, 0, (windowedDataLength),real(smooths),imag(smooths),20)
	Variable peaksFound= AutoFindPeaks(coefficientsum, 0, (windowedDataLength),real(smooths),smoothing,20)

	if( peaksFound > 0 )
		WAVE W_AutoPeakInfo
		// Remove too-small peaks
		peaksFound= TrimAmpAutoPeakInfo(W_AutoPeakInfo,percent/100)		
		//print num2str(peaksFound)+" Peaks Found"

	else
		//print "No Peaks Found"
		return -1
	endif
	
	if( !ParamIsDefault(TotalNumber))
		if(peaksFound<=TotalNumber)
		
		else
			do
				
				make/o/n=(dimsize(W_AutoPeakInfo,0)) Garbage
				Garbage=W_AutoPeakInfo[p][2]
				wavestats/q Garbage
				print v_minloc
				deletepoints/M=0 v_minloc,1, W_AutoPeakInfo
				PeaksFound-=1
			while(PeaksFound>TotalNumber)
		endif


	endif
	
	make/o/n=(peaksfound) Peaks
	Peaks=W_AutoPeakInfo[p][0]
	ReorderWave(Peaks)
	make/o/n=(peaksfound,2) PeakLoc
	make/o/n=(peaksfound,3) Ruptures
	//PeakLoc=w2[W_AutoPeakInfo[p][0]*stepWindow]
	PeakLoc[][0]=w1[Peaks[p]*stepWindow]
	if(cmpstr(nameofwave(w1),nameofwave(w2))==0)

	PeakLoc[][1]=Peaks[p]*dimdelta(w1,0)+dimoffset(w1,0)
	else
	PeakLoc[][1]=w2[Peaks[p]*stepWindow]
	endif

	
end


Static function ReorderWave(w1)
	wave w1
	duplicate/o w1 Hold1
	Hold1=.11
	variable i
	wavestats/q w1
variable ends=v_npnts
	for(i=0;i<(ends);i+=1)
		wavestats/q w1
		Hold1[i]=w1[v_minloc]
		w1[v_minloc]=Nan
	endfor
	w1=Hold1
	killwaves hold1
	
end


Static Function DE_Scholl_FilCor(RawWave,FiltWave,Point,Decimation)
	wave RawWave,FiltWave
	variable Point,Decimation
	variable Result=FiltWave[Point]
	variable GuessPoint=Point*Decimation
	FindLevel/P/Q/R=[GuessPoint-Decimation,GuessPoint+Decimation] RawWave Result
	return (V_levelX)

end

Static Function DE_Scholl_MaxRange(ForceWave,Point)
	wave ForceWave
	variable Point
	wavestats/q/R=[point-5,point] ForceWave
	return x2pnt(ForceWave,v_minLoc)
end

Static Function/C LineFits(ForceWave,SepWave,FiltSep,decimation,point)
	wave ForceWave,SepWave,FiltSep
	variable decimation,point

	variable CorPoint=DE_Scholl_FilCor(SepWave,FiltSep,Point,Decimation)
	variable maxpoint=DE_Scholl_MaxRange(ForceWave,CorPoint)
	wave Parms
	make/o/n=3 Parms
	variable Offset= str2num(stringbykey("NumPtsPerSec",note(forcewave),":","\r"))*.003/str2num(stringbykey("RetractVelocity",note(forcewave),":","\r"))
	CurveFit/Q/NTHR=0 line  ForceWave[maxpoint-Offset,maxpoint]
	wave w_coef
	variable maxslope=w_coef[1]
	variable maxintercept=w_coef[0]
	CurveFit/Q/NTHR=0 line  ForceWave[CorPoint-1,CorPoint+1] 
	variable Rupslope=w_coef[1]
	variable Rupintercept=w_coef[0]
	variable magic=((Rupintercept-MaxIntercept)/(MaxSlope-RupSlope))
	variable magic2=MaxSlope*((Rupintercept-MaxIntercept)/(MaxSlope-RupSlope))+maxintercept
	return cmplx(magic,magic2)
 
end

Static Function Loading(ForceWave,SepWave,FiltSep,decimation,point)
	wave ForceWave,SepWave,FiltSep
	variable decimation,point

	variable CorPoint=DE_Scholl_FilCor(SepWave,FiltSep,Point,Decimation)
	variable maxpoint=DE_Scholl_MaxRange(ForceWave,CorPoint)
	variable Offset= str2num(stringbykey("NumPtsPerSec",note(forcewave),":","\r"))*.003/str2num(stringbykey("RetractVelocity",note(forcewave),":","\r"))
	CurveFit/Q/NTHR=0 line  ForceWave[maxpoint-Offset,maxpoint] 
	wave w_coef
	variable maxslope=w_coef[1]

	return maxslope
 
end

Static function Convert(number)
	variable Number
	string NumberFull
	sprintf  NumberFull, "%04g",Number

	string Origin="Image"+NUmberFull+"Defl_Ret"
	wave W_Origin=$Origin
	string DefLV="Image"+NUmberFull+"DeflV_Away"
	string ZSnsr="Image"+NUmberFull+"Zsnsr_Away"
	string Defl="Image"+NUmberFull+"Defl_Away"


	if(Waveexists($DefLV)==0)
		wave HBSaveDef
		rename HBSaveDef $DefLV
		wave HBSaveZsn
		rename HBSaveZsn $ZSnsr
	else
	endif
	wave W_DeflV=$DeflV
	wave W_ZSnsr=$ZSnsr
	duplicate/o W_DeflV $Defl
	wave W_Defl=$Defl 

	variable ZSens=str2num(stringbykey("ZLVDTSENS",note(W_origin),":","\r"))
	variable Invols=str2num(stringbykey("Invols",note(W_origin),":","\r"))

	W_Defl=Invols*W_DeflV
	W_ZSnsr*=ZSens
end

