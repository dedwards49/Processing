#pragma rtGlobals=1	// Use modern global access method.
#pragma ModuleName = DE_FECWiggle 
//#include "D:\Devin\Documents\Software\Boulder\Devin'sIgor\ProcessingFEC\FastCorrection"
#include ":\Misc_FECProcessing\FastCorrection"
function CorrFEC(WDef_fit,WZsen_fit,WDef_corr,WZsen_corr,rstar,rend,[NewName,ForceName,SepName,FitType,ResName])
	wave WDef_fit,WZsen_fit,WDef_corr,WZsen_corr
	variable rstar,rend
	string NewName,ForceName,SepName,FitType,ResName
	if( ParamIsDefault(FitType))
		Prompt FitType,"Fit Type",popup,"LinSin;Poly;Spline"
		DoPrompt "Color Trace",FitType
		if (V_Flag)
			return 0									// user canceled
		endif


	endif
	
	//If we don't recognize what the fit type is (or it isn't provided), we ask
	if(cmpstr(FitType,"LinSin")==0&&cmpstr(FitType,"Poly")==0&&cmpstr(FitType,"MeanSub")==0&&cmpstr(FitType,"Spline")==0&&cmpstr(FitType,"None")==0)
		Prompt FitType,"Fit Type",popup,"LinSin;Poly"
		DoPrompt "Color Trace",FitType
		if (V_Flag)
			return 0									// user canceled
		endif
	
	endif
	
	string Name=nameofwave(WDef_corr)
	
	if( ParamIsDefault(ResName))
		ResName=Name+"_Fit"

	endif
	
	duplicate/o WDef_fit $ResName
	
	//Stetched sine-wave fitting
	if(cmpstr(FitType,"LinSin")==0)
		Make/D/N=6/O W_coef
		
		W_coef[0] = {(WDef_fit[rstar]),-00e-9,.3e-9,.2e8,0,0}
		//W_coef[0] = {-1.2074e-006,0.32438,-6.2937e-008,2.9686e+007,2.995,0.020003}
		FuncFit/Q/NTHR=0 linearsin2 W_coef  WDef_fit[rstar,rend] /X=WZsen_fit/R=$ResName
	elseif(cmpstr(FitType,"MeanSub")==0)
		//wavestats/q/r=[rstart,rend[
	

	
	elseif(cmpstr(FitType,"Poly")==0)

		CurveFit/Q/NTHR=0 poly 20, WDef_fit[rstar,rend] /X=WZsen_fit/R=$ResName
	
	elseif(cmpstr(FitType,"Spline")==0)
		make/o/n=0 $ResName
		make/free/n=0 HoldFit

		wave Res=$ResName
		DE_FECWiggle#MakeSmoothingSpline(WDef_fit,WZsen_fit,WZsen_fit,rstar,rend,.9,HoldFit)
		duplicate/o HoldFit Res
		Res=WDef_fit-Res
	elseif(cmpstr(FitType,"None")==0)



	endif


	

	if( ParamIsDefault(NewName))
		NewName=Name+"_Corr"
	endif
	if( ParamIsDefault(ForceName))
		ForceName=Name+"_Force"
	endif
	if( ParamIsDefault(SepName))
		SepName=Name+"_Sep"

	endif
	


	duplicate/o WDef_corr $NewName
	duplicate/o WZsen_corr $sepName

	wave w1=$NewName
	wave w2=$SepName
	duplicate/free w1 SubHold,SubHold2
	if(cmpstr(FitType,"LinSin")==0)
		SubHold=linearsin2calc(w_coef,WZsen_corr(x),WZsen_fit[rend])
	elseif(cmpstr(FitType,"MeanSub")==0)
		

	elseif(cmpstr(FitType,"Poly")==0)

		SubHold=PolyCalc(w_coef,WZsen_corr(x),WZsen_fit[rend])
	elseif(cmpstr(FitType,"Spline")==0)
	//WDef_fit,WZsen_fit,WDef_corr,WZsen_corr
		DE_FECWiggle#MakeSmoothingSpline(WDef_fit,WZsen_fit,WZsen_corr,rstar,rend,.9,SubHold)
		

	elseif(cmpstr(FitType,"None")==0)
		SubHold=0
	endif

	 
	FastOP w1=SubHold2-SubHold
	FastOP w2=(-1)*Wzsen_corr+w1

	duplicate/o w1 $ForceName
	wave w3=$ForceName
	variable spring=str2num(stringbykey("SpringConstant", note(w1),":","\r"))

	if(numtype(spring)==2)
		spring=str2num(stringbykey("Spring Constant", note(w1),":","\r"))
	endif
	if(numtype(spring)==2)
		spring=str2num(stringbykey("K", note(w1),"=","\r"))
	endif
	w3=spring*w3
	wave w_coef,w_sigma
	
	if(cmpstr(FitType,"LinSin")==0)
	killwaves w_coef, w_Sigma

	elseif(cmpstr(FitType,"Poly")==0)
		killwaves w_Sigma,subhold,subhold2,w_paramconfidenceinterval
	elseif(cmpstr(FitType,"None")==0)
		killwaves  w_Sigma,subhold,subhold2
	endif

end


Static function SingleProcess(ParmWave,WZsen_fit,WDef_corr,WZsen_corr,rstar,rend,[NewName,ForceName,SepName,FitType])
	wave ParmWave,WZsen_fit,WDef_corr,WZsen_corr
	variable rstar,rend
	string NewName,ForceName,SepName,FitType

	if( ParamIsDefault(FitType))
		Prompt FitType,"Fit Type",popup,"LinSin;Poly"
		DoPrompt "Color Trace",FitType
		if (V_Flag)
			return 0									// user canceled
		endif


	endif
	
	//If we don't recognize what the fit type is (or it isn't provided), we ask
	if(cmpstr(FitType,"LinSin")==0&&cmpstr(FitType,"Poly")==0&&cmpstr(FitType,"MeanSub")==0)
		Prompt FitType,"Fit Type",popup,"LinSin;Poly"
		DoPrompt "Color Trace",FitType
		if (V_Flag)
			return 0									// user canceled
		endif
	
	endif
	
	string Name=nameofwave(WDef_corr)
	
	
	
	if( ParamIsDefault(NewName))
		NewName=Name+"_Corr"
	endif
	if( ParamIsDefault(ForceName))
		ForceName=Name+"_Force"
	endif
	if( ParamIsDefault(SepName))
		SepName=Name+"_Sep"

	endif
	


	duplicate/o WDef_corr $NewName
	duplicate/o WZsen_corr $sepName
	
	wave w1=$NewName
	wave w2=$SepName

	duplicate/free w1 SubHold,SubHold2

	if(cmpstr(FitType,"LinSin")==0)
		SubHold=linearsin2calc(ParmWave,WZsen_corr(x),WZsen_fit[rend])
		
	elseif(cmpstr(FitType,"MeanSub")==0)
		

	elseif(cmpstr(FitType,"Poly")==0)

		SubHold=PolyCalc(ParmWave,WZsen_corr(x),WZsen_fit[rend])

	endif
		
	FastOP w1=SubHold2-SubHold
	FastOP w2=(-1)*Wzsen_corr+w1

	duplicate/o w1 $ForceName
	wave w3=$ForceName
	variable spring=str2num(stringbykey("SpringConstant", note(w1),":","\r"))

	if(numtype(spring)==2)
		spring=str2num(stringbykey("Spring Constant", note(w1),":","\r"))
	endif
		if(numtype(spring)==2)
		spring=str2num(stringbykey("K", note(w1),"=","\r"))
	endif
	w3=spring*w3
	
	//killwaves subhold,subhold2

end

Static function FindWiggleParms(WDef_fit,WZsen_fit, ResultWave, rstar,rend,[FitType,Resname])
	wave WDef_fit,WZsen_fit,ResultWave
	variable rstar,rend
	string FitType,Resname
	
	if( ParamIsDefault(FitType))
		Prompt FitType,"Fit Type",popup,"LinSin;Poly"
		DoPrompt "Color Trace",FitType
		if (V_Flag)
			return 0									// user canceled
		endif


	endif
	
	//If we don't recognize what the fit type is (or it isn't provided), we ask
	if(cmpstr(FitType,"LinSin")==0&&cmpstr(FitType,"Poly")==0&&cmpstr(FitType,"MeanSub")==0)
		Prompt FitType,"Fit Type",popup,"LinSin;Poly'"
		DoPrompt "Color Trace",FitType
		if (V_Flag)
			return 0									// user canceled
		endif
	
	endif
	
	string Name=nameofwave(WDef_corr)
	if( ParamIsDefault(ResName))
		ResName=Name+"_Fit"

	endif
	
	duplicate/o WDef_fit $ResName

	
	
	//Stetched sine-wave fitting
	if(cmpstr(FitType,"LinSin")==0)
		Make/D/N=6/O W_coef
		
		W_coef[0] = {(WDef_fit[rstar]),-00e-9,.3e-9,.2e8,0,0}
		//W_coef[0] = {-1.2074e-006,0.32438,-6.2937e-008,2.9686e+007,2.995,0.020003}
		FuncFit/Q/NTHR=0 linearsin2 W_coef  WDef_fit[rstar,rend] /X=WZsen_fit/R=$ResName
	elseif(cmpstr(FitType,"MeanSub")==0)
	
	elseif(cmpstr(FitType,"Poly")==0)

		CurveFit/Q/NTHR=0 poly 20, WDef_fit[rstar,rend] /X=WZsen_fit/R=$ResName

	endif

	wave w_coef
	duplicate/o w_coef ResultWave

	

	

	if(cmpstr(FitType,"LinSin")==0)
		killwaves w_coef, w_Sigma

	elseif(cmpstr(FitType,"Poly")==0)
		killwaves w_coef,w_Sigma,w_paramconfidenceinterval
	endif
	
end


function linearsin2(w,x):fitfunc
	wave w
	variable x
	return w[0]+w[1]*x+(w[2]+w[5]*x)*Sin(w[3]*x+w[4])
end

function linearsin2calc(w,x,cutoff)
	wave w
	variable x,cutoff
	variable y=w[0]+w[1]*x+(w[2]+w[5]*x)*Sin(w[3]*x+w[4])
	if(x>cutoff)
		y=w[0]+w[1]*cutoff+(w[2]+w[5]*cutoff)*Sin(w[3]*cutoff+w[4])
	endif
	return y
end

function PolyCalc(w,x,cutoff)
	wave w
	variable x,cutoff
	wavestats/q w
	variable n=v_npnts
	variable i=0
	variable y=0
	for(i=0;i<n;i+=1)
		y+=w[i]*(x)^(i)
	endfor
	if(x>cutoff)
		y=0
		for(i=0;i<n;i+=1)
			y+=w[i]*(cutoff)^(i)
		endfor
	endif
	return y
end

function GenFastX(slowz,fastdeflV,offset)
	wave slowz,fastdeflV
	variable offset
	variable spacing=dimdelta(fastdeflV,0)

	string DeflName=ReplaceString("DeflV",nameofwave(fastDeflv),"Defl",0)
	string ZName=ReplaceString("DeflV",nameofwave(fastDeflv),"Zsnsr",0)

	duplicate/o fastdeflV $DeflName
	variable xstart=x2pnt(fastdeflV,(pnt2x(slowz,0)+offset))
	wave w1=$DeflName
	deletepoints 0, xstart,w1
	SetScale/P x pnt2x(slowz,0),spacing,"s", w1
	wavestats/q slowz
	variable lastpoint=x2pnt(w1,pnt2x(slowz,(v_npnts-1)))
	wavestats/q w1
	if(lastpoint>=(v_npnts))
		duplicate slowz holdz
		lastpoint=x2pnt(slowz,pnt2x(w1,(v_npnts-1)))
		deletepoints lastpoint, 1e8,holdz

	else
		duplicate slowz holdz

		deletepoints lastpoint, 1e8,w1
	endif
	
	
	variable Invols=str2num(stringbykey("Invols",note(fastdeflV),":","\r"))
	FastOp w1=(Invols)*w1
	duplicate/o w1  $ZName
	wave w2=$ZName
	wavestats/q w2

	Interpolate2/T=1/N=(v_npnts)/Y=w2 holdz
	killwaves holdz

end

Static Function MakeSmoothingSpline(YWavetoFit,XWaveToFit,XOutput,startpnt,endpnt,smoothing,YOut)
	wave YWavetoFit,XWaveToFit,XOutput,YOut
	variable startpnt,endpnt,smoothing
	variable numberofpnts=endpnt-startpnt
	duplicate/free/R=[startpnt,endpnt] YWavetoFit YFreeFit,YResult,YResidual
	duplicate/free/R=[startpnt,endpnt] XWavetoFit XFreeFit
	duplicate/o	XOutput XBlah
	
	CurveFit/Q/NTHR=0 poly 20, YFreeFit /R=YResidual
	wavestats/Q YResidual
	make/free/n=0 YBlah
	make/free/n=0 IIndex,SortedX,YSorted
	SortWithIndex(XBlah,SortedX,IIndex)
	Interpolate2/S=(V_Sdev)/T=3/N=(numberofpnts)/F=(smoothing)/Y=YBLAH/X=SortedX/I=3 XFreeFit,YFreeFit
	ReSortWithIndex(IIndex,YBLAH,YSorted)
	duplicate/o YSorted YOut
	wave W_coef,W_paramConfidenceInterval,W_Sigma
	killwaves XBlah,YBlah,W_coef,W_paramConfidenceInterval,W_Sigma

end

Static Function SortWithIndex(Xin,XOut,IOut)

	wave Xin,XOut,IOut
	duplicate/free Xin IndexWave,XFree
	IndexWave=p
	Sort Xin,XFree,IndexWave
	duplicate/o  XFree XOut
	duplicate/o  IndexWave IOut

end

Static Function ReSortWithIndex(IIn,Xin,XOut)

	wave IIn,Xin,XOut
	duplicate/free Xin XFree
	Sort IIn,XFree
	duplicate/o  XFree XOut

end