#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma modulename=DE_RateMapProcessing
#include ":DE_BellEvans"
#include ":\Misc_PanelPrograms\Panel Progs"
#include ":\Misc_PanelPrograms\AsylumNaming"
Static Function PointToFolder(DRef)

	DFREF DRef
	DFREF dfSav= GetDataFolderDFR()
	setdatafolder DRef
	String AllFoldHist=wavelist("*AccFold_Hist",";","")
	String AllUnFoldHist=wavelist("*AccUnFold_Hist",";","")
	variable n,m


	make/free/n=0 MaxForces,MinForces
	for(n=0;n<itemsinlist(AllFoldHist);n+=1)
		wave Hist=$stringfromlist(n,AllFoldHist)
		insertpoints 0,1, MaxForces,MinForces
		MinForces[0]=pnt2x(Hist,numpnts(Hist)-1)
		MaxForces[0]=	pnt2x(Hist,0)
	endfor
	variable StartBin=wavemin(Maxforces)
	variable endbin=wavemax(MinForces)
	variable numbins=round((-startbin+endbin)/5e-13)
	make/free/n=(numbins+1) FreeFoldHist,FreeFoldRate
	FreeFoldHist=0
	FreeFoldRate=0
	SetScale/P x StartBin,5e-13,"", FreeFoldHist
	SetScale/P x StartBin,5e-13,"", FreeFoldRate

	for(n=0;n<itemsinlist(AllFoldHist);n+=1)
		wave Hist=$stringfromlist(n,AllFoldHist)
		wave Rates=$replacestring("Hist",stringfromlist(n,AllFoldHist),"Rate")

		duplicate/free Hist CurrHist
		duplicate/free Rates CurrRate

		PadHistoGrams(CurrHist,StartBin,EndBin,5e-13)
		PadHistoGrams(CurrRate,StartBin,EndBin,5e-13)
		for(m=0;m<numpnts(CurrHist);m+=1)
			if(numtype(CurrRate[m])!=0)
				CurrRate[m]=0
				CurrHist[m]=0

			endif
									if(CurrHist[m]<=2)
				CurrRate[m]=0
				CurrHist[m]=0

			endif
		endfor
		FreeFoldHist+=CurrHist
		FreeFoldRate+=CurrRate[p]*CurrHist
	endfor
	FreeFoldRate/=FreeFoldHist
	
	//////
	make/free/n=0 MaxForces,MinForces
	for(n=0;n<itemsinlist(AllUnFoldHist);n+=1)
		wave Hist=$stringfromlist(n,AllUnFoldHist)
		insertpoints 0,1, MaxForces,MinForces
		MinForces[0]=pnt2x(Hist,numpnts(Hist)-1)
		MaxForces[0]=	pnt2x(Hist,0)
	endfor

	StartBin=wavemin(Maxforces)
	endbin=wavemax(MinForces)
	numbins=round((-startbin+endbin)/10e-13)

	make/free/n=(numbins+1) FreeUnFoldHist,FreeUnFoldRate
	FreeUnFoldHist=0
	FreeUnFoldRate=0
	SetScale/P x StartBin,10e-13,"", FreeUnFoldHist
	SetScale/P x StartBin,10e-13,"", FreeUnFoldRate
	
	for(n=0;n<itemsinlist(AllUnFoldHist);n+=1)
		wave Hist=$stringfromlist(n,AllUnFoldHist)
		wave Rates=$replacestring("Hist",stringfromlist(n,AllUnFoldHist),"Rate")

		duplicate/free Hist CurrHist
		duplicate/free Rates CurrRate

		PadHistoGrams(CurrHist,StartBin,EndBin,10e-13)
		PadHistoGrams(CurrRate,StartBin,EndBin,10e-13)
		for(m=0;m<numpnts(CurrHist);m+=1)

			if(numtype(CurrRate[m])!=0)
				CurrRate[m]=0
				CurrHist[m]=0

			endif
						if(CurrHist[m]<=2)
				CurrRate[m]=0
				CurrHist[m]=0

			endif
		endfor

		FreeUnFoldHist+=CurrHist
		FreeUnFoldRate+=CurrRate[p]*CurrHist
	endfor
	FreeUnFoldRate/=FreeUnFoldHist
	string basename= stringfromlist(0,nameofwave(Hist),"_")
	
	//PareOut(FreeFoldRate,FreeFoldHist)
	//	PareOut(FreeUnFoldRate,FreeUnFoldHist)

	duplicate/o FreeFoldRate $(basename+"_FreeFoldRate")
	duplicate/o FreeFoldHist $(basename+"_FreeFoldHist")
	duplicate/o FreeUnFoldRate $(basename+"_FreeUnFoldRate")
	duplicate/o FreeUnFoldHist $(basename+"_FreeUnFoldHist")

	/////
	//killwaves CurrHist,CurrRate

	setdatafolder dfSav 

end

Static Function PareOut(RateIn,HistIn)

	wave RateIn,HistIn
	variable n,Threshold
	//Threshold=max(sum(HistIN)/200,2)
	Threshold=3
	for(n=0;n<numpnts(RateIn);n+=1)
	
		if(HistIn[n]<Threshold)
			Ratein[n]=NaN
		endif
	
	endfor

end

Static Function MakeNiceRatePlot(DRef)

	DFREF DRef
	DFREF dfSav= GetDataFolderDFR()
	setdatafolder DRef
	String AllFoldHist=wavelist("*AccFold_Hist",";","")
	String AllUnFoldHist=wavelist("*AccUnFold_Hist",";","")
	String AvgFoldRate=wavelist("*FreeFoldRate",";","")
	String AvgUnFoldRate=wavelist("*FreeUnFoldRate",";","")
	variable n
	String BaseName=stringfromlist(0,stringfromlist(0,AllFoldHist,";"),"_")
	DoWindow $(BaseName+"_Rate")

	if(V_flag==1)
		killwindow $(BaseName+"_Rate")
	endif

	display/N=$(BaseName+"_Rate")
	variable speed=0
	variable/c FitValues
	String NewFucker
	make/o/n=(0,5) $(BaseName+"_FitParms")
	wave FitParms=$(BaseName+"_FitParms")
	for(n=0;n<itemsinlist(AllFoldHist);n+=1)
		wave Rates=$replacestring("Hist",stringfromlist(n,AllFoldHist),"Rate")
		appendtograph/W=$(BaseName+"_Rate")  Rates
		ModifyGraph/W=$(BaseName+"_Rate") mode($nameofwave(Rates))=3,marker($nameofwave(Rates))=16
		duplicate/o Rates $(nameofwave(Rates)+"_Fit")
		wave Fit=$(nameofwave(Rates)+"_Fit")
		FitValues=DE_BellEvans#BellEvansFit(Rates,Fit)
		NewFucker=""
		NewFucker=replacestringbykey("Intercept",NewFucker,num2str(real(fitvalues)),":","\r")
		NewFucker=replacestringbykey("Slope",NewFucker,num2str(imag(fitvalues)),":","\r")
		note/k Fit, NewFucker
		Appendtograph/W=$(BaseName+"_Rate") Fit
		
		speed= str2num(stringfromlist(1,nameofwave(rates),"_"))
		
		switch(speed)
			case 20:
				ModifyGraph/W=$(BaseName+"_Rate") rgb($nameofwave(Rates))=(14135,32382,47288),rgb($nameofwave(Fit))=(14135,32382,47288)
				break
			case 50:
				ModifyGraph/W=$(BaseName+"_Rate") rgb($nameofwave(Rates))=(19789,45489,19018),rgb($nameofwave(Fit))=(19789,45489,19018)
				break
			case 100:
				ModifyGraph/W=$(BaseName+"_Rate") rgb($nameofwave(Rates))=(58596,6682,7196),rgb($nameofwave(Fit))=(58596,6682,7196)
				break
			default:
				break
		endswitch
		ModifyGraph/W=$(BaseName+"_Rate") mode($nameofwave(Rates))=3,marker($nameofwave(Rates))=16
	endfor
	
	for(n=0;n<itemsinlist(AllUnfoldHist);n+=1)
		wave Rates=$replacestring("Hist",stringfromlist(n,AllunFoldHist),"Rate")
		appendtograph  Rates
		ModifyGraph/W=$(BaseName+"_Rate")  mode($nameofwave(Rates))=3,marker($nameofwave(Rates))=19
		
		duplicate/o Rates $(nameofwave(Rates)+"_Fit")
		wave Fit=$(nameofwave(Rates)+"_Fit")
		FitValues=DE_BellEvans#BellEvansFit(Rates,Fit)
		NewFucker=""
		NewFucker=replacestringbykey("Intercept",NewFucker,num2str(real(fitvalues)),":","\r")
		NewFucker=replacestringbykey("Slope",NewFucker,num2str(imag(fitvalues)),":","\r")
		note/k Fit, NewFucker
		Appendtograph/W=$(BaseName+"_Rate") Fit

		speed= str2num(stringfromlist(1,nameofwave(rates),"_"))
		switch(speed)
			case 20:
				ModifyGraph/W=$(BaseName+"_Rate") rgb($nameofwave(Rates))=(14135,32382,47288),rgb($nameofwave(Fit))=(14135,32382,47288)
				break
			case 50:
				ModifyGraph/W=$(BaseName+"_Rate") rgb($nameofwave(Rates))=(19789,45489,19018),rgb($nameofwave(Fit))=(19789,45489,19018)
				break
			case 100:
				ModifyGraph/W=$(BaseName+"_Rate") rgb($nameofwave(Rates))=(58596,6682,7196),rgb($nameofwave(Fit))=(58596,6682,7196)
				break
			default:
				break
		endswitch
		wave UnfoldedFit=$replacestring("AccUnFold",Nameofwave(Fit),"AccFold")
		insertpoints/M=0 0,1, FitParms
		FitParms[0][0]=speed
		FitParms[0][1]=real(DE_BellEvans#CalcImportantParms(Fit,UnfoldedFit))
		FitParms[0][2]=imag(DE_BellEvans#CalcImportantParms(Fit,UnfoldedFit))
		//BellEvansDistance(WaveIn)
	FitParms[0][3]=DE_BEllEvans#BellEvansDistance(Fit)
		FitParms[0][4]=DE_BEllEvans#BellEvansDistance(UnfoldedFit)
	endfor
	
	appendtograph/W=$(BaseName+"_Rate") $stringfromlist(0,AvgFoldRate)
	wave Rates=$stringfromlist(0,AvgFoldRate)
	duplicate/o Rates $(nameofwave(Rates)+"_Fit")
	wave Fit=$(nameofwave(Rates)+"_Fit")
	wave FoldedFit=$(nameofwave(Rates)+"_Fit")

	FitValues=DE_BellEvans#BellEvansFit(Rates,Fit)
	NewFucker=""
	NewFucker=replacestringbykey("Intercept",NewFucker,num2str(real(fitvalues)),":","\r")
	NewFucker=replacestringbykey("Slope",NewFucker,num2str(imag(fitvalues)),":","\r")
	note/k Fit, NewFucker
	Appendtograph/W=$(BaseName+"_Rate") Fit
	ModifyGraph/W=$(BaseName+"_Rate") rgb($nameofwave(Fit))=(0,0,0)

	appendtograph/W=$(BaseName+"_Rate") $stringfromlist(0,AvgunFoldRate)
	wave Rates=$stringfromlist(0,AvgunFoldRate)
	duplicate/o Rates $(nameofwave(Rates)+"_Fit")
	wave Fit=$(nameofwave(Rates)+"_Fit")
	FitValues=DE_BellEvans#BellEvansFit(Rates,Fit)
	NewFucker=""
	NewFucker=replacestringbykey("Intercept",NewFucker,num2str(real(fitvalues)),":","\r")
	NewFucker=replacestringbykey("Slope",NewFucker,num2str(imag(fitvalues)),":","\r")
	note/k Fit, NewFucker
	Appendtograph/W=$(BaseName+"_Rate") Fit
	ModifyGraph/W=$(BaseName+"_Rate") rgb($nameofwave(Fit))=(0,0,0)

	ModifyGraph/W=$(BaseName+"_Rate")  mode($stringfromlist(0,AvgFoldRate))=3,marker($stringfromlist(0,AvgFoldRate))=16
	ModifyGraph/W=$(BaseName+"_Rate")  mode($stringfromlist(0,AvgunFoldRate))=3,marker($stringfromlist(0,AvgunFoldRate))=19

	ModifyGraph/W=$(BaseName+"_Rate") rgb($stringfromlist(0,AvgFoldRate))=(0,0,0),rgb($stringfromlist(0,AvgunFoldRate))=(0,0,0)

	ModifyGraph/W=$(BaseName+"_Rate") log(left)=1,useMrkStrokeRGB=1
	insertpoints/M=0 0,1, FitParms
	FitParms[0][0]=0
	FitParms[0][1]=real(DE_BellEvans#CalcImportantParms(FoldedFit,Fit))
	FitParms[0][2]=imag(DE_BellEvans#CalcImportantParms(FoldedFit,Fit))
	FitParms[0][3]=DE_BEllEvans#BellEvansDistance(Fit)
	FitParms[0][4]=DE_BEllEvans#BellEvansDistance(FoldedFit)
	DoWindow $(BaseName+"_Parms")

	if(V_flag==1)
		killwindow $(BaseName+"_Parms")
	endif
	
	edit/N=$(BaseName+"_Parms") FitParms
	setdatafolder dfSav 

end




Static Function MakeAverageRatePlots()
	String RelevantFolders= DE_PanelProgs#PrintAllFolders_String("*FreeFold*")
	DFREF dfSav= GetDataFolderDFR()
	variable FolderstoVisit=itemsinlist(RelevantFolders)
	variable n=0,one,two,three
	DoWindow AverageRates

	if(V_flag==1)
		killwindow AverageRates
	endif
	display/N=AverageRates
	for(n=0;n<FolderstoVisit;n+=1)
		string AvgFoldedRate,AvgUnFoldedRate
		DFREF FolderToVisit=$stringfromlist(n,RelevantFolders,";")
		setdatafolder FolderToVisit
		AvgFoldedRate= stringfromlist(0,wavelist("*FreeUnfold*Rate",";",""))
		AvgUnFoldedRate=stringfromlist(0,wavelist("*FreeFold*Rate",";",""))
		one=str2num(stringfromlist( 0,ColorList(n),","))
		two=str2num(stringfromlist( 1,ColorList(n),","))
		three=str2num(stringfromlist( 2,ColorList(n),","))

		Appendtograph/W=AverageRates $AvgFoldedRate
		ModifyGraph/W=AverageRates  mode($(AvgFoldedRate))=3,marker($(AvgFoldedRate))=16,rgb($(AvgFoldedRate))=(one,two,three)
		AppendtoGraph/W=AverageRates $AvgUnFoldedRate
		ModifyGraph/W=AverageRates  mode($(AvgUnFoldedRate))=3,marker($AvgUnFoldedRate)=19,rgb($(AvgunFoldedRate))=(one,two,three)

	endfor
	ModifyGraph/W=AverageRates log(left)=1,useMrkStrokeRGB=1

	setdatafolder dfSav 

end


Static Function PadHistoGrams(Histogram1,StartBin,EndBin,step)
	wave Histogram1
	variable StartBin,EndBin,step


	variable Start1=	pnt2x(Histogram1,0)
	variable End1=pnt2x(Histogram1,numpnts(Histogram1)-1)

	variable startdif=Start1-StartBin
	variable enddif=end1-endBin

	if(abs(startdif)<1e-13)
		
	elseif(StartBin<start1)
		insertpoints 0,round((-startbin+start1)/step), Histogram1

	else
		print "ERROR"
		return -1
	endif

	if(abs(enddif)<1e-13)

	elseif(endbin>end1)

		insertpoints numpnts(Histogram1),round((endbin-end1)/step),Histogram1
		SetScale/P x startbin,step,"", Histogram1

	else
		print "ERROR"
		return -1
	endif
	SetScale/P x startbin,step,"", Histogram1

end

Static Function/S ColorList(n)
	variable n
	String ResString
	switch (n)
		case 0:
		ResString="14135,32382,47288"
		break
		case 1:
	ResString="19789,45489,19018"
		break
		case 2:
	ResString="58596,6682,7196"
		break
		case 3:
	ResString="44253,29492,58982"
		break
		default:
	ResString="65535,32639,0"
	break
	endswitch
	return ResString

end

Static Function ButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	Strswitch(ba.ctrlname)
			string saveDF
			
		case "de_AvgRMaps_button0": 
			switch( ba.eventCode )
				case 2: // mouse up
				controlinfo/W=AvgRMaps de_AvgRMaps_popup0
				PointToFolder($S_Value)
				MakeNiceRatePlot($S_Value)
					break
				case -1: // control being killed
					break
			endswitch
			break
			case "de_AvgRMaps_button1": 
			switch( ba.eventCode )
				case 2: // mouse up
				MakeAverageRatePlots()
					break
				case -1: // control being killed
					break
			endswitch
			break
		
			
	endswitch
	return 0
End

Window AverageRateMapsPanel() : Panel

	PauseUpdate; Silent 1		// building window...
	NewPanel/N=AvgRMaps /W=(0,0,300,100)

	PopupMenu de_AvgRMaps_popup0,pos={2,2},size={125,21},Title="Target Folder"
	PopupMenu de_AvgRMaps_popup0,mode=1,popvalue="X",value= #"DE_PanelProgs#ListFolders()"

	Button de_AvgRMaps_button0,pos={2,30},size={50,21},proc=DE_RateMapProcessing#ButtonProc,title="GO!"
	Button de_AvgRMaps_button1,pos={2,55},size={150,21},proc=DE_RateMapProcessing#ButtonProc,title="MakeAveragePlot"

EndMacro

Static Function/S ListWaves(ControlStr,SearchStr)
	string ControlStr,SearchStr
	String saveDF

	saveDF = GetDataFolder(1)
	controlinfo $ControlStr
	SetDataFolder s_value
	String list = WaveList(SearchStr, ";", "")
	SetDataFolder saveDF
	return list

end

Menu "Ramp"
	"Open Rate map", AverageRateMapsPanel()


	
end

Static Function MakeNiceFits(DRef)

	DFREF DRef
	DFREF dfSav= GetDataFolderDFR()
	setdatafolder DRef
	
	
	
//	String AllFoldHist=wavelist("*AccFold_Hist",";","")
//	String AllUnFoldHist=wavelist("*AccUnFold_Hist",";","")
//	String AvgFoldRate=wavelist("*FreeFoldRate",";","")
//	String AvgUnFoldRate=wavelist("*FreeUnFoldRate",";","")
//	variable n
//	String BaseName=stringfromlist(0,stringfromlist(0,AllFoldHist,";"),"_")
//	DoWindow $(BaseName+"_Rate")
//
//	if(V_flag==1)
//		killwindow $(BaseName+"_Rate")
//	endif
//
//	display/N=$(BaseName+"_Rate")
//	variable speed=0
//
//	for(n=0;n<itemsinlist(AllFoldHist);n+=1)
//		wave Rates=$replacestring("Hist",stringfromlist(n,AllFoldHist),"Rate")
//		appendtograph  Rates
//		ModifyGraph/W=$(BaseName+"_Rate") mode($nameofwave(Rates))=3,marker($nameofwave(Rates))=16
//
//		speed= str2num(stringfromlist(1,nameofwave(rates),"_"))
//		
//		switch(speed)
//			case 20:
//				ModifyGraph/W=$(BaseName+"_Rate") rgb($nameofwave(Rates))=(14135,32382,47288)
//				break
//			case 50:
//				ModifyGraph/W=$(BaseName+"_Rate") rgb($nameofwave(Rates))=(19789,45489,19018)
//				break
//			case 100:
//				ModifyGraph/W=$(BaseName+"_Rate") rgb($nameofwave(Rates))=(58596,6682,7196)
//				break
//			default:
//				break
//		endswitch
//		ModifyGraph/W=$(BaseName+"_Rate") mode($nameofwave(Rates))=3,marker($nameofwave(Rates))=16
//	endfor
//	
//	for(n=0;n<itemsinlist(AllUnfoldHist);n+=1)
//		wave Rates=$replacestring("Hist",stringfromlist(n,AllunFoldHist),"Rate")
//		appendtograph  Rates
//		ModifyGraph/W=$(BaseName+"_Rate")  mode($nameofwave(Rates))=3,marker($nameofwave(Rates))=19
//		speed= str2num(stringfromlist(1,nameofwave(rates),"_"))
//		switch(speed)
//			case 20:
//				ModifyGraph/W=$(BaseName+"_Rate") rgb($nameofwave(Rates))=(14135,32382,47288)
//				break
//			case 50:
//				ModifyGraph/W=$(BaseName+"_Rate") rgb($nameofwave(Rates))=(19789,45489,19018)
//				break
//			case 100:
//				ModifyGraph/W=$(BaseName+"_Rate") rgb($nameofwave(Rates))=(58596,6682,7196)
//				break
//			default:
//				break
//		endswitch
//
//	endfor
//	
//	appendtograph/W=$(BaseName+"_Rate") $stringfromlist(0,AvgFoldRate)
//	appendtograph/W=$(BaseName+"_Rate") $stringfromlist(0,AvgunFoldRate)
//	ModifyGraph/W=$(BaseName+"_Rate")  mode($stringfromlist(0,AvgFoldRate))=3,marker($stringfromlist(0,AvgFoldRate))=19
//	ModifyGraph/W=$(BaseName+"_Rate")  mode($stringfromlist(0,AvgunFoldRate))=3,marker($stringfromlist(0,AvgunFoldRate))=19
//
//	ModifyGraph/W=$(BaseName+"_Rate") rgb($stringfromlist(0,AvgFoldRate))=(0,0,0),rgb($stringfromlist(0,AvgunFoldRate))=(0,0,0)
//	ModifyGraph/W=$(BaseName+"_Rate") log(left)=1,useMrkStrokeRGB=1


	setdatafolder dfSav 

end