#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma modulename=DE_MultiWLC
#include ":SimpleWLCPrograms"




//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

static function DeltaLC(Force,Sep,L0start,LStep,Lper,Temp,Crs,[Folder])
	wave Force,Sep
	variable L0start,LStep,Lper,Temp
	wave/t Crs
	DFREF Folder

	DFREF saveDFR = GetDataFolderDFR()	// Get reference to current data folder
	variable CStart,CEnd,totalsegs,n,totallength,totalnum,LC

	if(!ParamIsDefault(Folder))
		setdatafolder Folder
	endif

	make/o/n=(0,2) ComSep
	make/o/n=(0) ComForce,ResForce
	variable dims=dimsize(Crs,0)
	for(n=0;n<dims;n+=1)

		CStart=str2num(crs[n][0])
		CEnd=str2num(crs[n][1])

		LC=L0start+n*LStep
		print Lper
		print Temp
		print LC
		totalsegs+=1
		insertpoints (TotalLength), (CEnd-CStart+1), ComSep,ComForce,ResForce
		ComSep[TotalLength,(TotalLength+CEnd-CStart)][0]=Sep[p-TotalLength+CStart]
		ComSep[TotalLength,(TotalLength+CEnd-CStart)][1]=totalsegs-1
		ComForce[TotalLength,(TotalLength+CEnd-CStart)]=Force[p-TotalLength+CStart]
		ResForce[TotalLength,(TotalLength+CEnd-CStart)]=WLC(Sep[p-TotalLength+CStart],Lper,LC,Temp)
////		insertpoints 4,2,FitParms
////		FitParms[4]=TotalLength
		TotalLength+=(Cend-CStart)+1
////
////		FitParms[5]=TotalLength-1
//	
	endfor
//
////	FitParms[0]=TotSegs
////	FitParms[1]=TotalLength
////
//	//if(!ParamIsDefault(Folder))
	setdatafolder saveDFR
//	//endif
//
end


static function FinalLC(Force,Sep,L0Final,LStep,Lper,Temp,Crs,[Folder])
	wave Force,Sep
	variable L0Final,LStep,Lper,Temp
	wave/t Crs
	DFREF Folder

	DFREF saveDFR = GetDataFolderDFR()	// Get reference to current data folder
	variable CStart,CEnd,totalsegs,n,totallength,totalnum,LC

	if(!ParamIsDefault(Folder))
		setdatafolder Folder
	endif

	make/o/n=(0,2) ComSep
	make/o/n=(0) ComForce,ResForce
	variable dims=dimsize(Crs,0)
	for(n=0;n<dims;n+=1)

		CStart=str2num(crs[n][0])
		CEnd=str2num(crs[n][1])
		LC=L0final-(dims-n-1)*LStep
		totalsegs+=1
		insertpoints (TotalLength), (CEnd-CStart+1), ComSep,ComForce,ResForce
		ComSep[TotalLength,(TotalLength+CEnd-CStart)][0]=Sep[p-TotalLength+CStart]
		ComSep[TotalLength,(TotalLength+CEnd-CStart)][1]=totalsegs-1
		ComForce[TotalLength,(TotalLength+CEnd-CStart)]=Force[p-TotalLength+CStart]
		ResForce[TotalLength,(TotalLength+CEnd-CStart)]=WLC(Sep[p-TotalLength+CStart],Lper,LC,Temp)
////		insertpoints 4,2,FitParms
////		FitParms[4]=TotalLength
		TotalLength+=(Cend-CStart)+1
////
////		FitParms[5]=TotalLength-1
//	
	endfor
//
////	FitParms[0]=TotSegs
////	FitParms[1]=TotalLength
////
//	//if(!ParamIsDefault(Folder))
	setdatafolder saveDFR
//	//endif
//
end

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


static function ArbLcs(Force,Sep,Lper,Temp,Crs,Parms,[Folder])
	wave Force,Sep
	variable Lper,Temp
	wave/T Crs,parms
	DFREF Folder

	DFREF saveDFR = GetDataFolderDFR()	// Get reference to current data folder
	variable CStart,CEnd,totalsegs,n,totallength,totalnum,LC

	if(!ParamIsDefault(Folder))
		setdatafolder Folder
	endif

	make/o/n=(0,2) ComSep
	make/o/n=(0) ComForce,ResForce
	variable dims=dimsize(crs,0)
	for(n=0;n<dims;n+=1)

		CStart=str2num(crs[n][0])
		CEnd=str2num(crs[n][1])


		LC=str2num(Parms[n+2][1])*1e-9
		totalsegs+=1
		insertpoints (TotalLength), (CEnd-CStart+1), ComSep,ComForce,ResForce
		ComSep[TotalLength,(TotalLength+CEnd-CStart)][0]=Sep[p-TotalLength+CStart]
		ComSep[TotalLength,(TotalLength+CEnd-CStart)][1]=totalsegs-1
		ComForce[TotalLength,(TotalLength+CEnd-CStart)]=Force[p-TotalLength+CStart]
		ResForce[TotalLength,(TotalLength+CEnd-CStart)]=WLC(Sep[p-TotalLength+CStart],Lper,LC,Temp)
//		insertpoints 4,2,FitParms
//		FitParms[4]=TotalLength
		TotalLength+=(Cend-CStart)+1
//
//		FitParms[5]=TotalLength-1
	
	endfor

//	FitParms[0]=TotSegs
//	FitParms[1]=TotalLength
//
	//if(!ParamIsDefault(Folder))
	setdatafolder saveDFR
	//endif

end

function DE_Fit_MWC_DLC(w,x1,x2) : FitFunc
	wave w
	variable x1,x2


	variable y
	y=WLC(x1,w[0],w[1]+x2*w[2],w[3])
	return y

end
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

function DE_Fit_MWC_LCs(w,x1,x2) : FitFunc
	wave w
	variable x1,x2
	variable y
	y=WLC(x1,w[0],w[2+x2],w[1])
	return y

end
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

function DE_Fit_MWC_LCF(w,x1,x2) : FitFunc
	wave w
	variable x1,x2
	variable y
	variable pers=w[0]
	variable LC=w[1]-x2*w[2]
	variable Temp=w[3]
	y=WLC(x1,pers,Lc,Temp)
	return y

end
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------function DE_Fit_MWC_LCF(w,x1,x2) : FitFunc
function DE_Fit_MWC_LCsps(w,x1,x2) : FitFunc

	wave w
	variable x1,x2

	variable y
	y=WLC(x1,w[1+x2*2],w[2+x2*2],w[0])	
	return y

end
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------