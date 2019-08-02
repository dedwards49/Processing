#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma modulename=DE_CTFCMenuMulti
#include "DE_CTFCMenuProcess"

Static Function InsertASortedTimeWave(WaveToInsert,SortedTextWave)

	wave WaveToInsert
	wave/T SortedTextWave
	variable/D TimeValue=DE_TimeFuncs#ReturnTimeFromWave(WaveToInsert)+DE_TimeFuncs#ReturnDateFromWave(WaveToInsert)	
	STRING Datestring=Secs2Date(TimeValue,0)
	sTRING timestring=Secs2Time(tIMEvALUE,0)
	String FolderName= GetWavesDataFolder(WaveToInsert,1)
	if(dimsize(SortedTextWave,0)==0)
		make/T/free/n=(1,4) FreeSorted
		FreeSorted[0][0]=nameofwave(WaveToInsert)
		FreeSorted[0][1]=num2istr(TimeValue)
		FreeSorted[0][2]=Datestring+","+timestring
		FreeSorted[0][3]=FolderName

	else
		duplicate/T/free SortedTextWave FreeSorted
		insertpoints 0,1,FreeSorted
		FreeSorted[0][0]=nameofwave(WaveToInsert)
		FreeSorted[0][1]=num2istr(TimeValue)
		FreeSorted[0][2]=Datestring+","+timestring
		FreeSorted[0][3]=FolderName

		make/D/free/n=(dimsize(FreeSorted,0)) VariableTimes
		make/T/free/n=(dimsize(FreeSorted,0)) SingleName,SingleSec,SingleDates,SingleFolder
		SingleName=FreeSorted[p][0]
		SingleSec=FreeSorted[p][1]
		SingleDates=FreeSorted[p][2]
		SingleFolder=FreeSorted[p][3]
		VariableTimes=str2num(FreeSorted[p][1])
		Sort VariableTimes SingleName
		Sort VariableTimes SingleSec
		Sort VariableTimes SingleDates
		Sort VariableTimes SingleFolder

		FreeSorted[][0]=SingleName[p]
		FreeSorted[][1]=SingleSec[p]
		FreeSorted[][2]=SingleDates[p]
		FreeSorted[][3]=SingleFolder[p]

	endif
	duplicate/T/o FreeSorted SortedTextWave

end

Static Function ScanFolder(Folder,WaveIn)
	
	DFREF Folder
	wave/T wavein
	DFREF savedDF= GetDataFolderDFR()
	SetDataFolder Folder
	
	string WaveListForMe=wavelist("*D*",";","")
	variable num=itemsinlist(WaveListForMe)
	variable n
	for(n=0;n<num;n+=1)
		Wave TempWave=$stringfromlist(n,WaveListForMe)
		if(cmpstr(DetermineWaveType(TempWave),"Skip")==0)
		else
		InsertASortedTimeWave(TempWave,WaveIn)
		endif
	endfor
	
		 WaveListForMe=wavelist("CXRX*",";","")
	 num=itemsinlist(WaveListForMe)
	for(n=0;n<num;n+=1)
		Wave TempWave=$stringfromlist(n,WaveListForMe)
		if(cmpstr(DetermineWaveType(TempWave),"Skip")==0)
		else
		InsertASortedTimeWave(TempWave,WaveIn)
		endif
	endfor
	
	
	SetDataFolder savedDF

end

Static Function ScanAllFolders(WaveIn)
	Wave/T WaveIn
	variable n=0
	string FolderName
	do
	
		FolderName= GetIndexedObjNameDFR(root:, 4,n )

		if (strlen(FolderName) == 0)
			break
		endif
		DFREF Folder=$FolderName
		ScanFolder(Folder,WaveIn)
		n+=1
	while(1>0)

end


Static Function ProcessTheWaveListInOrder(WaveListIn,BaseString)
	Wave/T WaveListIn
	String BaseString
	variable totnum=dimsize(WaveListIn,0)
	variable n=0
	string WaveNameandPath,Type,ZWaveString,NumFour
	if(dimsize(WaveListIn,1)==4)
	print 4
	insertpoints/M=1 4,1, WaveListIn
	endif
	for(n=0;n<totnum;n+=1)
		WaveNameandPath=WaveListIn[n][3]+WaveListIn[n][0]
		wave DefV=$WaveNameAndPath
		Type=DetermineWaveType(DefV)
		sprintf NumFour "%04.4G",n
		WaveListIn[n][4]=BaseString+	NumFour

		Strswitch(Type)
		
		case "Initial":
				wave ZVolt=$replacestring("InitRamp_D",WaveNameandPath,"InitRamp_Z")
				DE_CTFCMenuProc#ProcessVoltWaves(Type,DefV,ZVolt,BaseString,n)
		break
		
		case "Centering":
				wave CXRX=DefV
				wave CXRZ=$replacestring("CXRX",WaveNameandPath,"CXRZ")
				wave CYRY=$replacestring("CXRX",WaveNameandPath,"CYRX")
				wave CYRZ=$replacestring("CXRX",WaveNameandPath,"CYRZ")
				DE_CTFCMenuProc#ProcessCentering(CXRX,CXRZ,CYRY,CYRZ,BaseString,n)
		break
		
		case "Single":
				wave ZVolt=$replacestring("Sramp_D",WaveNameandPath,"Sramp_Z")
				DE_CTFCMenuProc#ProcessVoltWaves(Type,DefV,ZVolt,BaseString,n)
		break
		
		case "Touch":
				wave ZVolt=$replacestring("Touch_D",WaveNameandPath,"Touch_Z")
				DE_CTFCMenuProc#ProcessVoltWaves(Type,DefV,ZVolt,BaseString,n)
		break
		
		case "Multi":
				wave ZVolt=$replacestring("Mramp_D",WaveNameandPath,"Mramp_Z")
				DE_CTFCMenuProc#ProcessVoltWaves(Type,DefV,ZVolt,BaseString,n)
		break
		case "EquilRamp":
				wave ZVolt=$replacestring("Mramp_D",WaveNameandPath,"Mramp_Z")
				DE_CTFCMenuProc#ProcessVoltWaves(Type,DefV,ZVolt,BaseString,n)
		break
		case "Equil":
				wave ZVolt=$replacestring("Mramp_D",WaveNameandPath,"Mramp_Z")
				DE_CTFCMenuProc#ProcessVoltWaves(Type,DefV,ZVolt,BaseString,n)
		break
		
		
		default:
		print nameofwave(DefV)
		break
		endswitch

	endfor
	

end

Static Function/S DetermineWaveType(WaveIn)
	wave WaveIn
	
	if(strsearch(nameofwave(WaveIn),"Init",0)!=-1)
		return "Initial"
	
	
	elseif(strsearch(nameofwave(WaveIn),"CXRX",0)!=-1)
		return "Centering"

	
	elseif(strsearch(nameofwave(WaveIn),"Sramp",0)!=-1)
		return "Single"

	elseif(strsearch(nameofwave(WaveIn),"Mramp",0)!=-1)
		return "Multi"
	
	elseif(strsearch(nameofwave(WaveIn),"Touch",0)!=-1)
		return "Touch"
	
	elseif(strsearch(nameofwave(WaveIn),"SqStep",0)!=-1)
		return "EquilRamp"
	elseif(strsearch(nameofwave(WaveIn),"SJEquil",0)!=-1)
		return "Equil"
	else
		return "Skip"
	
	endif

end