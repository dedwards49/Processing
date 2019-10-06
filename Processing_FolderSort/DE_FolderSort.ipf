#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
Function Test()
	NewPath/o/Q DataPath
	String AllWavesInFolder=SortList(IndexedFile(DataPath, -1, ".ibw"),";",16)
	variable tot=itemsinlist(AllWavesInFolder)
	variable n,namelen,CurNum
	 string CurrentName,BaseName,HistoricalBaseName
	 pathinfo Datapath
	 String BasePath= S_path,FileString
	 variable currentloop,NumPerFolder
	 NumPerFolder=500
	for(n=0;n<tot;n+=1)
		namelen=strlen(stringfromlist(n,AllWavesInFolder))
		FileString=stringfromlist(n,AllWavesInFolder)
		BaseName	= stringfromlist(n,AllWavesInFolder)[0,namelen-9]
		CurNum=str2num(stringfromlist(n,AllWavesInFolder)[namelen-8,namelen-5])
		if(n==0)
			NewPath/C/O CurrentPath BasePath+BaseName+"_0"
			HistoricalBaseName=BaseName
			currentloop=0
		else
			If(cmpstr(HistoricalBaseName,BaseName)!=0)
				NewPath/C/O CurrentPath BasePath+BaseName+"_0"
				currentloop=0
				HistoricalBaseName=BaseName

			elseIf(CurNum>(currentloop+1)*NumPerFolder)
				NewPath/C/O CurrentPath BasePath+BaseName+num2str(currentloop+1)
				currentloop+=1
			endif
		endif
			 pathinfo CurrentPath
			MoveFile/D/P=DataPath FileString as S_path
	endfor
	//killpath Datapath

end