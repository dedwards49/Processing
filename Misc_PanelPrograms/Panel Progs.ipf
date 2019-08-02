#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma ModuleName = DE_PanelProgs


Static Function/S ListFolders()

	string list=DE_PanelProgs#PrintAllFolders_String("*")
	return list
End

Static function ListWaves(Folder,CriteriaString,w1)
DFREF  Folder
string CriteriaString
wave/t w1


variable index,m,n
string ListName,DefName
	Redimension/n=0 w1
	w1=""
do

	Wave/Z w = WaveRefIndexedDFR(Folder, index)
	
	if (!WaveExists(w))
		break
	endif
	if (stringmatch(nameofwave(w),CriteriaString)==1)

		DefName=nameofwave(w)
		//sscanf Defname,"Image%04gForce_Ret",m
		//sprintf Listname, "Trace%04g",m
		Insertpoints n,1, w1
		//w1[n]=ListName
		w1[n]=DefName
		n+=1
	else

	endif

index+=1

while(index>-1)

end

Static Function/S PrintAllFolders(w1,SearchParm)
	wave/t w1
	string SearchParm
	String objName
	variable tote,limits,bigind,WListNum,countup,n
	SetdataFolder root:
	Variable index = 0
	DFREF dfr = GetDataFolderDFR() // Reference to current data folder
	DFREF odfr=dfr	
	String NameListInitial,NameListFinal
	WListNum=itemsinlist(WaveList(SearchParm, ";",""))

	NameListInitial=GetDataFolder(1)
	NameListFinal=""
	
	do
		Tote=CountObjectsDFR(dfr, 4)
		if (tote==0)
		
		else
			limits+=Tote
			do
				SetdataFolder GetIndexedObjNameDFR(dfr, 4, index)
				NameListInitial=AddListItem(GetDataFolder(1),NameListInitial,";",inf)
				
				SetdataFolder dfr
				index+=1
			while(index<tote)
		endif
		index=0
		bigind+=1
		SetdataFolder StringFromList(bigind,NameListInitial)
		DFREF dfr = GetDataFolderDFR() 
	while(bigind<=limits)	
	
	SetdataFolder odfr
	Redimension/n=0 w1
	w1=""
	index=0
	//make/o/T/n=0 $nameofwave(w1)
	countup=itemsinlist(NameListInitial)
	
	
	
	do
	
	
		SetdataFolder StringFromList(index,NameListInitial)
		WListNum=itemsinlist(WaveList(SearchParm, ";",""))

		if(WListNum==0)
	
		else
			NameListFinal=AddListItem(GetDataFolder(1),NameListFinal,";",inf)
	
	
	

			Insertpoints n,1, w1
			w1[n]=GetDataFolder(1)
			n+=1
	
	
		endif
		index+=1
	while(countup>index)
	SetdataFolder odfr
	return NameListFinal

End //PrintAllFolders


Static Function/S PrintAllFolders_String(SearchParm)
	string SearchParm
	String objName
	variable tote,limits,bigind,WListNum,countup,n
	SetdataFolder root:
	Variable index = 0
	DFREF dfr = GetDataFolderDFR() // Reference to current data folder
	DFREF odfr=dfr	
	String NameListInitial,NameListFinal
	WListNum=itemsinlist(WaveList(SearchParm, ";",""))

	NameListInitial=GetDataFolder(1)
	NameListFinal=""
	
	do
		Tote=CountObjectsDFR(dfr, 4)
		if (tote==0)
		
		else
			limits+=Tote
			do
				SetdataFolder GetIndexedObjNameDFR(dfr, 4, index)
				NameListInitial=AddListItem(GetDataFolder(1),NameListInitial,";",inf)
				
				SetdataFolder dfr
				index+=1
			while(index<tote)
		endif
		index=0
		bigind+=1
		SetdataFolder StringFromList(bigind,NameListInitial)
		DFREF dfr = GetDataFolderDFR() 
	while(bigind<=limits)	
	
	SetdataFolder odfr
	//Redimension/n=0 w1
	index=0
	//make/o/T/n=0 $nameofwave(w1)
	countup=itemsinlist(NameListInitial)
	
	
	
	do
	
	
		SetdataFolder StringFromList(index,NameListInitial)
		WListNum=itemsinlist(WaveList(SearchParm, ";",""))

		if(WListNum==0)
	
		else
			NameListFinal=AddListItem(GetDataFolder(1),NameListFinal,";",inf)
	
	
	

			//Insertpoints n,1, w1
			n+=1
	
	
		endif
		index+=1
	while(countup>index)
	SetdataFolder odfr
	return NameListFinal

End //PrintAllFolders



Static function ReplaceNote(w1,s1,k1)
	wave w1
	string s1, k1
	string Total=note(w1)
	Total=replacestringbykey(k1,Total,s1,":","\r")
	note/k w1 Total

end