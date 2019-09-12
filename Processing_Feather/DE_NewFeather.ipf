#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
//#include "C:\Users\dedwards\src_prh\IgorUtil\PythonApplications\FEATHER\Example\MainFeather"
#include "D:\Devin\Documents\Software\AppFEATHER\AppIgor\Example\MainFeather"
#pragma modulename=DE_NewFeather
Static Function OutportForce(ForceWave,SepWave)
	wave ForceWave,SepWave

	duplicate/o ForceWave $(replaceString("Force",nameofwave(ForceWave),"Time"))
	wave TimeWave=$(replaceString("Force",nameofwave(ForceWave),"Time"))
	TimeWave=pnt2x(ForceWave,p)
	display/N=TMP_D Forcewave vs SepWave 
	Appendtograph TimeWave
	String Path="D:\Data\Feather\Hold.pxp"
	SaveGraphCopy/o as Path
	KillWindow TMP_D
	killwaves Timewave

end

Static Function RunFeatheronOutput(OptionsWave)
	wave OptionsWave
	String Location = "D:\Data\Feather\Hold.pxp"
	
	///ModMainFEATHER#Main(base="C:/Users/dedwards/src_prh/",Input_file=Location,OptionsWave=OptionsWave)

ModMainFEATHER#Main(base="D:/Devin/Documents/Software/AppFEATHER/",Input_file=Location,OptionsWave=OptionsWave)
end


Static Function MakeSingleWaves(SearchString)
	STRING SearchString
	string AllForceRet= wavelist(SearchString,";","")
	String ForceWaveList="",SepWaveList=""
	variable n
	string ForceWave, SepWave
	for(n=0;n<itemsinlist(AllForceRet);n+=1)
	print "N:"+num2str(n)
		//for(n=0;n<itemsinlist(AllFOrceRet);n+=1)
		wave ForceRetWave=$stringfromlist(n,AllForceRet)
		wave ForceExtWave=$replacestring("Ret",nameofwave(ForceRetWave),"Ext")
		wave SepRetWave=$replacestring("Force",nameofwave(ForceRetWave),"Sep")
		wave SepExtWave=$replacestring("Force",nameofwave(ForceExtWave),"Sep")
		
		make/free/n=0 ForceAll,SepAll
		Concatenate/o {ForceExtWave,ForceRetWave},ForceAll
		Concatenate/o {SepExtWave,SepRetWave},SepAll
		duplicate/o ForceAll $replacestring("Force_Ret",nameofwave(ForceRetWave),"Force")
		duplicate/o SepAll $replacestring("Force_ret",nameofwave(ForceRetWave),"Sep")

	endfor
	killwaves ForceAll,SepAll
end


Static Function SaveOutAllWaves(SearchString)
	String SearchString
	string AllForce= wavelist(SearchString,";","")
	
	variable n,tot=itemsinlist(AllForce)
	display/N=TMP_D
	for(n=0;n<tot;n+=1)
		Wave ForceWave=$stringfromlist(n,AllForce)
		Wave SepWave=$ReplaceString("Force",nameofwave(ForceWave),"Sep")
		duplicate/o ForceWave $(replaceString("Force",nameofwave(ForceWave),"Time"))
		wave TimeWave=$(replaceString("Force",nameofwave(ForceWave),"Time"))
		Appendtograph/W=TMP_D TimeWave
		Appendtograph/W=TMP_D ForceWave vs SepWave
		TimeWave=pnt2x(ForceWave,p)

	endfor
	
	
	

	String Path="D:\Data\Feather\Hold.pxp"
	SaveGraphCopy/o as Path
	KillWindow TMP_D
	killwaves Timewave

end

Static Function LoadTheWaves(SearchString)
	String SearchString
	string AllForce= wavelist(SearchString,";","")

	NewPath/o/Q DataPath "D:Data:Feather:"
	string AllFiles= IndexedFile(DataPath, -1, ".txt")
	string CurrentFile,CurrentName
	variable n,tot=itemsInList(AllFiles)

	for(n=0;n<tot;n+=1)
		CurrentFile=stringfromlist(n,AllFiles)
		CurrentName=CurrentFile[0,strlen(CurrentFile)-5]
		if(FindListItem(Currentname,AllForce)!=-1)
			CurrentName=replaceString("Force",CurrentName,"Starts")
			LoadWave/Q/G/P=DataPath/N=GARBAGE CurrentFile
			wave garbage0
			duplicate/o GARBAGE0 $CurrentName
		else
		endif
	endfor

end