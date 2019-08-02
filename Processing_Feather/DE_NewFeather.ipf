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