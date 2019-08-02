#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma ModuleName = DE_BatchProcessing
#include ":\Misc_PanelPrograms\AsylumNaming"
Function MaxForceofTrace()

	NewPath/o/q Pa
	PathINfo Pa
	string fname
	variable maxforce, maxSep
	variable n
	make/o/n=0 MForce, MSep
	do

		fName= IndexedFile(Pa,n,".ibw")
		if (strlen(fName) == 0)
			Break
		else
		endif
		fName=DE_SingleLoadIBW(S_path,n)
		wave DefRet=$(fName+"Defl_Ret")
		wave ZRet=$(fName+"Zsnsr_Ret")
		wave DefExt=$(fName+"Defl_Ext")
		wave ZExt=$(fName+"Zsnsr_Ext")
		GenerateForceandSep(DefRet,ZRet)
		struct ForceWave Name1
		DE_Naming#WavetoStruc(nameofwave(DefRet),Name1)
		string ForceName=DE_Naming#StringCreate(Name1.Name,Name1.VNum,"Force",Name1.SDirec)
		string SepName=DE_Naming#StringCreate(Name1.Name,Name1.VNum,"Sep",Name1.SDirec)
		wave FWave=$ForceName
		wave SWave=$SepName
		
		wavestats/q Fwave
		 maxforce=V_min-Fwave[numpnts(Fwave)-1]
		maxsep=SWave(v_minloc)-SWave[0]
		insertpoints n,1, MForce, MSep
		MForce[n]=maxforce
		MSep[n]=maxsep
		killwaves DefRet,ZRet,DefExt,ZExt,FWave,SWave
	n+=1
	while(1<3)

end

Function GenerateForceandSep(DefWave,ZWave)
	wave DefWave,ZWave

	variable Invols,Stiff

	if(cmpstr(stringbykey("Invols",note(DefWave),":","\r"),"NaN")==0)
		print "Bad Invols"
		return -1
	elseif(cmpstr(stringbykey("SpringConstant",note(DefWave),":","\r"),"NaN")==0)
		print "Bad Spring Constant"
		return -1
	endif


	Invols=str2num(stringbykey("Invols",note(DefWave),":","\r"))
	Stiff=str2num(stringbykey("SpringConstant",note(DefWave),":","\r"))

	struct ForceWave Name1
	DE_Naming#WavetoStruc(nameofwave(DefWave),Name1)
	string ForceName=DE_Naming#StringCreate(Name1.Name,Name1.VNum,"Force",Name1.SDirec)
	string SepName=DE_Naming#StringCreate(Name1.Name,Name1.VNum,"Sep",Name1.SDirec)


	
	duplicate/o ZWave $SepName
	duplicate/o DefWave $ForceName
	wave FWave=$ForceName
	wave SWave=$SepName
	FWave*=Stiff
	SWave=-ZWave+DefWave
	
	
end