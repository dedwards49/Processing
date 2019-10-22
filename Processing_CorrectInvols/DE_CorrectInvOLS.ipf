#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma ModuleName=DE_CorrectFEC

Static Function RecalculateForceSep(ForceWave,SepWave,NewSpring)
	Wave ForceWave,SepWave
	variable NewSpring
	
	String WaveNote=Note(ForceWave)
	variable currentInvols=ReturnCurrentInvols(ForceWave)
	Variable currentSpring=ReturnCurrentSpringConstant(ForceWave)
	variable newInvols

	if(NewSpring==0)
		NewSpring=str2num(stringbykey("SpringConstant",WaveNote,":","\r"))
		newInvols=sqrt(currentInvols^2*currentspring/newSpring)
	else
		newInvols=sqrt(currentInvols^2*currentspring/newSpring)
	endif
	
	WaveNote=replacestringbykey("DE_NewSpring",WaveNote,num2str(NewSpring),":","\r")
	WaveNote=replacestringbykey("DE_NewInvols",WaveNote,num2str(NewInvols),":","\r")

	duplicate/free Forcewave Defl,DeflV,ZSensor
	Defl=ForceWave/currentSpring
	DeflV=Defl/currentInvols
	ZSensor=Defl-SepWave
	Defl=NewInvols*DeflV
	Forcewave=Defl*NewSpring
	SepWave=Defl-ZSensor
	note/K ForceWave,WaveNote
	note/K SepWave,WaveNote

end

Static Function ReturnCurrentSpringConstant(ForceWave)
	wave ForceWave
	string Wavenote=note(ForceWave)
	variable result	
	If(cmpstr("",Stringbykey("DE_NewSpring",Wavenote,":","\r"))==0)
		result=str2num(stringbykey("SpringConstant",WaveNote,":","\r"))
	else
		result=str2num(stringbykey("DE_NewSpring",WaveNote,":","\r"))

	endif
	
	return result
end

Static Function ReturnCurrentInvols(ForceWave)
	wave ForceWave
	string Wavenote=note(ForceWave)
	variable result	
	If(cmpstr("",Stringbykey("DE_NewInvols",Wavenote,":","\r"))==0)
		result=str2num(stringbykey("Invols",WaveNote,":","\r"))
	else
		result=str2num(stringbykey("DE_NewInvols",WaveNote,":","\r"))

	endif
	
	return result
end