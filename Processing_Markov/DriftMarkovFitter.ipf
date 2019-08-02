#pragma rtGlobals=3		// Use modern global access method and strict wave access.

Function DriftMarkovFitter( UseWave, stateCount, modeCount, timeStep, driftBound, sigmaBound, transitionBound, iterationCount, [RAM, Threads])//Variables demanded by MarkovFit Code
Wave UseWave//Input wave
Variable stateCount, modeCount, timeStep, driftBound, sigmaBound, iterationCount, RAM, Threads
Variable TransitionBound
killwaves /z HidMar0, HidMar1, HidMar2, HidMar3, HidMar4, usable//Getting rid of generated waves to generate new ones
RAM = paramIsDefault(RAM) ? 4:RAM
Threads = paramIsDefault(Threads) ? 1000:Threads
killwaves /z Used
duplicate /o UseWave Used
if(timeStep==0)
timestep = 1.0
endif
Variable hold
if(iterationCount==0)
iterationCount = 4
endif
if(RAM == 0)
RAM = 4
endif
if(modeCount ==0)
Variable i
for(i=0;i<numpnts(Used); i+=1)
Used[i] += -driftBound*i
endfor
endif
String InfoPass = "java -Xmx" + num2str(RAM) +"g -jar C:\MarkovFitter\DriftMarkov2.jar C:\MarkovFitter\UseWave.txt " + num2str(stateCount)+" 0 "//infopass exists to hold the command sent to DOS
InfoPass = InfoPass + num2str(modeCount)+" "+num2str(timeStep)+" "+num2str(driftBound)+" "+num2str(sigmaBound)+" "+num2str(transitionBound)+" "+num2str(iterationCount)+" "+num2str(Threads)
Save/J/W Used as "C:\MarkovFitter\UseWave.txt"//saving the wave that was given to  proper location
print(InfoPass)//gives view of command line in case anything is wrong
executescripttext InfoPass//sendng command to command line
LoadWave/A=HidMar/J/D/W/K=0 "C:MarkovFitter:DriftMarkovOut.txt"//getting waves from location jar tosses them to(waves have base name HidMar
Display UseWave//displaying wave given
variable Temp
duplicate/o $"HidMar1" usable//while wave1 is created through this code it cannot regonize it so it must be duplicated
Temp =dimoffset(UseWave,0)
print(temp)
setscale/P x dimoffset(UseWave,0), dimdelta(UseWave,0), "s", usable//ensuring scaling of input and output wave are the same
if(modeCount ==0)
for(i=0;i<numpnts(UseWave);i+=1)
usable[i] += driftBound*i
endfor
endif
AppendToGraph usable//putting on same graph
ModifyGraph rgb(usable)=(0,0,65280)//changing color so both waves are visible
display $"HidMar2"//displaying simple jump wave
executescripttext "java -jar C:\MarkovFitter\GetRidOfUseWave.jar"//Eliminates file created earlier to prevent problems on future runs
end