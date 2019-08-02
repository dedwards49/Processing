#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma modulename=DE_Slippage
#include "DE_Filtering"
Menu "Slippage"
	"PlotSlippage", SelectandPlot()
End
Function SelectandPlot()
	String WaveList1=WaveList("*Defl_Ret*", ";", "")
	String WaveList2=WaveList("*Raw*_Ret", ";", "")
variable Smoothing=51
	string DeflectionWave,RawWave
		Prompt Smoothing,"Smoothing"

	Prompt DeflectionWave,"Defl",popup,WaveList1
	Prompt RawWave,"Raw",popup,WaveList2
	
	DoPrompt "Pick",Smoothing,DeflectionWave,RawWave



	
	if (V_Flag)
		return -1 // User canceled
	endif
	
	wave w1=$DeflectionWave
	wave w2=$RawWave
	
	PlottheWayILike(w1,w2,Smoothing=Smoothing)
End

Static Function PlotTheWayILike(DeflWave,RawWave,[Smoothing])

	wave DeflWave,RawWave
	variable Smoothing
	DoWindow Slippage
	if(ParamisDefault(Smoothing))
	Smoothing=51
	endif

	if(V_Flag==1)
		killwindow Slippage
	endif
	duplicate/o DeflWave SmoothedDef
	duplicate/o RawWave,SmoothedRaw
	if(smoothing<1)
		DE_Filtering#TVD1D_denoise(DeflWave,smoothing,SmoothedDef)
		DE_Filtering#TVD1D_denoise(RawWave,smoothing,SmoothedRaw)

	elseif(smoothing>5)
		Smooth/S=2 Smoothing, SmoothedDef
		Smooth/S=2 Smoothing, SmoothedRaw
	else
	
	endif

	Display/N=Slippage DeflWave
		AppendToGraph/W=Slippage SmoothedDef
	AppendToGraph/W=Slippage/L=L1 RawWave
	AppendToGraph/W=Slippage/L=L1 SmoothedRaw

	ModifyGraph/W=Slippage axisEnab(left)={0,0.45},axisEnab(L1)={0.55,1},freePos(L1)={0,bottom}
	ModifyGraph/W=Slippage lblPosMode(left)=1,lblPosMode(L1)=1
	ModifyGraph/W=Slippage muloffset($nameofwave(DeflWave))={0,-1},muloffset($nameofwave(SmoothedDef))={0,-1}
	
	ModifyGraph/W=Slippage rgb($nameofwave(DeflWave))=(65280,48896,48896);DelayUpdate
ModifyGraph/W=Slippage rgb($nameofwave(RawWave))=(40448,49664,57344);DelayUpdate
ModifyGraph/W=Slippage rgb(SmoothedRaw)=(14848,32256,47104)
	
end

Static Function MakeSomeLines(NumberofLines,Spacing,Start,SepWave,Outwave)

	variable NumberofLines,Spacing,Start
	wave SepWave,Outwave
	make/free/n=(numberoflines*1000,2) FreeLines
	variable starttime=pnt2x(Sepwave,0)
	variable delta=dimdelta(Sepwave,0)

	variable n
	for(n=0;n<numberoflines;n+=1)
	FreeLines[1000*n,1000*n+998][0]=Start+n*Spacing
		FreeLines[1000*n+999][0]=NaN

		FreeLines[1000*n,1000*n+999][1]=starttime+(p-999*n)*delta*numpnts(SepWave)/999

	//FreeLinesX
	
	
	endfor
	duplicate/o FreeLines Outwave
end

Static Function CuttoRegion(WaveIn)
	wave WaveIn
	make/free/n=(dimsize(WaveIn,0)) ColA,ColB
	ColA=Wavein[p][0]
	ColB=Wavein[p][1]
	FindLevels/Q ColB hcsr(A)
	variable Number=(v_LevelsFound-1)/2
	wave W_FindLevels
	duplicate/free W_FindLevels Firstlevels

	FindLevels/Q ColB hcsr(B)
	duplicate/free W_FindLevels Secondlevels
	variable n,firstloc,secondloc
	firstloc=floor(Firstlevels[0])-2

	secondloc=(ceil(Secondlevels[0]))+2
	for(n=number;n>=0;n-=1)
	deletepoints (secondloc+n*1000),(999-secondloc), ColA,ColB
	deletepoints n*1000,(firstloc-1), ColA,ColB

	
	endfor

	make/free/n=(numpnts(ColA),2) FreeFinal
	
	FreeFinal[][0]=ColA[p]
	FreeFinal[][1]=ColB[p]
	duplicate/o FreeFinal WaveIn
	killwaves w_findlevels
end