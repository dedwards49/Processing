#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma modulename=DE_TwoWLCFit
#include "SimpleWLCPrograms"
#include ":\Misc_PanelPrograms\Panel Progs"
#include ":\Misc_PanelPrograms\AsylumNaming"

Static Function CombineCurves(Force1,Sep1,Force2,Sep2,Fout,SOut)
	wave Force1,Sep1,Force2,Sep2,Fout,SOut
	make/free/n=(numpnts(Sep1),2) tpSep1
	make/free/n=(numpnts(Sep2),2) tpSep2
	
	tpSep1[][0]=Sep1[p]
tpSep1[][1]=0
tpSep2[][0]=Sep2[p]
tpSep2[][1]=1
	concatenate/o/np=0 {tpSep1,tpSep2}, STemp
concatenate/o/np=0 {Force1,Force2}, FTemp
	duplicate/o STemp, SOut 
	duplicate/o FTemp, FOut
	killwaves STemp,FTemp 

end

Static Function MakeAMultiFit(SepIn,WLCParms,ForceOut)
	wave SepIn,WLCParms,ForceOut
	make/free/n=(dimsize(Sepin,0)) TForce
	TForce=DE_FitTwo(WLCParms,SepIn[p][0],SepIn[p][1])
	duplicate/o TForce ForceOut
end

Static Function FitForcePair(FoldedForce,FoldedSep,UnFoldedForce,UnFoldedSep,ParmOut,FitOut,[Constrained])
	wave FoldedForce,FoldedSep,UnFoldedForce,UnFoldedSep,ParmOut,FitOut,Constrained
	  variable 	timerRefNum = startMSTimer
	if(ParamisDefault(Constrained))
	make/free/n=0 Constrained
	Constrained={-4e-10,0,0,298,0,0}
	endif

	//note that WLCGuess should have the format: Lp,Lc1,Lc2,T,Xoff,Foff
	make/free/n=0 Fout,Sout
	CombineCurves(FoldedForce,FoldedSep,UnFoldedForce,UnFoldedSep,Fout,SOut)
	make/free/n=5 w_coef
	variable  fitstart=wavemin(Foldedsep)-50e-9
	variable n
	w_coef={-.4e-9,150e-9,175e-9,298,fitstart,0}
	string ConStr="100100"
	for(n=0;n<numpnts(W_coef);n+=1)
		if(Constrained[n]==0)
		
		else
			w_coef[n]=Constrained[n]
			ConStr[n,n]="1"
		endif
	endfor
	
	FuncFit/N/Q/W=2/H=ConStr/NTHR=0 DE_FitTwo W_coef  Fout /X=SOut 
	make/free/n=0 FWLCFit
	MakeAMultiFit(SOut,w_coef,FWLCFit)
	make/free/n=(Numpnts(FWLCFit),2) FitTest
	FitTest[][1]=FWLCFit[p]
	FitTest[][0]=SOut[p][0]
	duplicate/o w_coef ParmOut
	duplicate/o FitTest FitOut
	variable	microSeconds = stopMSTimer(timerRefNum)
	//Print microSeconds/1e6/60, "minutes per iteration"
	//print w_coef[2]-w_coef[1]
end

Static Function/C CalculateFandSOffsetfromWLC(WLCBase,WLCComp)

	wave WLCBase,WLCComp
	if(WLCComp[0]!=WLCBase[0]||WLCComp[1]!=WLCBase[1]||WLCComp[2]!=WLCBase[3]||WLCComp[2]!=WLCBase[3])
	endif
	
	return cmplx(WLCComp[4]-WLCBase[4],WLCComp[5]-WLCBase[5])


end

Static Function GenerateFitFromONeWave(ForceWave,SepWave,pnt1,pnt2,pnt3,pnt4,Waveout,FitReturn,[Constraint])
	wave ForceWave,SepWave,Waveout,FitReturn,Constraint
	variable pnt1,pnt2,pnt3,pnt4
	duplicate/free/r=[pnt1,pnt2] ForceWave, FoldedForce
	duplicate/free/r=[pnt3,pnt4] ForceWave, UnFoldedForce
	duplicate/free/r=[pnt1,pnt2] SepWave, FoldedSep
	duplicate/free/r=[pnt3,pnt4] SepWave, UnFoldedSep
	make/free/n=0 ParmOut,FitOut
	if(ParamisDefault(Constraint))
	FitForcePair(FoldedForce,FoldedSep,UnFoldedForce,UnFoldedSep,ParmOut,FitOut)
	else
		FitForcePair(FoldedForce,FoldedSep,UnFoldedForce,UnFoldedSep,ParmOut,FitOut,Constrained=Constraint)

	endif
	duplicate/o FitOut,FitReturn
	duplicate/o ParmOut,Waveout

end


Static Function/S TwoFitFromWave(ForceWave,SepWave,[Constraint,pnt1,pnt2,pnt3,pnt4])

	wave ForceWave,SepWave,Constraint
	variable pnt1,pnt2,pnt3,pnt4
	DoWindow TwoWLC
	if(V_flag==1)
	killwindow TwoWLC
	endif
	display/N=TwoWLC ForceWave vs SepWave
	ShowInfo
		Variable autoAbortSecs=20
	if(Paramisdefault(pnt1)||Paramisdefault(pnt2)||Paramisdefault(pnt3)||Paramisdefault(pnt4))
	if (UserCursorAdjust("TwoWLC",autoAbortSecs) != 0)
		return "-1"
	endif
	if (strlen(CsrWave(A))>0 && strlen(CsrWave(B))>0)	// Cursors are on trace?
		pnt1=pcsr(A)
		pnt2=pcsr(B)
	endif
	if (UserCursorAdjust("TwoWLC",autoAbortSecs) != 0)
		return "-1"
	endif
	if (strlen(CsrWave(A))>0 && strlen(CsrWave(B))>0)	// Cursors are on trace?
		pnt3=pcsr(A)
		pnt4=pcsr(B)
	endif
	else
	endif
	make/free/n=0 WaveOut,FitOut
	if(ParamisDefault(Constraint))
	
	GenerateFitFromONeWave(ForceWave,SepWave,pnt1,pnt2,pnt3,pnt4,Waveout,FitOut)
	else
	GenerateFitFromONeWave(ForceWave,SepWave,pnt1,pnt2,pnt3,pnt4,Waveout,FitOut,Constraint=Constraint)

	endif
	
	duplicate/o Waveout $(nameofwave(ForceWave)+"_WLCParms")
	duplicate/o FitOut $(nameofwave(ForceWave)+"_WLCFit")
	wave w1=$(nameofwave(ForceWave)+"_WLCFit")
	appendtograph/w=TwoWLC w1[][1] vs w1[][0]
	ModifyGraph rgb($nameofwave(w1))=(0,0,0)
	return "pnt1="+num2str(pnt1)+","+"pnt2="+num2str(pnt2)+","+"pnt3="+num2str(pnt3)+","+"pnt4="+num2str(pnt4)

end

Static Function/S TwoFitFromTwoWaves(FW1,Sw1,FW2,SW2,[Constraint])

	wave FW1,Sw1,FW2,SW2,Constraint
//	DoWindow TwoWLC
//	if(V_flag==1)
//	killwindow TwoWLC
//	endif
//	display/N=TwoWLC ForceWave vs SepWave
//	ShowInfo
//		Variable autoAbortSecs=20
//	if(Paramisdefault(pnt1)||Paramisdefault(pnt2)||Paramisdefault(pnt3)||Paramisdefault(pnt4))
//	if (UserCursorAdjust("TwoWLC",autoAbortSecs) != 0)
//		return "-1"
//	endif
//	if (strlen(CsrWave(A))>0 && strlen(CsrWave(B))>0)	// Cursors are on trace?
//		pnt1=pcsr(A)
//		pnt2=pcsr(B)
//	endif
//	if (UserCursorAdjust("TwoWLC",autoAbortSecs) != 0)
//		return "-1"
//	endif
//	if (strlen(CsrWave(A))>0 && strlen(CsrWave(B))>0)	// Cursors are on trace?
//		pnt3=pcsr(A)
//		pnt4=pcsr(B)
//	endif
//	else
//	endif
	make/free/n=0 WaveOut,FitOut
	if(ParamisDefault(Constraint))
	FitForcePair(FW1,SW1,FW2,SW2,WaveOut,FitOut)
	else
		FitForcePair(FW1,SW1,FW2,SW2,WaveOut,FitOut,Constrained=Constraint)

	endif
	
	duplicate/o Waveout $(nameofwave(FW1)+"_WLCParms")
	duplicate/o FitOut $(nameofwave(FW1)+"_WLCFit")
	//wave w1=$(nameofwave(FW1)+"_WLCFit")
	//appendtograph/w=TwoWLC w1[][1] vs w1[][0]
//	ModifyGraph rgb($nameofwave(w1))=(0,0,0)
	//return "pnt1="+num2str(pnt1)+","+"pnt2="+num2str(pnt2)+","+"pnt3="+num2str(pnt3)+","+"pnt4="+num2str(pnt4)

end


Function UserCursorAdjust(graphName,autoAbortSecs)
	String graphName
	Variable autoAbortSecs

	DoWindow/F $graphName							// Bring graph to front
	if (V_Flag == 0)									// Verify that graph exists
		Abort "UserCursorAdjust: No such graph."
		return -1
	endif

	NewPanel /K=2 /W=(187,368,437,531) as "Pause for Cursor"
	DoWindow/C tmp_PauseforCursor					// Set to an unlikely name
	AutoPositionWindow/E/M=1/R=$graphName			// Put panel near the graph

	DrawText 21,20,"Adjust the cursors and then"
	DrawText 21,40,"Click Continue."
	Button button0,pos={80,58},size={92,20},title="Continue"
	Button button0,proc=UserCursorAdjust_ContButtonProc
	Variable didAbort= 0
	if( autoAbortSecs == 0 )
		PauseForUser tmp_PauseforCursor,$graphName
	else
		SetDrawEnv textyjust= 1
		DrawText 162,103,"sec"
		SetVariable sv0,pos={48,97},size={107,15},title="Aborting in "
		SetVariable sv0,limits={-inf,inf,0},value= _NUM:10
		Variable td= 10,newTd
		Variable t0= ticks
		Do
			newTd= autoAbortSecs - round((ticks-t0)/60)
			if( td != newTd )
				td= newTd
				SetVariable sv0,value= _NUM:newTd,win=tmp_PauseforCursor
				if( td <= 10 )
					SetVariable sv0,valueColor= (65535,0,0),win=tmp_PauseforCursor
				endif
			endif
			if( td <= 0 )
				DoWindow/K tmp_PauseforCursor
				didAbort= 1
				break
			endif
				
			PauseForUser/C tmp_PauseforCursor,$graphName
		while(V_flag)
	endif
	return didAbort
End

Function UserCursorAdjust_ContButtonProc(ctrlName) : ButtonControl
	String ctrlName

	DoWindow/K tmp_PauseforCursor				// Kill self
End


End

function DE_FitTwo(w,x1,x2) : FitFunc
	wave w
	variable x1,x2
	

	variable y
	variable lc
	if(x2==0)
		LC=w[1]
	elseif(x2==1)
		LC=w[2]
	endif
	//w[0] Lp
	//w[1] Lc1
	//w[2] Lc2
	//w[3] Temp
	//w[4] xoff
	//w[5] Foff
	y=WLC(x1-w[4],w[0],LC,w[3])+w[5]
	return y

end

Function TestMacro(startnum,endnum)
	variable startnum,endnum
	variable n=startnum
	string ForceNameb,ForceNamea,Sepnameb,SepNameA,NumName,Points
	variable/C First,Second
	make/o/n=0 macroFalign,macroSalign,macroFalignZero
	make/o/n=0/T macroName
	for(n=startnum;n<=endnum;n+=1)

		sprintf NumName "%04.4G",n
		ForceNameb="pH6p2M_b"+NumName+"Force_Ret"
		if(waveexists($ForceNameb)==0)
		print "SKIPPING: "+NumName
		else
		ForceNamea="pH6p2M_a"+NumName+"Force_Ret"
		Sepnameb="pH6p2M_b"+NumName+"Sep_Ret"
		SepNameA="pH6p2M_a"+NumName+"Sep_Ret"
		Points=DE_TwoWLCFit#TwoFitFromWave($ForceNameb,$Sepnameb,Constraint=$"wave0")
		First= DE_TwoWLCFit#CalculateFandSOffsetfromWLC($(ForceNameb+"_WLCParms"),$"pH6p2M_b0002Force_Ret_WLCParms")
		Execute "DE_TwoWLCFit#TwoFitFromWave("+ForceNameb+","+Sepnameb+",Constraint=wave1,"+Points+")"
		Second= DE_TwoWLCFit#CalculateFandSOffsetfromWLC($(ForceNameb+"_WLCParms"),$"pH6p2M_b0002Force_Ret_WLCParms")
		print ForceNameb+":"+num2str(real(First))+":"+num2str(imag(First))
		print ForceNameb+":"+num2str(real(Second))+":"+num2str(imag(Second))
		insertpoints (2*(n-startnum)),1, macroName,macroFalign,macroSalign,macroFalignZero
		macroName[2*(n-startnum)]=ForceNameb
		macroFalign[2*(n-startnum)]=imag(First)
		macroSalign[2*(n-startnum)]=real(First)
		macroFalignZero[2*(n-startnum)]=imag(Second)
		Points=DE_TwoWLCFit#TwoFitFromWave($ForceNamea,$Sepnamea,Constraint=$"wave0")
		First= DE_TwoWLCFit#CalculateFandSOffsetfromWLC($(ForceNamea+"_WLCParms"),$"pH6p2M_b0002Force_Ret_WLCParms")
		Execute "DE_TwoWLCFit#TwoFitFromWave("+ForceNamea+","+Sepnamea+",Constraint=wave1,"+Points+")"
		Second= DE_TwoWLCFit#CalculateFandSOffsetfromWLC($(ForceNamea+"_WLCParms"),$"pH6p2M_b0002Force_Ret_WLCParms")
		print ForceNamea+":"+num2str(real(First))+":"+num2str(imag(First))
		print ForceNamea+":"+num2str(real(Second))+":"+num2str(imag(Second))
		insertpoints (2*(n-startnum)+1),1, macroName,macroFalign,macroSalign,macroFalignZero
				macroName[2*(n-startnum)+1]=ForceNamea
		macroFalign[2*(n-startnum)+1]=imag(First)
		macroSalign[2*(n-startnum)+1]=real(First)
		macroFalignZero[2*(n-startnum)+1]=imag(Second)
		endif
	endfor
end