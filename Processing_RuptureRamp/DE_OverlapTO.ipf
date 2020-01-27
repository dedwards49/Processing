#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma modulename=DE_OverlapTO

//Static Function MakeIntoForceSep(Deflection,ZSnsr,ForceOut,SepOut)
//	wave Deflection,ZSnsr,ForceOut,SepOut
//
//	String noteStr=note(Deflection)
//	variable SpringConstant= str2num(stringbykey("Spring Constant",noteStr,":","\r"))
//	duplicate/free Deflection FreeForce,FreeSep
//	FastOp FreeSep=ZSnsr-Deflection
//	FastOp FreeForce=(SpringConstant)*Deflection
//	duplicate/o FreeSep,SepOut
//	duplicate/o FreeForce ForceOut
//end

Function/C ReturnShifts(CtrlDef,CtrlZSnsr,CompDef,CompZSnsr)
	wave CtrlDef,CtrlZSnsr,CompDef,CompZSnsr
	variable/C FitParms=LineFit(CtrlDef,CtrlZSnsr)
	make/free/n=2 LineParms
	LineParms[0]=real(FitParms)
	LineParms[1]=imag(FitParms)
	variable/C Shifts=ShiftLineFit(CompDef,CompZSnsr,Lineparms)
	String WindowName="ShiftGraph"
//
	dowindow $WindowName
	if(v_flag==1)
		killwindow $WindowName
	endif
	
	variable ctrlx,ctrly,compx,compy
	ctrlx=CtrlZSnsr[numpnts(CtrlZSnsr)-1]
	compx=CompZSnsr[numpnts(CompZSnsr)-1]
	ctrly=imag(FitParms)*ctrlx+real(FitParms)

	compy=imag(FitParms)*(compx+real(shifts))+real(FitParms)

	variable shifty=ctrly-compy
	
	variable adaptx=(shifty)/imag(fitparms)

	display/N=$WindowName CtrlDef vs CtrlZSnsr
	appendtograph/W=$WindowName CompDef vs CompZSnsr
	ModifyGraph/W=$WindowName rgb($nameofwave(CTrlDef))=(0,52224,26368)
	ModifyGraph offset($nameofwave(CompDef))={real(shifts)+adaptx,shifty-imag(shifts)}
	return cmplx(real(shifts)+adaptx,shifty-imag(shifts))
end
Function DE_ShiftLine(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = k*(x+x0)+b+y0
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 4
	//CurveFitDialog/ w[0] = k
	//CurveFitDialog/ w[1] = b
	//CurveFitDialog/ w[2] = x0
	//CurveFitDialog/ w[3] = y0

	return w[0]*(x+w[2])+w[1]+w[3]
End

 Function/C LineFit(Yin,XIn)
	wave Yin,XIn
		String WindowName="LineFitWin"

	dowindow $WindowName
	if(v_flag==1)
	
	killwindow $WindowName
	endif
	
	display/N=$WindowName Yin vs XIn
	ModifyGraph/W=$WindowName rgb($nameofwave(YIn))=(0,0,0)

	ShowInfo

	Variable rval= UserCursorAdjust(WindowName)
	if (rval == -1)							// Graph name error?
		return -1;
	endif

	if (rval == 1)								// User canceled?
		DoAlert 0,"Canceled"
		return -1;
	endif

	CurveFit/Q/W=2 line YIn[pcsr(A),pcsr(B)] /X=XIn /D 
	wave w_coef
	variable/c Result= cmplx(w_coef[0],w_coef[1])

	killwaves w_coef
		killwindow $windowname

	return Result
end


Function/C ShiftLineFit(YIn,XIn,Lineparms)
	wave Yin,XIn,Lineparms
	String WindowName="ShiftFitLine"
	dowindow $WindowName
	if(v_flag==1)
	
	killwindow $WindowName
	endif
	
	
	
	display/N=$WindowName YIn vs Xin
	ModifyGraph/W=$WindowName rgb($nameofwave(YIn))=(0,0,0)
	ShowInfo

	Variable rval= UserCursorAdjust(WindowName)
	if (rval == -1)							// Graph name error?
		return -1;
	endif

	if (rval == 1)								// User canceled?
		DoAlert 0,"Canceled"
		return -1;
	endif
	Make/D/N=4/O W_coef
	W_coef[0] = {Lineparms[1],Lineparms[0],0,0}
	FuncFit/Q/W=2/H="1100" DE_ShiftLine W_coef YIn[pcsr(B),pcsr(A)] /X=XIn /D
	variable/c Result= cmplx(w_coef[2],w_coef[3])
	killwaves w_coef
	killwindow $windowname
	return Result
end


Function UserCursorAdjust(graphName)
	String graphName

	DoWindow/F $graphName			// Bring graph to front
	if (V_Flag == 0)					// Verify that graph exists
		Abort "UserCursorAdjust: No such graph."
		return -1
	endif

	NewDataFolder/O root:tmp_PauseforCursorDF
	Variable/G root:tmp_PauseforCursorDF:canceled= 0

	NewPanel/K=2 /W=(139,341,382,450) as "Pause for Cursor"
	DoWindow/C tmp_PauseforCursor					// Set to an unlikely name
	AutoPositionWindow/E/M=1/R=$graphName			// Put panel near the graph

	DrawText 21,20,"Adjust the cursors and then"
	DrawText 21,40,"press Continue."
	Button button0,pos={80,58},size={92,20},title="Continue"
	Button button0,proc=UserCursorAdjust_ContButtonProc
	Button button1,pos={80,80},size={92,20}
	Button button1,proc=UserCursorAdjust_CancelBProc,title="Cancel"

	PauseForUser tmp_PauseforCursor,$graphName

	NVAR gCaneled= root:tmp_PauseforCursorDF:canceled
	Variable canceled= gCaneled			// Copy from global to local before global is killed
	KillDataFolder root:tmp_PauseforCursorDF

	return canceled
End

Function UserCursorAdjust_ContButtonProc(ctrlName) : ButtonControl
	String ctrlName
	DoWindow/K tmp_PauseforCursor			// Kill self
End

Function Demo()
	Make/O jack;SetScale x,-5,5,jack
	jack= exp(-x^2)+gnoise(0.1)
	DoWindow Graph0

	if (V_Flag==0)
		Display jack
	endif


End
