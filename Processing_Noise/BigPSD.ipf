


function BigPSDFunc(w,seglen,window,Sqroot)		//now a function with a macro front end
	string w
	variable seglen, window, Sqroot

variable start = StopMSTimer(-2)
	
	variable npsd = 2^(7+seglen)				// number of points in group (resultant psd wave len= npsd/2+1)
	variable psdOffset = npsd/2					// offset each group by this amount
	variable psdFirst = 0							// start of current group
	variable nsrc = numpnts($w)
	variable nsegs,winNorm						// count of number of segements and window normalization factor
	string destw = w+"_psd",srctmp = w+"_tmp"
	string winw = w+"_psdWin"					// window goes here
	
	if( npsd > nsrc/2 )
		Abort "psd: source wave should be MUCH longer than the segment length"
	endif
	make/O/N=(npsd/2+1) $destw
	make/O/N=(npsd) $srctmp,$winw
	wave WinWave = $winw
	wave SourceTemp = $srctmp
	wave DestWave = $destw
	WinWave = 1

	switch (window)
		case 1:
			winNorm = 1
			break
		case 2:
			Hanning WinWave
			winNorm = 0.372				//  winNorm is avg squared value
			break
		case 3:
			winNorm = Parzen(WinWave)
			break
		case 4:
			winNorm = Welch(WinWave)
			break
		case 5:
			winNorm = Hamming(WinWave)
			break
		case 6:
			winNorm = BlackmanHarris3(WinWave)
			break
		case 7:
			winNorm = KaiserBessel(WinWave)
			break
		default:
			Abort "unknown window index"
	endswitch
		
	Duplicate/O/R=[0,npsd-1] $w SourceTemp
	SourceTemp *= WinWave
	fft SourceTemp
	CopyScales/P SourceTemp, DestWave
	DestWave = magsqr(SourceTemp)
	psdFirst = psdOffset
	nsegs = 1
	do
		Duplicate/O/R=[psdFirst,psdFirst+npsd-1] $w SourceTemp
		SourceTemp *= WinWave
		fft SourceTemp
		DestWave += magsqr(SourceTemp)
		psdFirst += psdOffset
		nsegs += 1
	while (psdFirst+npsd < nsrc)
	print nsegs
	
	winNorm = 2*deltax($w)/(winNorm*nsegs*npsd)
	DestWave *= winNorm
	DestWave[0] /= 2

	if (Sqroot == 1)
		DestWave = Sqrt(DestWave)
	endif
	KillWaves SourceTemp, WinWave

print (StopMSTimer(-2)-start)/1e6

end //BigPSDFunc
function GenerateNoiseFigures(w,AllanMax,PSDPnts)
wave w
variable AllanMax,PSDPnts
variable scaling, mult
variable stiffness,invols

	
Prompt Scaling,"What units do you want?",popup,"volts;nm,N"
	DoPrompt "Scaling",Scaling
	
	if (V_Flag)
		return 0									// user canceled
	endif

	if (Scaling == 1)
		mult=1
	
	elseif(Scaling == 2)
	Prompt Invols,"Invols?"
	DoPrompt "Invols",Invols
	mult=Invols
		elseif(Scaling == 1)
	Prompt Invols,"Invols?"
	Prompt stiffness,"stiffness?"
	DoPrompt "Parameters",Invols, stiffness
		mult=Invols*stiffness

	endif
string  AllanName=GenerateAllan(w,AllanMax)

string PSDName=GeneratePSD(w,PSDPnts)



string IFNName=GenerateIFN($PSDName,1e4)
print nameofwave($AllanName)
print nameofwave($PSDName)
print nameofwave($IFNName)
//wave wA=$AllanName
//FastOp wa=(mult)*wa
//wave wP=$PSDName
//FastOp wp=(mult^2)*wp
//wave wi=$IFNName
//FastOp wi=(mult)*wi


end


function/S GenerateAllan(w,nMax)


	wave w
	Variable nMax 
	string wname=nameofwave(w)
	Variable wd
	wd = WaveDims(w)
	If(wd !=1)
		Abort "Allan variance for 1-D waves only."
	endif
	String aaxis,avar,alertStr,adev
	aaxis = wName + "_avx"  // x-axis for output wave
	avar = wName + "_avar" // output wave
	adev = wName + "_adev" // output wave

	If((Exists(aaxis) == 1) %| (Exists(avar)==1))
	alertStr = "Wave(s) for Allan variance of \'"+wName
	alertStr += "\' already exist, overwrite them?"
		DoAlert 1, alertStr
		If(V_flag == 2)
			Abort "User aborted AllanVariance."
		endif
	endif	
	Variable npt,dX,n2,n3,n4
	npt = numpnts(w)
	dX = deltax(w)
	n2 = 2*floor(sqrt(npt))  // later limit n2 to nMax points
	Make/O/N=(n2) $aaxis
	WAVE aaxisW = $aaxis
	aaxisW = (p+1)*dX
	aaxisW[(n2/2+1),(n2-1)] = dX*round(npt/(n2-p+1))
	If(n2 > nMax)  // replace middle point by  exp spaced pts
		n3 = floor(nMax/3)
		n4 = n2-2*n3
		DeletePoints n3,n4,aaxisW
		InsertPoints n3,n3,aaxisW
		Variable j=0
		Variable h = (aaxisW[(2*n3)]/aaxisW[(n3-1)])^(1/(n3+1))
		Variable hh = h
		Do
			aaxisW[j+n3] = aaxisW[(n3-1)]*hh
			hh *=h
			j += 1
		While(j <= n3)
	endif
	Make/O/N=(numpnts(aaxisW)) $avar
	WAVE avarW = $avar
	avarW = FAllanVar(w, aaxisW[p])
	WAVE adevW = $adev
	duplicate/o avarW $adev
	adevW=sqrt(avarW)
	//Display/K=1  avarW vs aaxisW
	//AllanVarianceStyle()
	return adev
end





function/S GeneratePSD(w,npnts)
	wave w
	variable npnts
	string outwave=nameofwave(w)+"_psd"
	BigPSDFunc(nameofwave(w),npnts,2,0)
	return outwave
	
end

function/S GenerateIFN(w,tstart)
wave w
variable tstart
variable pstart=x2pnt(w,  tstart)
string outwave=nameofwave(w)+"_psd"
if(pstart==0)
print "Minimum Time is below resolution"
return ""
endif

duplicate/o w toInt
toInt[0,(pstart-1)]=0
wave w1=$outwave
Integrate/METH=1 toInt/D=$outwave
w1=sqrt(w1)
killwaves toInt
return outwave


end



function BigPSDFunc_Single(w,window,Sqroot)		//now a function with a macro front end
	string w
	variable window, Sqroot

variable start = StopMSTimer(-2)
	
//	variable npsd = 2^(7+seglen)				// number of points in group (resultant psd wave len= npsd/2+1)
//	variable psdOffset = npsd/2					// offset each group by this amount
//	variable psdFirst = 0	
						// start of current group
variable npsd= numpnts($w)
	variable psdOffset = npsd/2					// offset each group by this amount
	variable psdFirst = 0	
	variable nsrc = numpnts($w)
	variable nsegs,winNorm						// count of number of segements and window normalization factor
	string destw = w+"_psd",srctmp = w+"_tmp"
	string winw = w+"_psdWin"					// window goes here
//	
//	if( npsd > nsrc/2 )
//		Abort "psd: source wave should be MUCH longer than the segment length"
//	endif
	make/O/N=(npsd/2+1) $destw
	make/O/N=(npsd) $srctmp,$winw
	wave WinWave = $winw
	wave SourceTemp = $srctmp
	wave DestWave = $destw
	WinWave = 1

	switch (window)
		case 1:
			winNorm = 1
			break
		case 2:
			Hanning WinWave
			winNorm = 0.372				//  winNorm is avg squared value
			break
		case 3:
			winNorm = Parzen(WinWave)
			break
		case 4:
			winNorm = Welch(WinWave)
			break
		case 5:
			winNorm = Hamming(WinWave)
			break
		case 6:
			winNorm = BlackmanHarris3(WinWave)
			break
		case 7:
			winNorm = KaiserBessel(WinWave)
			break
		default:
			Abort "unknown window index"
	endswitch
		
	Duplicate/O/R=[0,npsd-1] $w SourceTemp
	SourceTemp *= WinWave
	fft SourceTemp
	CopyScales/P SourceTemp, DestWave
	DestWave = magsqr(SourceTemp)
	psdFirst = psdOffset
	nsegs = 1
//	do
	//	Duplicate/O/R=[psdFirst,psdFirst+npsd-1] $w SourceTemp
			Duplicate/O $w SourceTemp

		SourceTemp *= WinWave
		fft SourceTemp
		DestWave += magsqr(SourceTemp)
		psdFirst += psdOffset
		//nsegs += 1
//	while (psdFirst+npsd < nsrc)
//	print nsegs
	
	//winNorm =2* deltax($w)/(winNorm*nsegs*npsd)
		winNorm = deltax($w)/(winNorm*nsegs*npsd)

	DestWave *= winNorm
	//DestWave[0] /= 2

	if (Sqroot == 1)
		DestWave = Sqrt(DestWave)
	endif
	KillWaves SourceTemp, WinWave

print (StopMSTimer(-2)-start)/1e6

end //BigPSDFunc