﻿#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.


#pragma ModuleName=DE_heatMap

Static Function/C ReturnXYOffsetfromGraph(WaveIn,GraphNameString)
	wave WaveIn
	string GraphNameString
	String TraceNames=tracenamelist(GraphNameString,";",1)
	if(Strsearch(tracenames,nameofwave(wavein),0)!=0)
		
		variable itemnumber=WhichListItem(nameofwave(wavein),tracenames)
		string Thistracename=stringfromlist(itemnumber,traceNames)
		String offsets=stringbykey("offset(x)",traceinfo(GraphNameString,Thistracename,0),"=",";")
		Variable xOffset,yOffset
    	sscanf offsets, "{%g,%g}",xOffset,yOffset
    	variable/C Returns=cmplx(xoffset,yoffset)
    	return Returns
	else
	return -1
	endif
end

Static Function GenerateSingleWaveswithShifts(GraphNameString,YIn,Xin,Yout,Xout)
	wave YIn,Xin,Yout,Xout
String GraphNameString
	duplicate/free Yin YFree
	duplicate/free Xin Xfree
	variable/c Offsets=ReturnXYOffsetfromGraph(Yin,GraphNameString)
	YFree-=imag(Offsets)
	XFree+=real(Offsets)
	YFree*=-1
	duplicate/o YFree Yout
	duplicate/o Xfree Xout
end

Static Function MakeSingleHeatMap(GraphNameString,YIn,Xin,Xdims,Ydims,MapOut)
	string GraphNameString
	wave YIn,Xin,Xdims,Ydims,MapOut
	make/o/n=0 Yfree,XFree
	GenerateSingleWaveswithShifts(GraphNameString,YIn,Xin,Yfree,XFree)	
	JointHistogram/XBWV=Xdims/YBWV=Ydims XFree, Yfree
	wave M_JointHistogram
	duplicate/o M_JointHistogram MapOut
	variable tot=sum(MapOut)
	MapOut/=tot
	killwaves M_JointHistogram
end

Static Function AccumulateHeatMaps(graphnamestring)
	String GraphNameString
	String ForceWaveNames=tracenamelist(GraphNameString,";",1)
	variable tot=itemsinlist(ForceWaveNames)
	//ForceWaveNames=RemoveListItem(tot-1, ForceWaveNames)
	tot=itemsinlist(ForceWaveNames)
	variable n=0
	variable nx=200
	variable ny=100
	make/o/n=(nx) XBins;XBins=0e-9+200e-9*x/(nx-1)
	make/o/n=(ny) YBins;YBins=10e-12+500e-12*x/(ny-1)
	make/o/n=(nx-1,ny-1) Final=0
	for(n=0;n<tot;n+=1)
	print n
		wave FW=$stringfromlist(n,ForceWaveNames)
		wave SW=$replacestring("fORCE",nameofwave(FW),"sEP")
		make/free/n=(0,0) MapOut
		MakeSingleHeatMap(GraphNameString,FW,SW,XBins,YBins,MapOut)
		Final+=MapOut
	endfor







end