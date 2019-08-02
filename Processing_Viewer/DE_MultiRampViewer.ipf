#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma modulename=DE_MultiRampViewer
#include "DE_Filtering"
#include "SimpleWLCPrograms"
#include ":\Misc_PanelPrograms\Panel Progs"
#include ":\Misc_PanelPrograms\AsylumNaming"

Static Function SelectRamp()

	string saveDF = GetDataFolder(1)
	controlinfo de_Viewer_popup0
	SetDataFolder s_value
	
	if(waveExists(root:DE_Viewer:FullForceTrace)==0)
	
		return 0
	else
	
	endif
	
	wave ForceWave=root:DE_Viewer:FullForceTrace
	wave ForceWaveSm=root:DE_Viewer:FullForceSm

	wave SepWave=root:DE_Viewer:FullSepTrace
	wave SepWaveSm=root:DE_Viewer:FullSepSm
	wave ForceWaveLessSm=root:DE_Viewer:FullForceLessSm
	wave SepWaveLessSm=root:DE_Viewer:FullSepLessSm

	controlinfo de_Viewer_list0
	variable row=v_value
	variable startret,mid,endext
		controlinfo/W=MRViewer de_Viewer_check2

	variable dwell=v_value
	if(row==0)
				
		startRet=0
		Mid= str2num(StringfromList(0,stringbykey("DE_Ind",note(ForceWave),":","\r")	))
		if(dwell==1)
						EndExt=str2num(StringfromList(0,stringbykey("DE_PauseLoc",note(ForceWave),":","\r")	))

		else
				EndExt=str2num(StringfromList(1,stringbykey("DE_Ind",note(ForceWave),":","\r")	))

		endif
	else
		startRet=str2num(StringfromList(2*row-1,stringbykey("DE_Ind",note(ForceWave),":","\r")	))
		Mid= str2num(StringfromList(2*row,stringbykey("DE_Ind",note(ForceWave),":","\r")	))
		if(dwell==1)
					EndExt=str2num(StringfromList(2*row,stringbykey("DE_PauseLoc",note(ForceWave),":","\r")	))

		else
			EndExt=str2num(StringfromList(2*row+1,stringbykey("DE_Ind",note(ForceWave),":","\r")	))

		endif
				
	endif
	Removefromgraph/W=MRViewer#Data YDisp
	Removefromgraph/W=MRViewer#Data YDisp
	Removefromgraph/W=MRViewer#Data YDispSm
	Removefromgraph/W=MRViewer#Data YDispSm
	Removefromgraph/W=MRViewer#Data YDispLessSm
	Removefromgraph/W=MRViewer#Data YDispLessSm
	
	duplicate/o/r=[startRet,EndExt] ForceWave root:DE_Viewer:YDisp
	wave Ydisp=root:DE_Viewer:YDisp
	duplicate/o/r=[startRet,EndExt] SepWave root:DE_Viewer:XDisp
	wave Xdisp=root:DE_Viewer:XDisp
	
	duplicate/o/r=[startRet,EndExt] ForceWaveSm root:DE_Viewer:YDispSm
	wave YDispSm=root:DE_Viewer:YDispSm
	duplicate/o/r=[startRet,EndExt] SepWaveSm root:DE_Viewer:XDispSm
	wave XDispSm=root:DE_Viewer:XDispSm
	duplicate/o/r=[startRet,EndExt] ForceWaveLessSm root:DE_Viewer:YDispLessSm
	wave YDispLessSm=root:DE_Viewer:YDispLessSm
	duplicate/o/r=[startRet,EndExt] SepWaveLessSm root:DE_Viewer:XDispLessSm
	wave XDispLessSm=root:DE_Viewer:XDispLessSm

	Controlinfo de_Viewer_setvar0

	Appendtograph/W=MRViewer#Data YDisp[0,Mid-startRet] vs XDisp[0,Mid-startRet]
	Appendtograph/W=MRViewer#Data YDisp[Mid-startRet,EndExt-startRet] vs XDisp[Mid-startRet,EndExt-startRet] 
	Appendtograph/W=MRViewer#Data YDispLessSm[0,Mid-startRet] vs XDispLessSm[0,Mid-startRet]
	Appendtograph/W=MRViewer#Data YDispLessSm[Mid-startRet,EndExt-startRet] vs XDispLessSm[Mid-startRet,EndExt-startRet] 
	Appendtograph/W=MRViewer#Data YDispSm[0,Mid-startRet] vs XDispSm[0,Mid-startRet]
	Appendtograph/W=MRViewer#Data YDispSm[Mid-startRet,EndExt-startRet] vs XDispSm[Mid-startRet,EndExt-startRet]
	ModifyGraph/W=MRViewer#Data hidetrace(YDisp)=1,hidetrace(YDisp#1)=1
	ModifyGraph/W=MRViewer#Data rgb(YDisp)=(63232,45824,45824),rgb(YDisp#1)=(49152,54784,59648)
	ModifyGraph/W=MRViewer#Data rgb(YDispLessSm)=(61184,28672,28672),rgb(YDispLessSm#1)=(34304,46080,55552)
	ModifyGraph/W=MRViewer#Data rgb(YDispSm)=(58368,6656,7168),rgb(YDispSm#1)=(14848,32256,47104)
	SetDataFolder saveDF
	
end

Static Function LBP(ctrlName,row,col,event) : ListBoxControl

	String ctrlName
	Variable row
	Variable col
	Variable event	//1=mouse down, 2=up, 3=dbl click, 4=cell select with mouse or keys
	//5=cell select with shift key, 6=begin edit, 7=end
	//
	switch (event)
		case -1:
			break
		case 1:
			break
		case 3:
		
			DE_MultiRampViewer#SelectRamp()
			UpdateLocalPoints()
			ListBox de_Viewer_list2 selrow=0
			ListBox de_Viewer_list1 selrow=0
			wave ForceStart=root:DE_Viewer:FullForceSm

			wave Wave1=root:DE_Viewer:CurrUpF
			if(numpnts(Wave1)==0)

			else
			
				wave Wave2=root:DE_Viewer:CurrUpP
				make/o/n=1 root:DE_Viewer:SelUpY, root:DE_Viewer:SelUpX
				wave wave3=root:DE_Viewer:SelUpY
				wave wave4=root:DE_Viewer:SelUpX
				wave3=Wave1[0]
				wave4=pnt2x(ForceStart,Wave2[0])
			endif
			wave Wave1=root:DE_Viewer:CurrDownF
			if(numpnts(Wave1)==0)
			else
				wave Wave2=root:DE_Viewer:CurrDownP
				make/o/n=1 root:DE_Viewer:SelDownY, root:DE_Viewer:SelDownX
				wave wave3=root:DE_Viewer:SelDownY
				wave wave4=root:DE_Viewer:SelDownX
				wave3=Wave1[0]
				wave4 =pnt2x(ForceStart,Wave2[0])
			endif
			//DE_MultiRampViewer#MakeSlopeFits()

			break			
		case 4:	
			DE_MultiRampViewer#SelectRamp()
			UpdateLocalPoints()
			ListBox de_Viewer_list2 selrow=0
			ListBox de_Viewer_list1 selrow=0
			wave ForceStart=root:DE_Viewer:FullForceSm

			wave Wave1=root:DE_Viewer:CurrUpF
			if(numpnts(Wave1)==0)
							make/o/n=0 root:DE_Viewer:SelUpY, root:DE_Viewer:SelUpX

			else
				wave Wave2=root:DE_Viewer:CurrUpP
				make/o/n=1 root:DE_Viewer:SelUpY, root:DE_Viewer:SelUpX
				wave wave3=root:DE_Viewer:SelUpY
				wave wave4=root:DE_Viewer:SelUpX
				wave3=Wave1[0]
				wave4=pnt2x(ForceStart,Wave2[0])
			endif
			wave Wave1=root:DE_Viewer:CurrDownF
			if(numpnts(Wave1)==0)
							make/o/n=0 root:DE_Viewer:SelDownY, root:DE_Viewer:SelDownX

			else
				wave Wave2=root:DE_Viewer:CurrDownP
				make/o/n=1 root:DE_Viewer:SelDownY, root:DE_Viewer:SelDownX
				wave wave3=root:DE_Viewer:SelDownY
				wave wave4=root:DE_Viewer:SelDownX
				wave3=Wave1[0]
				wave4 =pnt2x(ForceStart,Wave2[0])
			endif
			//	DE_MultiRampViewer#MakeSlopeFits()
			break
	endswitch
	
	return 0
End

Static Function LBP1(ctrlName,row,col,event) : ListBoxControl

	String ctrlName
	Variable row
	Variable col
	Variable event	//1=mouse down, 2=up, 3=dbl click, 4=cell select with mouse or keys
	//5=cell select with shift key, 6=begin edit, 7=end

	//
	variable start,stop
	switch (event)
		case -1:
			break
		case 1:
			break
		case 3:


			break			
		case 4:
			
			wave ForceStart=root:DE_Viewer:FullForceSm

			strswitch(CtrlName)
				case "de_Viewer_list1":
					wave Wave1=root:DE_Viewer:CurrUpF
					if(numpnts(Wave1)==0)
					else
						wave Wave2=root:DE_Viewer:CurrUpP
						make/o/n=1 root:DE_Viewer:SelUpY, root:DE_Viewer:SelUpX
						wave wave3=root:DE_Viewer:SelUpY
						wave wave4=root:DE_Viewer:SelUpX
						wave3=Wave1[row]
						wave4=pnt2x(ForceStart,Wave2[row])
						start=pnt2x(ForceStart,Wave2[row])-20e-3
						stop=pnt2x(ForceStart,Wave2[row])+50e-3
					endif
					break
			
				case "de_Viewer_list2":
					wave Wave1=root:DE_Viewer:CurrDownF
					if(numpnts(Wave1)==0)
					else
						wave Wave2=root:DE_Viewer:CurrDownP
						make/o/n=1 root:DE_Viewer:SelDownY, root:DE_Viewer:SelDownX
						wave wave3=root:DE_Viewer:SelDownY
						wave wave4=root:DE_Viewer:SelDownX
						wave3=Wave1[row]
						wave4 =pnt2x(ForceStart,Wave2[row])
						start=pnt2x(ForceStart,Wave2[row])-20e-3
						stop=pnt2x(ForceStart,Wave2[row])+50e-3
					endif


					break
			
					break
			endswitch
			UpdateAxis(start,stop)
	endswitch
	
	return 0
End

Static Function LBP2(ctrlName,row,col,event) : ListBoxControl

	String ctrlName
	Variable row
	Variable col
	Variable event	//1=mouse down, 2=up, 3=dbl click, 4=cell select with mouse or keys
	//5=cell select with shift key, 6=begin edit, 7=end
	string saveDF = GetDataFolder(1)
	variable startRet,Mid,EndExt

	//
	switch (event)
		case -1:
			break
		case 1:
			break
		case 3:
			break			
		case 7:	
			wave CUpLength=root:DE_Viewer:CurrUpLength
			wave CDownLength=root:DE_Viewer:CurrDownLength
			wave/T TCUpLength=root:DE_Viewer:TCurrUpL
			wave/T TCDownLength=root:DE_Viewer:TCurrDownL
			CUpLength=str2num(TCUpLength)
			CDownLength=str2num(TCDownLength)
			//DE_MultiRampViewer#MakeSlopeFits()
			//UpdateLengthsandLoading()
			break
	endswitch
	
	return 0
End


Static Function LoadForce()

	string saveDF = GetDataFolder(1)
	controlinfo de_Viewer_popup0
	SetDataFolder s_value
	
	controlinfo de_Viewer_popup1
	wave ForceStart=$S_value
	wave SepStart=$ReplaceString("Force",S_value,"Sep")
	duplicate/o ForceStart root:DE_Viewer:FullForceTrace
	duplicate/o ForceStart root:DE_Viewer:FullForceSm,root:DE_Viewer:FullForceLessSm
	duplicate/o SepStart root:DE_Viewer:FullSepTrace
	duplicate/o SepStart root:DE_Viewer:FullSepSm, root:DE_Viewer:FullSepLessSm
	Controlinfo de_Viewer_setvar0
	if(V_value>5)
	De_Filtering#FilterForceSep(root:DE_Viewer:FullForceTrace,root:DE_Viewer:FullSepTrace,root:DE_Viewer:FullForceSm,root:DE_Viewer:FullSepSm,"SVG",v_value)
	
	else
		De_Filtering#FilterForceSep(root:DE_Viewer:FullForceTrace,root:DE_Viewer:FullSepTrace,root:DE_Viewer:FullForceSm,root:DE_Viewer:FullSepSm,"TVD",v_value)

	endif
	De_Filtering#FilterForceSep(root:DE_Viewer:FullForceTrace,root:DE_Viewer:FullSepTrace,root:DE_Viewer:FullForceLessSm,root:DE_Viewer:FullSepLessSm,"SVG",51)

	DE_MultiRampViewer#UpdateSegmentList()
end

Static Function UpdateSegmentList()	
	string saveDF = GetDataFolder(1)
	controlinfo de_Viewer_popup0
	SetDataFolder s_value
	make/T/o/n=0 root:DE_Viewer:ListWave1
	wave/T/z LW=root:DE_Viewer:ListWave1
	wave ForceWave=root:DE_Viewer:FullForceTrace

	wave SepWave=root:DE_Viewer:FullSepTrace
	String Me= stringbykey("DE_Ind",note(ForceWave),":","\r")
	
	String Response=""
	SetDataFolder saveDF
	variable Num=ItemsInList(Me)/2
	variable n
	For(n=0;n<num;n+=1)
		Insertpoints n,1,LW
		LW[n]=num2str(n+1)
	endfor
	make/o/n=(n) root:DE_Viewer:SelWave1
	wave SW=root:DE_Viewer:SelWave1
	SW=0
	SW[0]=1

end

Static Function UpdateHistograms()
	wave AUpp=root:DE_Viewer:AllUpPoints
	wave ADown=root:DE_Viewer:AllDownPoints
	
	wave UpForce=root:DE_Viewer:AllUpForce
	wave DownForce=root:DE_Viewer:AllDownForce
	Wave ForceSm=root:DE_Viewer:FullForceSm
	controlinfo de_Viewer_popup0
	SetDataFolder s_value
	controlinfo de_Viewer_popup1
	wave ForceWave=root:DE_Viewer:FullForceTrace
	wave SepWave=root:DE_Viewer:FullSepTrace
	duplicate/o AUpp UpForce
	duplicate/o ADown DownForce
	
	Controlinfo/W=MRViewer de_viewer_check1
	if(numpnts(UpForce)>0)
		if(V_value==0)
		UpForce=ForceSm[AUpp]
		else
		UpForce=ForceWave[AUpp]
	
		endif
	else
	endif


	if(numpnts(DownForce)>0)
		if(V_value==0)
		DownForce=ForceSm[ADown]
		else
		DownForce=ForceWave[ADown]
	
		endif
	else
	endif


	
	Make/N=10/O root:DE_Viewer:AllUpForceHist
	Make/N=10/O root:DE_Viewer:AllDownForceHist

	wave UpForceHist=root:DE_Viewer:AllUpForceHist
	wave DownForceHist=root:DE_Viewer:AllDownForceHist
	if(numpnts(DownForce)>0)
		Histogram/B=1 DownForce,DownForceHist

	else
	endif
	
	if(numpnts(UpForce)>0)
		Histogram/B=1 UpForce,UpForceHist

	else
	endif

end

Static Function LoadRuptures()

	string saveDF = GetDataFolder(1)
	controlinfo de_Viewer_popup2
	SetDataFolder s_value
	controlinfo de_Viewer_popup3
	Wave UpPoints=$S_value
	wave DownPoints=$ReplaceString("U", S_value, "D",1)

	duplicate/o UpPoints root:DE_Viewer:AllUpPoints	
	duplicate/o DownPoints root:DE_Viewer:AllDownPoints	
	duplicate/o UpPoints root:DE_Viewer:AllUpLength
	duplicate/o DownPoints root:DE_Viewer:AllDownLength	
	PareResults( root:DE_Viewer:AllUpPoints)
	PareResults( root:DE_Viewer:AllDownPoints)
	controlinfo de_Viewer_setvar1
	wave UpL= root:DE_Viewer:AllUpLength
	UpL= v_value
	wave DownL=root:DE_Viewer:AllDownLength	
	DownL=v_value
	De_MultiRampViewer#MakeAllLoadingSlopes()
	De_MultiRampViewer#UpdateHistograms()
	De_MultiRampViewer#UpdateLocalPoints()
		
	
	SetDataFolder saveDF

end

Static Function PopMenuProc(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa
	Variable popNum = pa.popNum
	String popStr = pa.popStr
	string saveDF
	strswitch(pa.ctrlName)
	
		case "de_Viewer_popup1":
			switch( pa.eventCode )
				case 2: // mouse up
					DE_MultiRampViewer#LoadForce()
					break
				case -1: // control being killed
					break
			endswitch
			break 
		case "de_Viewer_popup3":
			switch( pa.eventCode )
				case 2: // mouse up
					DE_MultiRampViewer#LoadRuptures()

					break
				case -1: // control being killed
					break
			endswitch
			break 
	endswitch
	return 0
End

Static Function UpdateLocalPoints()
	wave AUpp=root:DE_Viewer:AllUpPoints
	wave ADown=root:DE_Viewer:AllDownPoints
	wave AUpL=root:DE_Viewer:AllUpLength
	wave ADownL=root:DE_Viewer:AllDownLength
	wave AUpLoad=root:DE_Viewer:UpLoad
	wave ADownLoad=root:DE_Viewer:Download
	make/D/o/n=0 root:DE_Viewer:CurrUpF,root:DE_Viewer:CurrUpX,root:DE_Viewer:CurrDownF,root:DE_Viewer:CurrDownX,root:DE_Viewer:CurrDownP,root:DE_Viewer:CurrUpP,root:DE_Viewer:CurrDownLength,root:DE_Viewer:CurrUpLength
	make/D/o/n=0 root:DE_Viewer:CurrUpLoad,root:DE_Viewer:CurrDownLoad
	wave CUpF=	root:DE_Viewer:CurrUpF
	wave CUpX=	root:DE_Viewer:CurrUpX
	wave CDownF= root:DE_Viewer:CurrDownF
	wave CDownX= root:DE_Viewer:CurrDownX
	wave CurrDownP= root:DE_Viewer:CurrDownP
	wave CurrUpP= root:DE_Viewer:CurrUpP
	Wave CurrDownL=root:DE_Viewer:CurrDownLength
	Wave CurrUpL=root:DE_Viewer:CurrUpLength
	wave FW= root:DE_Viewer:YDispSm
	wave FWUnSm=root:DE_Viewer:YDisp
	wave CurrUpLoad=root:DE_Viewer:CurrUpLoad
	wave CurrDownLoad=root:DE_Viewer:CurrDownLoad
	controlinfo de_Viewer_list0
	variable row=v_value, startRet,Mid,EndExt
	
	wave ForceWave=root:DE_Viewer:FullForceTrace
	if(row==0)
			
		startRet=0
		Mid= str2num(StringfromList(0,stringbykey("DE_Ind",note(ForceWave),":","\r")	))
		EndExt=str2num(StringfromList(1,stringbykey("DE_Ind",note(ForceWave),":","\r")	))
	else
			
		startRet=str2num(StringfromList(2*row-1,stringbykey("DE_Ind",note(ForceWave),":","\r")	))
		Mid= str2num(StringfromList(2*row,stringbykey("DE_Ind",note(ForceWave),":","\r")	))
		EndExt=str2num(StringfromList(2*row+1,stringbykey("DE_Ind",note(ForceWave),":","\r")	))
			
	endif
	variable n=0
	variable m=0

	for(n=0;n<numpnts(AUpp);n+=1)
		if(AUpp[n]>=startRet&&AUpp[n]<=EndExt)
			insertpoints m,1,CUpF,CUpX,CurrUpP,CurrUpL,CurrUpLoad
			Controlinfo/W=MRViewer de_viewer_check1
			if(V_value==0)
			CUpX[m]=leftx(FW)+deltax(FW )*(AUpp[n]-startret)
				//CUpX[m]=pnt2x(FW,AUpp[n]-startret)
				CUpF[m]=FW[AUpp[n]-startret]
			else
				CUpX[m]=leftx(FWUnSm)+deltax(FWUnSm)*(AUpp[n]-startret)

			//	CUpX[m]=pnt2x(FWUnSm,AUpp[n]-startret)
				CUpF[m]=FWUnSm[AUpp[n]-startret]
	
			endif


			CurrUpP[m]=AUpp[n]
			CurrUpL[m]=AUpL[n]
			CurrUpLoad[m]=AUpLoad[n]
			m+=1

		endif
	endfor
	n=0
	m=0
	for(n=0;n<numpnts(ADown);n+=1)
		if(ADown[n]>=startRet&&ADown[n]<=EndExt)
		
			insertpoints m,1,CDownF,CDownX,CurrDownP,CurrDownL,CurrDownLoad
			Controlinfo/W=MRViewer de_viewer_check1
			if(V_value==0)
				CDownX[m]=leftx(FW)+deltax(FW )*(ADown[n]-startret)

				//CDownX[m]=pnt2x(FW,ADown[n]-startret)
				CDownF[m]=FW[ADown[n]-startret]
			else
									CDownX[m]=leftx(FWUnSm)+deltax(FWUnSm)*(ADown[n]-startret)

				//CDownX[m]=pnt2x(FWUnSm,ADown[n]-startret)
				CDownF[m]=FWUnSm[ADown[n]-startret]
	
			endif


			CurrDownP[m]=ADown[n]
			CurrDownL[m]=ADownL[n]
			CurrDownLoad[m]=ADownLoad[n]
			m+=1

		endif
	endfor

	if( numpnts(CupF)==0)
		make/t/o/n=(0) root:DE_Viewer:TCurrUpF,root:DE_Viewer:TCurrUpL
		make/o/n=(0) root:DE_Viewer:TCurrUpL_Sel

	else
		make/t/o/n=(numpnts(CUpF)) root:DE_Viewer:TCurrUpF,root:DE_Viewer:TCurrUpL
		make/o/n=(numpnts(CUpF)) root:DE_Viewer:TCurrUpL_Sel

		wave/T TCUp=root:DE_Viewer:TCurrUpF
		wave/T TCUL=root:DE_Viewer:TCurrUpL
		wave TCULSEL=root:DE_Viewer:TCurrUpL_Sel
		TCULSEL=3
		TCUp[]=num2str(CUpF[p])
		TCUL[]=num2str(CurrUpL[p])

	endif
	
	if(numpnts(CDownF)==0)
		make/t/o/n=(0) root:DE_Viewer:TCurrDownF,root:DE_Viewer:TCurrDownL
		make/o/n=(0) root:DE_Viewer:TCurrDownL_Sel
	else
		make/t/o/n=(numpnts(CDownF)) root:DE_Viewer:TCurrDownF,root:DE_Viewer:TCurrDownl
		make/o/n=(numpnts(CDownF)) root:DE_Viewer:TCurrDownl_Sel

		wave/T TCDown=root:DE_Viewer:TCurrDownF
		wave/T TCDownL=root:DE_Viewer:TCurrDownL
		wave TCDownLSel=root:DE_Viewer:TCurrDownL_Sel
		TCDownLSel=3
		TCDown[][0]=num2str(CDownF[p])
		TCDownL[][]=num2str(CurrDownL[p])

	endif

end


Static Function SVP(ctrlName,varNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable varNum
	String varStr
	String varName
	strswitch(ctrlName)
		case "de_Viewer_setvar0":
			Controlinfo de_Viewer_setvar0
			//Smooth/S=2 v_value, root:DE_Viewer:FullForceSlopeSm
			//De_Filtering#FilterForceSep(root:DE_Viewer:FullForceTrace,root:DE_Viewer:FullSepTrace,root:DE_Viewer:FullForceSm,root:DE_Viewer:FullSepSm,"TVD",15000e-12)
			
			if(V_value>5)
			De_Filtering#FilterForceSep(root:DE_Viewer:FullForceTrace,root:DE_Viewer:FullSepTrace,root:DE_Viewer:FullForceSm,root:DE_Viewer:FullSepSm,"SVG",varnum)
	
			else
			De_Filtering#FilterForceSep(root:DE_Viewer:FullForceTrace,root:DE_Viewer:FullSepTrace,root:DE_Viewer:FullForceSm,root:DE_Viewer:FullSepSm,"TVD",varnum)

			endif
			break
	endswitch
End		

Static Function RemoveThisTrace()
	string saveDF = GetDataFolder(1)
	controlinfo de_Viewer_popup0
	SetDataFolder s_value
	

	if(waveExists(root:DE_Viewer:FullForceTrace)==0)
	
		return 0
	else
	
	endif
	
	wave ForceWave=root:DE_Viewer:FullForceTrace
	wave SepWave=root:DE_Viewer:FullSepTrace
	wave ForceSm=root:DE_Viewer:FullForceSm
	wave SepSm=root:DE_Viewer:FullSepSm
	wave ForceWaveLessSm=root:DE_Viewer:FullForceLessSm
	wave SepWaveLessSm=root:DE_Viewer:FullSepLessSm
	
	
	String DE_Ind=stringbykey("DE_Ind",note(ForceWave),":","\r")
	String DE_Dir=stringbykey("DE_Dir",note(ForceWave),":","\r")
	String DE_PauseLoc=stringbykey("DE_PauseLoc",note(ForceWave),":","\r")
	String DE_PauseState=stringbykey("DE_PauseState",note(ForceWave),":","\r")
		
		
	controlinfo de_Viewer_list0
	variable row=v_value
	variable startret,mid,endext
	if(row==0)
				
		startRet=0
		Mid= str2num(StringfromList(0,stringbykey("DE_Ind",note(ForceWave),":","\r")	))
		EndExt=str2num(StringfromList(1,stringbykey("DE_Ind",note(ForceWave),":","\r")	))
	else
		startRet=str2num(StringfromList(2*row-1,stringbykey("DE_Ind",note(ForceWave),":","\r")	))
		Mid= str2num(StringfromList(2*row,stringbykey("DE_Ind",note(ForceWave),":","\r")	))
		EndExt=str2num(StringfromList(2*row+1,stringbykey("DE_Ind",note(ForceWave),":","\r")	))
				
	endif

	deletepoints startRet,(EndExt-startRet), ForceWave,SepWave,ForceSm,SepSm,ForceWaveLessSm,SepWaveLessSm
	string IndString=""
	string DirString=""
	string IDString=""
	variable n=0,Current
	
	variable IndexOne=row*2
	variable indexTwo=row*2+1
	string numbertochange
	For(n=IndexOne;n<(ItemsInList(DE_Ind));n+=1)
		numbertochange=num2istr(str2num(stringfromlist(n,DE_Ind))-(EndExt-startRet))
		DE_Ind=ReplaceStringbyIndex(DE_Ind,numbertochange,n,";")
		numbertochange=num2istr(str2num(stringfromlist(n,DE_PauseLoc))-(EndExt-startRet))
		DE_PauseLoc=ReplaceStringbyIndex(DE_PauseLoc,numbertochange,n,";")
	
	
	endfor
	
 	DE_Ind=RemoveListItem(IndexTwo, DE_Ind)
 	DE_Ind=RemoveListItem(IndexOne, DE_Ind)
 	DE_Dir=RemoveListItem(IndexTwo, DE_Dir)
 	DE_Dir=RemoveListItem(IndexOne, DE_Dir)
 	DE_PauseLoc=RemoveListItem(IndexTwo, DE_PauseLoc)
 	DE_PauseLoc=RemoveListItem(IndexOne, DE_PauseLoc)
 	DE_PauseState=RemoveListItem(IndexTwo, DE_PauseState)
 	DE_PauseState=RemoveListItem(IndexOne, DE_PauseState)

//	For(n=0;n<(ItemsInList(stringbykey("DE_Ind",note(ForceWave),":","\r")));n+=1)
//		Current=str2num(StringfromList(n,stringbykey("DE_Ind",note(ForceWave),":","\r")	))
//		if(StartRet==0)
//
//			if(Current<StartRet)
//
//				IndString+=StringfromList(n,stringbykey("DE_Ind",note(ForceWave),":","\r")	)+";"
//				DirString+=StringfromList(n,stringbykey("DE_Dir",note(ForceWave),":","\r")	)+";"
//				IDString+=StringfromList(n,stringbykey("DE_ID",note(ForceWave),":","\r")	)+";"
//			elseif(Current>EndExt)
//
//				Current-=(EndExt-startRet)
//				IndString+=num2istr(Current)+";"
//				DirString+=StringfromList(n,stringbykey("DE_Dir",note(ForceWave),":","\r")	)+";"
//				IDString+=StringfromList(n,stringbykey("DE_ID",note(ForceWave),":","\r")	)+";"
//			else
//			endif
//			//elseif(EndExt==
//		else
//
//			if(Current<=StartRet)
//
//				IndString+=StringfromList(n,stringbykey("DE_Ind",note(ForceWave),":","\r")	)+";"
//				DirString+=StringfromList(n,stringbykey("DE_Dir",note(ForceWave),":","\r")	)+";"
//				IDString+=StringfromList(n,stringbykey("DE_ID",note(ForceWave),":","\r")	)+";"
//			elseif(Current>EndExt)
//				Current-=(EndExt-startRet)
//				IndString+=num2istr(Current)+";"
//				DirString+=StringfromList(n,stringbykey("DE_Dir",note(ForceWave),":","\r")	)+";"
//				IDString+=StringfromList(n,stringbykey("DE_ID",note(ForceWave),":","\r")	)+";"
//			else
//			endif
//		
//		endif
//	endfor
	string TotalNote=note(ForceWave)
	//print stringbykey("DE_Ind",Totalnote,":","\r")
	
	string S1= ReplaceStringByKey( "DE_Ind",ReplaceStringByKey( "DisplayGraphParms", note(ForceWave),"Color",":","\r"),DE_Ind,":","\r")
	string S2= ReplaceStringByKey( "DE_Dir",S1,DE_Dir,":","\r")
	string S3= ReplaceStringByKey( "DE_PauseLoc",S2,DE_PauseLoc,":","\r")
	string S4= ReplaceStringByKey( "DE_PauseState",S3,DE_PauseState,":","\r")
	note/K  ForceWave,S4

	//print stringbykey("DE_Ind",S4,":","\r")

 	
 	
	DE_MultiRampViewer#UpdateSegmentList()
	De_MultiRampViewer#UpdateHistograms()		

	De_MultiRampViewer#UpdateLocalPoints()
	
	SetDataFolder saveDF

end
	Static Function/S ReplaceStringbyIndex(StringIn,Replacement,index,separator)
	
	string StringIn,Replacement,separator
	variable index
	String StringOut=AddListItem(Replacement, StringIn,separator,index)
	 StringOut=RemoveListItem(index+1, StringOut)
//	variable n
//	for(n=0;n<itemsinlist(StringIn);n+=1)
//		if(n==index)
//				StringOut+=Replacement+";"
//
//		else
//		StringOut+=Stringfromlist(n,StringIn)+";"
//		endif
//	
//	endfor
	return StringOut
	end
Static Function RemoveCurrentMarkers()
	string saveDF = GetDataFolder(1)
	controlinfo de_Viewer_popup0

	SetDataFolder s_value
	
	wave AUpp=root:DE_Viewer:AllUpPoints
	wave ADown=root:DE_Viewer:AllDownPoints
	wave AUpForce=root:DE_Viewer:AllUpForce
	wave ADownForce=root:DE_Viewer:AllDownForce
	wave AUpLoad=root:DE_Viewer:UpLoad
	wave ADownLoad=root:DE_Viewer:DownLoad
	wave CDownForce=root:DE_Viewer:CurrDownF
	wave CUpForce=root:DE_Viewer:CurrUpF
	variable row=v_value, startRet,Mid,EndExt

	variable n=0	
	if(numpnts(CUpforce)==0)
		controlinfo de_Viewer_list0
		row=v_value
		wave ForceWave=root:DE_Viewer:FullForceTrace
		if(row==0)
			
			//startRet=0
			//Mid= str2num(StringfromList(0,stringbykey("DE_Ind",note(ForceWave),":","\r")	))
			EndExt=str2num(StringfromList(1,stringbykey("DE_Ind",note(ForceWave),":","\r")	))
		else
			
			//startRet=str2num(StringfromList(2*row-1,stringbykey("DE_Ind",note(ForceWave),":","\r")	))
			//Mid= str2num(StringfromList(2*row,stringbykey("DE_Ind",note(ForceWave),":","\r")	))
			EndExt=str2num(StringfromList(2*row+1,stringbykey("DE_Ind",note(ForceWave),":","\r")	))
		endif
		
		if(EndExt<AUpp[0])
			AUpp-=+numpnts(root:DE_Viewer:YDisp)
		elseif(EndExt>AUpp[numpnts(Aupp)-1])
		else
			FindLevel/Q/P  AUpp,EndExt
			if(numpnts(AUpp)==V_levelx)
			else
				AUpp[V_levelx,]-=+numpnts(root:DE_Viewer:YDisp)
			endif
		endif
	else
		for(n=0;n<numpnts(CUpforce);n+=1)
			FindValue/T=1e-16 /V=(CUpforce[n]) AUpForce
			deletepoints V_value,1, AUpForce	
			deletepoints V_value,1, AUpp	
			deletepoints V_value,1, AUpLoad	
		endfor

		if(numpnts(AUpp)==V_value)
		else
			AUpp[V_value,]-=+numpnts(root:DE_Viewer:YDisp)
			
		endif
	endif
	n=0
	
	if(numpnts(CDownforce)==0)
		controlinfo de_Viewer_list0
		row=v_value
	
		wave ForceWave=root:DE_Viewer:FullForceTrace
		if(row==0)
			
			//	startRet=0
			//	Mid= str2num(StringfromList(0,stringbykey("DE_Ind",note(ForceWave),":","\r")	))
			EndExt=str2num(StringfromList(1,stringbykey("DE_Ind",note(ForceWave),":","\r")	))
		else
			
			//	startRet=str2num(StringfromList(2*row-1,stringbykey("DE_Ind",note(ForceWave),":","\r")	))
			//	Mid= str2num(StringfromList(2*row,stringbykey("DE_Ind",note(ForceWave),":","\r")	))
			EndExt=str2num(StringfromList(2*row+1,stringbykey("DE_Ind",note(ForceWave),":","\r")	))
		
		endif
		if(EndExt<ADown[0])
			ADown-=+numpnts(root:DE_Viewer:YDisp)
		
		elseif(EndExt>ADown[numpnts(Adown)-1])
		
		else
			FindLevel/Q/P  ADown,EndExt

			if(numpnts(ADown)==V_levelx)
			else
				ADown[(V_LevelX),]-=+numpnts(root:DE_Viewer:YDisp)
			endif
		endif
	else
		for(n=0;n<numpnts(CDownforce);n+=1)
			FindValue/T=1e-16 /V=(CDownforce[n]) ADownForce
			print v_value
			deletepoints V_value,1, ADownForce	
			deletepoints V_value,1, ADown	
			deletepoints V_value,1, ADownLoad	

		endfor
		if(numpnts(ADown)==V_value)
		else
			ADown[V_value,]-=+numpnts(root:DE_Viewer:YDisp)
		endif
	endif
	FIxWaveORder()

	SetDataFolder saveDF

end


Static Function RemoveSelectMarker(Direction,Num)
	String Direction
	Variable Num

	If(cmpstr(Direction,"Up")==0)
		wave AllP=root:DE_Viewer:AllUpPoints
		wave AllF=root:DE_Viewer:AllUpForce
		wave CForce=root:DE_Viewer:CurrUpF
		wave ALoad=root:DE_Viewer:UpLoad
	
	
	elseif(cmpstr(Direction,"Down")==0)
		wave AllP=root:DE_Viewer:AllDownPoints
		wave AllF=root:DE_Viewer:AllDownForce
		wave CForce=root:DE_Viewer:CurrDownF
		wave ALoad=root:DE_Viewer:DownLoad

	else
		print "Bad Direction Dummy"
		return -1
	endif
	
	if(Num>=numpnts(CForce))

		print "Not a valid point, dummy"
		return -1
	endif
	
	FindValue/T=1e-18 /V=(CForce[Num]) AllF
	deletepoints V_value,1, AllF	
	deletepoints V_value,1, AllP	
	deletepoints V_value,1, ALoad	

	FIxWaveORder()

end

Static Function FixWaveOrder()
	wave AllUp=root:DE_Viewer:AllUpPoints
	wave AllDown=root:DE_Viewer:AllDownPoints
	Sort AllUp AllUp
	Sort AllDown AllDown
end

Static Function AddMarker(Direction,Num)
	String Direction
	Variable Num
	
	If(cmpstr(Direction,"Up")==0)
		wave AllP=root:DE_Viewer:AllUpPoints
		wave AllF=root:DE_Viewer:AllUpForce
		wave CForce=root:DE_Viewer:CurrUpF
		wave ALoad=root:DE_Viewer:UpLoad
		Wave AllLength=root:DE_Viewer:AllUpLength
	
	
	elseif(cmpstr(Direction,"Down")==0)
		wave AllP=root:DE_Viewer:AllDownPoints
		wave AllF=root:DE_Viewer:AllDownForce
		wave CForce=root:DE_Viewer:CurrDownF
		wave ALoad=root:DE_Viewer:DownLoad
		Wave AllLength=root:DE_Viewer:AllDownLength

	else
		print "Bad Direction Dummy"
		return -1
	endif
	

	
	wave ForceSm=root:DE_Viewer:FullForceSm
	variable NewPoint=x2pnt(ForceSm,(xcsr(A,"MRViewer#TimeData")))

	Variable NewForce=ForceSm[NewPoint]
	variable location
	if(NewPoint<AllP[0])
	
		location=0
	elseif(NewPoint>AllP[numpnts(AllP)-1])
		location=numpnts(AllP)
	else
		FindLevel/Q/P AllP NewPoint
	
		location=floor(V_levelx)	
	endif
	insertpoints location, 1, AllP,AllLength,ALoad
	AllP[location]=NewPoint
	controlinfo de_Viewer_setvar1
	AllLength[location]=v_value
	ALoad[location]=0
	FIxWaveORder()

end



Static Function TweakSelectMarker(Direction,Num)
	String Direction
	Variable Num
	Variable aExists= strlen(CsrInfo(A,"MRViewer#TimeData")) > 0	// A is a name, not a string
	//print Direction
	
	if( aexists==0)
		print "Drop a Cursor, Dummy"
		return -1
	endif
		
	If(cmpstr(Direction,"Up")==0)
		wave AllP=root:DE_Viewer:AllUpPoints
		wave AllF=root:DE_Viewer:AllUpForce
		wave CForce=root:DE_Viewer:CurrUpF
		wave ALoad=root:DE_Viewer:UpLoad

	
	elseif(cmpstr(Direction,"Down")==0)
		wave AllP=root:DE_Viewer:AllDownPoints
		wave AllF=root:DE_Viewer:AllDownForce
		wave CForce=root:DE_Viewer:CurrDownF
		wave ALoad=root:DE_Viewer:DownLoad

	else
		print "Bad Direction Dummy"
		return -1
	endif

	if(Num>=numpnts(CForce))
		print "Not a valid point, dummy"
		return -1
	endif
	
	wave FullForcetrace=root:DE_Viewer:FullForcetrace
	wave FullSeptrace=root:DE_Viewer:FullSeptrace
	Controlinfo de_Viewer_setvar0
	wave ForceSm=root:DE_Viewer:FullForceSm
	FindValue/T=1e-15 /V=(CForce[Num]) AllF
	variable CurrentP=v_value
	if(v_value<0)
		print "That didn't go right"
		return -1
	endif
	
	variable location=x2pnt(ForceSm,(xcsr(A,"MRViewer#TimeData")))
	//FindValue/T=1e-18 /V=(vcsr(A,"MRViewer#TimeData")) ForceSm
	variable NewSpot
	if(location<AllP[0])
		NewSpot=0
	elseif(location>AllP[numpnts(AllP)-1])
		NewSpot=numpnts(AllP)-1
	
	else
		FindLevel/Q AllP location
		Newspot=floor(V_levelX)
	endif


	Deletepoints CurrentP,1,ALoad,AllP,AllF
	InsertPoints Newspot,1,ALoad,AllP,AllF
	
	AllP[ Newspot]=location
	AllF[Newspot]=vcsr(A,"MRViewer#TimeData")
	ALoad[Newspot]=0

	//	AllP[CurrentP]=location
	//	AllF[CurrentP]=vcsr(A,"MRViewer#TimeData")
	FIxWaveORder()
	UpdateLocalPoints()
	//MakeSlopeFits()
	//UpdateLengthsandLoading()
end

Static Function MakeAllLoadingSlopes()
	wave UpPoints=root:DE_Viewer:AllUpPoints
	wave DownPoints=root:DE_Viewer:AllDownPoints
	wave UpLength=root:DE_Viewer:AllUpLength
	wave Downlength=root:DE_Viewer:AllDownLength
	wave ForceWave=root:DE_Viewer:FullForceTrace
	wave SepWave=root:DE_Viewer:FullSepTrace
	make/free/n=0 ForceSm,SepSm
	Controlinfo de_Viewer_setvar0
	wave ForceSm=root:DE_Viewer:FullForceSm
	duplicate/o UpPoints root:DE_viewer:UpLoad
	duplicate/o DownPoints root:DE_Viewer:DownLoad
	wave UpLoad= root:DE_viewer:UpLoad
	wave DownLoad= root:DE_viewer:DownLoad



	variable n=0, StartP
	for(n=0;n<numpnts(UpPoints);n+=1)
		if((UpPoints[n]-UpLength[n])<0)
			StartP=0
		else
			StartP=UpPoints[n]-UpLength[n]
		endif
		//CurveFit/Q/NTHR=0 line ForceSm[StartP,UpPoints[n]] 
		CurveFit/Q/NTHR=0 line ForceWave[StartP,UpPoints[n]] 
		wave w_coef
		UpLoad[n]=w_coef[1]
	endfor

	n=0

	for(n=0;n<numpnts(DownPoints);n+=1)
		if((DownPoints[n]-DownLength[n])<0)
			StartP=0
		else
			StartP=DownPoints[n]-DownLength[n]
		endif
		//CurveFit/Q/NTHR=0 line ForceSm[StartP,DownPoints[n]] 
		CurveFit/Q/NTHR=0 line ForceWave[StartP,DownPoints[n]] 
		wave w_coef
		DownLoad[n]=w_coef[1]		
	endfor
		wave w_coef,W_Sigma
		killwaves w_coef,W_Sigma

end

Static Function DeleteFits()
	//	
	string saveDF = GetDataFolder(1)

	SetDataFolder root:

	string wavelists=WaveList("SlopeFit*", ";","" )
	variable numtoremove=ItemsinList(wavelists)
	variable n=0
	
	for(n=0;n<numtoremove;n+=1)
		killwaves/z $StringFromList(n, wavelists)
 
	endfor
	SetDataFolder saveDF

end

Static Function ClearFitsOff()
	string wavelists=TraceNameList("MRViewer#TimeData", ";", 1)
	string matches=ListMatch(wavelists,"SlopeFit*",";")

	variable numtoremove=ItemsinList(matches)
	variable n=0
	
	for(n=0;n<numtoremove;n+=1)
		RemoveFromGraph/W=MRViewer#TimeData $StringFromList(n, matches)
  
	endfor
end



Static Function ResettoBase()
	DE_MultiRampViewer#LoadForce()
	DE_MultiRampViewer#LoadRuptures()

end

Static Function SaveOut()

	controlinfo de_Viewer_popup1
	String x=S_Value+"_Adj"
	controlinfo de_Viewer_popup3
	String y=S_Value+"_Adj"
	Prompt x, "Enter X component: "		// Set prompt for x param
	Prompt y, "Enter Y component: "		// Set prompt for y param
	DoPrompt "Enter X and Y", x, y
	if (V_Flag)
		return -1								// User canceled
	endif
	
	duplicate/o root:DE_Viewer:FullForceTrace $X
	duplicate/o root:DE_Viewer:FullSepTrace $ReplaceString("Force",X,"Sep")
	duplicate/o root:DE_Viewer:AllUpPoints $Y
	duplicate/o root:DE_Viewer:AllDownPoints $ReplaceString("U",Y,"D",1)

end

	
Static Function ButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			strswitch(ba.ctrlName)
				case "de_Viewer_but0":
					DE_MultiRampViewer#RemoveCurrentMarkers()

					DE_MultiRampViewer#RemoveThisTrace()
					break
			
				case "de_Viewer_but1":
					controlinfo de_Viewer_list1
					DE_MultiRampViewer#TweakSelectMarker("Up",v_value)
					De_MultiRampViewer#UpdateHistograms()
					//MakeAllLoadingSlopes()
					DE_MultiRampViewer#SelectRamp()
					UpdateLocalPoints()

					break
			
				case "de_Viewer_but2":
					controlinfo de_Viewer_list2
					DE_MultiRampViewer#TweakSelectMarker("Down",v_value)	
					De_MultiRampViewer#UpdateHistograms()
					//MakeAllLoadingSlopes()
					DE_MultiRampViewer#SelectRamp()
					UpdateLocalPoints()
					break
			
			
				case "de_Viewer_but3":
					DE_MultiRampViewer#ResettoBase()
					break
			
				case "de_Viewer_but4":
					DE_MultiRampViewer#SaveOut()
					break
			
				case "de_Viewer_but5":
					//DE_MultiRampViewer#UpdateLengths()
					break
				case "de_Viewer_but6":
					controlinfo de_Viewer_list1
					RemoveSelectMarker("Up",v_value)
					De_MultiRampViewer#UpdateHistograms()
					//MakeAllLoadingSlopes()
					DE_MultiRampViewer#SelectRamp()
					UpdateLocalPoints()
					break
				case "de_Viewer_but7":
					controlinfo de_Viewer_list2
					RemoveSelectMarker("Down",v_value)
					De_MultiRampViewer#UpdateHistograms()
					//MakeAllLoadingSlopes()
					DE_MultiRampViewer#SelectRamp()
					UpdateLocalPoints()
					break			
				case "de_Viewer_but8":
					controlinfo de_Viewer_list1

					AddMarker("Up",v_value)
					De_MultiRampViewer#UpdateHistograms()
					//MakeAllLoadingSlopes()
					DE_MultiRampViewer#SelectRamp()
					UpdateLocalPoints()
					break
				case "de_Viewer_but9":
					controlinfo de_Viewer_list2

					AddMarker("Down",v_value)
					De_MultiRampViewer#UpdateHistograms()
					//MakeAllLoadingSlopes()

					DE_MultiRampViewer#SelectRamp()
					UpdateLocalPoints()
					break
				case "de_Viewer_but10":
					PlotOverlays()
					break				
				case "de_Viewer_but11":
					PlotStartandEnd()
					break
				case "de_Viewer_but12":
					PlotWLC()
					break
				case "de_Viewer_but13":
					SetAxis/A/W=MRViewer#TimeData
					break
			
			
			endswitch
	
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End
Static Function UpdateAxis(start,stop)
	variable start,stop
	SetAxis/W=MRViewer#TimeData bottom start,stop
	SetAxis/A=2/W=MRViewer#TimeData left

end
Static Function PlotOverlays()
	wave ForceSm=root:DE_Viewer:FullForceSm
	wave SepSm=root:DE_Viewer:FullSepSm
	
	DoWindow Overlays
	if(V_flag==1)
		killwindow Overlays
	endif
	display/n=Overlays ForceSm vs SepSm
end

Static Function NumberofRampsinTrace(WaveIn)
	wave Wavein
	String Me= stringbykey("DE_Ind",note(Wavein),":","\r")
	String Response=""
	//	SetDataFolder save
	variable Num=ItemsInList(Me)/2
	return Num


end


Static Function PlotStartandEnd()


	controlinfo de_Viewer_popup1
	String OriginalWave=S_value
	variable loc=strsearch( OriginalWave,"Force",0)
	String Result=OriginalWave[loc-4,loc-1]
	String ExtNUmber
	sprintf ExtNUmber "%04.0f", (str2num(Result)-1)
	variable a=strsearch(OriginalWave,"Force",0)
	String ExtWave=OriginalWave[0,a-1]+"Force"
	ExtWave=ReplaceString(Result,ExtWave,ExtNUmber)+"_ext"
	//controlinfo/W=MRViewer de_Viewer_list0
	variable numRamps=NumberofRampsinTrace($OriginalWave)
	String RetNUmber
	sprintf RetNUmber "%04.0f", (str2num(Result)+numRamps)
	a=strsearch(OriginalWave,"Force",0)
	String RetWave=OriginalWave[0,a-1]+"Force"
	 RetWave=ReplaceString(Result,RetWave,RetNUmber)+"_towd"
	DoWindow StartandEnd
	if(V_flag==1)
		killwindow StartandEnd
	endif
	if(DataFolderExists("HoldingWaves"))
		wave EWave=$("root:HoldingWaves:"+ExtWave)
		wave RWave=$("root:HoldingWaves:"+RetWave)
	else
		wave EWave=$ExtWave
		wave RWave=$RetWave
	endif
	Display/n=StartandEnd RWave
	appendtograph/w=StartandEnd  EWave
	appendtograph/w=StartandEnd  EWave[numpnts(EWave)*.02,numpnts(EWave)*.12]
	appendtograph/w=StartandEnd  RWave[numpnts(RWave)*.88,numpnts(RWave)*.98]
	
	ModifyGraph/w=StartandEnd offset($(nameofwave(RWave)))={0,-2e-11};DelayUpdate
	ModifyGraph/w=StartandEnd offset($(nameofwave(RWave)+"#1"))={0,-2e-11}

	ModifyGraph rgb($(nameofwave(RWave)+"#1"))=(0,0,0);DelayUpdate
	ModifyGraph rgb($(nameofwave(EWave)+"#1"))=(0,0,0)	
	wavestats/q/r=[numpnts(EWave)*.02,numpnts(EWave)*.12] EWave
	print v_avg
	wavestats/q/r=[numpnts(RWave)*.88,numpnts(RWave)*98] RWave
	print v_avg

end

Static Function PlotWLC()
	controlinfo de_Viewer_popup1
	String OriginalWave=S_value
	variable loc=strsearch( OriginalWave,"Force",0)
	String Result=OriginalWave[loc-4,loc-1]
	
	
	String FirstNumStr
	sprintf FirstNumStr "%04.0f", (str2num(Result)-1)
	variable a=strsearch(OriginalWave,"Force",0)
	String ForceFirst=OriginalWave[0,a-1]+"Force"

	ForceFirst=ReplaceString(Result,ForceFirst,FirstNumStr)+"_ret"
	String SepFirst=ReplaceString("Force",ForceFirst,"Sep")
	
	variable numRamps=NumberofRampsinTrace($OriginalWave)
	
	String SecondNumStr
	sprintf SecondNumStr "%04.0f", (str2num(Result)+numRamps)
	 a=strsearch(OriginalWave,"Force",0)
	String ForceSecond=OriginalWave[0,a-1]+"Force"

	 ForceSecond=ReplaceString(Result,ForceSecond,SecondNumStr)+"_towd"
	String SepSecond=ReplaceString("Force",ForceSecond,"Sep")

	DoWindow StartandEnd
	if(V_flag==1)
		killwindow StartandEnd
	endif
	make/o/n=0 root:DE_Viewer:WLCwaveForce,root:DE_Viewer:WLCwaveSep
	wave WLCForce=root:DE_Viewer:WLCwaveForce
	wave WLCSep=root:DE_Viewer:WLCwaveSep

	MakeAlpha3DWLC(WLCForce,WLCSep)
	WLCForce*=-1	
	
	controlinfo/W=MRViewer de_viewer_check0

	if(DataFolderExists("HoldingWaves"))
		wave FSecond=$("root:HoldingWaves:"+ForceSecond)
		wave SSecond=$("root:HoldingWaves:"+SepSecond)
		wave FFirst=$("root:HoldingWaves:"+ForceFirst)
		wave SFirst=$("root:HoldingWaves:"+SepFirst)
	else
		wave FSecond=$(ForceSecond)
		wave SSecond=$(SepSecond)
		wave FFirst=$(ForceFirst)
		wave SFirst=$(SepFirst)
	endif
	Display/n=StartandEnd  FSecond vs SSecond
	appendtograph/w=StartandEnd FFirst vs SFirst
	appendtograph/w=StartandEnd WLCForce vs WLCSep
	if(v_value==1)
		wave RampForce=root:DE_Viewer:FullForceSm
		wave RampSep=root:DE_Viewer:FullSepSm
		
	else
		wave RampForce=root:DE_Viewer:YDispSm
		wave RampSep=root:DE_Viewer:XDispSm

	endif
	appendtograph/w=StartandEnd RampForce vs RampSep
	ModifyGraph/w=StartandEnd rgb($nameofwave(FSecond))=(19712,44800,18944);DelayUpdate
	ModifyGraph/w=StartandEnd rgb($nameofwave(FFirst))=(58368,6656,7168)
	ModifyGraph/w=StartandEnd lstyle($nameofwave(WLCForce))=3,lsize($nameofwave(WLCForce))=2;DelayUpdate
	ModifyGraph/w=StartandEnd rgb($nameofwave(WLCForce))=(0,0,0)
	ModifyGraph/w=StartandEnd muloffset={0,-1}
	ModifyGraph rgb($nameofwave(RampForce))=(14848,32256,47104)

end

Static Function MakeAlpha3DWLC(OutForceWave,OutSepWave)

	wave OutForceWave,OutSepWave
	make/free/n=1000 Force,Ext
	variable LCend=150e-9
	variable LCddflN4=15e-9
	variable LCAlpha=25e-9
	make/free/n=(5,3) WLCParms
	
	WLCParms[][0]=.4e-9
	WLCParms[][1]=298


	WLCParms[0][2]=LCend-2*LCddflN4-LCAlpha
	WLCParms[1][2]=LCend-2*LCddflN4
	WLCParms[2][2]=LCend-LCddflN4
	WLCParms[3][2]=LCend
	DE_WLC#WLC_ArbWave(WLCParms,Force,Ext)

	duplicate/o Force OutForceWave
	duplicate/o Ext OutSepWave

end

Static Function CheckProc(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Macro MultiRampViewer() : Panel
	//
	PauseUpdate; Silent 1		// building window...
	NewPanel/N=MRViewer /W=(50,50,1300,750)
	NewDataFolder/o root:DE_Viewer
	make/T/o/n=0 root:DE_Viewer:ListWave1
	make/o/n=0 root:DE_Viewer:SelWave1
	PopupMenu de_Viewer_popup0,pos={250,2},size={129,21}
	PopupMenu de_Viewer_popup0,mode=1,popvalue="X",value= #"DE_PanelProgs#ListFolders()"
	PopupMenu de_Viewer_popup1,pos={250,40},size={129,21},proc=DE_MultiRampViewer#PopMenuProc
	PopupMenu de_Viewer_popup1,mode=1,popvalue="X",value= #"DE_MultiRampViewer#ListWaves(\"de_Viewer_popup0\",\"*Force\")"
	
	PopupMenu de_Viewer_popup2,pos={550,2},size={129,21}
	PopupMenu de_Viewer_popup2,mode=1,popvalue="X",value= #"DE_PanelProgs#ListFolders()"
	PopupMenu de_Viewer_popup3,pos={550,40},size={129,21},proc=DE_MultiRampViewer#PopMenuProc
	PopupMenu de_Viewer_popup3,mode=1,popvalue="X",value= #"DE_MultiRampViewer#ListWaves(\"de_Viewer_popup0\",\"*RupPntU*\")"
	
	make/o/n=0 root:DE_Viewer:YDisp, root:DE_Viewer:XDisp,root:DE_Viewer:YDispSm, root:DE_Viewer:XDispSm,root:DE_Viewer:YDispLessSm, root:DE_Viewer:XDispLessSm

	display/host=MRViewer/N=Data/W=(250,100,700,400) root:DE_Viewer:YDisp vs root:DE_Viewer:XDisp
	appendtograph/W=MRViewer#Data root:DE_Viewer:YDisp vs root:DE_Viewer:XDisp
	appendtograph/W=MRViewer#Data root:DE_Viewer:YDispLessSm vs root:DE_Viewer:XDispLessSm
	appendtograph/W=MRViewer#Data root:DE_Viewer:YDispLessSm vs root:DE_Viewer:XDispLessSm

	appendtograph/W=MRViewer#Data root:DE_Viewer:YDispSm vs root:DE_Viewer:XDispSm
	appendtograph/W=MRViewer#Data root:DE_Viewer:YDispSm vs root:DE_Viewer:XDispSm

	ModifyGraph/W=MRViewer#Data rgb(YDisp)=(63232,45824,45824),rgb(YDisp#1)=(49152,54784,59648);DelayUpdate
	ModifyGraph/W=MRViewer#Data rgb(YDispSm)=(58368,6656,7168),rgb(YDispSm#1)=(14848,32256,47104)
	
	
	display/host=MRViewer/N=TimeData/W=(250,400,700,700) root:DE_Viewer:YDisp 
	appendtograph/W=MRViewer#TimeData root:DE_Viewer:YDispLessSm
	appendtograph/W=MRViewer#TimeData root:DE_Viewer:YDispSm


	ModifyGraph/W=MRViewer#TimeData rgb(YDisp)=(63232,45824,45824)
	ModifyGraph/W=MRViewer#TimeData hidetrace(YDisp)=1
		ModifyGraph/W=MRViewer#TimeData rgb(YDispLessSm)=(65280,48896,48896)

	ModifyGraph/W=MRViewer#TimeData rgb(YDispSm)=(58368,6656,7168)

	ListBox de_Viewer_list0,pos={10,40},size={100,100},proc=DE_MultiRampViewer#LBP,listWave=root:DE_Viewer:ListWave1
	ListBox de_Viewer_list0,row= 0,mode=2,selRow= 0//,selWave=root:DE_Viewer:SelWave1
	
	SetVariable de_Viewer_setvar0,pos={310,10},size={100,20},proc=DE_MultiRampViewer#SVP,value= _NUM:0,title="Smoothing"
	SetVariable de_Viewer_setvar0 value= _NUM:5e-9
	SetVariable de_Viewer_setvar1,pos={710,10},size={180,20},proc=DE_MultiRampViewer#SVP,value= _NUM:0,title="Starting BackGround"
	SetVariable de_Viewer_setvar1 value= _NUM:1000,limits={1,inf,1}
	
	
	make/D/o/n=0 root:DE_Viewer:AllUpPoints,root:DE_Viewer:AllDownPoints,root:DE_Viewer:CurrUpF,root:DE_Viewer:CurrUpX,root:DE_Viewer:CurrDownF,root:DE_Viewer:CurrDownX
	make/D/o/n=0 root:DE_Viewer:SelUpX,root:DE_Viewer:SelUpY,root:DE_Viewer:SelDownX,root:DE_Viewer:SelDownY
	appendtograph/W=MRViewer#TimeData  root:DE_Viewer:CurrUpF vs root:DE_Viewer:CurrUpX
	appendtograph/W=MRViewer#TimeData root:DE_Viewer:CurrDownF vs root:DE_Viewer:CurrDownX
	appendtograph/W=MRViewer#TimeData root:DE_Viewer:SelUpY vs root:DE_Viewer:SelUpX
	appendtograph/W=MRViewer#TimeData root:DE_Viewer:SelDownY vs root:DE_Viewer:SelDownX
	ModifyGraph mode(SelUpY)=3,marker(SelUpY)=17,mode(SelDownY)=3,marker(SelDownY)=17
	ModifyGraph msize(SelUpY)=4,rgb(SelUpY)=(14848,32256,47104),msize(SelDownY)=4;DelayUpdate
	ModifyGraph rgb(SelDownY)=(19712,44800,18944)
	ModifyGraph/W=MRViewer#TimeData mode(CurrUpF)=3,mode(CurrDownF)=3
	ModifyGraph/W=MRViewer#TimeData marker(CurrUpF)=19,rgb(CurrUpF)=(2816,5632,8192),marker(CurrDownF)=19;DelayUpdate
	ModifyGraph/W=MRViewer#TimeData rgb(CurrDownF)=(29440,0,58880)
	Cursor/W=MRViewer#TimeData A, YDispSm, 0	// cursor A on first point of myWave

	make/t/o/n=(0) root:DE_Viewer:TCurrUpF,root:DE_Viewer:TCurrDownF,root:DE_Viewer:TCurrUpL,root:DE_Viewer:TCurrDownL
	make/o/n=0 root:DE_Viewer:TCurrUpL_Sel,root:DE_Viewer:TCurrDownL_Sel
	ListBox de_Viewer_list1,pos={10,200},size={100,100},proc=DE_MultiRampViewer#LBP1,listWave=root:DE_Viewer:TCurrUpF
	ListBox de_Viewer_list1,row= 0,mode=2,selRow= 0//selWave=root:DE_Viewer:SelWave1
	ListBox de_Viewer_list2,pos={10,400},size={100,100},proc=DE_MultiRampViewer#LBP1,listWave=root:DE_Viewer:TCurrDownF
	ListBox de_Viewer_list2,row= 0,mode=2,selRow= 0//selWave=root:DE_Viewer:SelWave1
	ListBox de_Viewer_list3,pos={110,200},size={100,100},proc=DE_MultiRampViewer#LBP2,listWave=root:DE_Viewer:TCurrUpL
	ListBox de_Viewer_list3,row= 0,mode=2,selRow= 0,selWave=root:DE_Viewer:TCurrUpL_Sel
	ListBox de_Viewer_list4,pos={110,400},size={100,100},proc=DE_MultiRampViewer#LBP2,listWave=root:DE_Viewer:TCurrDownL
	ListBox de_Viewer_list4,row= 0,mode=2,selRow= 0,selWave=root:DE_Viewer:TCurrDownL_Sel

	make/D/o/n=0 root:DE_Viewer:AllUpForceHist,root:DE_Viewer:AllDownForceHist,root:DE_Viewer:AllUpForce,root:DE_Viewer:AllDownForce
	display/host=MRViewer/N=Hist/W=(750,100,950,400) root:DE_Viewer:AllUpForceHist
	appendtograph/W=MRViewer#Hist root:DE_Viewer:AllDownForceHist
	ModifyGraph/W=MRViewer#Hist rgb(AllUpForceHist)=(19712,44800,18944)
	ModifyGraph/W=MRViewer#Hist rgb(AllDownForceHist)=(0,0,65280)
	ModifyGraph/W=MRViewer#Hist mode=5,hbFill=5
	
	Button de_viewer_but0,pos={750,500},size={150,25},proc=DE_MultiRampViewer#ButtonProc,title="Delete This Trace"
	Button de_viewer_but1,pos={910,540},size={150,25},proc=DE_MultiRampViewer#ButtonProc,title="Modify Selected Unfolding"
	Button de_viewer_but6,pos={1070,540},size={150,25},proc=DE_MultiRampViewer#ButtonProc,title="Delete Selected Unfolding"
	Button de_viewer_but8,pos={750,540},size={150,25},proc=DE_MultiRampViewer#ButtonProc,title="Add Unfolding"

	Button de_viewer_but2,pos={910,580},size={150,25},proc=DE_MultiRampViewer#ButtonProc,title="Modify Selected Folding"
	Button de_viewer_but7,pos={1070,580},size={150,25},proc=DE_MultiRampViewer#ButtonProc,title="Delete Selected Folding"
	Button de_viewer_but9,pos={750,580},size={150,25},proc=DE_MultiRampViewer#ButtonProc,title="Add Folding"

	Button de_viewer_but3,pos={750,620},size={150,25},proc=DE_MultiRampViewer#ButtonProc,title="Reset"
	Button de_viewer_but4,pos={750,660},size={150,25},proc=DE_MultiRampViewer#ButtonProc,title="Save Working"
	Button de_viewer_but10,pos={10,550},size={150,25},proc=DE_MultiRampViewer#ButtonProc,title="PlotOverlays!"
	Button de_viewer_but11,pos={10,585},size={150,25},proc=DE_MultiRampViewer#ButtonProc,title="PlotBeginningAndEnd!"
	Button de_viewer_but12,pos={10,645},size={150,25},proc=DE_MultiRampViewer#ButtonProc,title="PlotWLC"
	
	Button de_Viewer_but13,pos={700,450},size={50,25},proc=DE_MultiRampViewer#ButtonProc,title="Axis"
	
	
	CheckBox de_viewer_check0 title="Plot ALL?",pos={200,645},size={150,25},proc=DE_MultiRampViewer#CheckProc
	CheckBox de_viewer_check1 title="Use Raw Wave?",pos={450,10},size={150,25},proc=DE_MultiRampViewer#CheckProc
	CheckBox de_viewer_check2 title="hide Puase?",pos={450,40},size={150,25},proc=DE_MultiRampViewer#CheckProc

	ControlUpdate/A/W=MRViewer
	//Button de_viewer_but5,pos={100,750},size={150,25},proc=DE_MultiRampViewer#ButtonProc,title="Update Lengths"
	variable/G CurrentUps,CurrentDowns,Ordering
	CurrentUps:=numpnts(root:DE_Viewer:TCurrUpF)
	CurrentDowns:=numpnts(root:DE_Viewer:TCurrDownF)
	make/o/n=0 root:DE_Viewer:CurrUpP,root:DE_Viewer:CurrDownP
	Ordering:=DE_MultiRampViewer#CheckWaveForOrdering(root:DE_Viewer:CurrUpP,root:DE_Viewer:CurrDownP) 
	ValDisplay de_viewer_valdisp0 pos={50,300},size={150,25},value=CurrentUps,title="Ups"
	ValDisplay de_viewer_valdisp1 pos={50,510},size={150,25},value=CurrentDowns,title="Downs"
	ValDisplay de_viewer_valdisp2 pos={50,150},size={150,25},value=Ordering,title="Ordered?"

EndMacro

Static Function CheckWaveForOrdering(w0,w1)
	wave w0,w1

	if(numpnts(w1)==0||numpnts(w0)==0)
		return 0
	endif

	variable n,Errors
	if(w0[0]>w1[0])
		Errors+=1
	endif
	for(n=1;n<numpnts(w0);n+=1)
		//for(n=1;n<3;n+=1)
		if(n>=numpnts(w1))
		Errors+=1
	
		elseif(w0[n]<w1[n-1]||w0[n]>w1[n])
			Errors+=1
		endif
	
	endfor

	if( Errors==0)
		return 0
	else
		return 1
	endif
end


Static Function/S ListWaves(ControlStr,SrcStr)
	string ControlStr,SrcStr
	String saveDF

	saveDF = GetDataFolder(1)
	controlinfo $ControlStr
	SetDataFolder s_value
	String list = WaveList(SrcStr, ";", "")
	SetDataFolder saveDF
	return list

end

Static Function PareResults(Wave1)
	wave wave1

	variable n=0
	if(numpnts(wave1)==0||numpnts(wave1)==1)
	return 0
	endif
	
	do
	
		if(Wave1[n+1]==Wave1[n])
			deletepoints (n+1),1, Wave1
		else
			n+=1
		endif
		
	while(n<(numpnts(wave1)-1))


End


//
//Static Function UpdateLengthsandLoading()
//	wave CUpLength=root:DE_Viewer:CurrUpLength
//	wave CDownLength=root:DE_Viewer:CurrDownLength
//	wave CUpP=root:DE_Viewer:CurrUpP
//	wave CDownP=root:DE_Viewer:CurrDownP
//	wave CUpLoad=root:DE_Viewer:CurrUpLoad
//	wave CDownLoad=root:DE_Viewer:CurrDownLoad
//
//	
//	wave AUpLength=root:DE_Viewer:AllUpLength
//	wave ADownLength=root:DE_Viewer:AllDownLength
//	wave AUpP=root:DE_Viewer:AllUpPoints
//	wave ADownP=root:DE_Viewer:AllDownPoints
//	wave AUpLoad=root:DE_Viewer:UpLOad
//	wave ADownLoad=root:DE_Viewer:DownLoad
//	variable n=0
//	For(n=0;n<numpnts(CUpP);n+=1)
//		FindValue/T=1e-18 /V=(CUpP[n]) AUpP
//		AUpLength[v_value]=CUpLength[n]
//		AUpLoad[V_value]=CUpLoad[n]
//	endfor
//	n=0
//	For(n=0;n<numpnts(CDownP);n+=1)
//		FindValue/T=1e-18 /V=(CDownP[n]) ADownP
//		ADownLength[v_value]=CDownLength[n]
//		ADownLoad[V_value]=CDownLoad[n]
//	endfor
//	//MakeAllLoadingSlopes()
//end
//Static Function MakeSlopeFits()
//	DE_MultiRampViewer#ClearFitsOff()
//	DE_MultiRampViewer#DeleteFits()
//	
//	wave AUpp=root:DE_Viewer:AllUpPoints
//	wave ADown=root:DE_Viewer:AllDownPoints
//	variable Length
//	wave ForceWave=root:DE_Viewer:FullForceTrace
//	wave SepWave=root:DE_Viewer:FullSepTrace
//	wave ForceSm=root:DE_Viewer:FullForceSm
//	wave SepSm=root:DE_Viewer:FullSepSm
//	//make/free/n=0 ForceSm,SepSm
//	//Controlinfo de_Viewer_setvar0
//	//De_Filtering#FilterForceSep(ForceWave,SepWave,ForceSm,SepSm,"SVG",v_value)
//	//De_Filtering#FilterForceSep(ForceWave,SepWave,ForceSm,SepSm,"TVD",5000e-12)
//	wave CUpPoint=root:DE_Viewer:CurrUpP
//	wave CDownPoint=root:DE_Viewer:CurrDownP
//	wave CUpLength=root:DE_Viewer:CurrUpLength
//	wave CDownLength=root:DE_Viewer:CurrDownLength
//
//	wave CUpLoad=root:DE_Viewer:CurrUpLoad
//	wave CDownLoad=root:DE_Viewer:CurrDownLoad
//
//	variable n=0
//	variable m=0
//	variable SummedCount=0
//	string NameHere
//	for(n=0;n<numpnts(CUpPoint);n+=1)
//		length=CUpLength[n]
//		NameHere=("SlopeFit"+num2str(SummedCount))
//		make/o/n=(Length+1) $NameHere
//		SetScale/P x pnt2x(ForceWave,(CUpPoint[n]-Length)),DimDelta(ForceWave,0),"s", $NameHere
//		//duplicate/free/r=[CUpPoint[n]-Length,CUpPoint[n]] ForceSm Test
//		duplicate/free/r=[CUpPoint[n]-Length,CUpPoint[n]] ForceWave Test
//		CurveFit/Q/NTHR=0 line Test /D=$NameHere
//		wave w_coef
//		Appendtograph/W=MRViewer#TimeData $NameHere
//		ModifyGraph/W=MRViewer#TimeData rgb($NameHere)=(0,0,0)
//		CUpLoad[n]=w_coef[1]
//		m+=1
//		SummedCount+=1
//
//	endfor
//	//	
//	n=0
//	m=0
//	for(n=0;n<numpnts(CDownPoint);n+=1)
//		length=CDownLength[n]
//		NameHere=("SlopeFit"+num2str(SummedCount))
//		make/o/n=(Length+1) $NameHere
//		SetScale/P x pnt2x(ForceWave,(CDownPoint[n]-Length)),DimDelta(ForceWave, 0),"s", $NameHere
//		//	duplicate/free/r=[CDownPoint[n]-Length,CDownPoint[n]] ForceSm Test
//
//		duplicate/free/r=[CDownPoint[n]-Length,CDownPoint[n]] ForceWave Test
//		CurveFit/Q/NTHR=0 line Test /D=$NameHere
//		wave w_coef
//
//		Appendtograph/W=MRViewer#TimeData $NameHere
//		ModifyGraph/W=MRViewer#TimeData rgb($NameHere)=(29440,0,58880)
//		CDownLoad[n]=w_coef[1]
//
//		m+=1
//		SummedCount+=1
//
//	endfor
//	wave w_coef,W_Sigma
//	killwaves w_coef,W_Sigma
//end
//
//Static Function UpdateCurrentSlopes()
//	
//
//end