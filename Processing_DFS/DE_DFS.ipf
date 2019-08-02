#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma ModuleName = DE_DFS

#include ":\Processing_Smth\DE_Filtering"
#include ":\Misc_PanelPrograms\AsylumNaming"
#include ":\Misc_PanelPrograms\Panel Progs"
#include "SCholl_Panel"
Static Function LBP(ctrlName,row,col,event) : ListBoxControl

	String ctrlName
	Variable row
	Variable col
	Variable event	//1=mouse down, 2=up, 3=dbl click, 4=cell select with mouse or keys
//					//5=cell select with shift key, 6=begin edit, 7=end
	wave/t/z FL=root:DFSPanel:FolderList
	wave/t/z LW1=root:DFSPanel:ListWave1
	wave/t/z LW2=root:DFSPanel:ListWave2

	controlinfo de_DFS_setvar1
	string search1= S_Value
	controlinfo de_DFS_setvar2
	string search2= S_Value
	string/g CurrFold

	switch (event)
		case -1:
			break
		case 1:
			break
		case 3:
			DE_PanelProgs#ListWaves($Fl[row],Search1,LW1)
			DE_PanelProgs#ListWaves($Fl[row],Search2,LW2)

			CurrFold=FL[row]
			
		case 4:	
			DE_PanelProgs#ListWaves($Fl[row],Search1,LW1)
			DE_PanelProgs#ListWaves($Fl[row],Search2,LW2)
			CurrFold=FL[row]
	endswitch
	
	return 0
End


Static Function LBP_1(ctrlName,row,col,event) : ListBoxControl
//
	String ctrlName
	Variable row
	Variable col
	Variable event	//1=mouse down, 2=up, 3=dbl click, 4=cell select with mouse or keys
//	//5=cell select with shift key, 6=begin edit, 7=end
	if(event==10)
		return 0
	endif
	
	wave/t/z LW1=root:DFSPanel:ListWave1
	make/o/n=0  root:DFSpanel:DefV_1, root:DFSpanel:ZSnsr_1
	struct ForceWave Name1
	DE_Naming#WavetoStruc(LW1[row],Name1)
	controlinfo de_DFS_setvar6
	string ZsnsrName=DE_Naming#StringCreate(Name1.Name,Name1.VNum,S_Value,Name1.SDirec)
	wave w1=$LW1[row]
	wave w2=$ZsnsrName
	duplicate/o w1 root:DFSpanel:Force
	duplicate/o w2 root:DFSpanel:Sep
		
	DE_DFS#FilterUpdate()

	DoUpdate
	return 0
End

Static Function SVP(ctrlName,varNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable varNum
	String varStr
	String varName
	strswitch(ctrlName)
			//	case "de_scholl_setvar0":
			//	DE_PrintAllFolders(root:SchollPanel:FolderList,varstr)
			//break
		case "de_DFS_setvar3":
			DE_DFS#FilterUpdate()
			doupdate
			break
		case "de_DFS_setvar4":
			Cursor/W=DFSPanel#Plot1/P A,Force_fil,varnum
			doupdate
			break
		case "de_DFS_setvar5":
			Cursor/W=DFSPanel#Plot1/P B,Force_fil,varnum
			doupdate
			break
	endswitch
End		

Static Function BP_1(ctrlName) : ButtonControl
	String ctrlName
	wave w0=root:DFSPanel:Force
	wave w1=root:DFSPanel:Force_fil
	wave w2=root:DFSPanel:SmRup
	wave w3=root:DFSPanel:UnSmRup

	wave/t/z LW1=root:DFSPanel:ListWave1
	struct ForceWave Name1
	Controlinfo de_DFS_list1
	variable row=v_value
	DE_Naming#WavetoStruc(LW1[row],Name1)
	variable Direc=cmpstr(Name1.SDirec,"Ext")
	controlinfo de_DFS_setvar4
	variable c1=v_value
	controlinfo de_DFS_setvar5
	variable c2=v_value
	controlinfo de_DFS_setvar7
	variable l1=v_value
	controlinfo de_DFS_setvar8
	variable l2=v_value
	variable/c topline=DE_DFS#MakeSlopes(w1,c1,c2,l1,l2,0)
	variable/c botline=DE_DFS#MakeSlopes(w1,c1,c2,l1,l2,1)
	variable smoothed=DE_DFS#AnotherLastCrossing(w1,c1,c2,topline,botline,Direc)
	variable unsmoothed=DE_DFS#AnotherLastCrossing(w0,c1,c2,topline,botline,Direc)
	w2=w1[smoothed]
	SetScale/P x pnt2x(w1,smoothed),1,"m", w2
	w3=w0[unsmoothed]
	SetScale/P x pnt2x(w0,unsmoothed),1,"m", w3
	ValDisplay de_DFS_valdisp0, value=_NUM: (w1[smoothed]*-1e12)
	ValDisplay de_DFS_valdisp1, value=_NUM:  (w0[unsmoothed]*-1e12)
	ValDisplay de_DFS_valdisp2, value=_NUM:  (imag(topline)*-1e12)

End

Static Function BP_2(ctrlName) : ButtonControl
	String ctrlName
	wave Results
	wave/t/z LW1=root:DFSPanel:ListWave1

	struct ForceWave Name1
	Controlinfo de_DFS_list1
	variable row=v_value
	DE_Naming#WavetoStruc(LW1[row],Name1)
	variable Number=Name1.VNum
	
	if(waveexists(Results)==0)
		make/o/n=(0,4) Results
	endif	
	
	insertpoints/M=0 0,1, Results
	Results[0][0]=Number
	controlinfo de_DFS_valdisp0
	Results[0][1]=v_value
	controlinfo de_DFS_valdisp1
	Results[0][2]=v_value
		controlinfo de_DFS_valdisp2
	Results[0][3]=v_value

end

Static Function BP_3(ctrlName) : ButtonControl
	String ctrlName
	wave Results
	
	
	if(waveexists(Results)==0)
	return -1
	endif	
	
	make/o/n=(0,4) Results

end

Static Function BP_4(ctrlName) : ButtonControl
	String ctrlName
	wave Results
	wave w1=root:DFSPanel:Force_fil
	
	DE_SchollFitting#Scholl(w1,w1,1,1,10,99,TotalNumber=1)
end



Window DFSPanel() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel/N=DFSPanel /W=(10,10,700,700)
	SetWindow DFSPanel, hook(MyHook)=DE_DFS#WindowHook
	newdatafolder/s/o root:DFSPanel
	make/d/o/n=0 Force,Sep,Force_fil,Sep_fil,Fit1,Fit2//,CorForce,CorSep,FitResiduals,FitResidualsX,CorForceFil,CorSepFil,root:SchollPanel:DispRup,root:SchollPanel:DispXRup,root:SchollPanel:ZeroX,root:SchollPanel:ZeroY
	make/d/o/n=1 SmRup=Nan
	make/d/o/n=1 UnSmRup=Nan
	make/o/t/n=0 FolderList,ListWave1,Listwave2

	DE_PrintAllFolders(root:DFSPanel:FolderList,"Image*")

	Button DE_DFS_B1,pos={420,20},size={125,20},proc=DE_DFS#BP_1,title="Do it!"
	Button DE_DFS_B2,pos={420,50},size={125,20},proc=DE_DFS#BP_2,title="Copy"
	Button DE_DFS_B3,pos={420,80},size={125,20},proc=DE_DFS#BP_3,title="Clear"
	Button DE_DFS_B4,pos={420,110},size={125,20},proc=DE_DFS#BP_4,title="Ind. Scholl"


	ListBox de_DFS_list0,pos={10,20},size={150,150},proc=DE_DFS#LBP,listWave=root:DFSPanel:FolderList
	ListBox de_DFS_list0,row= 0,mode= 1,selRow= 0
	ListBox de_DFS_list1,pos={200,20},size={150,150},proc=DE_DFS#LBP_1,listWave=root:DFSPanel:ListWave1
	ListBox de_DFS_list1,row= 0,mode= 1,selRow= 0

	SetVariable de_DFS_setvar0,pos={10,0},size={150,16},proc=de_DFS#SVP,value= _STR:"Image*",title="Search String"
	SetVariable de_DFS_setvar1,pos={165,0},size={200,16},proc=de_DFS#SVP,value= _STR:"Image*Force_Ext",title="Search String 1"
	SetVariable de_DFS_setvar3,pos={50,620},size={100,30},proc=de_dfs#SVP,value= _NUM:51,title="Filtering"
	SetVariable de_DFS_setvar4,pos={200,620},size={200,30},proc=de_dfs#SVP,value=_NUM:100,title="Transition 1"
	SetVariable de_DFS_setvar5,pos={200,650},size={200,30},proc=de_dfs#SVP,value=_NUM:200,title="Transition 2"
	SetVariable de_DFS_setvar6,pos={165,180},size={200,16},proc=de_dfs#SVP,value= _STR:"Sep",title="Suffix Search"
	SetVariable de_DFS_setvar7,pos={400,620},size={200,30},proc=de_dfs#SVP,value=_NUM:200,title="Length 1"
	SetVariable de_DFS_setvar8,pos={400,650},size={200,30},proc=de_dfs#SVP,value=_NUM:200,title="Length 2"
	ValDisplay de_DFS_valdisp0,pos= {50,250},size={150,30}, title="Smoothed", value=_NUM:0
	ValDisplay de_DFS_valdisp1,pos= {250,250},size={150,30}, title="UnSmoothed", value=_NUM:0
	ValDisplay de_DFS_valdisp2,pos= {420,250},size={150,30}, title="Slope", value=_NUM:0

//	CheckBox de_scholl_check0,pos={1150,60},size={40,14},proc=DE_Scholl_CP,value= 0
	Display/W=(20,300,600,600)/HOST=DFSPanel/N=Plot1


	appendtograph/W=DFSPanel#Plot1 root:DFSPanel:Force //vs root:DFSPanel:Sep
	appendtograph/W=DFSPanel#Plot1 root:DFSPanel:Force_fil //vs root:DFSPanel:Sep_fil
	appendtograph/W=DFSPanel#Plot1 root:DFSPanel:Fit1
	appendtograph/W=DFSPanel#Plot1 root:DFSPanel:Fit2
	appendtograph/W=DFSPanel#Plot1 root:DFSPanel:SmRup
	appendtograph/W=DFSPanel#Plot1 root:DFSPanel:UnSmRup

	ModifyGraph/W=DFSPanel#Plot1  rgb(Force)=(63232,45824,45824),rgb(Force_fil)=(58368,6656,7168)
	ModifyGraph rgb(Fit1)=(0,0,0),rgb(Fit2)=(0,0,0)
	ModifyGraph mode(SmRup)=3,marker(SmRup)=19,rgb(SmRup)=(0,0,0),mode(UnSmRup)=3;DelayUpdate
	ModifyGraph marker(UnSmRup)=19,rgb(UnSmRup)=(0,0,65280)

EndMacro


Static Function WindowHook(s)
	STRUCT WMWinHookStruct &s
	GetWindow $s.winName activeSW
	String activeSubwindow = S_value
	if (CmpStr(activeSubwindow,"DFSPanel#Plot1") != 0)
		return 0
	endif
	Variable hookResult = 0
	switch(s.eventCode)
	case 7: // Deactivate
	SetVariable de_DFS_setvar4,value= _NUM:str2num(stringbykey("POINT",CsrInfo(A),":",";"))
	SetVariable de_DFS_setvar5,value= _NUM:str2num(stringbykey("POINT",CsrInfo(B),":",";"))

	break
	//// And so on . . .
	endswitch
	return hookResult // 0 if nothing done, else 1
End


Static function FilterUpdate()
//
	wave/t/z LW1=root:DFSPanel:ListWave1
	Controlinfo de_DFS_list1
	variable row=V_Value
	struct ForceWave Name1
	DE_Naming#WavetoStruc(LW1[row],Name1)
	controlinfo de_DFS_setvar6
	string ZsnsrName=DE_Naming#StringCreate(Name1.Name,Name1.VNum,S_Value,Name1.SDirec)
	wave w1=$LW1[row]
	wave w2=$ZsnsrName
	
	
	
	ControlInfo de_DFS_setvar3
	String DefWaveName=DE_Naming#StringCreate(Name1.Name,Name1.VNum,"DeflCor",Name1.SDirec)
	String ZSnsrWaveName=DE_Naming#StringCreate(Name1.Name,Name1.VNum,"ZSnsr",Name1.SDirec)
	
	wave DefWave=$DefWaveName

		if( waveexists(DefWave)==0)
		Print "Error, something invalid in the Deflection wave naming"
	
	endif
	wave ZWave=$ZSnsrWaveName
	if(waveexists(ZWave)==0)
		Print "Error, something invalid in the ZSnsr wave naming"
	
	endif

	if(waveexists(ZWave)==0||waveexists(DefWave)==0)
		return -1
	endif
	
	controlinfo de_scholl_setvar3
	if(V_Value<5)
	V_Value=5
	endif
	V_value=2*floor(V_value/2)+1
	DE_Filtering#FilterFECs(DefWave,ZWave,"SVG",V_Value)
	wave Force_Smth,Sep_Smth,Defl_Smth
	


	duplicate/o Force_Smth root:DFSPanel:Force_fil
	duplicate/o Sep_Smth root:DFSPanel:Sep_fil
	killwaves Defl_Smth,Sep_Smth,Force_Smth
//
end
//
Static function/C MakeSlopes(ForceWave,c1,c2,startpoints,endpoints,which)
	wave ForceWave
	variable c1, c2,startpoints,endpoints,which
	variable/c Res
	make/o/n=0 StartSlope,EndSlope
	CurveFit/Q/NTHR=0 line  ForceWave[c1-startpoints,c1]/D
	wave w1=$("fit_"+nameofwave(ForceWave))
	wave w_coef
	if(which==0)
	Res=cmplx(W_coef[0],w_coef[1])
	endif
	duplicate/o w1 root:DFSPanel:Fit1
	CurveFit/Q/NTHR=0 line  ForceWave[c2,c2+endpoints]/D
	if(which==1)
	Res=cmplx(W_coef[0],w_coef[1])
	endif
	duplicate/o w1 root:DFSPanel:Fit2
	removefromgraph/Z/W=DFSPanel#Plot1 $("fit_"+nameofwave(ForceWave))
	killwaves w1
	return res


end



Static function AnotherLastCrossing(w1,c1,c2,topline,botline,Direc)
	wave w1
	variable c1,c2,Direc
	variable/C topline,botline
	variable trigger
	variable n
	variable over
	//variable gap=7
	duplicate/o w1 CorOne,SubOne,CorTwo,SubTwo
	SubOne=real(topline)+imag(topline)*pnt2x(w1,p)
	SubTwo=real(botline)+imag(botline)*pnt2x(w1,p)
	FastOp CorOne=CorOne-SubOne
	FastOp CorTwo=CorTwo-SubTwo


	if(Direc==1)
		FindLevels/R=[c1,c2]/Q/Edge=1/P CorOne, 0
		wave w_findLevels
		duplicate/o w_findLevels FirstLevels
		FindLevels/R=[c1,c2]/Q/Edge=1/P CorTwo, 0
		duplicate/o w_findLevels SecondLevels
		killwaves w_findLevels
		n=DE_DFS#CompareCrossing(Firstlevels,SecondLevels)
	
	else
		
		FindLevels/R=[c1,c2]/Q/Edge=2/P CorOne, 0
		wave w_findLevels
		duplicate/o w_findLevels FirstLevels
		FindLevels/R=[c1,c2]/Q/Edge=2/P CorTwo, 0
		duplicate/o w_findLevels SecondLevels
		killwaves w_findLevels
		n=DE_DFS#CompareCrossing(Firstlevels,SecondLevels)
		
	endif
//
	killwaves CorOne,SubOne,CorTwo,SubTwo
	return n
end

static function CompareCrossing(w1,w2)
	wave w1,w2

	if (numpnts(w1)==1)
		return w1[0]
	
	endif
	
	duplicate/o w1 First
	duplicate/o w2 Second
	variable n,j
	variable result 
	For(n=0;n<=numpnts(Second)-1;n+=1)
		If(First[numpnts(First)-1]<Second[n])
			result= First[numpnts(First)-1]
			killwaves First,Second
			return Result
		endif
		If(First[0]<Second[n]&&First[numpnts(First)-1]>Second[n])
			For(j=0;j<=numpnts(First)-1;j+=1)
				If(First[j]<Second[n]&&First[j+1]>Second[n])
					result= First[j]
					killwaves First,Second
					return Result
				endif
			endfor
		endif
	endfor
	print "ERROR"
	killwaves First,Second

	return -1

end


//Legacy Crap
//
//Static function LastCrossing(w1,c1,c2,topline,Direc)
//	wave w1
//	variable c1,c2,Direc
//	variable/C topline
//	variable trigger
//	variable n
//	variable over
//	if(Direc==1)
//	for(n=c2;n>=c1;n-=1)
//		trigger=real(topline)+imag(topline)*pnt2x(w1,n)
//		if(w1[n]<trigger)
//			over=n
//			break
//		endif
//
//	endfor
//	
//	else
//		for(n=c1;n<=c2;n+=1)
//
//		trigger=real(topline)+imag(topline)*pnt2x(w1,n)
//		if(w1[n]<trigger)
//			over=n
//			break
//		endif
//
//	endfor
//	
//	endif
//	
//	
//	return n
//end
//
//
//Static function ModLastCrossing(w1,c1,c2,topline,Direc)
//	wave w1
//	variable c1,c2,Direc
//	variable/C topline
//	variable trigger
//	variable n
//	variable over
//	variable gap=7
//	duplicate/o w1 Test,Sub
//	Sub=real(topline)+imag(topline)*pnt2x(w1,p)
//	FastOp Test=w1-Sub
//	if(Direc==1)
//		FindLevels/Q/Edge=1/P/R=[c1,c2] Test, 0
//		wave w_findLevels
//		duplicate/o w_findLevels RiseLevels
//		FindLevels/Q/Edge=2/P/R=[c1,c2] Test, 0
//		duplicate/o w_findLevels FallLevels
//		killwaves w_findLevels
//		n=FindGap(RiseLevels,FallLevels,gap)
//
//	
//	else
//		
//		FindLevels/Q/Edge=1/P/R=[c1,c2] Test, 0
//		wave w_findLevels
//
//		duplicate/o w_findLevels RiseLevels
//		FindLevels/Q/Edge=2/P/R=[c1,c2] Test, 0
//		duplicate/o w_findLevels FallLevels
//		killwaves w_findLevels		
//
//		n=FindGap(FallLevels,RiseLevels,gap)
//	endif
//
//	//killwaves Test,Sub,w_findLevels
//	return n
//end
//
//static function FindGap(w1,w2,gap)
//	wave w1,w2
//	variable gap
//	if (numpnts(w1)==1)
//		return w1[0]
//	
//	endif
//	
//	duplicate/o w1 First
//	duplicate/o w2 Second
//	variable Final
//	variable n
//	variable q=0
//	if(w1[0]<w2[0])
//	else
//		deletepoints 0,1, Second
//	endif
//	variable last=min(numpnts(First),numpnts(Second))
//
//	for(n=q;n<=last-1;n+=1)
//		if((Second[n]-First[n])>=gap)
//			Final=First[n]
//			break
//		endif
//
//	endfor
//	//killwaves First,Second
//	return Final
//
//
//end