#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma ModuleName = DE_Naming


Static function/S StringCreate(Name,Num,Value,Direc,[Modifier])
	String Name
	variable num
	string Value, Direc
		string Modifier

	string numhold
	sprintf numhold "%04.4G",num
	string ReturnString
	if(cmpstr(Direc,"All")==0)
		if(cmpstr(Value,"All")==0)
			ReturnString= Name+numhold

		else
			ReturnString= Name+numhold+Value
		endif
	else
		ReturnString= Name+numhold+Value+"_"+Direc
	endif
	if(ParamisDefault(Modifier))
	else
	ReturnString+="_"+Modifier
	endif
	return ReturnString
end

Structure ForceWave
variable VNum //Curve Number
string SNum 
string Name
string SMeas
variable VMeas
DFREF Folder //Number associated with the measurement type.
//"DeflV"=1
//"Defl"=2
//"Force"=3
// "Raw"=4
// "Zsnsr"=5
// "Sep":=6
// "RawV":=7
// "DeflCor":=8
//
//Not recognized=-inf
string SDirec
variable VDirec//Number associated with the direction of the measurement
//"Ext"=1
//"Ret"=2
//"Away"=3
// "Towd"=4
//"All"=5
//Not recognized=-inf
	
EndStructure

Static function WavetoStruc(InWave,NameStruc,[FOlder])
	string InWave
	struct ForceWave &NameStruc
	DFREF Folder
	if(ParamisDefault(FOlder))
		wave WaveInWave=$InWave
	else
		wave WaveInWave=Folder:$InWave

	endif
	string Suff= stringbykey("BaseSuffix",note(WaveInWave),":","\r")
	variable len=strlen(Suff)
	Suff=Suff[len-5,len]
	variable n
	Suff=num2str(str2num(suff))
	for(n=strlen(suff);n<4;n+=1)
		Suff="0"+Suff

	endfor
	variable SuffPos=strsearch (Inwave,Suff,strlen(Inwave),1)
	string NameStr=Inwave[0,SuffPos-1]
	string SuffStr=Inwave[SuffPos,SuffPos+strlen(Suff)-1]
	string EndInfo=Inwave[SuffPos+strlen(Suff),strlen(Inwave)-1]
	variable divider=DE_FindLast(EndInfo,"_")
	
		string Meas=EndInfo[0,divider-1]

	string Dirr=EndInfo[divider+1,strlen(endinfo)-1]
	if(divider==-1)
			 Meas=EndInfo[0,99]
			 Dirr=""
	
	endif
	NameStruc.SNum=Suff
	NameStruc.VNum=str2num(NameStruc.SNum)
	NameStruc. Name=NameStr
	NameStruc.SMeas=CheckMeasS(Meas)

	NameStruc.VMeas=CheckMeasV(Meas)
	
	NameStruc.SDirec=CheckDirS(Dirr)
	NameStruc.VDirec=CheckDirV(Dirr)
end


Static function TrueWavetoStruc(WaveInWave,NameStruc)
	wave WaveInWave
	struct ForceWave &NameStruc
	string Inwave=nameofwave(WaveInWave)
	string Suff= stringbykey("BaseSuffix",note(WaveInWave),":","\r")
	variable len=strlen(Suff)
	Suff=Suff[len-5,len]
	variable n
	Suff=num2str(str2num(suff))
	for(n=strlen(suff);n<4;n+=1)
		Suff="0"+Suff

	endfor
	variable SuffPos=strsearch (Inwave,Suff,strlen(Inwave),1)
	string NameStr=Inwave[0,SuffPos-1]
	string SuffStr=Inwave[SuffPos,SuffPos+strlen(Suff)-1]
	string EndInfo=Inwave[SuffPos+strlen(Suff),strlen(Inwave)-1]
	variable divider=DE_FindLast(EndInfo,"_")
	
		string Meas=EndInfo[0,divider-1]

	string Dirr=EndInfo[divider+1,strlen(endinfo)-1]
	if(divider==-1)
			 Meas=EndInfo[0,99]
			 Dirr=""
	
	endif
	NameStruc.SNum=Suff
	NameStruc.VNum=str2num(NameStruc.SNum)
	NameStruc. Name=NameStr
	NameStruc.SMeas=CheckMeasS(Meas)

	NameStruc.VMeas=CheckMeasV(Meas)
	
	NameStruc.SDirec=CheckDirS(Dirr)
	NameStruc.VDirec=CheckDirV(Dirr)
	NameStruc.Folder=GetWavesDataFolderDFR(WaveInWave)

end

Static function RealWavetoStruc(InputWave,NameStruc)
	Wave InputWave
	struct ForceWave &NameStruc
		string InWave=nameofwave(InputWave)

	string Suff= stringbykey("BaseSuffix",note(InputWave),":","\r")

	variable len=strlen(Suff)

	Suff=Suff[len-5,len]

	variable n
	Suff=num2str(str2num(suff))
	for(n=strlen(suff);n<4;n+=1)
		Suff="0"+Suff

	endfor
	variable SuffPos=strsearch (Inwave,Suff,0)

	string NameStr=Inwave[0,SuffPos-1]
	string SuffStr=Inwave[SuffPos,SuffPos+strlen(Suff)-1]
	string EndInfo=Inwave[SuffPos+strlen(Suff),strlen(Inwave)-1]
	variable divider=DE_FindLast(EndInfo,"_")
	
		string Meas=EndInfo[0,divider-1]

	string Dirr=EndInfo[divider+1,strlen(endinfo)-1]
	if(divider==-1)
			 Meas=EndInfo[0,99]
			 Dirr=""
	
	endif
	NameStruc.SNum=Suff
	NameStruc.VNum=str2num(NameStruc.SNum)
	NameStruc. Name=NameStr
	NameStruc.SMeas=CheckMeasS(Meas)

	NameStruc.VMeas=CheckMeasV(Meas)
	
	NameStruc.SDirec=CheckDirS(Dirr)
	NameStruc.VDirec=CheckDirV(Dirr)
end

function CheckDirV(Str)
	String Str
	variable res

	strswitch(Str)
		case "Ext":
			res=1
			break
		case "Ret":
			res=2
			break
		case "Away":
			res=3
			break
		case "Towd":
			res= 4
			break
		case "All":
			res= 5
			break
		case "":
			res= 5
			break		
		case "Equil":
			res= 6
			break
					case "Fast":
			res= 7
			break
		default:
			res=-inf
		

	endswitch
	
	return res

end

function/S CheckDirS(Str)
	String Str
	string res

	strswitch(Str)
		case "Ext":
			res= "Ext"
			break
		case "Ret":
			res="Ret"
			break
		case "Away":
			res="Away"
			break
		case "Towd":
			res= "Towd"
			break
		case "All":
			res= "All"
			break
		case "Fast":
			res="Fast"
			break
		case "":
			res= "All"
			break
					case "Equil":
			res= "Equil"
			break
								case "Fast":
			res= "Fast"
			break
		default:
			res="Unk"
		

	endswitch
	return res


end

function CheckMeasV(Str)
	String Str
	variable res

	strswitch(Str)
		case "DeflV":
			res=1
			break
		case "Defl":
			res=2
			break
		case "Force":
			res=3
			break
		case "Raw":
			res= 4
			break		
		case "Zsnsr":
			res= 5
			break
		case "Sep":
			res= 6
			break
		
		case "RawV":
			res= 7
			break
	
		case "DeflCor":
			res= 8
			break
		case "":
			res= 9
			break
		default:
			res=-inf
	endswitch

	
	return res

end

function/S CheckMeasS(Str)
	String Str
	string res

	strswitch(Str)
		case "DeflV":
			res="DeflV"
			break
		case "Defl":
			res="Defl"
			break
		case "Force":
			res="Force"
			break
		case "Raw":
			res= "Raw"
			break		
		case "Zsnsr":
			res= "Zsnsr"
			break
		case "Sep":
			res= "Sep"
			break
			
		case "RawV":
			res= "RawV"
			break
		
		case "DeflCor":
			res= "DeflCor"
			break
		case "":
			res= "All"
			break
			
		default:
			res="Unk"
	endswitch


	
	return res

end


Function DE_FindLast(Str,Str2find)
	String Str, Str2find
	
	Return(StrSearch(Str,Str2Find,Strlen(Str)-1,3))
	
	
	Variable ind, cnt
	//get the last occurance of the str2find in String str
	//New Version of Find_last, it is 2X Faster, and does the same thing
	cnt = strlen(Str)-strlen(Str2find)
	ind = -1
	do
		ind = StrSearch(str,str2find,cnt)
		cnt -= 1
	while ((ind == -1) && (cnt > -1))
	return(ind)
End //FindLast	
