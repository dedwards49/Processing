#pragma rtGlobals=1		// Use modern global access method.

#pragma ModuleName = DE_LoadIBW


function DE_LoadIBW()
	NewPath/o/q Pa
	string fname,wname,wname_defl_ret,wname_defl_ext,wname_zsnsr_ret,wname_zsnsr_ext
	variable l,m,n,t1,t2
	do


		fName= IndexedFile(Pa,n,".ibw")
		if (strlen(fName) == 0)
			Break
		else
		endif

		if (t1!=t2)
			print "File Naming Error"
			Break
		else
		endif
	
		loadwave/o/A/Q/P=Pa fname
		wname=fname[0,strlen(fname)-5]
		wname_defl_ret=wname+"Defl_Ret"
		wname_defl_ext=wname+"Defl_Ext"
		wname_zsnsr_ret=wname+"Zsnsr_Ret"
		wname_zsnsr_ext=wname+"Zsnsr_Ext"
		wave w1=$wname

		//make/o/n=(dimsize(w1,0)) $wname_defl_ret,$wname_defl_ext,$wname_zsnsr_ret,$wname_zsnsr_ext

		strswitch (stringbykey("Direction",note(w1),":","\r"))
			case "Inf,1,0,-1,0,":
				duplicate/o/R=[0,str2num(StringFromList(1, stringbykey("Indexes",note(w1),":","\r"),","))][1]  w1,$wname_defl_ext
				duplicate/o/R=[0,str2num(StringFromList(1, stringbykey("Indexes",note(w1),":","\r"),","))][2]  w1,$wname_zsnsr_ext

				duplicate/o/R=[str2num(StringFromList(2, stringbykey("Indexes",note(w1),":","\r"),",")),str2num(StringFromList(3, stringbykey("Indexes",note(w1),":","\r"),","))][1]  w1,$wname_defl_ret
				duplicate/o/R=[str2num(StringFromList(2, stringbykey("Indexes",note(w1),":","\r"),",")),str2num(StringFromList(3, stringbykey("Indexes",note(w1),":","\r"),","))][2]  w1,$wname_zsnsr_ret
				Redimension/N=-1 $wname_defl_ext
				Redimension/N=-1 $wname_zsnsr_ext
				Redimension/N=-1 $wname_defl_ret

				Redimension/N=-1 $wname_zsnsr_ret

				break

			case " NaN,1,0,-1":
				duplicate/o/R=[0,str2num(StringFromList(1, stringbykey("Indexes",note(w1),":","\r"),","))][1]  w1,$wname_defl_ext
				duplicate/o/R=[0,str2num(StringFromList(1, stringbykey("Indexes",note(w1),":","\r"),","))][2]  w1,$wname_zsnsr_ext

				duplicate/o/R=[str2num(StringFromList(2, stringbykey("Indexes",note(w1),":","\r"),",")),str2num(StringFromList(3, stringbykey("Indexes",note(w1),":","\r"),","))][1]  w1,$wname_defl_ret
				duplicate/o/R=[str2num(StringFromList(2, stringbykey("Indexes",note(w1),":","\r"),",")),str2num(StringFromList(3, stringbykey("Indexes",note(w1),":","\r"),","))][2]  w1,$wname_zsnsr_ret
				Redimension/N=-1 $wname_defl_ext
				Redimension/N=-1 $wname_zsnsr_ext
				Redimension/N=-1 $wname_defl_ret

				Redimension/N=-1 $wname_zsnsr_ret
				break
	
			case  "0,-1,0,":
		
				duplicate/o/R=[0,str2num(StringFromList(1, stringbykey("Indexes",note(w1),":","\r"),","))][1]  w1,$wname_defl_ret
				duplicate/o/R=[0,str2num(StringFromList(1, stringbykey("Indexes",note(w1),":","\r"),","))][2]  w1,$wname_zsnsr_ret

				Redimension/N=-1 $wname_defl_ret

				Redimension/N=-1 $wname_zsnsr_ret


				break
			default:
				print stringbykey("Direction",note(w1),":","\r")
				print n
				return-1
		

	
		endswitch
		killwaves w1
	
	if(mod(n,10)==0)
	print n
	endif
		n+=1
	While(n>-1)

end

function/S DE_SingleLoadIBW(Path,number)
	string Path
	variable number
	NewPath/o/q Pa, Path
	string fname,wname,wname_defl_ret,wname_defl_ext,wname_zsnsr_ret,wname_zsnsr_ext
	variable l,m,n,t1,t2


		fName= IndexedFile(Pa,number,".ibw")
		if (strlen(fName) == 0)
			Print "Not Found"
			return "-1"
		else
		endif

	
		loadwave/o/A/Q/P=Pa fname
		wname=fname[0,strlen(fname)-5]
		wname_defl_ret=wname+"Defl_Ret"
		wname_defl_ext=wname+"Defl_Ext"
		wname_zsnsr_ret=wname+"Zsnsr_Ret"
		wname_zsnsr_ext=wname+"Zsnsr_Ext"
		wave w1=$wname


		strswitch (stringbykey("Direction",note(w1),":","\r"))
			case "Inf,1,0,-1,0,":
				duplicate/o/R=[0,str2num(StringFromList(1, stringbykey("Indexes",note(w1),":","\r"),","))][1]  w1,$wname_defl_ext
				duplicate/o/R=[0,str2num(StringFromList(1, stringbykey("Indexes",note(w1),":","\r"),","))][2]  w1,$wname_zsnsr_ext

				duplicate/o/R=[str2num(StringFromList(2, stringbykey("Indexes",note(w1),":","\r"),",")),str2num(StringFromList(3, stringbykey("Indexes",note(w1),":","\r"),","))][1]  w1,$wname_defl_ret
				duplicate/o/R=[str2num(StringFromList(2, stringbykey("Indexes",note(w1),":","\r"),",")),str2num(StringFromList(3, stringbykey("Indexes",note(w1),":","\r"),","))][2]  w1,$wname_zsnsr_ret
				Redimension/N=-1 $wname_defl_ext
				Redimension/N=-1 $wname_zsnsr_ext
				Redimension/N=-1 $wname_defl_ret

				Redimension/N=-1 $wname_zsnsr_ret

				break

			case " NaN,1,0,-1":
				duplicate/o/R=[0,str2num(StringFromList(1, stringbykey("Indexes",note(w1),":","\r"),","))][1]  w1,$wname_defl_ext
				duplicate/o/R=[0,str2num(StringFromList(1, stringbykey("Indexes",note(w1),":","\r"),","))][2]  w1,$wname_zsnsr_ext

				duplicate/o/R=[str2num(StringFromList(2, stringbykey("Indexes",note(w1),":","\r"),",")),str2num(StringFromList(3, stringbykey("Indexes",note(w1),":","\r"),","))][1]  w1,$wname_defl_ret
				duplicate/o/R=[str2num(StringFromList(2, stringbykey("Indexes",note(w1),":","\r"),",")),str2num(StringFromList(3, stringbykey("Indexes",note(w1),":","\r"),","))][2]  w1,$wname_zsnsr_ret
				Redimension/N=-1 $wname_defl_ext
				Redimension/N=-1 $wname_zsnsr_ext
				Redimension/N=-1 $wname_defl_ret

				Redimension/N=-1 $wname_zsnsr_ret
				break
	
			case  "0,-1,0,":
		
				duplicate/o/R=[0,str2num(StringFromList(1, stringbykey("Indexes",note(w1),":","\r"),","))][1]  w1,$wname_defl_ret
				duplicate/o/R=[0,str2num(StringFromList(1, stringbykey("Indexes",note(w1),":","\r"),","))][2]  w1,$wname_zsnsr_ret

				Redimension/N=-1 $wname_defl_ret

				Redimension/N=-1 $wname_zsnsr_ret


				break
			default:
				print stringbykey("Direction",note(w1),":","\r")
				print n
				return "-1"
		

	
		endswitch
		string nameofw1=nameofwave(w1)
		killwaves w1
		return nameofw1
end