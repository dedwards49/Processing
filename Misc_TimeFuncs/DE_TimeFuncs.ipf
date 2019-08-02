#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma modulename=DE_TimeFuncs
Static Function ReturnTimeFromWave(InputWave)
	wave InputWave

	String TimeString= stringbykey("Time",note(InputWave),":","\r")
	variable Hour=str2num(stringfromlist(0,TimeString,":"))
	variable Minute=str2num( stringfromlist(1,TimeString,":"))
	variable Seconds=str2num(Stringfromlist(0,stringfromlist(2,TimeString,":")," "))
	string AMPM=Stringfromlist(1,stringfromlist(2,TimeString,":")," ")
	if(stringmatch(AMPM,"AM")==1)
	elseif(stringmatch(AMPM,"PM")==1)
		Hour+=12
	endif
	variable Result=3600*Hour+60*Minute+Seconds

	return Result

end
Static Function ReturnDateFromWave(InputWave)
	wave InputWave

	String DateString= stringbykey("Date",note(InputWave),":","\r")
	String Month= Stringfromlist(1,stringfromlist(1,DateString,",")," ")
	variable MonthVar=MonthtoNum(Month)
	variable Day= str2num(Stringfromlist(2,stringfromlist(1,DateString,",")," "))
	variable Year= str2num(stringfromlist(2,DateString,", "))
	return date2secs(Year,MonthVar,Day)

end



Static Function MonthtoNum(StringMonth)
	String StringMonth
	variable MonthVar
	strswitch(StringMonth)
		case "Jan":
			MonthVar=1
			break
		case "Feb":
			MonthVar=2
			break
		case "Mar":
			MonthVar=3
			break
		case "Apr":
			MonthVar=4
			break
		case "May":
			MonthVar=5
			break
		case "Jun":
			MonthVar=6
			break
		case "Jul":
			MonthVar=7
			break
		case "Aug":
			MonthVar=8
			break
		case "Sep":
			MonthVar=9
			break
		case "Oct":
			MonthVar=10
		case "Nov":
			MonthVar=11
		case "Dec":
			MonthVar=12
			break
			
	endswitch
	
	return MonthVar
end