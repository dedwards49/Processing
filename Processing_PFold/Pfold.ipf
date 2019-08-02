#pragma rtGlobals=1		// Use modern global access method.

//------------------------------------------------------------------------------------------------------------
// Calculate Pfold from trajectory, scanning through the whole range between 
//two boundaries representing the folded and unfolded states.
//------------------------------------------------------------------------------------------------------------
Function EmpiricalPfoldScan(inwave,x_f,x_u,stepsize)
Wave inwave
Variable x_f		//unfolded state position
Variable x_u 		//folded state position
Variable stepsize	//resolution for pfold calculation

Wave QCrossing,U_crossing, F_crossing
Variable nsteps
Variable i,j,k,m
Variable nu,nf,nq		//u = no. of crossings at x_u, f at x_f, q at q
Variable nfold, nunfold, n_unfinished
Variable q, current

Variable flag_u,flag_f
Wave Pfold

nsteps = trunc((x_u - x_f)/stepsize)

Make/O/n=0 U_crossing
FindLevels/Q/D=U_crossing inwave, x_u
nu = V_LevelsFound
InsertPoints nu,1,U_crossing
U_crossing[nu] = inf

Make/O/n=0 F_crossing
FindLevels/Q/D=F_crossing inwave, x_f
nf = V_LevelsFound
InsertPoints nf,1,F_crossing
F_crossing[nf] = inf

//Print "nsteps",nsteps,"nu",nu,"nf",nf

Make/O/N=(nsteps+1) Pfold
Pfold[0] = 1		//Don't need to calculate for when it's at the boundaries!
Pfold[nsteps] = 0
m = 0 	//counter for running through q positions

Do		//iterate over each q position
	m += 1		//Start counting at 1, since Pfold[0] = 1 by definition
	q = x_f + m*stepsize
	Make/O/n=0 Qcrossing
	FindLevels/Q/D=QCrossing inwave, q
	nq =V_LevelsFound
//	Print "step",m,"q",q,"nq",nq

	i=0		//index for the q crossings
	j=0		//index for the unfolded state crossings
	k=0		//index for the folded state crossings
	Nfold = 0
	Nunfold = 0
	N_unfinished = 0
	Do
		current = Qcrossing[i]
		flag_u = 0
		Do						//Scroll through unfolded boundary crossings to find first one after q crossing
			If (j == (nu))			//Prevent running off end of wave:
				flag_u = 1		//flag to stop scrolling through unfolded state
			Else				//Find index for first boundary crossings after q:
				If (current > u_crossing[j])	//If crossing at q is after last crossing of unfolded boundary, keep
					j += 1		//scrolling through boundary crossings till get to first one after q crossing
				Else
					flag_u = 1	//flag to stop scrolling through unfolded state
				EndIf
			EndIf
		while (flag_u < 0.5)

		flag_f = 0
		Do						//Scroll through folded boundary crossings to find first one after q crossing
			If (k == (nf))			//Prevent running off end of wave:
				flag_f = 1		//flag to stop scrolling through folded state
			Else	
				If (current > f_crossing[k])	//If crossing at q is after last crossing of folded boundary, keep
					k += 1		//scrolling through boundary crossings till get to first one after q crossing
				Else
					flag_f = 1	//flag to stop scrolling through folded state
				EndIf
			EndIf
		while (flag_f < 0.5)

		//See if F reached before U:
		If ((j == nu) && (k == nf))
			 N_unfinished += 1	//Neither boundary crossed after last q crossing
		Else 					//Find which boundary crossed first
			If (f_crossing[k] < u_crossing[j])
				Nfold += 1
			ElseIf (f_crossing[k] > u_crossing[j])
				Nunfold += 1
			Endif
		EndIf

		//Move to next q crossing
		i += 1
	while (i < nq)
	
	//Calculate pfold for this q value: 
	Pfold[m] = Nfold/(Nfold+Nunfold)
	//Print "q",q,"Nfold",nfold,"Nunfold",nunfold,"pfold",pfold[m]
//Move to next q value
while (m < (nsteps-1))
SetScale/P x x_f,stepsize,"", Pfold
end


//--------------------------------------------------------------------------------------------
//Generate landscape from Pfold(x). Written for boxcar smoothing, 
//so use odd smoothing factor
//--------------------------------------------------------------------------------------------
Function PfoldLandscape(inwave,smthfactor)
Wave inwave
Variable smthfactor

Wave temp,temp_dif,landscape
Variable kT = 2.45	//kT in kJ/mol
Duplicate/O inwave temp,landscape
Smooth/E=3/B smthfactor, temp
Differentiate temp/D=temp_dif
Landscape = kT*ln(-1*temp_dif)
Duplicate/O Landscape $(nameofwave(inwave)+"_landscape")
//Display $(nameofwave(inwave)+"_landscape")

End


//--------------------------------------------------------------------------
//Version of Pfold calculation to use for single q value
//--------------------------------------------------------------------------
Function EmpiricalPfoldSingleQ(inwave,x_f,x_u,q)
Wave inwave
Variable x_f			//unfolded state position
Variable x_u 			//folded state position
Variable q			//extension value for calculation

Wave temp
Wave QCrossing,U_crossing, F_crossing
Variable i,j,k
Variable nu,nf,nq		//u = no. of crossings at x_u, f at x_f, q at q
Variable nfold, nunfold,n_unfinished
//Wave W_FindLevels
Variable current

Variable flag_u,flag_f
Variable Pfold_q

Duplicate/O inwave, temp
Make/O/n=0 U_crossing
FindLevels/Q/D=U_crossing inwave, x_u
nu = V_LevelsFound
InsertPoints nu,1,U_crossing
U_crossing[nu] = inf

Make/O/n=0 F_crossing
FindLevels/Q/D=F_crossing inwave, x_f
nf = V_LevelsFound
InsertPoints nf,1,F_crossing
F_crossing[nf] = inf

Make/O/n=0 Qcrossing
FindLevels/Q/D=QCrossing inwave, q
nq =V_LevelsFound
//Print "nu",nu,"nf",nf,"nq",nq

i=0		//index for the q crossings
j=0		//index for the unfolded state crossings
k=0		//index for the folded state crossings
Nfold = 0
Nunfold = 0
N_unfinished = 0
Do							//Scroll through each q crossing too see if folds or unfolds first
	current = Qcrossing[i]
	flag_u = 0
	Do						//Scroll through unfolded boundary crossings to find first one after q crossing
		If (j == (nu))			//Prevent running off end of wave:
			flag_u = 1		//flag to stop scrolling through unfolded state
		Else				//Find index for first boundary crossings after q:
			If (current > u_crossing[j])	//If crossing at q is after last crossing of unfolded boundary, keep
				j += 1		//scrolling through boundary crossings till get to first one after q crossing
			Else
				flag_u = 1	//flag to stop scrolling through unfolded state
			EndIf
		EndIf
	while (flag_u < 0.5)

	flag_f = 0
	Do						//Scroll through folded boundary crossings to find first one after q crossing
		If (k == (nf))			//Prevent running off end of wave:
			flag_f = 1		//flag to stop scrolling through folded state
		Else	
			If (current > f_crossing[k])	//If crossing at q is after last crossing of folded boundary, keep
				k += 1		//scrolling through boundary crossings till get to first one after q crossing
			Else
				flag_f = 1	//flag to stop scrolling through folded state
			EndIf
		EndIf
	while (flag_f < 0.5)

	//See if F reached before U:
	If ((j == (nu)) && (k == (nf)))
		 N_unfinished += 1	//Neither boundary crossed after last q crossing
	Else 					//Find which boundary crossed first
		If (f_crossing[k] < u_crossing[j])
			Nfold += 1
		ElseIf (f_crossing[k] > u_crossing[j])
			Nunfold += 1
		Endif
	EndIf
//	Print "q crossing",qcrossing[i],"next u crossing",u_crossing[j], "next f crossing",f_crossing[k],"Nfold",Nfold,"Nunfold",Nunfold

	//Move to next q crossing
	i += 1
while (i < nq)
	
//Calculate pfold for this q value: 
Pfold_q = Nfold/(Nfold+Nunfold)
//Return Nfold/(Nfold+Nunfold)
Print "q",q, "nu",nu,"nf",nf,"nq",nq, "Nfold",nfold,"Nunfold",nunfold,"Unclear",N_unfinished,"Pfold",pfold_q
end

