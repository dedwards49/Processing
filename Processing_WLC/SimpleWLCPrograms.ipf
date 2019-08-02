#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma ModuleName=DE_WLC
//This is the original WLC function which is not as accurate. Lp is in meters, while z and Lo should be as well, but in principle can be in anything so long as their
//units agree.


function oldWLC(z,Lp,L0,T)
	variable z,Lp,L0,T

	variable pre,su,res

	if (z>=L0)
		res=NaN
	else

		pre=1.3806488*1e-23*T/Lp
		su=(1/4/(1-z/L0)^2-1/4+z/L0)
		res=pre*su
	endif

	return res
end

Static Function FindLC_Woodside	(F,z,Lp,K,T)
	variable F,z,Lp,T,K
	variable start= stopmstimer(-2)
	variable pre,su,res
	make/o/n=1000 xs,Resid
		xs=1*p/999
		pre=1.3806488*1e-23*T/Lp
		Resid=pre*(1/4*(1-xs)^(-2)-1/4+xs)-F
	FindLevel/Q Resid,0
	
	return z/xs[v_levelx]

end
//
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//
// WLC simply calculates the WLC force at an extension z. Lp is the persistance length (in m), T
//is temperature (in K), while L0 is the contour length. L0 and z can be in any units, so long as they
//are the same units, but working in m is probably smartest. This is the WLC refined from 
//"Estimating the Persistence Length of a Worm-Like Chain Molecule from Force-Extension Measurements"
// C. Bouchiat, M.D. Wang, J.-F. Allemand, T. Strick, S.M. Block, V. Croquette
//Biophysical Journal Volume 76, Issue 1, January 1999, Pages 409–413

function WLC(z,Lp,L0,T)
	variable z,Lp,L0,T

	variable pre,su,res
	variable a2,a3,a4,a5,a6,a7
	a2=-.5164228
	a3=-2.737418
	a4=16.07497
	a5=-38.87607
	a6=39.49949
	a7=-14.17718

	if (z>=L0)//Does not true to compute anything if z is larger than the contour length.
		res=NaN
	else

		pre=1.3806488*1e-23*T/Lp
		su=1/(4*(1-z/L0)^2)-1/4+z/L0+a2*(z/L0)^2+a3*(z/L0)^3+a4*(z/L0)^4+a5*(z/L0)^5+a6*(z/L0)^6+a7*(z/L0)^7
		res=pre*su
	endif

	if(abs(res)>2e-9) //This eliminates values that are above 2nN, which is reasonable for our work.
		res=nan
	endif
	return res
end

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Function qmWLC(Force,Pers,T,LC)
	variable Force,Pers,T,LC
	variable y1=27.4
	variable y2=109.8

	variable f=Force*Pers/(1.3806488*1e-23*T)
	variable b=exp((900/f)^(1/4))


	Variable Lcorr=LC/2/y2*(sqrt(4*Force*y2+y1^2)-y1+2*y2)
	// LCorr
	return Lcorr*(4/3-4/3*1/sqrt(f+1)-10*b/sqrt(f)/(b-1)^2+f^1.62/(3.55+3.8*f^2.2))

end

function WLCSlope(z,Lp,L0,T)
    variable z,Lp,L0,T

    variable pre,su,res
    variable a2,a3,a4,a5,a6,a7
    a2=-.5164228
    a3=-2.737418
    a4=16.07497
    a5=-38.87607
    a6=39.49949
    a7=-14.17718

    if (z>=L0)//Does not true to compute anything if z is larger than the contour length.
        res=NaN
    else

        pre=1.3806488*1e-23*T/Lp/L0
	  su=1/2*(1-z/L0)^(-3)+1+2*a2*(z/L0)^1+3*a3*(z/L0)^2+4*a4*(z/L0)^3+5*a5*(z/L0)^4+6*a6*(z/L0)^5+7*a7*(z/L0)^6
        res=pre*su
    endif

//    if(abs(res)>2e-9) //This eliminates values that are above 2nN, which is reasonable for our work.
//        res=nan
//    endif
    return res
end 

Static Function ReturnExtentionatForce(F,Lp,L0,T)
	variable F,Lp,L0,T
	make/free/n=1000 Ext,Force
	Ext=L0/999*x
	Force=WLC(Ext,Lp,L0,T) 
	FindLevel/Q Force F
	return Ext[V_levelx]

end

Static Function ContourTransform(Force,x,Pers,T)
	variable Force,x,Pers,T
	if (Force<2e-12)
	 Force=2e-12
		endif
	
	variable f=Force*Pers/(1.3806488*1e-23*T)
	variable b=exp((900/f)^(1/4))
	variable y1=27.4
	variable y2=109.8

	variable res= x/(4/3-4/3*1/sqrt(f+1)-10*b/sqrt(f)/(b-1)^2+f^1.62/(3.55+3.8*f^2.2))

	if(res<1e-9)
	 f=2e-12*Pers/(1.3806488*1e-23*T)
	 b=exp((900/f)^(1/4))
	 y1=27.4
	 y2=109.8

	 res= x/(4/3-4/3*1/sqrt(f+1)-10*b/sqrt(f)/(b-1)^2+f^1.62/(3.55+3.8*f^2.2))
	endif
	return res

end

Function ContourTransformqm(Force,x,Pers,T)
	variable Force,x,Pers,T
	variable f=Force*Pers/(1.3806488*1e-23*T)
	variable b=exp((900/f)^(1/4))


	return x/(4/3-4/3*1/sqrt(f+1)-10*b/sqrt(f)/(b-1)^2+f^1.62/(3.55+3.8*f^2.2))


end



//Simple WLC fitting with an x and y offset (xoff and off respectively).
Function WLC_FIT(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = WLC(x-xoff,Lp,L0,T)-off
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 5
	//CurveFitDialog/ w[0] = Lp
	//CurveFitDialog/ w[1] = L0
	//CurveFitDialog/ w[2] = T
	//CurveFitDialog/ w[3] = off
	//CurveFitDialog/ w[4] = xoff

	return WLC(x-w[4],w[0],w[1],w[2])-w[3]
End



//Simple WLC fitting with an x and y offset (xoff and off respectively).
Function WLCOLD_FIT(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = WLC(x-xoff,Lp,L0,T)-off
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 5
	//CurveFitDialog/ w[0] = Lp
	//CurveFitDialog/ w[1] = L0
	//CurveFitDialog/ w[2] = T
	//CurveFitDialog/ w[3] = off
	//CurveFitDialog/ w[4] = xoff

	return WLC(x-w[4],w[0],w[1],w[2])-w[3]
End

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//Creates a series of WLC curves for a homo-polyprotein based on the input extension wave Ext (should be in m). n is the number of subunits, while L0 is the contour length of the shortest unit, and DeltaL0 is the 
//step size in LC. Both of these should be in m as well, but could be in another unit if Ext is also given in the matching units. Lp is the persistence length and must be in m. T is Temperature in K.
function WLCSeries(Ext,n,L0Init,DeltaL0,Lp,T)
	wave Ext
	variable n,L0Init,DeltaL0,Lp,T
	variable split,i
	if (n==0)
		return 1
	endif
	wavestats/q Ext
	//make/o/n=(v_npnts*n+n) WLC_Force,WLC_Ext,WLC_Color
	make/o/n=(v_npnts*n+n) WLC_Force,WLC_Ext
	WLC_Ext=Ext[mod(p,v_npnts)]
	WLC_Force=WLC(WLC_Ext[p],Lp,L0init+floor(p/v_npnts)*DeltaL0,T)

	wavestats/q WLC_Force
	for(i=0;i<(v_npnts+V_numNans);i+=1)
		if(i<=3000)
			if(abs(WLC_Force[i])>.12e-9)//Discoards forces above 500pN
				WLC_FOrce[i]=NaN
			endif
		elseif(i<=4000)
			if(abs(WLC_Force[i])>.18e-9)//Discoards forces above 500pN
				WLC_FOrce[i]=NaN
			endif
		else
			if(abs(WLC_Force[i])>.22e-9)//Discoards forces above 500pN
				WLC_FOrce[i]=NaN
			endif
		endif
		
	endfor

	//for(i=0;i<n;i+=1)
	//WLC_Color[v_npnts*i,(i+1)*v_npnts-1]=i-5
	//WLC_Force[(i+1)*v_npnts]=NaN
	//Endfor

end

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//As above but for the polyprotein construct of Nug2 and Cal with another foldering intermediate...well can be the combination of any 2 polyproteins but this assumes that polyprotein 1 (all of them) unfold first.
function WLC_Arb(Parms)
	wave Parms


	variable num=dimsize(Parms,0)
	variable v_npnts=1000
	make/o/n=(1000*num) WLC_Force,WLC_Ext,WLC_Color
	variable i
	for(i=0;i<=num-1;i+=1)
		WLC_Ext[i*v_npnts,(i+1)*v_npnts-1]=500e-9*(p-i*1000)/999

		WLC_Force[i*v_npnts,(i+1)*v_npnts-1]=WLC(WLC_Ext[p],Parms[i][0],Parms[i][2],Parms[i][1])
		WLC_Color[i*v_npnts,(i+1)*v_npnts-1]=0
	endfor

	
	wavestats/q WLC_Force
	for(i=0;i<(v_npnts+V_numNans);i+=1)
		if(i<=2000)
			if(abs(WLC_Force[i])>.1e-9)//Discoards forces above 500pN
				WLC_FOrce[i]=NaN
			endif
		elseif(i<=3000)
			if(abs(WLC_Force[i])>.15e-9)//Discoards forces above 500pN
				WLC_FOrce[i]=NaN
			endif
					elseif(i<=5000)
			if(abs(WLC_Force[i])>.3e-9)//Discoards forces above 500pN
				WLC_FOrce[i]=NaN
			endif
		else
			if(abs(WLC_Force[i])>.5e-9)//Discoards forces above 500pN
				WLC_FOrce[i]=NaN
			endif
		endif
		
	endfor



end

Static function WLC_ArbWave(Parms,Fout,XOut)
	wave Parms,Fout,XOut


	variable num=dimsize(Parms,0)
	variable v_npnts=1000
	make/free/n=(1000*num) WLC_Force,WLC_Ext,WLC_Color
	variable i
	for(i=0;i<=num-1;i+=1)
		WLC_Ext[i*v_npnts,(i+1)*v_npnts-1]=500e-9*(p-i*1000)/999

		WLC_Force[i*v_npnts,(i+1)*v_npnts-1]=WLC(WLC_Ext[p],Parms[i][0],Parms[i][2],Parms[i][1])
		WLC_Color[i*v_npnts,(i+1)*v_npnts-1]=0
	endfor

	
	wavestats/q WLC_Force
	for(i=0;i<(v_npnts+V_numNans);i+=1)
		if(i<=1000)
			if(abs(WLC_Force[i])>.1e-9)//Discoards forces above 500pN
				WLC_FOrce[i]=NaN
			endif
		elseif(i<=5000)
			if(abs(WLC_Force[i])>.2e-9)//Discoards forces above 500pN
				WLC_FOrce[i]=NaN
			endif
		else
			if(abs(WLC_Force[i])>.35e-9)//Discoards forces above 500pN
				WLC_FOrce[i]=NaN
			endif
		endif
		
	endfor

	duplicate/o WLC_Force Fout
	duplicate/o WLC_Ext XOut
end


function TestMakeParms()

	make/o/n=(3,3) WLCParms
	WLCParms[][0]=.4e-9
	WLCParms[][1]=298
	variable last=225e-9
		WLCParms[2][2]=last
		WLCParms[1][2]=WLCParms[2][2]-16.6e-9
		WLCParms[0][2]=WLCParms[1][2]-14e-9
		//WLCParms[0][2]=WLCParms[1][2]-25e-9


	





	WLC_Arb(WLCParms)
	


end

function TestMakeParmsNuG2(last)
variable last
	make/o/n=(3,3) WLCParms
	WLCParms[][0]=.4e-9
	WLCParms[][1]=298
		WLCParms[2][2]=last
		WLCParms[1][2]=WLCParms[2][2]-17e-9
		WLCParms[0][2]=WLCParms[1][2]-17e-9
	//	WLCParms[0][2]=WLCParms[1][2]-23e-9


	





	WLC_Arb(WLCParms)
	wave WLC_Force,WLC_Ext
	duplicate/o WLC_Force NuG2_WLC_Force
	duplicate/o WLC_Ext NUG2_WLC_Ext



end

function TestMakeParmsNewalpha(alpha,marker,last)
variable alpha,marker,last
	make/o/n=(4,3) WLCParms
	WLCParms[][0]=.4e-9
	WLCParms[][1]=298
	WLCParms[3][2]=last
	WLCParms[2][2]=WLCParms[3][2]-marker
	WLCParms[1][2]=WLCParms[2][2]-marker
	WLCParms[0][2]=WLCParms[1][2]-alpha

	WLC_Arb(WLCParms)
	wave WLC_Force,WLC_Ext
	duplicate/o WLC_Force alpha_WLC_Force
	duplicate/o WLC_Ext alpha_WLC_Ext

end

function TestMakeParmsELC(last,Nug2Step,ELCstep)
variable last,Nug2Step,ELCstep
	make/o/n=(6,3) WLCParms
	WLCParms[][0]=.4e-9
	WLCParms[][1]=298
	WLCParms[5][2]=last
	WLCParms[4][2]=last-Nug2Step
	WLCParms[3][2]=last-2*Nug2Step
	WLCParms[2][2]=last-3*Nug2Step
	WLCParms[1][2]=last-4*Nug2Step
	WLCParms[0][2]=last-4*Nug2Step-ELCStep

	WLC_Arb(WLCParms)
	wave WLC_Force,WLC_Ext
	duplicate/o WLC_Force ELC_WLC_Force
	duplicate/o WLC_Ext ELC_WLC_Ext



end

function TestMakeParmsNickMCD(lastMCD,Nug2Step)
variable lastMCD,Nug2Step
	make/o/n=(11,3) WLCParms
	variable step1=70e-9
	variable step2=26.2e-9
	variable step3=37.92e-9
	variable step4=14.12e-9
	variable step5=61e-9
	variable step6=74e-9
	WLCParms[][0]=.4e-9
	WLCParms[][1]=298

	WLCParms[6][2]=lastMCD
	WLCParms[7][2]=WLCParms[6][2]+Nug2Step
	WLCParms[8][2]=WLCParms[7][2]+Nug2Step
	WLCParms[9][2]=WLCParms[8][2]+Nug2Step
	WLCParms[10][2]=WLCParms[9][2]+Nug2Step
	
	WLCParms[5][2]=	WLCParms[6][2]-step6
	WLCParms[4][2]=	WLCParms[5][2]-step5
	WLCParms[3][2]=	WLCParms[4][2]-step4
	WLCParms[2][2]=	WLCParms[3][2]-step3
	WLCParms[1][2]=	WLCParms[2][2]-step2

	WLCParms[0][2]=	WLCParms[1][2]-step1

	WLC_Arb(WLCParms)
	wave WLC_Force,WLC_Ext
	duplicate/o WLC_Force MCD_WLC_Force
	duplicate/o WLC_Ext MCD_WLC_Ext



end

function TestMakeParmsMCD(last,Nug2Step,ELCstep)
variable last,Nug2Step,ELCstep
	make/o/n=(6,3) WLCParms
	WLCParms[][0]=.4e-9
	WLCParms[][1]=298
	WLCParms[5][2]=last
	WLCParms[4][2]=last-Nug2Step
	WLCParms[3][2]=last-2*Nug2Step
	WLCParms[2][2]=last-3*Nug2Step
	WLCParms[1][2]=last-4*Nug2Step
	WLCParms[0][2]=last-4*Nug2Step-ELCStep

	WLC_Arb(WLCParms)
	wave WLC_Force,WLC_Ext
	duplicate/o WLC_Force MCD_WLC_Force
	duplicate/o WLC_Ext MCD_WLC_Ext



end
function TestMakeParmsRLCBD(last,Nug2Step,RLCstep1,RLCStep2)
variable last,Nug2Step,RLCstep1,RLCStep2
	make/o/n=(7,3) WLCParms
	WLCParms[][0]=.4e-9
	WLCParms[][1]=298
	WLCParms[6][2]=last
	WLCParms[5][2]=last-Nug2Step
	WLCParms[4][2]=last-4*Nug2Step
	WLCParms[3][2]=last-3*Nug2Step
	WLCParms[2][2]=last-4*Nug2Step
	WLCParms[1][2]=last-4*Nug2Step-RLCstep1
	WLCParms[0][2]=last-4*Nug2Step-RLCstep1-RLCStep2

	WLC_Arb(WLCParms)
	wave WLC_Force,WLC_Ext
	duplicate/o WLC_Force RLC_WLC_Force
	duplicate/o WLC_Ext RLC_WLC_Ext
	killwaves WLC_Force,WLC_Ext


end

function TestMakeParmsRLCBDCoh(last,ddFln1,ddFln2,RLCstep1,RLCStep2)
variable last,ddFln1,ddFln2,RLCstep1,RLCStep2
	make/o/n=(5,3) WLCParms
	WLCParms[][0]=.4e-9
	WLCParms[][1]=298
	WLCParms[4][2]=last
	WLCParms[3][2]=last-ddFln2
	WLCParms[2][2]=last-ddFln2-ddFln1
	WLCParms[1][2]=last-ddFln2-ddFln1-RLCstep1
	WLCParms[0][2]=last-ddFln2-ddFln1-RLCstep1-RLCStep2

	WLC_Arb(WLCParms)
	wave WLC_Force,WLC_Ext
	duplicate/o WLC_Force RLC_WLC_Force
	duplicate/o WLC_Ext RLC_WLC_Ext
	killwaves WLC_Force,WLC_Ext


end
function TestMakeParmsRLCBDddFLN4(last,ddFLN4Step1,ddFLN4Step2,RLCstep1,RLCStep2)
variable last,ddFLN4Step1,ddFLN4Step2,RLCstep1,RLCStep2
	make/o/n=(5,3) WLCParms
	WLCParms[][0]=.4e-9
	WLCParms[][1]=298
	WLCParms[4][2]=last
	WLCParms[3][2]=last-ddFLN4Step1
	WLCParms[2][2]=last-ddFLN4Step1-ddFLN4Step2
	WLCParms[1][2]=last-ddFLN4Step1-ddFLN4Step2-RLCstep1
	WLCParms[0][2]=last-ddFLN4Step1-ddFLN4Step2-RLCstep1-RLCstep2

	WLC_Arb(WLCParms)
	wave WLC_Force,WLC_Ext
	duplicate/o WLC_Force RLC_WLC_Force
	duplicate/o WLC_Ext RLC_WLC_Ext
	killwaves WLC_Force,WLC_Ext


end
function TestMakeParms_alpha()

	make/o/n=(4,3) WLCParms
	WLCParms[][0]=.4e-9
	WLCParms[][1]=298
	variable Split1=16.8e-9
	variable Split2=25e-9
	variable FirstNUg2=167e-9
//	WLCParms[5][2]=FirstNUg2+Split1*4
//	WLCParms[4][2]=FirstNUg2+Split1*3
	WLCParms[3][2]=FirstNUg2+Split1*2
	WLCParms[2][2]=FirstNUg2+Split1*1
	WLCParms[1][2]=FirstNUg2
	WLCParms[0][2]=FirstNUg2-Split2
//	WLCParms[5][2]=133.28e-9
//	WLCParms[4][2]=117.55e-9
//	WLCParms[3][2]=100.18e-9
//	WLCParms[2][2]=83.7e-9
//	WLCParms[1][2]=66.05e-9
//	WLCParms[0][2]=41.2e-9

	


	//Parms[][2]-=1e-8



	WLC_Arb(WLCParms)



end

function TestMakeParms_Titin(Final,Nug2,Titin)
	variable Final,Nug2,Titin
	make/o/n=(4,3) WLCParms
	WLCParms[][0]=.4e-9
	WLCParms[][1]=298

	WLCParms[3][2]=Final
	WLCParms[2][2]=Final-Titin
	WLCParms[1][2]=Final-Nug2-Titin
	WLCParms[0][2]=Final-Nug2*2-Titin

	WLC_Arb(WLCParms)
	wave WLC_Force,WLC_Ext
	duplicate/o WLC_Force Titin_WLC_Force
	duplicate/o WLC_Ext Titin_WLC_Ext
	
	make/o/n=(3,3) WLCParms
	WLCParms[][0]=.4e-9
	WLCParms[][1]=298

	WLCParms[2][2]=Final
	WLCParms[1][2]=Final-Nug2
	WLCParms[0][2]=Final-Nug2*2

	WLC_Arb(WLCParms)
	duplicate/o WLC_Force MT_WLC_Force
	duplicate/o WLC_Ext MT_WLC_Ext
	
	
end

function TestMakeParms_Titin2(Final,Nug2,Titin)
	variable Final,Nug2,Titin
	make/o/n=(4,3) WLCParms
	WLCParms[][0]=.4e-9
	WLCParms[][1]=298

	WLCParms[3][2]=Final
	WLCParms[2][2]=Final-Nug2
	WLCParms[1][2]=Final-Nug2*2
	WLCParms[0][2]=Final-Nug2*2-Titin

	WLC_Arb(WLCParms)
	wave WLC_Force,WLC_Ext
	duplicate/o WLC_Force Titin2_WLC_Force
	duplicate/o WLC_Ext Titin2_WLC_Ext

	
	
end

function WLCSeries_Cal(Ext,n1,n2,L0Init,DeltaL0_n1,DeltaL0_n2,Lp,T)
	wave Ext
	variable n1,n2,L0Init,DeltaL0_n1,DeltaL0_n2,Lp,T
	variable split
	if (n1==0||n2==0)
		return 1
	endif
	wavestats/q Ext
	make/o/n=(v_npnts*(n1+n2+1)+(n1+n2)) WLC_Force,WLC_Ext,WLC_Color
	variable i
	for(i=0;i<=n1;i+=1)
		WLC_Ext[i*v_npnts,(i+1)*v_npnts-1]=Ext[p-i*v_npnts]
		WLC_Force[i*v_npnts,(i+1)*v_npnts-1]=WLC(WLC_Ext[p],Lp,L0init+i*DeltaL0_n1,T)
		WLC_Color[i*v_npnts,(i+1)*v_npnts-1]=0
	endfor

	for(i=n1+1;i<=(n2+n1);i+=1)
		WLC_Ext[i*v_npnts,(i+1)*v_npnts-1]=Ext[p-(i)*v_npnts]
		WLC_Force[i*v_npnts,(i+1)*v_npnts-1]=WLC(WLC_Ext[p],Lp,L0init+n1*DeltaL0_n1+(i-n1)*DeltaL0_n2,T)
		WLC_Color[i*v_npnts,(i+1)*v_npnts-1]=6
	endfor
	wavestats/q WLC_Force
	for(i=0;i<(v_npnts+V_numNans);i+=1)
		if(i<=3000)
			if(abs(WLC_Force[i])>.05e-9)//Discoards forces above 500pN
				WLC_FOrce[i]=NaN
			endif
		elseif(i<=7000)
			if(abs(WLC_Force[i])>.15e-9)//Discoards forces above 500pN
				WLC_FOrce[i]=NaN
			endif
		else
			if(abs(WLC_Force[i])>.15e-9)//Discoards forces above 500pN
				WLC_FOrce[i]=NaN
			endif
		endif
		
	endfor


end

function FindLC(F,Lp,ex,T)
	variable F,Lp,ex,T
	make/o/n=100000 LoTest,FTest,RTest
	LoTest=1e-6*x/99999
	FTest=WLC(ex,Lp,LoTest,T)
	FastOp RTest=FTest+(-F)

	FindLevel/Q RTest,0
	variable ret= LoTest[V_LevelX]

	killwaves RTest,LoTest,FTest
	return ret
end


//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//FindX() returns the extension value that yields a particular force based on the WLC model. This is a really crude program based on zero-finding so it's quite slow.
//F is the force (Newtons), Lp is the persistence length (m), Lo is the contour length (m) and T is the temperature (K)


function FindX(F,Lp,Lo,T)
	variable F,Lp,Lo,T
	make/o/n=10000 ExtTest,FTest,RTest
	ExtTest=Lo*x/9999
	FTest=WLC(ExtTest,Lp,Lo,T)
	FastOp RTest=FTest+(-F)

	FindLevel/Q RTest,0
	variable ret= ExtTest[V_LevelX]

	killwaves RTest,Exttest,FTest
	return ret
end

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//GenerateFX() is a crude program to generate the extension curve for a series of forces of two tether chains of different stiffnesses. Each has a contour length (lo1/lo2) and a persistence length (lp1/lp2) which
//should be given in meters. The temperature (T, in K) is shared. The script generates the extensions corresponding to a linear list of applied forces (F, in n) from Fstart in steps of FStep for Fnum. The resulting
//extension for each chain is then added together. This thing is slow as the calculation of the extension is based on FindX(), which is a zero-finder.

function GenerateFX(Fstart,FStep,FNum,Lo1,Lp1,Lo2,Lp2,T)

	variable  Fstart,FStep,FNum,Lo1,Lp1,Lo2,Lp2,T
	variable timerRefNum = startMSTimer
	make/o/n=(FNum) Forces,Ext1,Ext2,TotExt
	Forces=Fstart+FStep*x
	Ext1=FindX(Forces,Lp1,Lo1,T)

	Ext2=FindX(Forces,Lp2,Lo2,T)
	TotExt=Ext1+Ext2
	killwaves Ext1,Ext2
	variable microSeconds = stopMSTimer(timerRefNum)
	Print microSeconds/10000, "microseconds per iteration"

end

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Extensible WLC model with K. This is done as a zero-finder to handle solving the analytic problem.

function WLC_Exten(z,Lp,L0,T,K)
	
	variable z,Lp,L0,T,K

	variable pre,res
	variable a2,a3,a4,a5,a6,a7
	a2=-.5164228
	a3=-2.737418
	a4=16.07497
	a5=-38.87607
	a6=39.49949
	a7=-14.17718
	make/o/n=1000 ForceTest
	ForceTest=50e-10*x/999
	duplicate/o Forcetest LTest,SU,TotalRHS,ZeroFun
	LTest=z/L0-ForceTest/K
	if (z>=1.1*L0)
		res=NaN
	else

		//pre=4.11434*(T/298)/(Lp)
		pre=1.3806488*1e-23*T/Lp
		su=1/(4*(1-ltest)^2)-1/4+ltest+a2*ltest^2+a3*ltest^3+a4*ltest^4+a5*ltest^5+a6*ltest^6+a7*ltest^7
		TotalRHS=pre*su
		ZeroFun=ForceTest-TotalRHS
		//su=1/(4*(1-z/L0)^2)-1/4+z/L0+a2*(z/L0)^2+a3*(z/L0)^3+a4*(z/L0)^4+a5*(z/L0)^5+a6*(z/L0)^6+a7*(z/L0)^7
		//res=pre*su
		Findlevel/q zerofun, 0
		res=Forcetest[v_levelx]
	endif 
	killwaves ZeroFun,Forcetest,LTest,SU,TotalRHS
	if(abs(res)>1e-9)
		res=nan
	endif
	return res


end


Function FitEXWLC(w,x,y) : FitFunc
Wave w
Variable x
Variable y
//CurveFitDialog/
//CurveFitDialog/ Coefficients 6
//CurveFitDialog/ w[0] = T
//CurveFitDialog/ w[1] = L0
//CurveFitDialog/ w[2] = K
//CurveFitDialog/ w[3] = Lp
//CurveFitDialog/ w[4] = x0
//CurveFitDialog/ w[5] = y0

 		variable pre=1.3806488*1e-23*w[0]/w[3]
 		variable z=(x-w[4])/w[1]
 		variable F=(y-w[5] )/w[2]
return pre*(1/4*(1-z+F)^(-2)-1/4+z-F)-F*w[2]
End
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//A simle fitting program for the extensible WLC model. This is quite slow.

Function WLC_EXT_FIT(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = WLC(x-xoff,Lp,L0,T)-off
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 6
	//CurveFitDialog/ w[0] = Lp
	//CurveFitDialog/ w[1] = L0
	//CurveFitDialog/ w[2] = T
	//CurveFitDialog/ w[3] = off
	//CurveFitDialog/ w[4] = xoff
	//CurveFitDialog/ w[5] =K


	return WLC_Exten(x-w[4],w[0],w[1],w[2],w[5])-w[3]
End

function WLC_Energy(zstart,z,Lp,L0,T)
	variable z,zstart,Lp,L0,T
	make/free/n=1000 Ext,Force
	Ext=zstart+x*(z-zstart)/999
	Force=WLC(Ext,Lp,L0,T)
	Integrate/METH=1 Force/X=Ext/D=Force_int
	variable Total=Force_int[999]*1e9*1e12*2.479/4.114
	killwaves Force_int

	return Total //Returns in kJ/M

end

function WLC_Energy_Force(zmin,ForceTop,Lp,L0,T)
	variable zmin,ForceTop,Lp,L0,T
	make/free/n=1000 Ext,Force
	Ext=zmin+x*(L0-zmin)/999
	Force=WLC(Ext,Lp,L0,T)
	FindLevel/q Force, ForceTop
	Integrate/METH=1 Force/X=Ext/D=Force_int
	variable Total=Force_int[v_levelx]*1e9*1e12/4.11

	killwaves Force_int
	return Total //Returns in kJ/M

end


function DE_TiltLandscape(originallandscape,forcetilt,zmin,Cstart,Cend)
	wave originallandscape
	variable forcetilt,zmin,Cstart,Cend

	string newname=nameofwave(originallandscape)+"_tilt"
	duplicate/o originallandscape WLCTilt,SpringTilt
	duplicate/o originallandscape $newname
	wave outwave=$newname
	SpringTilt=1/2*((1020.7-36.6*x))^2/36.6*2.48/4.11
	make/o/N=(numpnts(WLCTilt)) Lcs
	Lcs=Cstart+(Cend-Cstart)/(numpnts(WLCTilt)-1)*p
	print Lcs[0]
	WLCTilt[]=WLC_Energy_Force(zmin,forcetilt,.4e-9,Lcs[p],300)
	outwave=originallandscape+WLCTilt+SpringTilt
end

Function  low_liva_c(x,kT,a,L)
	variable x,kT,a,L
	return ((x*3.0*kT)/(a*L))

end

Function  mid_liva_c(x,kT,l_pers,L)
	variable x,kT,l_pers,L
	return (kT/(4.0* l_pers))*(1.0/(1.0-(x/L) )^2.0 )
end

Function  lhigh_liva_c(x,kT,b,L)
	variable x,kT,b,L
	return (kT/(2.0*b)) * ( 1.0 / (1.0-(x/L)) )
end

 
Static Function/S  CalculateDeltaLCs(wave2)
	wave wave2
	wave wave1=root:WLCFIT:W_coef
	print nameofwave(wave2)
	variable length=numpnts(wave1)
	variable n=2
	for(n=3;n<length;n+=1)
		print (wave1[n]-wave1[n-1])/1e-9
	
	
	endfor
	print stringbykey("DE_Xoff",note(wave2),":","\r")
	print stringbykey("DE_Yoff",note(wave2),":","\r")

end

Static Function PetrosyanWLC(z,Lp,L0,T)

	variable z,Lp,L0,T

variable res,pre,x
	x=z/L0
	

	if (z>=L0)//Does not true to compute anything if z is larger than the contour length.
		res=NaN
	else
		pre=1.3806488*1e-23*T/Lp

		res=(z-.8*x^(2.15)+0.25/(1-x)^2-.25)*pre
	endif

	if(abs(res)>2e-9) //This eliminates values that are above 2nN, which is reasonable for our work.
		res=nan
	endif
	return res
end

Static Function InvPetrosyanWLC(F,Lp,L0,T)

	variable F,Lp,L0,T

	variable res,pre,fred
			pre=1.3806488*1e-23*T/Lp

	fred=F/pre
	

//	if (z>=L0)//Does not true to compute anything if z is larger than the contour length.
	//	res=NaN
//	else

		res=4/3-4/(3*sqrt(fred+1))-10*exp((900/fred)^(1/4))/(sqrt(fred)*(exp((900/fred)^(1/4))-1)^2)+(fred^1.62)/(3.55+3.8*(fred^2.2))
	res*=L0
//	endif

	//if(abs(res)>2e-9) //This eliminates values that are above 2nN, which is reasonable for our work.
	//	res=nan
//	endif
	return res
end

