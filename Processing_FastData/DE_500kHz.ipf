#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma modulename=DE_500kHz
#include "EquilibriumJumpData"
#include "DE_Filtering"

Static Function Process500kHz(FastDeflV,FastZSnsrV,SlowDeflection,SlowZSnsr)

	wave FastDeflV,FastZSnsrV,SlowDeflection,SlowZSnsr
	variable InVols=str2num(stringbykey("Invols",note(FastDeflV),":","\r"))
	variable k=str2num(stringbykey("Spring Constant",note(FastDeflV),":","\r"))
	variable LVDTSens=str2num(stringbykey("ZLVDTSens",note(SlowDeflection),":","\r"))
	duplicate/free FastDeflV TestDefl,TestZSnsr
	FastOp TestDefl=(InVols)*TestDefl
	FastOp TestZSnsr=(LVDTSens)*FastZSnsrV
	make/free/n=0 SmoothedZSNSR
	DE_Filtering#TVD1D_denoise(TestZSnsr,20e-9,SmoothedZSNSR)
	DE_EquilJump#GenerateRegions(SmoothedZSNSR)
	wave Midpoints
	duplicate/o MidPoints MidPointsFast
	killwaves Midpoints
	DE_EquilJump#GenerateRegions(SlowZSnsr)
	wave Midpoints
	duplicate/free MidPoints MidPointsSlow	
	killwaves MidPoints
	deletepoints (numpnts(MidPointsSlow)-1),1, MidPointsSlow,MidPointsFast
	deletepoints 0,1, MidPointsSlow,MidPointsFast
	make/D/free/n=2 W_Coefs
	W_Coefs={0,1}
	CurveFit/W=2/Q/H="01"/NTHR=0 line kwCWave=W_Coefs  MidPointsSlow /X=MidPointsFast
	wave W_sigma
	killwaves W_SIGMA,MidPointsSlow,MidPointsFast
	deletepoints 0, x2pnt(TestDefl,-w_coefs[0]), TestDefl,TestZSNSR
	duplicate/free SlowDeflection TestSlowDefl
	Resample/DOWN=10 TestDefl, 
	Deletepoints (numpnts(TestDefl)),1e8,TestSlowDefl
	W_Coefs={0,1}
	CurveFit/W=2/Q/H="01"/NTHR=0 line kwCWave=W_Coefs  TestSlowDefl /X=TestDefl

	duplicate/o FastDeflV, $ReplaceString("Ret", nameofwave(SlowDeflection), "Fast")
	wave w1=$ReplaceString("Ret", nameofwave(SlowDeflection), "Fast")
	FastOp w1=(InVols)*w1+(w_coefs[0])
	
	duplicate/o FastZSnsrV $ReplaceString("Ret", nameofwave(SlowZSnsr), "Fast")
	wave w2=$ReplaceString("Ret", nameofwave(SlowZSnsr), "Fast")
	FastOp w2=(LVDTSens)*w2
	wave W_SIGMA
	killwaves W_Sigma
end

