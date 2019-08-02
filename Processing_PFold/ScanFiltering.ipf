#pragma rtGlobals=3		// Use modern global access method and strict wave access.

#pragma modulename=LandscapeFiltering
#include "pfold"

static function GeneralizedIterativePFold(BaseWave,folded,unfolded, edge,StartFilter,EndFilter,FilterSteps,PFoldStep,PFoldSmoothing,Dq)
	wave basewave
	variable folded,unfolded, edge,StartFilter,EndFilter,FilterSteps,PFoldStep,PFoldSmoothing,Dq

	variable n
	variable k=0
	string BName,PName,Notes
	variable BoltzmannPeak,BoltzmannValley,PFoldPeak
	NewDataFolder/O root:FilterScan		// Make sure this exists.
	duplicate/o BaseWave root:FilterScan:DataWave
	
	SetDataFolder root:FilterScan	
	wave DataWave
	make/o/n=0 Smoothing,PFold_Height,PFold_Distance,Boltzmann_Height,Boltzmann_Distance,CWell,CBar,Ka
	
	for(n=StartFilter;n<EndFilter;n+=FilterSteps)
		duplicate/o DataWave smoothed
		smooth/S=2 (n), smoothed
		
		if(mean(smoothed,0,100)<=1e-6)
		smoothed*=1e9
		endif


		EmpiricalPfoldScan(smoothed,(folded-edge),(unfolded+edge),PFoldStep)  
		wave pfold
		PfoldLandscape(Pfold,PFoldSmoothing)
		wave pfold_landscape
		pfold_landscape/=2.45
		insertpoints 0,1, Smoothing,PFold_Height,PFold_Distance,Boltzmann_Height,Boltzmann_Distance,CWell,CBar,Ka
		Smoothing[0]=n
		
		
		wavestats/q/r=(folded,unfolded) pfold_landscape
		PFoldPeak=v_maxloc
		wavestats/q/r=(folded,PFoldPeak) pfold_landscape
		PFold_Distance[0]= PFoldPeak-v_minloc
		PFold_Height[0]= v_max-pfold_landscape(v_minloc)
		Notes=AddListItem(num2str(x2pnt(pfold_landscape,PFoldPeak)),"")
		Notes=AddListItem(num2str(x2pnt(pfold_landscape,v_minloc)),Notes)
		Notes=AddListItem(num2str(x2pnt(pfold_landscape,unfolded)),Notes)



		//Here is a different way to calculate parameters based on where we input the folded and unfolded states are.
//		wavestats/q/r=(folded,unfolded) pfold_landscape
//
//		PFold_Distance[0]= folded-v_maxloc
//		PFold_Height[0]= v_max-pfold_landscape(folded)
//
//		Notes=AddListItem(num2str(x2pnt(pfold_landscape,v_maxloc)),"")
//		Notes=AddListItem(num2str(x2pnt(pfold_landscape,folded)),Notes)
//		Notes=AddListItem(num2str(x2pnt(pfold_landscape,unfolded)),Notes)
////		
		
		
		note/k pfold_landscape Notes
		PName="PFold_"+num2str(n)
		duplicate/o pfold_landscape $Pname
		
		
		Make/N=40/O smoothed_Hist;DelayUpdate
		Histogram/C/B=1 smoothed,smoothed_Hist
		duplicate/o smoothed_Hist Boltzmann
		Boltzmann=-ln(smoothed_Hist)
		wavestats/q/r=(folded,unfolded) Boltzmann
		Boltzmann_Height[0]=v_max-Boltzmann(folded)
		Boltzmann_Distance[0]= v_maxloc-folded


		Make/o/D/N=3/O W_coef
		W_coef[0] = {-2.8,Boltzmann(v_maxloc),v_maxloc}
		FuncFit/Q/NTHR=0 Stiffness W_coef  Boltzmann(v_maxloc-2*edge,v_maxloc+2*edge) /D 
		wave w_coef
		CBar[0]=w_coef[0]*4.1
		Notes=AddListItem(num2str(x2pnt(Boltzmann,v_maxloc)),"")

		wavestats/Q/r=(folded-2*edge,folded+2*edge) Boltzmann

		Make/o/D/N=3/O W_coef
		W_coef[0] = {1,Boltzmann(v_minloc),v_minloc}
		Notes=AddListItem(num2str(x2pnt(Boltzmann,v_minloc)),Notes)

		FuncFit/Q/NTHR=0 Stiffness W_coef  Boltzmann(v_minloc-2*edge,v_minloc+2*edge) /D 
		Cwell[0]=w_coef[0]*4.1
	
		note/k Boltzmann Notes
		BName="Boltz_"+num2str(n)
		duplicate/o Boltzmann $Bname
	
		
		Ka[0]=1/2/pi*sqrt(abs(Cwell[0]*Cbar[0]))/4.11*Dq*exp(-Boltzmann_Height[0])

	endfor
	SetDataFolder root:
end


static function IterativePFold(BaseWave,folded,unfolded, edge)
	wave basewave
	variable folded,unfolded, edge
	variable n=5
	variable k=0
	string BName,PName,Notes
	variable BoltzmannPeak,BoltzmannValley
	make/o/n=0 Smoothing,PFold_Height,PFold_Distance,Boltzmann_Height,Boltzmann_Distance,CWell,CBar,Ka
	
	for(n=5;n<202;n+=2)
		duplicate/o basewave smoothed
		smooth/S=2 (n), smoothed
		smoothed*=1e9


		EmpiricalPfoldScan(smoothed,(folded-edge),(unfolded+edge),.01)  
		wave pfold
		PfoldLandscape(Pfold,5)
		wave pfold_landscape
		pfold_landscape/=2.45
		insertpoints 0,1, Smoothing,PFold_Height,PFold_Distance,Boltzmann_Height,Boltzmann_Distance,CWell,CBar,Ka
		Smoothing[0]=n
		wavestats/q/r=(folded,unfolded) pfold_landscape
		PFold_Distance[0]= v_maxloc-folded
		PFold_Height[0]= v_max-pfold_landscape(folded)
		Notes=AddListItem(num2str(x2pnt(pfold_landscape,v_maxloc)),"")
		Notes=AddListItem(num2str(x2pnt(pfold_landscape,folded)),Notes)
		Notes=AddListItem(num2str(x2pnt(pfold_landscape,unfolded)),Notes)

		note/k pfold_landscape Notes
		PName="PFold_"+num2str(n)
		duplicate/o pfold_landscape $Pname
		
		
		Make/N=40/O smoothed_Hist;DelayUpdate
		Histogram/C/B=1 smoothed,smoothed_Hist
		duplicate/o smoothed_Hist Boltzmann
		Boltzmann=-ln(smoothed_Hist)
		wavestats/q/r=(folded,unfolded) Boltzmann
		Boltzmann_Height[0]=v_max-Boltzmann(folded)
		Boltzmann_Distance[0]= v_maxloc-folded


		Make/o/D/N=3/O W_coef
		W_coef[0] = {-2.8,Boltzmann(v_maxloc),v_maxloc}
		FuncFit/Q/NTHR=0 Stiffness W_coef  Boltzmann(v_maxloc-2*edge,v_maxloc+2*edge) /D 
		wave w_coef
		CBar[0]=w_coef[0]
		Notes=AddListItem(num2str(x2pnt(Boltzmann,v_maxloc)),"")

		wavestats/Q/r=(folded-2*edge,folded+2*edge) Boltzmann

		Make/o/D/N=3/O W_coef
		W_coef[0] = {1,Boltzmann(v_minloc),v_minloc}
		Notes=AddListItem(num2str(x2pnt(Boltzmann,v_minloc)),Notes)

		FuncFit/Q/NTHR=0 Stiffness W_coef  Boltzmann(v_minloc-2*edge,v_minloc+2*edge) /D 
		Cwell[0]=w_coef[0]
	
		note/k Boltzmann Notes
		BName="Boltz_"+num2str(n)
		duplicate/o Boltzmann $Bname


		Ka[0]=1/2/pi*sqrt(abs(Cwell[0]*Cbar[0]))/4.11*1e4*exp(-Boltzmann_Height[0])
	endfor
end

Function Stiffness(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = y0+k*(x-x0)^2
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 3
	//CurveFitDialog/ w[0] = k
	//CurveFitDialog/ w[1] = y0
	//CurveFitDialog/ w[2] = x0

	return w[1]+w[0]*(x-w[2])^2
End
