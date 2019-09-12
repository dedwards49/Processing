
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

Function DE_RUN(waven)
wave waven
wave  wave1,wave2,wave3,wave4,wave5,wave6,wave7,wave8
InsertPoints (dimsize(wave1,0)), 1, wave1,wave2,wave3,wave4,wave5,wave6,wave7,wave8
wave1[(dimsize(wave1,0))-1]=DE_Centroid(waven,"X")
wave2[(dimsize(wave1,0))-1]=DE_Centroid(waven,"Y")

wave3[(dimsize(wave1,0))-1]=DE_D4sigma(waven,"X")
wave4[(dimsize(wave1,0))-1]=DE_D4sigma(waven,"Y")
K6=0
CurveFit/Q/H="0000001"/NTHR=0 Gauss2D  waven 
wave w_coef
wave5[(dimsize(wave1,0))-1]=w_coef[2]
wave6[(dimsize(wave1,0))-1]=w_coef[4]

wave7[(dimsize(wave1,0))-1]=4*w_coef[3]
wave8[(dimsize(wave1,0))-1]=4*W_COEF[5]
end

function DE_D4sigma(W1,axis)
wave w1
string axis

variable result

duplicate/free w1 Data,DX

variable bkgd_point = adjust_bkgd(Data, axis)
variable Cent=DE_Centroid(w1,axis)

if(cmpstr(axis,"x")==0)
DX=Data*(x-Cent)^2
Integrate/dim=0 Data /D=Data_Int
Integrate/dim=0 DX/D=DX_Int
make/free/n=(dimsize(Data,1)) Data1,DX1
Data1=Data_Int[dimsize(Data,0)-1][p]
DX1=Dx_Int[dimsize(Data,0)-1][p]
result=4*sqrt(area(DX1)/area(Data1))
wave/z vertical_boundaries,horizontal_boundaries
killwaves/z DX_int,Data_int,vertical_boundaries,horizontal_boundaries
return result



elseif(cmpstr(axis,"y")==0)

DX=Data*(y-Cent)^2
Integrate/dim=1 Data /D=Data_Int
Integrate/dim=1 DX/D=DX_Int
make/free/n=(dimsize(Data,0)) Data1,DX1
Data1=Data_Int[p][dimsize(Data,1)-1]
DX1=Dx_Int[p][dimsize(Data,1)-1]
result=4*sqrt(area(DX1)/area(Data1))
killwaves DX_int,Data_int,vertical_boundaries,horizontal_boundaries
return result


else
return 0
endif



end

function PullOne(W1,axis)
wave W1
variable axis
make/o/n=(dimsize(W1,0),dimsize(W1,1)) NewWave
NewWave=W1[p][q][axis]
SetScale/P x 0,dimdelta(W1,0),"", NewWave
SetScale/P y 0,dimdelta(W1,1),"", NewWave

end

	
function DE_Centroid(W1,axis)
wave w1
string axis
variable result
duplicate/free w1 Data,DX
//Data[][]=w1[p][q][2]
Data[][]=w1[p][q]
variable bkgd_point = adjust_bkgd(Data, axis)
if(cmpstr(axis,"x")==0)
DX=Data*x
Integrate/dim=0 Data /D=Data_Int
Integrate/dim=0 DX/D=DX_Int
make/free/n=(dimsize(Data,1)) Data1,DX1

Data1=Data_Int[dimsize(Data,0)-1][p]
DX1=Dx_Int[dimsize(Data,0)-1][p]
result=area(DX1)/area(Data1)
wave/z vertical_boundaries,horizontal_boundaries
//killwaves/z DX_int,Data_int,vertical_boundaries,horizontal_boundaries
return result



elseif(cmpstr(axis,"y")==0)

DX=Data*y
Integrate/dim=1 Data /D=Data_Int
Integrate/dim=1 DX/D=DX_Int
make/free/n=(dimsize(Data,0)) Data1,DX1
Data1=Data_Int[p][dimsize(Data,1)-1]
DX1=Dx_Int[p][dimsize(Data,1)-1]
result=area(DX1)/area(Data1)
wave/z vertical_boundaries,horizontal_boundaries

killwaves/z DX_int,Data_int,vertical_boundaries,horizontal_boundaries
return result


else
return 0
endif
end

function/c get_centroids(data, x_ax, y_ax)
	wave data, x_ax, y_ax
	duplicate/free data, temp_data
	adjust_data(temp_data)
	dfref df = getwavesdatafolderdfr(data)
	string expr = "(.*)_(.*)", wave_id, rest_of_wavename
	splitstring/e=(expr) nameofwave(data), rest_of_wavename, wave_id
	variable x0 = get_centroid(temp_data, x_ax, "x")
	variable y0 = get_centroid(temp_data, y_ax, "y")
	variable/c centroids = cmplx(x0, y0)
	variable/g df:$(wave_id + "_x0") = x0
	variable/g df:$(wave_id + "_y0") = y0
	return centroids
end

function adjust_data(data)
	wave data
	variable n = dimsize(data, 0), m = dimsize(data, 1)
	// normalise data
	variable wmax = wavemax(data)
	data /= wmax
	// invert data if necessary
	if (mean(get_outer_ring(data)) > data[n/2][m/2])
		print "inverting", nameofwave(data)
		data = 1 - data
	endif
end

function/wave get_outer_ring(data)
	wave data
	variable n = dimsize(data, 0), m = dimsize(data, 1)
	make/free/n=(m) first_row = data[p][0]
	make/free/n=(m) last_row = data[p][m-1]
	make/free/n=(n-2) first_column = data[0][p+1]
	make/free/n=(n-2) last_column = data[n-1][p+1]
	make/free/n=0 outer_ring
	concatenate/np/kill/o {first_row, last_row, first_column, last_column}, outer_ring
	return outer_ring
end

function get_centroid(data, ax, axis)
	wave data, ax
	string axis
	// create a free wave duplicate of the data to manipulate
	make/free/n=(dimsize(data,0), dimsize(data,1)) temp_data
	temp_data = data
	waveclear data
	variable bkgd_point = adjust_bkgd(temp_data, axis)
	if (bkgd_point != 0)
		return bkgd_point
	endif
	variable centroid = 0, num = 0, denom = 0
	variable i, j
	if (stringmatch(axis, "x"))
		for (i = 0; i < dimsize(temp_data, 0); i += 1)
			for (j = 0; j < dimsize(temp_data, 1); j += 1)
				num += ax[i] * temp_data[i][j]
				denom += temp_data[i][j]
			endfor
		endfor

	elseif (stringmatch(axis, "y"))
		for (i = 0; i < dimsize(temp_data, 0); i += 1)
			for (j = 0; j < dimsize(temp_data, 1); j += 1)
				num += ax[j] * temp_data[i][j]
				denom += temp_data[i][j]
			endfor
		endfor
	endif
	centroid = num/denom
	return centroid
end

function adjust_bkgd(data, axis)
	wave data
	string axis
	variable n = dimsize(data, 0), m = dimsize(data, 1)
	// subtract/adjust background
	if (stringmatch(axis, "x"))
		wave boundaries = get_vertical_boundaries(data)
	elseif (stringmatch(axis, "y"))
		wave boundaries = get_horizontal_boundaries(data)
	endif
	variable bkgd = wavemax(boundaries)
	variable bkgd_point
	if (bkgd < wavemax(data))
		data -= bkgd
		variable i, j
		for (i = 0; i < dimsize(data, 0); i += 1)
			for (j = 0; j < dimsize(data, 1); j += 1)
				if (data[i][j] < 0)
					data[i][j] = 0
				endif
			endfor
		endfor
		bkgd_point = 0
	elseif (bkgd == wavemax(data))
		imagestats/q data
		if (stringmatch(axis,"x"))
			bkgd_point = V_maxRowLoc
		elseif (stringmatch(axis,"y"))
			bkgd_point = V_maxColLoc
		endif
	endif
	return bkgd_point
end

function/wave get_vertical_boundaries(data)
	// required for x-axis centroid
	// columns are displayed vertically in Igor images
	// vertical edges in an image correspond to the columns in the Igor tables
	wave data
	// sum over first and last columns: rows 0 and n-1 held constant
	variable n = dimsize(data, 0)
	variable q=dimsize(data, 1)
	make/free/n=(q) first_column = data[0][p]
	make/free/n=(q) last_column = data[n-1][p]
	make/free/n=0 vertical_boundaries
	concatenate/np/kill/o {first_column, last_column}, vertical_boundaries
	return vertical_boundaries
end

function/wave get_horizontal_boundaries(data)
	// required for y-axis centroid
	// rows are displayed horizontally in Igor images
	// horizontal edges in an image correspond to the rows in the Igor tables
	wave data
	// sum over first and last rows: columns 0 and m-1 held constant
	variable m = dimsize(data, 1)
		variable q=dimsize(data, 0)
	make/free/n=(q) first_row = data[p][0]
	make/free/n=(q) last_row = data[p][m-1]
	make/free/n=0 horizontal_boundaries
	concatenate/np/kill/o {first_row, last_row}, horizontal_boundaries
	return horizontal_boundaries
end

// Test Functions //
// gaussian profiles have {bkgd, amp, x0, x_width, y0, y_width, correlation}
function create_gaussian()
	variable nx=256, ny=256
	make/o/n=(nx,ny) test_gaussian
	make/o/n=7 w_coef = {0, 5, nx/2, 010, ny/2, 10, 0}
	make/o/n=256 x_ax=x, y_ax=x
	wave test_gaussian, x_ax, y_ax, w_coef
	//test_gaussian[][] = gauss2d(w_coef, x_ax[p], y_ax[q]) + enoise(2e-3)
	test_gaussian[][] = gauss2d(w_coef, x_ax[p], y_ax[q])
end

function create_gaussian2()
	variable nx=6, ny=6
	make/o/n=(nx,ny) test_gaussian
	make/o/n=7 w_coef = {0, 5, nx/2, 50, ny/2, 50, 0}
	make/o/n=10 x_ax=x, y_ax=x
	wave test_gaussian, x_ax, y_ax, w_coef
	test_gaussian[][] = gauss2d(w_coef, x_ax[p], y_ax[q]) + enoise(2e-3)
	//test_gaussian[][] = gauss2d(w_coef, x_ax[p], y_ax[q]) 
end

function test1()
	create_gaussian()
	print get_centroids(root:test_gaussian, root:x_ax, root:y_ax)
	create_gaussian2()
	print get_centroids(root:test_gaussian, root:x_ax, root:y_ax)
end

Function ConvertData()
	string WaveListString=WaveList("SLD*",";","")
	variable tot=itemsinlist(WaveListstring)
	
	variable HighDist=50e-6
	variable LowDist=-50e-6
	variable steps=50
	variable HighDistCount=HighDist/2e-7
	variable LowDistCount=LowDist/2e-7
	variable stepsCount=(HighDist-LowDist)/steps
	string Number
	variable n
	for(n=0;n<tot;n+=1)
		wave w1=$stringfromlist(n,WaveListString)
			
		Number=num2str(round(1e6*( LowDist+stepsCount*n)))
		if(str2num(number)<0)
			Number="NewCypherStart_m"+num2str(abs(str2num(Number)))

		else
			Number="NewCypherStart_p"+num2str(abs(str2num(Number)))

			
		endif
		make/o/n=(dimsize(w1,0),dimsize(w1,1)) $Number
		wave w2=$Number
		w2=w1[p][q][2]
		SetScale/P x 0,dimdelta(w1,0),"", w2
		SetScale/P y 0,dimdelta(w1,1),"", w2
		
	endfor


end

Function RunAll()
	string WaveListString=WaveList("NewCypherStart_p*",";","")

	variable tot=itemsinlist(WaveListstring)
 wave wave0
 variable zfocus,start
	variable n
	for(n=0;n<tot;n+=1)
		wave w1=$stringfromlist(n,WaveListString)
		start=strsearch(nameofwave(w1),"_p",0)
		zfocus=str2num(nameofwave(w1)[start+2,start+5])
		InsertPoints (dimsize(wave0,0)), 1, wave0
		wave0[(dimsize(wave0,0))-1]=zfocus*1e-6
		DE_RUN(w1)
		
	endfor
	
		 WaveListString=WaveList("NewCypherStart_m*",";","")

	 tot=itemsinlist(WaveListstring)
	for(n=0;n<tot;n+=1)
		wave w1=$stringfromlist(n,WaveListString)
		start=strsearch(nameofwave(w1),"_m",0)
		zfocus=str2num(nameofwave(w1)[start+2,start+5])
		InsertPoints (dimsize(wave0,0)), 1, wave0
		wave0[(dimsize(wave0,0))-1]=zfocus*-1e-6
		DE_RUN(w1)
		
	endfor


end