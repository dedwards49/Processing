#pragma rtGlobals=3		// Use modern global access method and strict wave access.

function PareDownData(w1,start,decimate,smoothing,SmoothAll)
wave w1
variable start,decimate,smoothing,SmoothAll

variable startpnt=x2pnt(w1,start)
Duplicate/O w1,Hold1, Ox,HOld,Hold2
Ox=(deltax(w1))*p+(pnt2x(w1,0))
if(smoothing==0)
else
if(SMoothAll==1)
Smooth/M=0 smoothing, Hold1
Duplicate/O Hold1, Hold2

else
Smooth/M=0 smoothing, Hold2
endif
endif

//Duplicate/O Hold1, Hold2, Nx
Duplicate/O Hold1, Nx
 Nx=(deltax(Hold1))*p+(pnt2x(Hold1,0))

Resample/DOWN=(decimate)/n=1/WINF=None Hold2
Resample/DOWN=(decimate)/n=1/WINF=None Nx
variable startpnt2=x2pnt(Hold2,start)
 
deletepoints (startpnt+1), 1e8, Hold1,Ox

deletepoints 0,startpnt2, Hold2,Nx
variable tot=numpnts(Hold2)+numpnts(hold1)-1
make/o/n=(tot+1) finaly,FinalX
Finaly[0,numpnts(Hold1)-1]=Hold1
Finalx[0,numpnts(Hold1)-1]=Ox

Finaly[numpnts(Hold1),numpnts(Hold1)+numpnts(Hold2)-1]=Hold2[p-numpnts(Hold1)]
Finalx[numpnts(Hold1),numpnts(Hold1)+numpnts(Hold2)-1]=NX[p-numpnts(Hold1)]
 
string FNameX, FnameY
FNameX=nameofwave(w1)+"_redx"
FNamey=nameofwave(w1)+"_redy"
duplicate/o Finalx $FNameX

duplicate/o Finaly $FNamey
killwaves Hold,Hold1,Hold2,Nx,OX,Finalx,Finaly

end