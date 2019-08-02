#pragma rtGlobals=3		// Use modern global access method and strict wave access.
function Sader_Re(p,w,b,n)
variable  p,w,b,n

return p*w*b^2/4/n

end

function/C Sader_GammaCirc(Re)
variable  Re

return 1+cmplx(0,4)*BesselK(1,cmplx(sqrt(Re)/sqrt(2),-sqrt(Re)/sqrt(2)))/BesselK(0,cmplx(sqrt(Re)/sqrt(2),-sqrt(Re)/sqrt(2)))/cmplx(sqrt(Re/2),sqrt(Re/2))

end


function/C Sader_Omegar(Re)
variable  Re
variable tau=log(Re)
variable real1=(.91324-.48274*tau+.46842*tau^2-.12886*tau^3+.044055*tau^4-.0035117*tau^5+.00069085*tau^6)/(1-.56964*tau+.48690*tau^2-.13444*tau^3+.045155*tau^4-.0035862*tau^5+0.00069085*tau^6)
variable imag1=(-.024134-.029256*tau+.016294*tau^2-0.00010961*tau^3+0.000064577*tau^4-.000044510*tau^5)/(1-.59702*tau+.55182*tau^2-.18357*tau^3+.079156*tau^4-.014369*tau^5+.0028361*tau^6)

return cmplx(real1,imag1)
end


function/C Sader_GammaRect(Re)
variable  Re
return Sader_Omegar(Re)*Sader_GammaCirc(Re)

end


function Sader_k(f,b,L,Q)
variable f,b,L,Q
variable w=2*pi*f
variable pf=1.18
variable n=1.86e-5
variable Re=Sader_RE(pf,w,b,n)
return 0.1906*pf*b^2*imag(Sader_GammaRect(Re))*L*Q*w^2
end
