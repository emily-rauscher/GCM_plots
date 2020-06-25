pro solve,alpha

Tdeep=1700.
Ttop=800.
Tirr=1000.
Tint=500.

Ptop=1.e-5
Pdeep=1.e2
Piso=1.

g=2.2e3
f=0.5

kth=g*(alpha+1)/(pdeep*1.e6*pdeep^alpha)*(4./3*(tdeep/tint)^4-2./3)
kv=kth*ptop^alpha*2./sqrt(3)*(2./f*(ttop/tirr)^4-1./f*(tint/tirr)^4-1.)

if piso eq 1 then begin
   talpha=g/sqrt(3.)/(piso*1.e6)*kv/(kth^2-kv^2)
endif

if n_elements(alpha) eq 1 then begin
   print,'kth, kv, talpha:',kth,kv,talpha
endif else begin
   window,0
   plot,alpha,kth
   window,1
   plot,alpha,kv
   window,2
   plot,alpha,talpha
endelse

end
