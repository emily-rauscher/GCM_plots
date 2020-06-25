function tmag,rho,temp,Bfield,simp=simp,eta=eta,metals=metals
; in cgs
; simp uses potassium only
; rho and T can be 1D arrays, of the same size

gascon=3523.

asize=n_elements(rho)
if n_elements(temp) ne asize then stop

if keyword_set(simp) then begin
   xe=6.47e-13 * sqrt(1.e-7/1.e-7) *(temp/1.e3)^0.75 * sqrt(2.4e15/rho*1.6726e-24) * exp(-25188./temp)/1.15e-11
   eta=230.*sqrt(temp)/xe
endif else begin
   f = [2.66e10, 1.8e9, 60, 1.2, 9, 1.11e7, 2.31e6, 1.84e7, $
      780, 2.6e6, 6.0e4, 1.06e6, 8.5e4, 1.0e6, 6500, 5.0e5, 4740, $
      1.06e5, 3500, 6.25e4, 31, 2400, 254, 1.27e4, 9300, 9.0e5, $
      2200, 4.78e4]
   if keyword_set(metals) then begin
      f*=3.
      f[0:1]=[2.66e10, 1.8e9]
   endif

   c = [13.6, 24.6, 5.4, 9.32, 8.30, 11.26, 14.54, 13.61, $
        17.42, 21.56, 5.14, 7.64, 5.98, 8.15, 10.55, 10.36, 13.01, $
        15.76, 4.34, 6.11, 6.56, 6.83, 6.74, 6.76, 7.43, 7.90, $
        7.86, 7.63]

;   mu = 1.67e-24
   mu = 8.3145/gascon*1.e3/6.022e23  ; cgs mass of particle
   k = 1.38e-16
   h = 6.62e-27
   me = 9.11e-28
;   mp = 1.67e-24
   e = 4.8e-10
   Cl = 3.e10  ;speed o light 

   ndens=Rho/mu ;asize
   xe=fltarr(asize)

   for i=27,0,-1 do begin
      f(i)=f(i)/f(0)  ; abundance relative to H 
      ni=ndens*f(i)
      Ki=(1.0/(ni*k*temp))*((2.0*!pi*me)^(3.0/2.0)/h)* $
         ((k*temp)^(5.0/2.0)/h) *(exp(-c(i)*11605./temp)/h)
      xi=(Ki/(1.+Ki))^0.5 
      xe=xe+f(i)*xi
   endfor

   eta=230.*sqrt(temp)/xe

endelse

tdrag=4.*!pi*rho*eta/Bfield^2

return,tdrag

end
