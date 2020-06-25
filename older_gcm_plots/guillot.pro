function guillot, press,cosa,f=f,nightside=nightside,tint=tint,kth=kth,alpha=alpha,$
                  tirr=tirr,g=g,kv=kv,plott=plott
; input press: pressure in bars

if not keyword_set(tint) then tint=500. ; in K
if not keyword_set(g) then g=8.e2 ; in cgs
if not keyword_set(tirr) then tirr=2077.5 ; in K, =2077.5 for fig 2, =2013 for solc=9.31e5
if not keyword_set(kth) then kth=1.e-2 ; in cgs
if not keyword_set(kv) then kv=4.e-3 ; in cgs
if keyword_set(nightside) then tirr=0.

;cosa=1. ; cos of angle from substellar
;f=1. ;1 at substellar, 0.5 for dayside ave, 0.25 for total ave
; f=0.375 seems like a good average for the initial profile

; tau= kth/g * P, convert from bar to dyne/cm2 with 1.e6

gamma=kv/kth
if keyword_set(alpha) then gamma=gamma/(press^alpha)
; (kth = k_o * (press/1 bar)^alpha, so gamma unitless

nlev=n_elements(press)

if not keyword_set(f) then begin
   if not keyword_set(alpha) then begin
   ; eqn 27
      temp=[ 0.75*Tint^4*(2./3. + (kth/g)*press*1.e6) $
             + 0.75*Tirr^4*cosa*(2./3. + cosa/gamma + $
                                 (gamma/3./cosa - cosa/gamma) $
                                 *exp(-gamma*kth/g*press*1.e6/cosa))]^0.25
   endif else begin
      ; tau= kth/g *P * (P/P_ref)^alpha
      tau=(kth/g)*1.e6*press^(alpha+1)  ;unitless
      intgrl=fltarr(nlev)
      func=[0,(press^alpha)*exp(-gamma*tau/cosa)]  ;unitless
      funcx=[0,press*1.e6]
;      plot,funcx[1:36],func[1:36],/xlog
      for i=0,nlev-1 do begin
         ;integration wrong if func decreases too quickly
         test=int_tabulated(funcx[0:i+1],func[0:i+1])
         if ((i gt 0) and (test lt intgrl[i-1])) then $
            intgrl[i]=intgrl[i-1]+func[i+1]*(funcx[i+1]-funcx[i]) $
         else intgrl[i]=test
;         oplot,funcx[0:i],func[0:i],thick=3
;         print,'',test
;         print,'',intgrl[i-1]+func[i+1]*(funcx[i+1]-funcx[i])
;         print,intgrl[i]
;         print,''
;         wait,3
      endfor
      temp=[0.75*Tint^4*(2./3.+tau/(alpha+1)) $
            + 0.75*Tirr^4*cosa*(2./3. + [gamma/3./cosa]*exp(-gamma*tau/cosa) $
                                + (kth/g)*intgrl)]^0.25
   endelse
endif else begin
   if not keyword_set(alpha) then begin
   ; eqn 29
      temp=[ 0.75*Tint^4*(2./3. + (kth/g)*press*1.e6) $
             + 0.75*Tirr^4*f*(2./3. + 1./(gamma*sqrt(3.)) + $
                              (gamma/sqrt(3.) - 1./(gamma*sqrt(3.))) $
                              *exp(-gamma*kth/g*press*1.e6*sqrt(3.)))]^0.25
   endif else begin
      tau=(kth/g)*1.e6*press^(alpha+1)
      intgrl=fltarr(nlev)
      func=[0,(press^alpha)*exp(-sqrt(3.)*gamma*tau)]
      funcx=[0,press*1.e6]
      for i=0,nlev-1 do begin
         ;integration wrong if func decreases too quickly
         test=int_tabulated(funcx[0:i+1],func[0:i+1])
         if ((i gt 0) and (test lt intgrl[i-1])) then $
            intgrl[i]=intgrl[i-1]+func[i+1]*(funcx[i+1]-funcx[i]) $
         else intgrl[i]=test
      endfor
      temp=[0.75*Tint^4*(2./3.+tau/(alpha+1)) $
            + 0.75*Tirr^4*f*(2./3. + [gamma/sqrt(3.)]*exp(-sqrt(3.)*gamma*tau) $
                                + (kth/g)*intgrl)]^0.25
   endelse
endelse

if keyword_set(plott) then begin
   plot,temp,press,/ylog,yr=[max(press),min(press)]
   if plott eq 1 then unstable=0.286 else unstable=plott
   convec=fltarr(nlev)
   test=deriv(alog(press),alog(temp))
   convi=where(test ge unstable)
   if convi[0] ne -1 then convec[convi]=temp[convi]
   oplot, convec, press, psym = 2
endif

return,temp

end
