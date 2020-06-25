pro planck_table

;from 1.5 to 100 micron, 1.5e4-1e6 A
wave=findgen(500)*9.85e5/499.+1.5e4

;from 100 to 3500 K
temp=findgen(35)*100.+100.

fac=fltarr(35)

openw,17,'planck_table.txt'
printf,17,'Temperature [K], Bbody flux integrated from 1.5 to 100 micron [W/m^2], sigma*T^4 [W/m^2], int/sig'

for i=0,34 do begin
   bbody=planck(wave,temp[i])  ; in cgs flux units: erg/cm^2/s/A
   intgrl=int_tabulated(wave,bbody)*1.e-3  ; erg/cm^2/s *1.e-3 = J/m^2/s
   sigt=5.6704e-8*(temp[i])^4.  ; J/m^2/s/K^4
   fac[i]=intgrl/sigt
   printf,17,temp[i],intgrl,sigt,fac[i]
endfor

close,17

;psopen,'plot.eps',/enc,/landscape
;!p.font=0
plot,temp,fac,/ystyle,xtit='Temperature [K]',ytit='Correction to blackbody flux';,xthick=5,ythick=5,thick=5,charsize=1.5
;psclose

end
