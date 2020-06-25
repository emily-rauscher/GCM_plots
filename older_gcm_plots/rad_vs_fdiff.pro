pro rad_vs_fdiff,nlat,nlon,nlev,mbar=mbar,fbase=fbase,psfile=psfile,fdiff=fdiff

; plots the total LW (and SW) fluxes at each flux level (or maybe heating rate
; at full levels?), as well as what radiative transfer alone would have given
; reads fort.63
;
; only plots cirrad fluxes if fdiff keyword set
;
; pressures in bar, unless mbar keyword used (Jun 2011 change to
; fort.63 output format)

; format of fort.63 is:
; LATITUDE, LONGITUDE:  -81.6505907502981       0.000000000000000E+000
;           UPWARD FLUX (Wm-2)     DOWNWARD FLUX (Wm-2)     NET FLUX (Wm-2)           SW NET FLUX (Wm-2)
;         0.616      0.21033E+06           0.00000E+00           0.21033E+06           0.15074E+06
;       .... (nlev+1 entries total, last being surface) ........
;    101645.531      0.33613E+07           0.33578E+07           0.35000E+04           0.00000E+00

entry = create_struct( 'lat', 0.0, $
                       'lon', 0.0, $
                       'pres', fltarr(nlev), $
                       'sflux', fltarr(nlev), $  ; SW flux
                       'flux', fltarr(nlev), $   ; actual IR flux
                       'rflux', fltarr(nlev))    ; cirrad flux
nent=(nlat-1)*nlon
rates=replicate(entry,nent)

read1=strarr(2)
data=fltarr(5,nlev+1)
read4=''
toss=fltarr(3*(nlev+1))
read5=''

openr,17,'fort.63'

for i=0,nent-1 do begin
   readf,17,read1
;   print,read1[0]
   rates[i].lat=float(strmid(read1[0],23,6))
;   print,rates[i].lat
   rates[i].lon=float(strmid(read1[0],47,6))
;   print,rates[i].lon
   readf,17,data
   rates[i].pres=reform(data[0,0:nlev-1])
   rates[i].rflux=reform(data[1,0:nlev-1]-data[2,0:nlev-1])
   rates[i].flux=reform(data[3,0:nlev-1])
   rates[i].sflux=reform(data[4,0:nlev-1])
   readf,17,read4
;   print,read4
   readf,17,toss
   readf,17,read5
;   print,read5
; fix to deal with rates too small the exponent gets incorrectly dropped (making them big again):
;   test=where(data[1,*] lt min(abs(data[2,0:nlev-1])))
;   if test[0] ne -1 then rates[i].sw[test[0]:nlev-1]=replicate(0.,nlev-test[0])
endfor

close,17

if keyword_set(psfile) then begin
   if (size(psfile))[1] eq 2 then filename='plot.eps' else filename=psfile
   psopen,filename,/enc,/color,xsize=20,ysize=15;,bits_per_pixel=24
   !p.thick=5
   !x.thick=5
   !y.thick=5
   !p.font=0
endif else begin
  device,true_color=24,decomposed=0,retain=2
endelse
  !p.charsize=1.5

umin=min(rates.flux,max=umax)
if keyword_set(fdiff) then begin
   rmin=min(rates.rflux,max=rmax)
   if rmin lt umin then umin=rmin
   if rmax gt umax then umax=rmax
endif
if keyword_set(fbase) then umin=0.
pmin=min(rates.pres,max=pmax)
plot,rates[0].flux,rates[0].pres,yr=[pmax,pmin],/ylog,xr=[umin,umax],$
     xtit='Net flux [W/m^2]',ytit='Pressure [bar]',/xlog,xmargin=[12,3]
if keyword_set(fdiff) then oplot,rates[0].rflux,rates[0].pres,linestyle=2
oplot,rates[0].sflux,rates[0].pres,linestyle=1
;for i=1,nent-1 do begin
;   oplot,rates[i].flux,rates[i].pres
;   oplot,rates[i].rflux,rates[i].pres,linestyle=2
;endfor
if keyword_set(fbase) then begin
   oplot,[fbase,fbase],[pmax,pmin],linestyle=3
   al_legend,['Corrected IR flux','Uncorrected IR flux','SW flux','Interior flux'],linestyle=[0,2,1,3],charsize=1.25
endif else al_legend,['Corrected IR flux','Uncorrected IR flux','SW flux'],linestyle=[0,2,1],charsize=1.25

;stop

if keyword_set(psfile) then begin
   !p.thick=1
   !x.thick=1
   !y.thick=1
   psclose
endif

end
