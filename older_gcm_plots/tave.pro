pro tave,p0,oom,tint=tint,kth=kth,alpha=alpha,tirr=tirr,g=g,kv=kv,psfile=psfile,xlog=xlog

;grab fort.26 info
openr, 17, 'fort.26'
readf,17,nlat,nlon,nlev
xy=fltarr(6,nlat,nlon,nlev)  ;0: lon, 1: lat, 2: level, 3: u, 4: v, 5: T
readf,17,xy
close,17

pres=p0*get_sigma(oom,nlev)

tmin=min(xy[5,0,*,*],max=tmax)

if not keyword_set(psfile) then begin
; Screen output
   device,true_color=24,decomposed=0,retain=2
   !p.font=-1
   !p.charsize=1.5
endif else begin
   if (size(psfile))[1] eq 2 then filename='plot.eps' else filename=psfile
   psopen,filename,/enc,/color,xsize=25,ysize=20 ;,bits_per_pixel=24
   !p.font=0
   !p.charsize=2
endelse

plot,xy[5,0,0,*],pres,/ylog,yr=[max(pres),min(pres)],xr=[tmin,tmax],$
     xtit='Temperature [K]',ytit='Pressure [bar]',xthick=6,ythick=6,xlog=xlog,xmargin=[12,3]
for ilat=0,nlat-1 do begin
   for ilon=0,nlon-1 do begin
      oplot,xy[5,ilat,ilon,*],pres,color=120
   endfor
endfor

tg=guillot(pres,1.,f=0.25,tint=tint,kth=kth,alpha=alpha,tirr=tirr,g=g,kv=kv)
avet=fltarr(nlev)
for ilev=0,nlev-1 do begin
   for ilat=0,nlat-1 do begin
      clat=cos(!pi/180.*xy[1,ilat,0,0])
      avet[ilev]+=mean(xy[5,ilat,*,ilev])*clat
   endfor
   avet[ilev]/=total(cos(!pi/180.*xy[1,*,0,0]))
endfor

;loadct,5,/silent
oplot,tg,pres,thick=6,linestyle=2;,color=100
oplot,avet,pres,thick=6;,color=100
al_legend,['Average profile','Guillot profile'],linestyle=[0,2],/right,charsize=1.5,thick=6;,color=100

loadct,0,/silent

if keyword_set(psfile) then begin
   psclose
;   spawn,'gs -r300 -dEPSCrop -dTextAlphaBits=4 -sDEVICE=png16m -sOutputFile=plot.png -dBATCH -dNOPAUSE plot.eps'
;   spawn,'convert plot.png eps3:plot.eps'
endif

;ga=8.  ;m/s^2
;radea=1.e8  ;m

;dP=fltarr(nlev)
;for ilev=1,nlev-2 do dP[ilev]=0.5*(pres[ilev+1]-pres[ilev-1]) ; to match swdp
;dP[0]=0.5*(pres[1])
;dP[nlev-1]=0.5*(p0-pres[nlev-2])

;dT=avet-tg  ; K
;set=(3523./0.286)*dT  ; J/kg
;et=set*dP/ga*radea^2*4.*!pi*1.e5  ;J
;trad=pres*(3523./0.286)/ga/4./5.67e-8/avet^3*1.e5  ;s
;let=et/trad  ;W

end
