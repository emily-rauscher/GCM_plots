pro Ri,psfile=psfile
;plot of zonal Ri along equator, as function of longitude and pressure

oom=5.34363
p0=220.612
kappa=0.321
gascon=4593.
ga=9.42


  OpenR, 17, 'fort.26'
  ReadF, 17, nlat,nlon,nlev	  
;	  print,nlat,nlon,nlev 
	  xy=fltarr(6,nlat,nlon,nlev)   
  ReadF, 17, xy
  Close, 17
    lon=reform(xy[0,*,*,*]) & lat=reform(xy[1,*,*,*])
    u=reform(xy[3,*,*,*]) & v=reform(xy[4,*,*,*]) & temp=reform(xy[5,*,*,*])
 
  rn=fltarr(nlon,nlev)
  n2=fltarr(nlon,nlev)
  flon=reform(lon[0,*,*])
  flev=fltarr(nlon,nlev)

   sig=get_sigma(oom,nlev)
   for i=0,nlev-1 do flev(*,i)=replicate(sig[i],nlon)
   flev*=p0
   pres=sig*p0

dtdp=fltarr(nlon,nlev)
dudp=fltarr(nlon,nlev)

ilat=nlat/2
dtdp[*,0]=(temp[ilat,*,0]-temp[ilat,*,1])/(pres[0]-pres[1])
dudp[*,0]=(u[ilat,*,0]-u[ilat,*,1])/(pres[0]-pres[1])
dtdp[*,nlev-1]=(temp[ilat,*,nlev-2]-temp[ilat,*,nlev-1])/(pres[nlev-2]-pres[nlev-1])
dudp[*,nlev-1]=(u[ilat,*,nlev-2]-u[ilat,*,nlev-1])/(pres[nlev-2]-pres[nlev-1])
for i=1,nlev-2 do begin
   dtdp[*,i]=0.5*((temp[ilat,*,i]-temp[ilat,*,i+1])/(pres[i]-pres[i+1]) $
                                  +(temp[ilat,*,i-1]-temp[ilat,*,i])/(pres[i-1]-pres[i]))
   dudp[*,i]=0.5*((u[ilat,*,i]-u[ilat,*,i+1])/(pres[i]-pres[i+1]) $
                                  +(u[ilat,*,i-1]-u[ilat,*,i])/(pres[i-1]-pres[i]))
endfor

test=fltarr(nlon,nlev)

for i=0,nlev-1 do begin
   n2[*,i]=(ga^2/gascon/temp[ilat,*,i])$
      *(replicate(kappa,nlon)-(pres[i]/temp[ilat,*,i])*dtdp[*,i])
   ;rn[*,i]=(gascon*temp[ilat,*,i]/(pres[i])^2.) $
   ;        *(replicate(kappa,nlon)-(pres[i]/temp[ilat,*,i])*dtdp[*,i]) $
   ;        /(dudp[*,i])^2.
   rn[*,i]=n2[*,i]*(gascon*temp[ilat,*,i]/pres[i]/ga)^2/(dudp[*,i])^2
   ;oops=where(test[*,i] ne rn[*,i])
   ;if oops[0] ne -1 then print, 'oops!'
   oops=where(n2[*,i] lt 0.)
   if oops[0] ne -1 then rn[oops,i]=0.
endfor

quant=dudp

umax=max(quant,min=umin)
  loadct, 4
;test=where(rn ne 0.)
;umin=min(rn[test])
;print,umin
umin=-1.e5
umax=1.e5

nlevels=25
step=(umax-umin)/nlevels
ulevels=indgen(nlevels)*step+umin

if not keyword_set(psfile) then begin
; Screen output
  device,true_color=24,decomposed=0,retain=2
;  set_plot, 'x'
  window, 0, xsize=x_size_window, ysize=y_size_window, retain=2
  !p.font=-1
endif else begin
   if (size(psfile))[1] eq 2 then filename='plot.eps' else filename=psfile
   psopen,filename,/enc,/color,xsize=25,ysize=20;,bits_per_pixel=24
  !p.font=0
endelse

  !x.style=1
  !p.charsize=1.5

   contour, quant,flon,flev, /cell_fill ,LEVELS=ulevels,/ystyle,ymargin=[4,4],$
            yr=[max(flev),min(flev)],/ylog,xtitle='Longitude [degrees]',$
            ytitle='Pressure [bar]',max_value=umax,min_value=umin,ytickname=ylabels

contour,quant,flon,flev,/overplot,levels=[0.25],color=255,thick=5,/closed

colorbar,position=[0.23,0.94,0.9,0.98],range=[umin,umax],charsize=1.5,color=0,format='(e8.1)'

   print, ' Min and Max Zonal Richardson numbers:      ',min(quant),max(quant)

if not keyword_set(psfile) then begin
  set_plot, 'x'
endif else begin
   psclose
endelse

loadct,0


end
