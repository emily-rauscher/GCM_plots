pro magdrag,allplots=allplots,etaplot=etaplot,tmagplot=tmagplot,usesimp=usesimp

  p0=1.e2  ;bar
  oom=4
  gascon=3523.  ;SI
  radea=4.68e7  ;SI
  ga=20.9  ;SI
  ww=2.53e-5

  BFIELD=1.0  ; in G
  TDRAG_MIN=0.005  ; in units of a planet day
  tdrag_min*=(2.*!pi/ww)  ; in seconds

openr, 17, 'fort.26'
readf,17,nlat,nlon,nlev
xy=fltarr(6,nlat,nlon,nlev)  ;0: lon, 1: lat, 2: level, 3: u, 4: v, 5: T
readf,17,xy
close,17

sigma=get_sigma(oom,nlev)

openr, 18, 'fort.50'
readf, 18, nlat,nlon
ab=fltarr(3,nlat,nlon)
readf, 18, ab
close, 18
sp=reform(ab[2,*,*])
sp = (sp+1.)*p0*1.e5  ;SI

rhocgs=fltarr(nlat,nlon,nlev)
tdrag=fltarr(nlat,nlon,nlev)
udrag=fltarr(nlat,nlon,nlev)
heating=fltarr(nlat,nlon,nlev)

for ilat=0,nlat-1 do begin
   ang_f=abs(cos(!pi/2.-xy[1,ilat,0,0]*!pi/180.))
   for ilev=0,nlev-1 do begin
      rhocgs[ilat,*,ilev]=reform(sp[ilat,*]*sigma[ilev]/gascon/xy[5,ilat,*,ilev])*1.e-3 ;cgs
      tdrag[ilat,*,ilev]=tmag(reform(rhocgs[ilat,*,ilev]),reform(xy[5,ilat,*,ilev]),$
                              bfield,simp=usesimp)/ang_f ; in sec
   endfor
endfor

tdmin=min(tdrag,max=tdmax,/nan)
print,'Tdrag_min, min, max tdrag [s]:',tdrag_min,tdmin,tdmax
;plot,tdrag,/ylog,/ystyle,/xstyle,ytit='Tdrag [s]'
;oplot,[0,nlat*nlon*nlev],[tdrag_min,tdrag_min]

if keyword_set(allplots) then begin
   window,0
   plot,xy[5,*,*,*],tdrag,/ylog,/ystyle,xtit='Temperature [s]',ytit='Tdrag [s]',psym=3
   oplot,[0,max(xy[5,*,*,*])],[tdrag_min,tdrag_min]
   window,2
   minrho=min(rhocgs,max=maxrho)
   plot,rhocgs,tdrag,xr=[minrho,maxrho],yr=[tdmin,tdmax],ytit='Tdrag [s]',xtit='Density [cgs]',psym=3,/xlog,/ylog
   oplot,[minrho,maxrho],[tdrag_min,tdrag_min]
endif

udrag=reform(-xy[3,*,*,*])/tdrag   ; in m/s^2
heating=-udrag*reform(xy[3,*,*,*]) ; in J/s/kg
pres=p0*get_sigma(oom,nlev)
dP=fltarr(nlev)     
for ilev=1,nlev-2 do dP[ilev]=0.5*(pres[ilev+1]-pres[ilev-1]) ; to match swdp
dP[0]=0.5*(pres[1])
dP[nlev-1]=0.5*(p0-pres[nlev-2])
for ilev=0,nlev-1 do heating[*,*,ilev]=heating[*,*,ilev]*(dP[ilev]*1.e5/ga) $
   *cos(xy[1,*,*,ilev]*!pi/180.)*radea^2*(2.*!pi/nlon)*(!pi/nlat) ; W

if keyword_set(allplots) then begin
   window,3
   umin=min(udrag,max=umax)
   plot,rhocgs,udrag,xr=[minrho,maxrho],yr=[umin,umax],ytit='Udrag [m/s^2]',xtit='Density [cgs]',psym=3,/xlog,/ylog
endif

tlimit=where(tdrag lt tdrag_min)
if tlimit[0] ne -1 then tdrag[tlimit]=tdrag_min
udrag=-reform(xy[3,*,*,*])/tdrag   ; in m/s^2

if keyword_set(allplots) then begin
   oplot,rhocgs,udrag,psym=4

   window,4
   hmin=min(heating,max=hmax)
   plot,rhocgs,heating,xr=[minrho,maxrho],yr=[hmin,hmax],ytit='Heating [W]',xtit='Density [cgs]',psym=3,/xlog,/ylog
endif
heating=-udrag*reform(xy[3,*,*,*]) ; in J/s/kg
for ilev=0,nlev-1 do heating[*,*,ilev]=heating[*,*,ilev]*(dP[ilev]*1.e5/ga) $
   *cos(xy[1,*,*,ilev]*!pi/180.)*radea^2*(2.*!pi/nlon)*(!pi/nlat) ; W
if keyword_set(allplots) then oplot,rhocgs,heating,psym=4

if keyword_set(etaplot) or keyword_set(tmagplot) then begin
   rho=reform(rhocgs,nlat*nlev*nlon)
   temp=reform(xy[5,*,*,*],nlat*nlev*nlon)
   tdrag=tmag(rho,temp,bfield,eta=eta,simp=usesimp)/abs(sin(!pi/180.*reform(xy[1,*,*,*])))
   psopen,'plot.eps',/enc,/color,/landscape
   !p.font=0
   plot,rho,temp,/xlog,psym=3,/ynoz,xtit='Density [cgs]',ytit='Temperature [K]',charsize=2,xthick=4,ythick=4
   if keyword_set(etaplot) then begin
      loadct,4,/silent
      leta=alog10(eta)
      lmax=max(leta,min=lmin)
      for i=0,n_elements(rho)-1 do oplot,[rho[i],rho[i]],[temp[i],temp[i]],color=255-(leta[i]-lmin)*255./(lmax-lmin),$
                                         symsize=2
      colorbar,range=[lmin,lmax],format='(f5.2)',position=[0.2,0.88,0.84,0.92],charsize=1.5,/reverse
   endif else begin
      loadct,13,/silent
      ltdrag=alog10(tdrag)
      tmax=max(ltdrag,min=tmin)
      for i=0,n_elements(rho)-1 do oplot,[rho[i],rho[i]],[temp[i],temp[i]],$
                                         color=255-(ltdrag[i]-tmin)*255./(tmax-tmin),symsize=2
      colorbar,range=[tmin,tmax],format='(f5.2)',position=[0.2,0.88,0.84,0.92],charsize=1.5,/reverse
   endelse      
   loadct,0,/silent
   !p.font=-1
   psclose
;   spawn,'gs -r300 -dEPSCrop -dTextAlphaBits=4 -sDEVICE=png16m -sOutputFile=plot.png -dBATCH -dNOPAUSE plot.eps'
;   spawn,'convert plot.png eps3:plot.eps'
endif

end
