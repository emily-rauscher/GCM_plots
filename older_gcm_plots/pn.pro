pro pn,oom=oom,p0=p0,ga=ga,radea=radea,oldheader=oldheader,missinglat=missinglat,psfile=psfile

; calculate advective heat transport between the hemispheres
; try a few ways of calculating this, compare terms
; include check of vertical heat transport

; uses: fort.26, fort.50, fort.63

if not keyword_set(oom) then oom=5.
if not keyword_set(p0) then p0=1.e2  ; bar
gascon=3523.  ; J/kg/K
akap=0.286
cp=gascon/akap  ; J/kg/K
if not keyword_set(ga) then ga=8. ; m/s^2
if not keyword_set(radea) then radea=1.e8  ; m

print, 'Are these variables correct? gascon, akap, p0, oom, ga, radea:'
print,gascon,akap,p0,oom,ga,radea

if keyword_set(psfile) then begin
   if (size(psfile))[1] eq 2 then filename='plot.eps' else filename=psfile
   !x.thick=5
   !y.thick=5
   !p.thick=5
   psopen,filename,/enc,/color,xsize=6,ysize=9,/inches
   !p.font=0
   !p.charsize=1.25
;   cthick=6
endif else begin
;   device,true_color=24,decomposed=0;,retain=2
;   cthick=3
endelse

  OpenR, 17, 'fort.26'
  ReadF, 17, nlat,nlon,nlev	  
	  xy=fltarr(6,nlat,nlon,nlev)   
  ReadF, 17, xy
  Close, 17
  lon=reform(xy[0,*,*,*]) & lat=reform(xy[1,*,*,*]) 
  u=reform(xy[3,*,*,*]) & v=reform(xy[4,*,*,*]) 
  temp=reform(xy[5,*,*,*])

  openr, 18, 'fort.50'
  readf, 18, nlat,nlon
  ab=fltarr(3,nlat,nlon)
  readf, 18, ab
  close, 18
  sp=reform(ab[2,*,*])
  sp = (sp+1.)*p0  ; bar

sigma=get_sigma(oom,nlev)
pres=p0*sigma   ; bar
dP=fltarr(nlev)   ; in bar                                                     
for ilev=1,nlev-2 do dP[ilev]=0.5*(pres[ilev+1]-pres[ilev-1])   ; to match swdp
dP[0]=0.5*(pres[1])
dP[nlev-1]=0.5*(p0-pres[nlev-2])

rho=fltarr(nlat,nlon,nlev)
for ilev=0,nlev-1 do begin
   rho[*,*,ilev]=(sigma[ilev]*sp*1.e5)/gascon/temp[*,*,ilev]  ; kg/m^3
endfor
   
vertvel,ga,gascon,radea,oom,p0,/psurf,/noplot,om=om  ; om[nlat,nlon,nlev] in bar/s

; horizontal advective heat transport, over terminator (from day to night)
terme=nlon/4   ; 90 degrees
termw=3*nlon/4 ; -90 degrees

dnadv=fltarr(3,nlev) ; 0:east term, 1:west, 2:total  (POSITIVE = DAYSIDE HEATING)
dlat=abs(deriv(lat[*,0,0]*!pi/180.))  ;dlat is latitude dependent (slightly, but might as well include)
for ilev=0,nlev-1 do begin
   dnadv[0,ilev]=-1./(4.*!pi*radea)*cp*total(reform(u[*,terme,ilev])*reform(temp[*,terme,ilev])*dlat) ;W/kg
   dnadv[1,ilev]=1./(4.*!pi*radea)*cp*total(reform(u[*,termw,ilev])*reform(temp[*,termw,ilev])*dlat) ;W/kg
endfor
dnadv[2,*]=dnadv[0,*]+dnadv[1,*]  ; W/kg

; vertical advection and adiabatic (expansion) heating
vadv=fltarr(3,nlev)  ; 0:over day side, 1:over night, 2:total
adheat=fltarr(3,nlev)  ; 0:over day side, 1:over night, 2:total

for ilat=0,nlat-1 do begin
   for ilon=0,terme-1 do begin
      dtdp=deriv(sp[ilat,ilon]*sigma,cp*reform(temp[ilat,ilon,*]))  ; J/kg/bar
      vadv[0,*]-=reform(om[ilat,ilon,*])*dtdp*cos(lat[ilat,0,0]*!pi/180.)*dlat[ilat]
      adheat[0,*]+=1.e5*reform(om[ilat,ilon,*])/reform(rho[ilat,ilon,*])$  ; W/kg
                   *cos(lat[ilat,0,0]*!pi/180.)*dlat[ilat]
   endfor
   for ilon=terme,termw-1 do begin
      dtdp=deriv(sp[ilat,ilon]*sigma,cp*reform(temp[ilat,ilon,*]))  ; J/kg/bar
      vadv[1,*]-=reform(om[ilat,ilon,*])*dtdp*cos(lat[ilat,0,0]*!pi/180.)*dlat[ilat]
      adheat[1,*]+=1.e5*reform(om[ilat,ilon,*])/reform(rho[ilat,ilon,*])$  ; W/kg
                   *cos(lat[ilat,0,0]*!pi/180.)*dlat[ilat]
   endfor
   for ilon=termw,nlon-1 do begin
      dtdp=deriv(sp[ilat,ilon]*sigma,cp*reform(temp[ilat,ilon,*]))  ; J/kg/bar
      vadv[0,*]-=reform(om[ilat,ilon,*])*dtdp*cos(lat[ilat,0,0]*!pi/180.)*dlat[ilat]
      adheat[0,*]+=1.e5*reform(om[ilat,ilon,*])/reform(rho[ilat,ilon,*])$  ; W/kg
                   *cos(lat[ilat,0,0]*!pi/180.)*dlat[ilat]
   endfor
endfor
dlon=2.*!pi/nlon

vadv[2,*]=vadv[0,*]+vadv[1,*]
vadv*=(dlon/4./!pi)  ; W/kg
adheat[2,*]=adheat[0,*]+adheat[1,*]
adheat*=(dlon/4./!pi)  ; W/kg

;print,'Global total vertical advection:',total(vadv[2,*]),' W'
;print,'(Should be zero?)'
;print,'Global total adiabatic heating:',total(adheat[2,*]),' W'

; radiative heating
entry = create_struct( 'lat', 0.0, $
                       'lon', 0.0, $
                       'p0',0.0, $
;                       'fpres', fltarr(nlev), $
;                       'sflux', fltarr(nlev), $  ; SW flux
;                       'flux', fltarr(nlev), $   ; actual IR flux
;                       'rflux', fltarr(nlev), $    ; cirrad flux
                       'pres', fltarr(nlev), $
                       'sw', fltarr(nlev), $  ;sw heating
                       'lw', fltarr(nlev))  ;lw heating
if keyword_set(missinglat) then nent=(nlat-1)*nlon else nent=nlat*nlon
rates=replicate(entry,nent)

read1=strarr(2)
data1=fltarr(5,nlev+1)
read4=''
data2=fltarr(3,nlev+1)
read5=''

openr,17,'fort.63'
if not keyword_set(oldheader) then begin
   toss=strarr(2)
   readf,17,toss
   print,toss[0]
endif
for i=0,nent-1 do begin
   readf,17,read1
   rates[i].lat=float(strmid(read1[0],23,6))
   rates[i].lon=float(strmid(read1[0],47,6))
   readf,17,data1
   readf,17,read4
   readf,17,data2
   readf,17,read5
   rates[i].pres=reform(data2[0,0:nlev-1])  ; in bar (or mbar, if keyword set)
   rates[i].p0=data2[0,nlev]
   rates[i].sw=reform(data2[1,0:nlev-1])/24./3600.*cp       ; in J/s/kg
   rates[i].lw=reform(data2[2,0:nlev-1])/24./3600.*cp       ; in J/s/kg
   ; fix to deal with rates too small the exponent gets incorrectly dropped (making them big again):
   test=where(data2[1,*] lt min(abs(data2[2,0:nlev-1])))
   if test[0] ne -1 then rates[i].sw[test[0]:nlev-1]=replicate(0.,nlev-test[0])
endfor
close,17

qrad=fltarr(3,nlev)  ; 0:day, 1:night, 2:total

night=where((rates.lon ge 90.) and (rates.lon lt 270.),complement=day)

; i'm just going to use the average dlat, instead of sorting the entries (this is an error of less than 0.1%)
adlat=mean(dlat)
for ilev=0,nlev-1 do begin
   qrad[0,ilev]=total(((rates[day].sw)[ilev,*]+(rates[day].lw)[ilev,*])*cos(rates[day].lat*!pi/180.))
;   qrad[1,ilev]=total(((rates[night].sw)[ilev,*]+(rates[night].lw)[ilev,*])*cos(rates[night].lat*!pi/180.))
   qrad[1,ilev]=total(((rates[night].lw)[ilev,*])*cos(rates[night].lat*!pi/180.))
endfor
if keyword_set(missinglat) then begin
   ; copy entry from southernmost lat to northernmost (same as just adding to total)
   mld=where(rates[day].lat eq min(rates.lat))
   mln=where(rates[night].lat eq min(rates.lat))
   for ilev=0,nlev-1 do begin
      qrad[0,ilev]+=total(((rates[day[mld]].lw)[ilev,*]+(rates[day[mld]].sw)[ilev,*])*cos(rates[day[mld]].lat*!pi/180.))
      qrad[1,ilev]+=total(((rates[night[mln]].lw)[ilev,*])*cos(rates[night[mln]].lat*!pi/180.))
   endfor
endif

qrad[2,*]=qrad[0,*]+qrad[1,*]
qrad*=adlat*dlon/4./!pi  ; W/kg
;qrad*=(radea^2/ga)*dlon*adlat  ; W/(kg/m s^2)

; vertical advection: vadv(3,nlev), 0:over day side, 1:over night, 2:total
; adiabatic (expansion) heating adheat(3,nlev), 0:over day side, 1:over night, 2:total
; radiative heating: qrad(3,nlev), 0:day, 1:night, 2:total
; horizontal advection over terminator: dnadv(3,nlev), 0:east term, 1:west, 2:total  (POSITIVE = DAYSIDE HEATING)

trate=fltarr(3,nlev) ; 0:day, 1:night, 2:total
trate[0,*]=vadv[0,*]+adheat[0,*]+qrad[0,*]+dnadv[2,*]
trate[1,*]=vadv[1,*]+adheat[1,*]+qrad[1,*]-dnadv[2,*]
trate[2,*]=vadv[2,*]+adheat[2,*]+qrad[2,*]
mp=dp*1.e5*4.*!pi*radea^2/ga  ; kg
print,'Total heating [W]:',total(trate[2,*]*mp)

xmin=min(vadv,max=xmax)
tmin=min(adheat,max=tmax)
if tmin lt xmin then xmin=tmin
if tmax gt xmax then xmax=tmax
tmin=min(qrad,max=tmax)
if tmin lt xmin then xmin=tmin
if tmax gt xmax then xmax=tmax
tmax=max(abs(dnadv[2,*]))
if tmax gt xmax then xmax=tmax
tmin=-tmax
if tmin lt xmin then xmin=tmin
tmin=min(trate,max=tmax)
if tmin lt xmin then xmin=tmin
if tmax gt xmax then xmax=tmax

if not keyword_set(psfile) then window, 0 else !P.MULTI=[0,1,2,0,0]

plot,[0,0],[p0,p0],/ylog,yr=[p0,min(pres)],xr=[xmin,xmax],$
     xtit='Heating rate [W/kg]',ytit='Pressure [bar]'
loadct,24,/silent
oplot,adheat[0,*],pres,color=250
oplot,adheat[1,*],pres,color=250,linestyle=2
oplot,adheat[2,*],pres,color=250,linestyle=1
oplot,vadv[0,*],pres,color=90
oplot,vadv[1,*],pres,color=90,linestyle=2
oplot,vadv[2,*],pres,color=90,linestyle=1
oplot,dnadv[2,*],pres,color=140
oplot,-dnadv[2,*],pres,color=140,linestyle=2
oplot,qrad[0,*],pres,color=125
oplot,qrad[1,*],pres,color=125,linestyle=2
oplot,qrad[2,*],pres,color=125,linestyle=1
loadct,0,/silent
oplot,trate[0,*],pres
oplot,trate[1,*],pres,linestyle=2
oplot,trate[2,*],pres,linestyle=1
al_legend,['On the dayside','On the nightside','Global total'],$
          linestyle=[0,2,1],/bottom,charsize=0.75
loadct,24,/silent
al_legend,['Total rate','Radiative heating','Horizontal advection','Vertical advection','Adiabatic expansion'],$
          textcolors=[190,125,140,90,250],/bottom,/right,charsize=0.75,box=0
loadct,0,/silent

if not keyword_set(psfile) then window,1
xmin=min(vadv[0,*]*mp,max=xmax)
tmin=min(adheat[0,*]*mp,max=tmax)
if tmin lt xmin then xmin=tmin
if tmax gt xmax then xmax=tmax
tmin=min(qrad[0,*]*mp,max=tmax)
if tmin lt xmin then xmin=tmin
if tmax gt xmax then xmax=tmax
tmax=max(abs(dnadv[2,*]*mp))
if tmax gt xmax then xmax=tmax
tmin=-tmax
if tmin lt xmin then xmin=tmin
tmin=min(trate[0,*]*mp,max=tmax)
if tmin lt xmin then xmin=tmin
if tmax gt xmax then xmax=tmax
plot,[0,0],[p0,p0],/ylog,yr=[p0,min(pres)],xr=[xmin,xmax],$
     xtit='Heating rate [W]',ytit='Pressure [bar]'
loadct,24,/silent
oplot,mp*adheat[0,*],pres,color=250
oplot,mp*adheat[1,*],pres,color=250,linestyle=2
oplot,mp*adheat[2,*],pres,color=250,linestyle=1
oplot,mp*vadv[0,*],pres,color=90
oplot,mp*vadv[1,*],pres,color=90,linestyle=2
oplot,mp*vadv[2,*],pres,color=90,linestyle=1
oplot,mp*dnadv[2,*],pres,color=140
oplot,-mp*dnadv[2,*],pres,color=140,linestyle=2
oplot,mp*qrad[0,*],pres,color=125
oplot,mp*qrad[1,*],pres,color=125,linestyle=2
oplot,mp*qrad[2,*],pres,color=125,linestyle=1
loadct,0,/silent
oplot,mp*trate[0,*],pres
oplot,mp*trate[1,*],pres,linestyle=2
oplot,mp*trate[2,*],pres,linestyle=1
al_legend,['On the dayside','On the nightside','Global total'],$
          linestyle=[0,2,1],charsize=0.75
loadct,24,/silent
al_legend,['Total rate','Radiative heating','Horizontal advection','Vertical advection','Adiabatic expansion'],$
          textcolors=[190,125,140,90,250],/right,charsize=0.75,box=0
loadct,0,/silent

if keyword_set(psfile) then begin
  !x.thick=1
  !y.thick=1
  !p.thick=1
  !p.charsize=1
  !p.multi=0
  psclose
endif

end
