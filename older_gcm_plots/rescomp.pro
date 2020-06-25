function rescomp,bfield,tdrag_min,p0=p0,oom=oom,radea=radea,ga=ga,metals=metals

gascon=3523.  ; J/kg/K
if not keyword_set(p0) then p0=1.e2
if not keyword_set(oom) then oom=5.
if not keyword_set(radea) then radea = 1.e8  ; m
if not keyword_set(ga) then ga=8.  ; m/s
print,'         CHECK: gascon, radea, ga, p0, oom:',gascon,radea,ga,p0,oom
print,'  DOUBLE-CHECK: bfield is',bfield
;tdrag_min = 0.005 * (2.*!pi)/ww

  OpenR, 17, 'fort.26'
  ReadF, 17, nlat,nlon,nlev	  
	  xy=fltarr(6,nlat,nlon,nlev)   
  ReadF, 17, xy
  Close, 17

  openr,18,'fort.50'
  readf,18,nlat,nlon
  ab=fltarr(3,nlat,nlon)
  readf,18,ab
  close,18

  lon=reform(xy[0,*,*,*]) & lat=reform(xy[1,*,*,*]) 
  u=reform(xy[3,*,*,*]) & v=reform(xy[4,*,*,*]) & temp=reform(xy[5,*,*,*])
  sp=reform(ab[2,*,*])
  sp=(sp+1.)*p0

  sigma=get_sigma(oom,nlev)
  pres=p0*sigma
  dP=fltarr(nlev)
  for ilev=1,nlev-2 do dP[ilev]=0.5*(pres[ilev+1]-pres[ilev-1]) ; to match swdp
  dP[0]=0.5*(pres[1])
  dP[nlev-1]=0.5*(p0-pres[nlev-2])

  ; total (heating, mom change) per radian
  latvals=fltarr(4,nlat)  ; 0: heating/rad, 1: u/tmag/rad, 2: latitudes, 3: rms u
  latvals[2,*]=reform(lat[*,0,0])
  latvals[3,*]=sqrt(total(total(u^2,3),2)/nlev/nlon)

  ke=u^2
  for ilev=0,nlev-1 do begin
     for ilat=0,nlat-1 do begin
        ang_f=abs(cos(!pi/2.-xy[1,ilat,0,0]*!pi/180.))
        rho=reform(sp[ilat,*]*sigma[ilev]/gascon/temp[ilat,*,ilev]*1.e2) ; g/cm^3
        tdrag=tmag(rho,reform(temp[ilat,*,ilev]),bfield,metals=metals)/ang_f   ; s
        for ilon=0,nlon-1 do tdrag[ilon]=max([tdrag[ilon],tdrag_min])
        latvals[0,ilat] += total(reform(ke[ilat,*,ilev])/tdrag) $ ; m^2/s^3  (or W/kg)
                        *cos(lat[ilat,0,0]*!pi/180.)*(!pi/nlat)*(2.*!pi/nlon) $
                        *(radea^2/ga)*dP[ilev]*1.e5 ; W
        latvals[1,ilat] += total(reform(u[ilat,*,ilev])/tdrag) $ ; m/s^2
                        *cos(lat[ilat,0,0]*!pi/180.)*(!pi/nlat)*(2.*!pi/nlon) $
                        *(radea^2/ga)*dP[ilev]*1.e5 ; kg m/s^2
     endfor
  endfor

  midlats=reform(lat[0:nlat-2,0,0]+lat[1:nlat-1,0,0])/2.
  deltlats=[90.,midlats[0:nlat-2]]-[midlats[0:nlat-2],-90.]
  deltlats*=!pi/180.
  latvals[0,*]/=deltlats
  latvals[1,*]/=deltlats

  plot,latvals[2,*],latvals[3,*],/xstyle,ytit='RMS zonal wind [m/s]',$
       position=[0.15,0.67,0.95,0.97],xticknam=replicate(' ',3),/ynoz
  plot,latvals[2,*],latvals[1,*],/xstyle,ytit='Momentum change [kg*m/s^2/rad]',$
       position=[0.15,0.37,0.95,0.67],/noerase,xticknam=replicate(' ',3)
  plot,latvals[2,*],latvals[0,*],/xstyle,ytit='Power [W/rad]',$
       position=[0.15,0.07,0.95,0.37],/noerase,xtit='Latitude [degrees]'
;  oplot,[-90,90],[0,0],linestyle=1

  return,latvals

end
