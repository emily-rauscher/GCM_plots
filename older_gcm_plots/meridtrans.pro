pro meridtrans,oom=oom,p0=p0,grav=grav,gascon=gascon

; calculates the zonally averaged meridional transport, as a function
; of pressure and latitude, from fort.26
; transport = < v (c_p T + gz) >

  if not keyword_set(oom) then stop
  if not keyword_set(p0) then stop ; p0 in bar
  if not keyword_set(grav) then stop ; grav in m/s^2
  if not keyword_set(gascon) then stop ; gascon in J/kg/K
  cp=0.286*gascon

  OpenR, 17, 'fort.26'
    ReadF, 17, nlat,nlon,nlev	  
    xy=fltarr(6,nlat,nlon,nlev)   
    ReadF, 17, xy
  Close, 17
  lon=reform(xy[0,*,*,*]) & lat=reform(xy[1,*,*,*]) 
  u=reform(xy[3,*,*,*]) & v=reform(xy[4,*,*,*]) 
  temp=reform(xy[5,*,*,*])

  sigma=get_sigma(oom,nlev)
  pres=p0*sigma
  dP=fltarr(nlev)
  for ilev=1,nlev-2 do dP[ilev]=0.5*(pres[ilev+1]-pres[ilev-1]) ; to match swdp
  dP[0]=0.5*(pres[1])
  dP[nlev-1]=0.5*(p0-pres[nlev-2])

  ; could also read in fort.50 to include surface pressure variation

  tgr=max(temp[*,*,nlev-1])
  ;create array to hold heights (z=0 at p=p0)
  z=fltarr(nlat,nlon,nlev)      ; units: m
  ;set altitude of first level (up from base=p0, where T=TGR)
  z[*,*,nlev-1]=(gascon/grav)*0.5*(temp[*,*,nlev-1]+TGR)*alog(1./sigma[nlev-1])
  ;integrate hydrostatic to solve for higher levels
  for i=nlev-2,0,-1 do begin
     z[*,*,i]=z[*,*,i+1]+(gascon/grav)*0.5*(temp[*,*,i]+temp[*,*,i+1])*alog(sigma[i+1]/sigma[i])
  endfor
  
  trans = fltarr(nlat,nlev)
  itrans = fltarr(nlat,nlev)
  flat=fltarr(nlat,nlev)
  for ilat=0,nlat-1 do begin
     flat[ilat,*]=lat[ilat,0,0]
     for ilev=0,nlev-1 do begin
        trans[ilat,ilev] = total(v[ilat,*,ilev] $
                                 *(cp*temp[ilat,*,ilev] + grav*z[ilat,*,ilev]) $
                                )/nlon ; units: m^3/s^3
        if ilev ne 0 then $
           itrans[ilat,ilev] = itrans[ilat,ilev-1]+trans[ilat,ilev]*dP[ilev]*1.e5/grav $
        else itrans[ilat,ilev]=trans[ilat,ilev]*pres[ilev]*1.e5/grav ; units: Rp*W/m^2
     endfor
  endfor
  flev=fltarr(nlat,nlev)
  for ilev=0,nlev-1 do flev[*,ilev]=pres[ilev]

  window,0
  ymin=min(flev,max=ymax)
;  qmin=min(trans,max=qmax)
  qmax=max(abs(trans))
  qmin=-qmax
  nlevels=25
  step=(qmax-qmin)/nlevels
  qlevels=indgen(nlevels)*step+qmin
  loadct,10,/silent
  contour, trans,flat,flev,/cell_fill,LEVELS=qlevels,/ystyle,/xstyle,charsize=2,$
           yr=[ymax,ymin],/ylog,xtitle='Latitude [degrees]',ytitle='Pressure [bar]',$
           max_value=qmax,min_value=qmin,xr=[-90,90],title='v (c_p T + g z)'
  colorbar,range=[qmin,qmax],charsize=1.5,format='(e9.2)',position=[0.165,0.15,0.9,0.2],/top
  loadct,0,/silent

  window,1
  qmax=max(abs(itrans))
  qmin=-qmax
  nlevels=25
  step=(qmax-qmin)/nlevels
  qlevels=indgen(nlevels)*step+qmin
  loadct,10,/silent
  contour, itrans,flat,flev,/cell_fill,LEVELS=qlevels,/ystyle,/xstyle,charsize=2,$
           yr=[ymax,ymin],/ylog,xtitle='Latitude [degrees]',ytitle='Pressure [bar]',$
           max_value=qmax,min_value=qmin,xr=[-90,90],title='int v(cpT+gz) dP/g'
  colorbar,range=[qmin,qmax],charsize=1.5,format='(e9.2)',position=[0.165,0.875,0.9,0.925]
  loadct,0,/silent

end
