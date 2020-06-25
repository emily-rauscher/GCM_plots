pro losvel,nl,psfile=psfile,velrange=velrange,vfrac=vfrac

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

  OpenR, 17, 'fort.26'
  ReadF, 17, nlat,nlon,nlev	  
;	  print,nlat,nlon,nlev 
	  xy=fltarr(6,nlat,nlon,nlev)   
  ReadF, 17, xy
  Close, 17

  lon=reform(xy[0,*,*,nl]) & lat=reform(xy[1,*,*,nl])
  u=reform(xy[3,*,*,nl]) & v=reform(xy[4,*,*,nl]) & temp=reform(xy[5,*,*,nl])

  losv=fltarr(nlat,nlon)  ;(not including rotation)
  for ilat=0,nlat-1 do begin
     for ilon=0,nlon-1 do begin
        losv[ilat,ilon]=u[ilat,ilon]*sin(lon[ilat,ilon]*!pi/180.)*cos(lat[ilat,ilon]*!pi/180.) $
                        +v[ilat,ilon]*cos(lon[ilat,ilon]*!pi/180.)*sin(lat[ilat,ilon]*!pi/180.)
     endfor
  endfor

flon=fltarr(nlat,nlon+1)
flon[*,0:nlon-1]=lon
flon[*,nlon]=reform(lon[*,0]+360.)
flat=fltarr(nlat,nlon+1)
flat[*,0:nlon-1]=lat
flat[*,nlon]=reform(lat[*,0])

qtp=fltarr(nlat,nlon+1)
qtp[*,0:nlon-1]=losv
qtp[*,nlon]=reform(losv[*,0])

if not keyword_set(velrange) then begin
   qmin=min(qtp,max=qmax)
endif else begin
   qmin=velrange[0]
   qmax=velrange[1]
endelse

absmax=max([abs(qmax),abs(qmin)])
range=qmax-qmin
maxrange=2.*absmax  ;from -absmax to +absmax
rat=range/maxrange
cbottom=(absmax+qmin)/maxrange*255.
nlevels=44.*rat+1.  ;so that delta color always = 255./44.
step=(qmax-qmin)/nlevels
mylevels=indgen(nlevels)*step+qmin

qtemp=fltarr(nlat,nlon+1)
qtemp[*,0:nlon-1]=temp
qtemp[*,nlon]=reform(temp[*,0])

if qmin ne qmax then begin
   loadct,6,/silent
   MAP_SET, /Miller_CYLINDRICAL,0,0,/ISOTROPIC
   contour, qtp, flon,flat ,/overplot ,/cell_fill,/closed ,NLEVELS=45  , min_value=qmin, max_value=qmax,$
            c_colors=findgen(nlevels)*255.*rat/(nlevels-1)+cbottom,levels=mylevels
;   on=where(abs(flat[*,0:nlon-1]) ne max(flat))
;   partvelvec,(u[*,0:nlon-1])[on],(v[*,0:nlon-1])[on],(flon[*,0:nlon-1])[on],(flat[*,0:nlon-1])[on], /over,fraction=vfrac,color=0
   contour,qtemp,flon,flat,/overplot,/closed,nlevels=20,/follow
   map_grid, /label,charsize=1  ;,color=0
   colorbar,  position=[0.1,0.07,0.90,0.10], range=[qmin,qmax], format='(e9.2)', charsize=1.5,divisions=5,$
              bottom=cbottom,ncolors=256*rat
   loadct,0,/silent
endif

if keyword_set(psfile) then psclose

end
