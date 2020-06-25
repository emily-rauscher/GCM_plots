pro hyperdiss

  OpenR, 17, 'fort.26'
  ReadF, 17, nlat,nlon,nlev	  
;	  print,nlat,nlon,nlev 
	  xy=fltarr(6,nlat,nlon,nlev)   
  ReadF, 17, xy
  Close, 17

nl=3
radea=9.44e7   ;m
WW = 2.06e-5  ;rad/sec
GASCON = 4593. ;J/kg/K
p0=220.612  ;bar
OOM = 5.34363
GA = 9.42   ;m/s^2
nu=8.54d47  ;m^8 s^-1

;    lon=reform(xy[0,*,*,nl]) & lat=reform(xy[1,*,*,nl]) 
;    u=reform(xy[3,*,*,nl]) & v=reform(xy[4,*,*,nl]) 
;    temp=reform(xy[5,*,*,nl])
    lon=reform(xy[0,*,*,*]) & lat=reform(xy[1,*,*,*]) & u=reform(xy[3,*,*,*]) & v=reform(xy[4,*,*,*]) & temp=reform(xy[5,*,*,*])

sigma=get_sigma(oom,nlev)

;mass elements
dm=fltarr(nlat,nlon,nlev)
for ilev=0,nlev-1 do begin
   if ilev gt 0 and ilev lt nlev-1 then dA=0.5*p0*(sigma[ilev+1]-sigma[ilev-1])*(2.*!pi/nlon)*(!pi/nlat)*(radea^2/ga)
   if ilev eq 0 then dA=p0*(sigma[0])*(2.*!pi/nlon)*(!pi/nlat)*(radea^2/ga)
   if ilev eq nlev-1 then dA=p0*(1.-sigma[nlev-1])*(2.*!pi/nlon)*(!pi/nlat)*(radea^2/ga)
   for ilat=0,nlat-1 do begin
      dm[ilat,*,ilev]=replicate(cos(lat[ilat,0,0]*!pi/180.)*dA,nlon)
   endfor
endfor
dm*=1.e5  ;for units conversion to kg

dtemp=dblarr(nlat,nlon,nlev)
dxtemp=dblarr(nlat,nlon,nlev)
dytemp=dblarr(nlat,nlon,nlev)

for l=0,nlev-1 do begin
   for i=0,nlat-1 do begin
      x=(radea*cos(lat[i,0,0]*!pi/180.)*lon[i,*,0]*!pi/180.)
      dxtemp[i,*,l]=deriv(x,temp[i,*,l])
   endfor
   ;print,minmax(dxtemp),strcompress('K/m^1')
   for o=2,8 do begin
      for i=0,nlat-1 do begin
         x=reform(radea*cos(lat[i,0,0]*!pi/180.)*lon[i,*,0]*!pi/180.)
         dxtemp[i,*,l]=deriv(x,dxtemp[i,*,l])
      endfor
   ;print,minmax(dxtemp),strcompress('K/m^'+string(o))
   endfor

   for j=0,nlon-1 do begin
      y=radea*lat[*,j,0]*!pi/180.
      dytemp[*,j,l]=deriv(y,temp[*,j,l])
   endfor
   ;print,minmax(dytemp),strcompress('K/m^1')
   for o=2,8 do begin
      for j=0,nlon-1 do begin
         y=radea*lat[*,j,0]*!pi/180.
         dytemp[*,j,l]=deriv(y,dytemp[*,j,l])
      endfor
   ;print,minmax(dytemp),strcompress('K/m^'+string(o))
   endfor
endfor


dtemp=(dxtemp+dytemp)*nu

print,'range: ',minmax(dtemp),'K/s'
print,minmax(dtemp*2.*!pi/ww),'K/day'
print,'total: ',total(dtemp),' K/s'
print,total(dtemp)*2.*!pi/ww,' K/day'

dtemp=dtemp*dm*gascon  ;W
print,'range: ',minmax(dtemp),'W'
print,'total: ',total(dtemp),' W'

stop

  flon=reform(lon[*,*,0])*360./lon[0,nlon-1,0]
  flat=reform(lat[*,*,0])*90./lat[nlat-1,0,0]
vert=fltarr(nlat,nlon,nlev)
for i=0,nlat-1 do begin
   for j=0,nlon-1 do begin
      vert[i,j,*]=p0*sigma
   endfor
endfor

loadct,4

window,0

umax=max(dtemp[*,*,nl],min=umin)
;absmax=max([abs(umax),abs(umin)])
range=umax-umin
;maxrange=2.*absmax  ;from -absmax to +absmax
;rat=range/maxrange
;cbottom=(absmax+umin)/maxrange*255.
;nlevels=44.*rat+1.  ;so that delta color always = 255./44.
nlevels=55
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
MAP_SET, /Miller_CYLINDRICAL,0.,180.,/ISOTROPIC;,position=[0.1,0.5,0.9,0.90]
contour,reform(dtemp[*,*,nl]),flon,flat,/cell_fill,min_value=umin,max_value=umax,levels=mylevels,/overplot,/closed;,c_colors=findgen(nlevels)*255.*rat/(nlevels-1)+cbottom
map_grid,/label,charsize=1,latlab=355.;,color=255
colorbar,range=[umin,umax],format='(e8.1)',divisions=4,charsize=1.25;,position=[0.15,0.52,0.85,0.54],bottom=cbottom,ncolors=256*rat

window,1

umax=max(dtemp,min=umin)
;absmax=max([abs(umax),abs(umin)])
;range=umax-umin
;maxrange=2.*absmax  ;from -absmax to +absmax
;rat=range/maxrange
;cbottom=(absmax+umin)/maxrange*255.
;nlevels=24.*rat+1.  ;so that delta color always = 255./24.
nlevels=25
step=(umax-umin)/nlevels
mylevels=indgen(nlevels)*step+umin
contour,reform(dtemp[24,*,*]),reform(lon[24,*,*]),reform(vert[24,*,*]),/cell_fill,yr=[max(p0*sigma),min(p0*sigma)],/ystyle,/ylog,levels=mylevels,ytitle='Pressure [bar]',min_value=umin,max_value=umax,/xstyle,xtickformat="(A1)";,c_colors=findgen(nlevels)*255.*rat/(nlevels-1)+cbottom;,position=[0.16,0.5,0.99,0.9]
colorbar,range=[umin,umax],format='(e8.1)',divisions=4,charsize=1.25,color=255;,position=[0.23,0.55,0.93,0.58],bottom=cbottom,ncolors=256*rat


loadct,0
end
