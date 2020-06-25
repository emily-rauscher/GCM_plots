pro KE_map,p0,levplot,latoff,lonoff

; plots m*v^2/2 at levplot level
; reads in KE map from fort.54
; needs p0 to normalize pressure (p0=pressure at z=0)
;          (enter p0 in bar, will convert to Pa)

p0*=1.e5

openr,17,'fort.54'
readf,17,nlat,nlon,nlev
if levplot ge nlev then print,'ERROR: levplot > nlev'
xy=fltarr(4,nlat,nlon,nlev)
readf,17,xy
close,17

lon=reform(xy[0,*,*,levplot]) & lat=reform(xy[1,*,*,levplot])
ke=reform(xy[3,*,*,levplot])

ke*=p0

flon=lon*360./lon(0,nlon-1)
flat=lat*90./lat(0,nlat-1)
if not keyword_set(latoff) then latoff=0.
if not keyword_set(lonoff) then lonoff=0.

!x.style=1
!p.font=-1
!p.charsize=1.5

print,minmax(ke)

loadct,4
MAP_SET, /CYLINDRICAL,latoff,lonoff,/ISOTROPIC, TITLE='Kinetic Energy Map'

if min(ke) ne max(ke) then colorbar,position=[0.1,0.08,0.90,0.12],range=[min(ke),max(ke)],$
                                      format='(E8.1)',charsize=1.5

contour, ke,flon,flat,/overplot,/cell_fill, nlevels=55

;temperature contours
;contour,temp,flon,flat,/overplot,nlevels=5,color=255

;partvelvec,u,v,flon,flat,/over,fraction=0.9,color=0

map_grid,/label

print,'L',levplot+1

end
