pro plot_inversion

;  device,true_color=24,decomposed=0,retain=2

  OpenR, 17, 'fort.26'
  ReadF, 17, nlat,nlon,nlev	  
;	  print,nlat,nlon,nlev 
	  xy=fltarr(6,nlat,nlon,nlev)   
  ReadF, 17, xy
  Close, 17

  lon=reform(xy[0,*,*,0]) & lat=reform(xy[1,*,*,0])

  temp=reform(xy[5,*,*,*])

  flon=lon*360./lon(0,nlon-1)
  flat=lat*90./lat(0,nlat-1)

  ;find inversions, calc as max(above 100 mbar)-min(100 mbar - 300 mbar)
  ;                            (levs=0-11)         (levs=12-14)
  inv=fltarr(nlat,nlon)

  for i=0,nlon-1 do begin
     for j=0,nlat-1 do begin
        inv[j,i]=max(temp[j,i,0:11])-min(temp[j,i,12:14])
     endfor
  endfor

MAP_SET, /Miller_CYLINDRICAL,0,0,/ISOTROPIC;, TITLE=' Temperature and Velocity Map'
loadct,4
prange=[min(inv),max(inv)]
contour,inv,flon,flat,/overplot,/cell_fill,/closed,nlevels=45;,max_value=prange[1],min_value=0
contour,inv,flon,flat,/overplot,levels=[0.],color=5,thick=4
contour,inv,flon,flat,/overplot,levels=[300.],color=255,thick=4
colorbar,position=[0.1,0.07,0.9,0.1],format='(i5)',charsize=2,range=[prange[0],prange[1]]
map_grid,/label,color=255

loadct,0

end
