$
$
$
PRO tmap

$ Screen output
  device,true_color=24,decomposed=0,retain=2
  loadct, 4 
  set_plot, 'x'

$ Postscript output
;   device,true_color=24,decomposed=0,retain=2
;  loadct,4
;  set_plot,'ps', /COPY
;  device, filename='tt.ps',bits_per_pixel=24,/color

;;;;;!P.BACKGROUND=255


    x_size_window=800
    y_size_window=600
 window, 0, xsize=x_size_window, ysize=y_size_window, retain=2

$       ntot=32*64*5
       xy = fltarr(6,32,64,5)
;       uvel 
;       vvel
;       Temp
  OpenR, 17, 'fort.26'
  ReadF, 17, nlat,nlon,nlev	  
	  print,nlat,nlon,nlev 
  ReadF, 17, xy
  Close, 17

$  lev=fltarr(ntot)
$  lon=fltarr(ntot)
$  lat=fltarr(ntot)
$  uvel=fltarr(ntot)
$  vvel=fltarr(ntot)
$  Temp=fltarr(ntot)


  nl=4

	  lon=reform(xy[0,*,*,nl]) & lat=reform(xy[1,*,*,nl]) & u=reform(xy[3,*,*,nl]) & v=reform(xy[4,*,*,nl]) & temp=reform(xy[5,*,*,nl])
 
  qtp=temp
$  qtp=u
$  qtp=v

$$ correct for grid not fitting 90 and 180 degrees exactly 
  lon=lon*360./lon(0,nlon-1)
  lat=lat*90./lat(0,nlat-1)


   MAP_SET, /STEREO, 0, 0, /ISOTROPIC, TITLE = 'DAY-SIDE TEMPERATURE' ;,color=0 

$$ This one fills the horizon with black, to avoid hole
$;  MAP_SET, /STEREO, 90, 0, /ISOTROPIC,/HORIZON, E_HORIZON={FILL:1}, TITLE = 'Temperature Map', color=0 


    loadct, 4
  contour, qtp,lon,lat, /overplot,/cell_fill ,NLEVELS=55  ;, min_value=1000, max_value=2000
    ;loadct,2
      contour, qtp,lon,lat,/overplot, NLEVELS=10, /follow,color=0

      map_grid, /label ;,color=0

	print, ' Min and Max Temp: ',min(temp),max(temp)


  set_plot, 'x'

END

