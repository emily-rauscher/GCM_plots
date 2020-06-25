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


  xy = fltarr(3,96,192)
;  OpenR, 17, 'hmass_f16.dat'
  OpenR, 17, 'hmass_f16_wf.dat'
;  OpenR, 17, 'pv_f16.dat'
  ReadF, 17, xy
  Close, 17
  lon=reform(xy[0,*,*]) & lat=reform(xy[1,*,*]) & temp=reform(xy[2,*,*])
 
$ Convert radian to degrees and 
$ correct for grid not fitting 90 and 180 degrees exactly 

  lonfac=2.0*3.1415926536/lon(0,191)*57.2957795132
  lon=lon*lonfac
  latfac=3.1415926536/2.0/lat(0,0)*57.2957795132
  lat=lat*latfac

   MAP_SET, /STEREO, 0, 0, /ISOTROPIC, TITLE = 'DAY-SIDE TEMPERATURE' ;,color=0 

$ This one fills the horizon with black, to avoid hole
;  MAP_SET, /STEREO, 90, 0, /ISOTROPIC,/HORIZON, E_HORIZON={FILL:1}, TITLE = 'Temperature Map', color=0 


    loadct, 4
  contour, temp,lon,lat, /overplot,/cell_fill ,NLEVELS=55  ;, min_value=1000, max_value=2000
    ;loadct,2
      contour, temp,lon,lat,/overplot, NLEVELS=10, /follow,color=0

      map_grid, /label ;,color=0

	print, ' Min and Max Temp: ',min(temp),max(temp)


  set_plot, 'x'

END

