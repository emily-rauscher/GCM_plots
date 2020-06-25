PRO moviemaker,p0,oom
;  To use the output from igcm to make frames for movies in lat, lon, level.
;  WARNING: Plots each frame in an IDL window; can't use if can't use
;  X windows.
;
;  To turn frames into an animated gif, use imagemagick: 
;    convert -delay 30 frame_lon0900.00_L*.png lon90movie.gif
;    convert -delay 10 frame_lon*_L14.png L14movie.gif
;    A big movie of EVERYTHING (will be BIG and computer will run SLOW): 
;      convert -delay 10 frame*.png bigmovie.gif  ;NOPE, makes them in
;      L, lon order (rather than rotating and then going down a level)
;  To make a quicktime movie, use Sequimago:
;    either drag and drop a list of files onto the icon, or
;    open the script and choose a list of files
;    (2 frames/sec is good for example the lon90movie.mov)
;    (4 frames/sec is good for example the L13movie.mov)


  device,true_color=24,decomposed=0,retain=2
  window, 1, xsize=800, ysize=900, retain=2

  loadct, 4
  vfrac=0.8
;!p.font=0
;!p.charsize=1.5
            
;---------orthographic projections, lon range, pressure range, one snapshot

  OpenR, 17, 'fort.26'
  ReadF, 17, nlat,nlon,nlev	  
;	  print,nlat,nlon,nlev 
	  xy=fltarr(6,nlat,nlon,nlev)   
  ReadF, 17, xy
  Close, 17

  lon=reform(xy[0,*,*,*]) & lat=reform(xy[1,*,*,*]) 
  u=reform(xy[3,*,*,*]) & v=reform(xy[4,*,*,*]) & temp=reform(xy[5,*,*,*])
  sigma=get_sigma(oom,nlev)
  plevels=p0*sigma

  flon=lon*360./lon(0,nlon-1)
  flat=lat*90./lat(0,nlat-1)

;begin looping over: pressure, lon
  for j=0,nlev-1 do begin

     ;title for each level
     ;if not keyword_set(plevels) then mtitle=strcompress('Level '+string(j+1,format='(I2.2)')) $
     mtitle=strcompress(string(plevels[j],format='(g7.2)')+' bar level')

     ;normalize colorscale for each level
     crange=[min(temp[*,*,j]),max(temp[*,*,j])]
  
     for i=0,nlon-1 do begin

        if lon[0,i,j] lt 10. then begin
           fpart=strcompress('00'+string(lon[0,i,j],format='(f6.2)'),/remove_all)
        endif else begin
           if lon[0,i,j] lt 100. then fpart=strcompress('0'+string(lon[0,i,j],format='(f6.2)'),/remove_all) $
           else fpart=string(lon[0,i,j],format='(f6.2)')
        endelse
                            
        filename=strcompress('~/research/images/movies/frames/frame_lon'+fpart+$
                             '_L'+string(j+1,format='(I2.2)')+'.png',/remove_all)
                           
        MAP_SET, /ortho,0,lon[0,i,j],/ISOTROPIC, TITLE=mtitle,charsize=2.5,position=[.05,.125,.95,.925]
                ;,/advance (for p.multi)
        contour,temp[*,*,j],flon[*,*,j],flat[*,*,j],$
                /overplot,/cell_fill,/closed,nlevels=45,min_value=crange[0],max_value=crange[1]
        partvelvec,u[*,*,j],v[*,*,j],flon[*,*,j],flat[*,*,j],/over,fraction=vfrac,color=0
        map_grid, /label ,color=255,charsize=2

        colorbar,range=crange,format='(i5)',position=[.1,.05,.9,.075],charsize=2,title='Temperature [K]'
           
        write_png,filename,tvrd(true=1)

        ;stop

     endfor
  endfor


;  MAP_SET, /STEREO, latoff, lonoff, /ISOTROPIC, TITLE = 'T/U/V MAP' ;,color=0 
;MAP_SET, /Miller_CYLINDRICAL,latoff,lonoff,/ISOTROPIC;, TITLE=' Temperature and Velocity Map'
;; This one fills the horizon with black, to avoid hole
;  MAP_SET, /STEREO, 90, 0, /ISOTROPIC,/HORIZON, E_HORIZON={FILL:1}, TITLE = 'Temperature Map', color=0 

;  if keyword_set(prange) then contour, qtp,flon,flat, /overplot,/cell_fill,/closed ,NLEVELS=45, min_value=prange[0], max_value=prange[1] $
;     else contour, qtp,flon,flat, /overplot,/cell_fill,/closed ,NLEVELS=45 ;, min_value=1000, max_value=2000
;  contour, qtp,flon,flat,/overplot, NLEVELS=10, /follow,color=0

;  if min(qtp) ne max(qtp) then begin
;     if not keyword_set(prange) then $
;        colorbar,position=[0.1,0.07,0.90,0.1],range=[min(qtp),max(qtp)],format='(i5)',charsize=2 $
;     else $
;        colorbar,  position=[0.1,0.07,0.90,0.1], range=prange, format='(i5)', charsize=2
;  endif

; Plot velocity vector     
;  if not keyword_set(vfrac) then vfrac=0.9
;  partvelvec,u,v,flon,flat, /over,fraction=vfrac, color=0


;   print, ' Min and Max Temperatures (K):            ',min(temp),max(temp)
;   print, ' Min and Max Zonal Wind Speed (m/s):      ',min(u),max(u)
;   print, ' Min and Max Meridional Wind Speed (m/s): ',min(v),max(v)
;   print, ' Min and Max Sound Speed (m/s):           ',72.74*sqrt(min(temp)),72.74*sqrt(max(temp))
;   print, strcompress('L'+string(levplot+1),/remove_all)


;stop

loadct,0

end
