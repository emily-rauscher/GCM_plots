pro orthoproj,nl,latoff,lonoff,rot=rot,vfrac=vfrac,grid=grid,folr=folr,$
              psfile=psfile,frange=frange,moviefiles=moviefiles,maxwind=maxwind
;
; plots temperature or flux map in orthographic projection
;
; keywords:
;   nl: level to plot (for temperature map, value doesn't matter
;       if flux plotted)
;   latoff: planet latitude at center of projection
;   lonoff: planet longitude at center of projection
; optional keywords:
;   rot: angle through which to rotate north axis
;   vfrac: fraction of velocity vectors to plot
;   grid: set to overplot grid lines
;   folr: set to plot flux instead of temperature
;   psfile: set to produce plot.eps, or set as string to produce "string".eps
;   frange: set to [min,max] values for contours
;   moviefiles: array of file names to create frames for movie
;               will be saved in subdirectory: frames/ortho_movie_'+*+'.png
;               Can then make movie by:  convert -delay 10 ortho_*.png movie.gif  
;   maxwind: wind speed of longest vector


if not keyword_set(rot) then rot=0

if not keyword_set(moviefiles) then begin

   if not keyword_set(folr) then begin
      loadct,4,/silent
      
      OpenR, 17, 'fort.26'
      ReadF, 17, nlat,nlon,nlev	  
;	  print,nlat,nlon,nlev 
      xy=fltarr(6,nlat,nlon,nlev)   
      ReadF, 17, xy
      Close, 17

      lon=reform(xy[0,*,*,nl]) & lat=reform(xy[1,*,*,nl])
      u=reform(xy[3,*,*,nl]) & v=reform(xy[4,*,*,nl]) & temp=reform(xy[5,*,*,nl])

      qtp=fltarr(nlat,nlon+1)
      qtp[*,0:nlon-1]=temp
      qtp[*,nlon]=reform(temp[*,0])

      psx=23.5
      cpos=[0.865, 0.10, 0.895, 0.90]
      cbx=0.735
      cby=0.93
      cbartit='Temperature [K]'

   endif else begin
      loadct,5,/silent
      
      OpenR, 17, 'fort.64'
      ReadF, 17, nlat,nlon,nlev	  
      xy=fltarr(3,nlat,nlon)   
      ReadF, 17, xy
      Close, 17
      
      lon=reform(xy[0,*,*]) & lat=reform(xy[1,*,*]) & olr=reform(xy[2,*,*])
      
      qtp=fltarr(nlat,nlon+1)
      qtp[*,0:nlon-1]=olr
      qtp[*,nlon]=reform(olr[*,0])

      psx=24
      cbx=0.79
      cby=0.93
      cbartit='Flux [W/m^2]'
      cpos=[0.855, 0.10, 0.885, 0.90]

   endelse

   if not keyword_set(psfile) then begin
; Screen output
      device,true_color=24,decomposed=0,retain=2
;  set_plot, 'x'
      window, 0, xsize=x_size_window, ysize=y_size_window, retain=2
      !p.font=-1
   endif else begin
      if (size(psfile))[1] eq 2 then filename='plot.eps' else filename=psfile
      psopen,filename,/enc,/color,xsize=psx,ysize=20,bcolor=1 ;,bits_per_pixel=24
      !p.font=0
   endelse

   flon=fltarr(nlat,nlon+1)
   flon[*,0:nlon-1]=lon
   flon[*,nlon]=reform(lon[*,0]+360.)
   flat=fltarr(nlat,nlon+1)
   flat[*,0:nlon-1]=lat
   flat[*,nlon]=reform(lat[*,0])

   if keyword_set(frange) then begin
      qmin=frange[0]
      qmax=frange[1]
   endif else qmin=min(qtp,max=qmax)
   print,qmin,qmax

   if keyword_set(psfile) then begin
      Polyfill, [1,1,0,0,1], [1,0,0,1,1], /NORMAL, COLOR=1
      MAP_SET,latoff,lonoff,rot,/ortho,/ISOTROPIC,/noerase,position=[0.0,0.0,0.85,1.]
   endif else MAP_SET,latoff,lonoff,rot,/ortho,/ISOTROPIC ;, TITLE=' Temperature and Velocity Map'

   nlevels=45.
   cbottom=55
   step=(qmax-qmin)/nlevels
   mylevels=indgen(nlevels)*step+qmin
   ccolors=findgen(nlevels)*(255.-cbottom)/(nlevels-1)+cbottom
   
   contour,qtp,flon,flat,/overplot,/cell_fill,/closed,c_colors=ccolors,levels=mylevels

   colorbar,range=[qmin,qmax],format='(i5)',charsize=2,color=255,$
            /vertical,/right,position=cpos,bottom=cbottom,$
            ncolors=255-cbottom ;,title='Temperature [K]'
   xyouts,cbx,cby,cbartit,charsize=2,/normal,color=255

   if not keyword_set(folr) then partvelvec,u,v,lon,lat,/over,fraction=vfrac,color=0

   if keyword_set(grid) then map_grid,/label,color=0

   loadct,0,/silent

   if keyword_set(psfile) then begin
      psclose
      !p.font=-1
;   spawn,'gs -r300 -dEPSCrop -dTextAlphaBits=4 -sDEVICE=png16m -sOutputFile=plot.png -dBATCH -dNOPAUSE plot.eps'
;   spawn,'convert plot.png eps3:plot.eps'
   endif

endif else begin

   if not keyword_set(folr) then begin
      loadct,4,/silent
      psx=847
      cpos=[0.88, 0.12, 0.91, 0.92]
      cbx=0.8
      cby=0.95
      cbartit='Temperature [K]'
      
      OpenR, 17, moviefiles[0]
      ReadF, 17, nlat,nlon,nlev	  
      xy=fltarr(6,nlat,nlon,nlev)   
      ReadF, 17, xy
      Close, 17

      lon=reform(xy[0,*,*,nl]) & lat=reform(xy[1,*,*,nl])
      u=reform(xy[3,*,*,nl]) & v=reform(xy[4,*,*,nl]) & temp=reform(xy[5,*,*,nl])

      qtp=fltarr(nlat,nlon+1)
      qtp[*,0:nlon-1]=temp
      qtp[*,nlon]=reform(temp[*,0])
   endif else begin
      loadct,5,/silent
      psx=864
      cbx=0.79
      cby=0.93
      cbartit='Flux [W/m^2]'
      cpos=[0.855, 0.10, 0.885, 0.90]
      
      OpenR, 17, moviefiles[0]
      ReadF, 17, nlat,nlon,nlev	  
      xy=fltarr(3,nlat,nlon)   
      ReadF, 17, xy
      Close, 17
      
      lon=reform(xy[0,*,*]) & lat=reform(xy[1,*,*]) & olr=reform(xy[2,*,*])
      
      qtp=fltarr(nlat,nlon+1)
      qtp[*,0:nlon-1]=olr
      qtp[*,nlon]=reform(olr[*,0])
   endelse

   device,true_color=24,decomposed=0,retain=2
   window, 0, xsize=psx, ysize=720, retain=2

   flon=fltarr(nlat,nlon+1)
   flon[*,0:nlon-1]=lon
   flon[*,nlon]=reform(lon[*,0]+360.)
   flat=fltarr(nlat,nlon+1)
   flat[*,0:nlon-1]=lat
   flat[*,nlon]=reform(lat[*,0])

   if keyword_set(frange) then begin
      qmin=frange[0]
      qmax=frange[1]
   endif else qmin=min(qtp,max=qmax)

   nlevels=45.
   cbottom=55
   step=(qmax-qmin)/nlevels
   mylevels=indgen(nlevels)*step+qmin
   ccolors=findgen(nlevels)*(255.-cbottom)/(nlevels-1)+cbottom

   MAP_SET,latoff,lonoff,rot,/ortho,/ISOTROPIC,position=[0.0,0.0,0.85,1.],/noborder
   colorbar,range=[qmin,qmax],format='(i5)',charsize=2,color=255,$
            /vertical,/right,position=cpos,bottom=cbottom,ncolors=255-cbottom
   xyouts,cbx,cby,cbartit,charsize=2,/normal,color=255
   if not keyword_set(folr) then begin
      loadct,0,/silent
      ; calculate arrow length, coordinate conversion from map to normal
      astart=convert_coord(lonoff,latoff,/to_normal)
      toss = 0.05 * ( max(lon)-min(lon) > max(lat)-min(lat) )  ; length of vector in partvelvec
      aend = convert_coord(lonoff+toss,latoff,/to_normal) ; (some slight loss due to curving away from subobserver)
      arrow,cbx,0.02,cbx+abs(astart[1]-aend[1]),0.02,/normal,hsize=20
      xyouts,cbx,0.07,strcompress('Max wind speed:'),/normal,charsize=2
      xyouts,cbx,0.04,strcompress(string(maxwind)+' m/s'),/normal,charsize=2
      loadct,4,/silent
   endif

   for iframe=0,n_elements(moviefiles)-1 do begin
      OpenR, 17, moviefiles[iframe]
      ReadF, 17, nlat,nlon,nlev	  
      if keyword_set(folr) then begin
         xy=fltarr(3,nlat,nlon,nlev)   
         ReadF, 17, xy
         Close, 17
         qtp[*,0:nlon-1]=reform(xy[3,*,*])
         qtp[*,nlon]=reform(xy[3,*,0])
      endif else begin
         xy=fltarr(6,nlat,nlon,nlev)   
         ReadF, 17, xy
         Close, 17
         u=reform(xy[3,*,*,nl]) & v=reform(xy[4,*,*,nl])
         if not keyword_set(maxwind) then vlength=0.05 else begin
            toss=max(abs(u)) > max(abs(v))  ; maximum wind speed
            if toss gt maxwind then print,'Warning: wind > maxwind:',toss
            vlength=0.05*toss/maxwind
         endelse
         qtp[*,0:nlon-1]=reform(xy[5,*,*,nl])
         qtp[*,nlon]=reform(xy[5,*,0,nl])
      endelse
      print,min(qtp),max(qtp)
      contour,qtp,flon,flat,/overplot,/cell_fill,/closed,c_colors=ccolors,levels=mylevels
      if keyword_set(grid) then map_grid,color=0,charsize=2,glinethick=2,glinestyle=2,$
                                         londel=30,latdel=30;,/label
      if not keyword_set(folr) then partvelvec,u,v,lon,lat,/over,fraction=vfrac,length=vlength
      if iframe lt 10 then mn=strcompress('00'+string(long(iframe)),/remov) else $
         mn=strcompress('0'+string(long(iframe)),/remov)
      filename=strcompress('frames/ortho_movie_'+mn+'.png',/remove_all)
      write_png,filename,tvrd(True=1)
   endfor
   loadct,0,/silent

endelse

end
