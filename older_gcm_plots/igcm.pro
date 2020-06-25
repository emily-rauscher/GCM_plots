PRO igcm, qtplot,levplot,latoff,lonoff,prange=prange,vfrac=vfrac,psfile=psfile,$
          zave=zave,zpres=zpres,moviefiles=moviefiles,label=label
;   qtplot = 1, 2 or 3 for U vel, V vel or Temp
;   levplot = level number (TOP:0 to BOTTOM:NL-1, i.e has to be known)
;   latoff = latitude offset in deg for map center
;   lonoff = longitude offset in deg for map center 
;   prange = forced [min,max] for contour plot
;   psfile = filename for plot, if set as /psfile, then defaults to
;   plot.eps
;   label = strarr(2), printed in top left corner
;
; MOVIEFILES gives the filenames (including full path) for a series of
; fort.50 files that will be read in sequence and plotted, to create a
; series of png files, which can then be converted into a movie via:
;    convert -delay 10 frame_*.png movie.gif
; WARNING: will dump frame files into local directory's subdirectory: frames

  if not keyword_set(latoff) then latoff=0
  if not keyword_set(lonoff) then lonoff=0

if not keyword_set(moviefiles) then begin
   if not keyword_set(psfile) then begin
;  device,true_color=24,decomposed=0,retain=2
   endif else begin
      if (size(psfile))[1] eq 2 then filename='plot.eps' else filename=psfile
      if keyword_set(zave) then xsz=40 else xsz=25
      psopen,filename,/enc,/color,xsize=xsz,ysize=20 ;,bits_per_pixel=24
      !p.font=0
      !x.thick=6
      !y.thick=6
      !p.thick=3
      !p.charsize=2
   endelse

; Read File with 3D fields
   OpenR, 17, 'fort.26'
   ReadF, 17, nlat,nlon,nlev	  
;	  print,nlat,nlon,nlev 
   xy=fltarr(6,nlat,nlon,nlev)   
   ReadF, 17, xy
   Close, 17

; Pick vertical level (lower layer is nlev-1)
;  nl=nlev   ; lower layer
;  nl=1        ; top layer
   nl=levplot                   ; from function command line
    
   lon=reform(xy[0,*,*,nl]) & lat=reform(xy[1,*,*,nl]) 
   u=reform(xy[3,*,*,nl]) & v=reform(xy[4,*,*,nl]) & temp=reform(xy[5,*,*,nl])
 
   if keyword_set(zave) then begin
      uz=fltarr(nlat)
      vz=fltarr(nlat)
      tz=fltarr(nlat)
      for i=0,nlat-1 do begin
         uz[i]=mean(u[i,*])
         vz[i]=mean(v[i,*])
         tz[i]=mean(temp[i,*])
      endfor
   endif

;; Poor fix on LON/LAT for grid not fitting 90 and 180 degrees exactly 
;; Stretches the LON/LAT scales to fill contour, incorrect but OK a hires
;  flon=lon*360./lon(0,nlon-1)
;  flat=lat*90./lat(0,nlat-1)
   flon=fltarr(nlat,nlon+1)
   flon[*,0:nlon-1]=lon
   flon[*,nlon]=reform(lon[*,0]+360.)
   flat=fltarr(nlat,nlon+1)
   flat[*,0:nlon-1]=lat
   flat[*,nlon]=reform(lat[*,0])

; PLOTTED Qty from command line
   qtp=fltarr(nlat,nlon+1)
   if (qtplot EQ 1) THEN begin
      qtp[*,0:nlon-1]=u
      qtp[*,nlon]=reform(u[*,0])
      if keyword_set(zave) then for i=0,nlat-1 do qtp[i,*]=qtp[i,*]-uz[i]
      ctn=25
   endif
   if (qtplot EQ 2) THEN begin
      qtp[*,0:nlon-1]=v
      qtp[*,nlon]=reform(v[*,0])
      if keyword_set(zave) then for i=0,nlat-1 do qtp[i,*]=qtp[i,*]-vz[i]
      ctn=10
   endif
   if (qtplot EQ 3) THEN begin
      qtp[*,0:nlon-1]=temp
      qtp[*,nlon]=reform(temp[*,0])
      if keyword_set(zave) then for i=0,nlat-1 do qtp[i,*]=qtp[i,*]-tz[i]
      ctn=4
   endif
   qmin=min(qtp,max=qmax)
   if keyword_set(prange) then begin
      qmin=prange[0]
      qmax=prange[1]
   endif

; For plotting velocity vector     
;  on=where(abs(flat) ne 90.)
   on=where(abs(flat[*,0:nlon-1]) ne max(flat))
   up=u[*,0:nlon-1]
   vp=v[*,0:nlon-1]
   if keyword_set(zave) then begin
      for i=0,nlat-1 do begin
         up[i,*]=up[i,*]-uz[i]
         vp[i,*]=vp[i,*]-vz[i]
      endfor
   endif

   if keyword_set(zave) then begin
      if keyword_set(zpres) then mtit=strcompress('Pressure level:'+string(zpres)+' bar')
      MAP_SET, /Miller_CYLINDRICAL,latoff,lonoff,/ISOTROPIC,xmargin=[37,1],title=mtit
      info=strcompress('Max. zonal, meridional wind deviations: '+ $
                       string(round(max(abs(up))))+','+string(round(max(abs(vp))))+' m/s')
      xyouts,0.48,0.85,info,/normal
   endif else $
      MAP_SET, /Miller_CYLINDRICAL,latoff,lonoff,/ISOTROPIC

   loadct, ctn,/silent
   cbottom=55.
   nlevels=45.
   step=(qmax-qmin)/nlevels
   mylevels=indgen(nlevels)*step+qmin
   ccolors=findgen(nlevels)*(255.-cbottom)/(nlevels-1)+cbottom

   contour, qtp,flon,flat, /overplot,/cell_fill,/closed,c_colors=ccolors,levels=mylevels

   if qmin ne qmax then begin
      if keyword_set(zave) then cbpos=[0.5,0.1,0.9,0.13] else cbpos=[0.1,0.07,0.90,0.1]
      colorbar,position=cbpos,range=[qmin,qmax],format='(i5)',charsize=2,bottom=cbottom,$
               ncolors=255-cbottom
   endif

   loadct,5,/silent
   partvelvec,up[on],vp[on],(flon[*,0:nlon-1])[on],(flat[*,0:nlon-1])[on],/over,$
              fraction=vfrac,veccolors=65
   loadct,0,/silent
   map_grid, /label ,color=255,charsize=1.
   if keyword_set(label) then begin
      xyouts,-175+lonoff,70,label[0],color=255
      xyouts,-175+lonoff,62,label[1],color=255
   endif

   if keyword_set(zave) then begin
      lmin=min(lat[*,0],max=lmax)
      loadct,5,/silent
      plot,tz,reform(lat[*,0]),xtit='Zonally averaged temperature [K]',$
           ytit='Latitude [degrees]',thick=6,xstyle=8,ystyle=1,/noerase,$
           position=[0.07,0.15,0.37,0.82],yr=[lmin,lmax];,ytickn=replicate(' ',3)
      oplot,tz,reform(lat[*,0]),thick=6,color=100
      wmin=min([uz,vz],max=wmax)
      axis,xaxis=1,xr=[wmin,wmax],/save,xtit='Zonally averaged wind speed [m/s]'
      oplot,uz,reform(lat[*,0]),color=50,thick=6
      oplot,vz,reform(lat[*,0]),color=50,thick=6,linestyle=2
      oplot,[0,0],[-90,90],linestyle=1
      al_legend,['Zonal wind','Meridional wind','Temperature'],linestyle=[0,2,0],$
                color=[50,50,100],/right,charsize=1.5,thick=5
      loadct,0,/silent
   endif

   print, strcompress('L'+string(levplot+1),/remove_all)
   print, ' Min and Max Temperatures (K):            ',min(temp),max(temp)
   print, ' Min and Max Zonal Wind Speed (m/s):      ',min(u),max(u)
   print, ' Min and Max Meridional Wind Speed (m/s): ',min(v),max(v)
   print, ' (For R=3523 and akap=0.286,)'
   print, '         Min and Max Sound Speed (m/s):   ',70.24*sqrt(minmax(temp))

   if keyword_set(psfile) then begin
      !p.font=-1
      !x.thick=1
      !y.thick=1
      !p.thick=1
      !p.charsize=1
      psclose
;      spawn,'convert plot.eps -crop 708x486+0+80 +repage plot.png'
      spawn,'gs -r300 -dEPSCrop -dTextAlphaBits=4 -sDEVICE=png16m -sOutputFile=plot.png -dBATCH -dNOPAUSE plot.eps'
      spawn,'convert plot.png eps3:plot.eps'
   endif

endif else begin

   nfiles=n_elements(moviefiles)
   for i=0,nfiles-1 do begin

      if i lt 10. then begin
         fpart=strcompress('00'+string(i,format='(i4)'),/remove_all)
      endif else begin
         if i lt 100. then $
            fpart=strcompress('0'+string(i,format='(i4)'),/remove_all) $
         else $
            fpart=string(i,format='(i4)')
      endelse
      pngname=strcompress('frames/frame_'+fpart+'.png',/remove_all)

      OpenR, 17, moviefiles[i]
      ReadF, 17, nlat,nlon,nlev	  
      xy=fltarr(6,nlat,nlon,nlev)   
      ReadF, 17, xy
      Close, 17

      nl=levplot                ; from function command line
      lon=reform(xy[0,*,*,nl]) & lat=reform(xy[1,*,*,nl]) 
      u=reform(xy[3,*,*,nl]) & v=reform(xy[4,*,*,nl]) & temp=reform(xy[5,*,*,nl])

      flon=fltarr(nlat,nlon+1)
      flon[*,0:nlon-1]=lon
      flon[*,nlon]=reform(lon[*,0]+360.)
      flat=fltarr(nlat,nlon+1)
      flat[*,0:nlon-1]=lat
      flat[*,nlon]=reform(lat[*,0])

      qtp=fltarr(nlat,nlon+1)
      if (qtplot EQ 1) THEN begin
         qtp[*,0:nlon-1]=u
         qtp[*,nlon]=reform(u[*,0])
      endif
      if (qtplot EQ 2) THEN begin
         qtp[*,0:nlon-1]=v
         qtp[*,nlon]=reform(v[*,0])
      endif
      if (qtplot EQ 3) THEN begin
         qtp[*,0:nlon-1]=temp
         qtp[*,nlon]=reform(temp[*,0])
      endif
      qmin=min(qtp,max=qmax)
      if keyword_set(prange) then begin
         qmin=prange[0]
         qmax=prange[1]
      endif

      on=where(abs(flat[*,0:nlon-1]) ne max(flat))
      up=u[*,0:nlon-1]
      vp=v[*,0:nlon-1]

      MAP_SET, /Miller_CYLINDRICAL,latoff,lonoff,/ISOTROPIC

      loadct, 4,/silent
      cbottom=55.
      nlevels=45.
      step=(qmax-qmin)/nlevels
      mylevels=indgen(nlevels)*step+qmin
      ccolors=findgen(nlevels)*(255.-cbottom)/(nlevels-1)+cbottom

      contour, qtp,flon,flat, /overplot,/cell_fill,/closed,c_colors=ccolors,levels=mylevels

      if qmin ne qmax then begin
         cbpos=[0.1,0.07,0.90,0.1]
         colorbar,position=cbpos,range=[qmin,qmax],format='(i5)',charsize=2,$
                  bottom=cbottom,ncolors=255-cbottom
      endif
      loadct,5,/silent
      partvelvec,up[on],vp[on],(flon[*,0:nlon-1])[on],(flat[*,0:nlon-1])[on],/over,$
                 fraction=vfrac,veccolors=65
      loadct,0,/silent
      map_grid, /label ,color=255,charsize=1.
      write_png,pngname,tvrd(true=1)
      print, ' Min and Max Temperatures (K):            ',min(temp),max(temp)
   endfor

endelse

END
