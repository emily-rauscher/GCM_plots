pro surfp,p0=p0,psfile=psfile,moviefiles=moviefiles,cmin=cmin,cmax=cmax,sp=sp,lon=lon,lat=lat
;
; MOVIEFILES gives the filenames (including full path) for a series of
; fort.50 files that will be read in sequence and plotted, to create a
; series of png files, which can then be converted into a movie via:
;    convert -delay 10 frame_*.png movie.gif
; WARNING: will dump frame files into local directory's subdirectory: frames
;
; can output surface pressure field as sp and coordinates as lon,lat  

  if not keyword_set(p0) then p0=1.

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

  if not keyword_set(moviefiles) then begin
     openr, 17, 'fort.50'
     readf, 17, nlat,nlon
     ab=fltarr(3,nlat,nlon)
     readf, 17, ab
     close, 17
     lon=reform(ab[0,*,*]) & lat=reform(ab[1,*,*]) & sp=reform(ab[2,*,*])
     sp = (sp+1.)*p0
     minsp = min(sp)
     maxsp = max(sp)
     print,minsp,maxsp

     flon=fltarr(nlat,nlon+1)
     flon[*,0:nlon-1]=lon
     flon[*,nlon]=reform(lon[*,0]+360.)
     flat=fltarr(nlat,nlon+1)
     flat[*,0:nlon-1]=lat
     flat[*,nlon]=reform(lat[*,0])
     psp=fltarr(nlat,nlon+1)
     psp[*,0:nlon-1]=sp
     psp[*,nlon]=reform(sp[*,0])
     
     map_set,/cylindrical,0.,0.,/isotropic ;,title='Surface pressure',charsize=1.5

     ctload,1,/silent,/reverse
     contour,psp,flon,flat,/overplot,/cell_fill,nlevels=55,min_value=cmin,max_value=cmax
     map_grid,/label,charsize=1

     if minsp ne maxsp then begin
        loadct,1,/silent
        if keyword_set(cmin) then minsp=cmin
        if keyword_set(cmax) then maxsp=cmax
        colorbar, position=[0.1,0.05,0.90,0.09],range=[minsp,maxsp],format='(f8.2)',$
                  charsize=1.5,/reverse
     endif else xyouts,10,10,'Constant',charsize=1.5
     loadct,0,/silent

     if keyword_set(psfile) then begin
        !p.font=-1
        !x.thick=1
        !y.thick=1
        !p.thick=1
        !p.charsize=1
        psclose
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

        openr, 17, moviefiles[i]
        readf, 17, nlat,nlon
        ab=fltarr(3,nlat,nlon)
        readf, 17, ab
        close, 17
        lon=reform(ab[0,*,*]) & lat=reform(ab[1,*,*]) & sp=reform(ab[2,*,*])
        sp = (sp+1.)*p0
        minsp = min(sp)
        maxsp = max(sp)
        print,minsp,maxsp

        flon=fltarr(nlat,nlon+1)
        flon[*,0:nlon-1]=lon
        flon[*,nlon]=reform(lon[*,0]+360.)
        flat=fltarr(nlat,nlon+1)
        flat[*,0:nlon-1]=lat
        flat[*,nlon]=reform(lat[*,0])
        psp=fltarr(nlat,nlon+1)
        psp[*,0:nlon-1]=sp
        psp[*,nlon]=reform(sp[*,0])
     
        map_set,/cylindrical,0.,0.,/isotropic ;,title='Surface pressure',charsize=1.5

        ctload,1,/silent,/reverse
        contour,psp,flon,flat,/overplot,/cell_fill,nlevels=55,min_value=cmin,max_value=cmax
        map_grid,/label,charsize=1

        if minsp ne maxsp then begin
           loadct,1,/silent
           if keyword_set(cmin) then minsp=cmin
           if keyword_set(cmax) then maxsp=cmax
           colorbar, position=[0.1,0.05,0.90,0.09],range=[minsp,maxsp],format='(f8.2)',$
                     charsize=1.5,/reverse
        endif else xyouts,10,10,'Constant',charsize=1.5
        loadct,0,/silent

        write_png,pngname,tvrd(true=1)

     endfor

  endelse
end
