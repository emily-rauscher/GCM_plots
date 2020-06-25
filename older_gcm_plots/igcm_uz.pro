PRO igcm_uz,oom=oom,p0=p0,files=files,psfile=psfile,urange=urange,ylabels=ylabels,$
            grayout=grayout,label=label,lchar=lchar,aveu=aveu,uz=uz,$
            ncb=ncb,nyl=nyl,cutvert=cutvert
;
; Routine to plot zonal [U] countour, reads from fort.26
;
; Keywords:
;  oom: orders of magnitude covered in pressure
;  p0: pressure at bottom boundary, in bar
;  files: a string array of filenames, for a time series of fort.26 files
;         to average and plot
;  psfile: saves plot as an eps file, either "plot.eps" or whatever
;          string psfile is set to
;  urange: force min,max values for contour plot
;  ylabels: a string array for yticknames (a crude way to force wanted
;           format)
;  grayout: pressure above which to overplot a rectangle of hatched lines,
;           indicating that the zonal ave shouldn't be trusted
;  label: strarr(2) to print in bottom left corner, label[0] above label[1]
; Optional output:  
;  aveu: returns nlev array with average u values at each level
;        (weighted by area, so cos(lat))
;  uz: nlat,nlev array holding [u] values (what is plotted)
;  ncb: if set will not plot colorbar
;  nyl: if set will not plot labels or title on yaxis (for grouping plots)
;  cutvert: if set (with ncb) will cut vertical extent of psfile to
;           eliminate whitespace leftover when colorbar is not plotted


if not keyword_set(psfile) then begin
; Screen output
      device,true_color=24,decomposed=0,retain=2
      !p.font=-1
      !p.charsize=1.5
endif else begin
   if (size(psfile))[1] eq 2 then filename='plot.eps' else filename=psfile
;   if keyword_set(ncb) then ysz=22.78 else ysz=25
   if keyword_set(nyl) then xsz=21.43 else xsz=25
   if keyword_set(cutvert) then ysz=22.22 else ysz=25
   psopen,filename,/enc,/color,xsize=xsz,ysize=ysz;,bits_per_pixel=24
  !p.font=0
  !p.charsize=2
endelse
!x.style=1

if not keyword_set(files) then begin
  
; Read File with 3D fields
  OpenR, 17, 'fort.26'
  ReadF, 17, nlat,nlon,nlev	  
	  print,nlat,nlon,nlev 
	  xy=fltarr(6,nlat,nlon,nlev)   
  ReadF, 17, xy
  Close, 17
;  lon=reform(xy[0,*,*,*]) & lat=reform(xy[1,*,*,*])
;  u=reform(xy[3,*,*,*]) & v=reform(xy[4,*,*,*]) & temp=reform(xy[5,*,*,*])
 
; Do zonal average and store LAT data for contour
  uz=total(reform(xy[3,*,*,*]),2)/nlon
  flat=reform(xy[1,*,0,*])

endif else begin

   nfiles=n_elements(files)
   
   OpenR, 17, files[0]
   ReadF, 17, nlat,nlon,nlev	  
   xy=fltarr(6,nlat,nlon,nlev)   
   ReadF, 17, xy
   Close, 17
   uz=total(reform(xy[3,*,*,*]),2)/nlon
   flat=reform(xy[1,*,0,*])
   
   for ifile=1,nfiles-1 do begin
      OpenR, 17, files[ifile]
      ReadF, 17, nlat,nlon,nlev	  
      xy=fltarr(6,nlat,nlon,nlev)   
      ReadF, 17, xy
      Close, 17
      uz+=total(reform(xy[3,*,*,*]),2)/nlon
   endfor

   uz/=nfiles
   
endelse
      
aveu=fltarr(nlev)
for i=0,nlev-1 do aveu[i]=total(uz[*,i]*cos(flat[*,0]*!pi/180.))/total(cos(flat[*,0]*!pi/180.))

;set flev to hold sigma levels (calculated based on nlev and whether oom set)
flev=fltarr(nlat,nlev)
if keyword_set(oom) then begin
   stp=-1.*oom/nlev
   flev(*,nlev-1)=replicate(10.^(stp/2.),nlat)
   for i=nlev-2,0,-1 do flev(*,i)=replicate(flev(0,i+1)*10.^(stp),nlat)
endif else begin
   for i=0,nlev-1 do flev(*,i)=replicate((0.5+i)/nlev,nlat)
endelse

if keyword_set(p0) then begin
   flev*=p0
   ytit='Pressure [bar]'
endif else begin
   ytit='sigma (normalized pressure)'
endelse
ymin=min(flev,max=ymax)

if keyword_set(urange) then begin
   umax=urange[1]
   umin=urange[0]
endif else begin
   umax=max(uz)
   umin=min(uz)
endelse

nlevels=25
step=(umax-umin)/nlevels
ulevels=indgen(nlevels)*step+umin

if umin eq umax then begin
   print,'Winds are uniform, at (m/s):',umin
endif else begin

   if keyword_set(cutvert) then cpos=[0.15,0.101,0.99,0.99] $
   else cpos=[0.15,0.09,0.99,0.88]
   if keyword_set(nyl) then begin
      yps=5
      ytit=''
      cpos[0]=0.01
   endif else yps=1

   if keyword_set(oom) then begin
;      if oom gt 5 then xmrg=[12,3] else xmrg=[10,3]
      ;plot in b&w to get axes in black, then overplot in color
      contour, uz,flat,flev, /cell_fill ,LEVELS=ulevels,ystyle=yps,$
               yr=[ymax,ymin],/ylog,xtitle='Latitude [degrees]',$
               position=cpos,$;ymargin=[4,4],xmargin=xmrg,$
               ytitle=ytit,max_value=umax,min_value=umin,ytickname=ylabels ;,/closed
;      loadct, 17,/silent
      ctload, 17,/silent,/reverse
      contour, uz,flat,flev,/overplot,/cell_fill ,LEVELS=ulevels,$
               max_value=umax,min_value=umin,position=cpos
   endif else begin
      contour, uz,flat,flev, /cell_fill ,LEVELS=ulevels,yr=[ymax,ymin],xmargin=[10,3],$
               /ystyle, xtitle='Latitude [degrees]',ytitle=ytit,ymargin=[4,4],$
               max_value=umax,min_value=umin,ytickname=ylabels
;      loadct, 17,/silent
      ctload, 17,/silent,/reverse
      contour, uz,flat,flev,/overplot,/cell_fill ,LEVELS=ulevels,$
               max_value=umax,min_value=umin,position=cpos
   endelse
   loadct,0,/silent
   contour,uz,flat,flev,/overplot,levels=[0.],color=0,thick=5 ;,/closed

   if keyword_set(grayout) then begin
      xmin=min(flat,max=xmax)
      xx=[xmin,xmin,xmax,xmax]
      yy=[grayout,ymin,ymin,grayout]
      polyfill,xx,yy,color=150,/data,/line_fill,orientation=45,thick=5
      polyfill,xx,yy,color=150,/data,/line_fill,orientation=135,thick=5
   endif

;   loadct,17,/silent
   ctload, 17,/silent,/reverse
   if not keyword_set(ncb) then colorbar,position=[0.23,0.935,0.9,0.95],$
                                title='Mean Zonal Wind Speed [m/s]',$
                                range=[umin,umax],format='(i5)',charsize=2,color=255
   
   print, ' Min and Max Zonal [U] Wind Speed (m/s):      ',min(uz),max(uz)

   loadct,0,/silent
   if keyword_set(label) then begin
      xyouts,min(flat)*.95,ymax/3.,label[0],color=0,charsize=lchar
      xyouts,min(flat)*.95,ymax/1.5,label[1],color=0,charsize=lchar
   endif

endelse

if keyword_set(psfile) then begin
   psclose
   spawn,'gs -r300 -dEPSCrop -dTextAlphaBits=4 -sDEVICE=png16m -sOutputFile=plot.png -dBATCH -dNOPAUSE plot.eps'
   spawn,'convert plot.png eps3:plot.eps'
endif
  
!x.style=0

END

