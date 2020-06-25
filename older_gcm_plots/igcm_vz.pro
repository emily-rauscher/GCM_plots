PRO igcm_vz,oom=oom,psfile=psfile,p0=p0,vrange=vrange

if not keyword_set(psfile) then begin
; Screen output
  device,true_color=24,decomposed=0,retain=2
  !p.font=-1
  !p.charsize=1.5
endif else begin
   if (size(psfile))[1] eq 2 then filename='plot.eps' else filename=psfile
   psopen,filename,/enc,/color,xsize=25,ysize=25;,bits_per_pixel=24
  !p.font=0
  !p.charsize=2
endelse

  !x.style=1

$ Read File with 3D fields

  OpenR, 17, 'fort.26'
  ReadF, 17, nlat,nlon,nlev	  
	  print,nlat,nlon,nlev 
	  xy=fltarr(6,nlat,nlon,nlev)   
  ReadF, 17, xy
  Close, 17


    lon=reform(xy[0,*,*,*]) & lat=reform(xy[1,*,*,*])
    u=reform(xy[3,*,*,*]) & v=reform(xy[4,*,*,*]) & temp=reform(xy[5,*,*,*])

  vz=fltarr(nlat,nlev)
  flat=fltarr(nlat,nlev)
  flev=fltarr(nlat,nlev)

$ Do zonal average and store LAT and LEV data for contour
  FOR K=0,nlev-1 DO BEGIN
  FOR I=0,nlat-1 DO BEGIN
       FOR J=0,nlon-1 DO BEGIN
       vz(I,K)=vz(I,K)+v(I,J,K)
       ENDFOR
       vz(I,K)=vz(I,K)/nlon
       flat(I,K)=lat(I,0,0)
  ENDFOR
  ENDFOR

;set flev to hold sigma levels (calculated based on nlev and whether oom set)
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

if keyword_set(vrange) then begin
   vmax=vrange[1]
   vmin=vrange[0]
endif else begin
   vmax=max(vz)
   vmin=min(vz)
endelse


nlevels=25
step=(vmax-vmin)/nlevels
vlevels=indgen(nlevels)*step+vmin

  loadct, 4,/silent
   if keyword_set(oom) then begin
      contour, vz,flat,flev, /cell_fill ,LEVELS=vlevels,/ystyle,ymargin=[4,4],xmargin=xmrg,$
               yr=[max(flev),min(flev)],/ylog,xtitle='Latitude [degrees]',$
               ytitle=ytit,max_value=vmax,min_value=vmin,ytickname=ylabels
   endif else begin
      contour, vz,flat,flev, /cell_fill ,LEVELS=vlevels,yr=[max(flev),min(flev)],xmargin=xmrg,$
               /ystyle, xtitle='Latitude [degrees]',ytitle=ytit,ymargin=[4,4],$
               max_value=vmax,min_value=vmin,ytickname=ylabels
   endelse
   contour,vz,flat,flev,/overplot,levels=[0.],color=255,thick=5,/closed
   colorbar,position=[0.23,0.965,0.9,0.98],range=[vmin,vmax],format='(i5)',charsize=2

   print, ' Min and Max Meridional [V] Wind Speed (m/s):      ',min(vz),max(vz)

loadct,0,/silent

if keyword_set(psfile) then begin
   psclose
   spawn,'gs -r300 -dEPSCrop -dTextAlphaBits=4 -sDEVICE=png16m -sOutputFile=plot.png -dBATCH -dNOPAUSE plot.eps'
   spawn,'convert plot.png eps3:plot.eps'
endif

!x.style=0

END

