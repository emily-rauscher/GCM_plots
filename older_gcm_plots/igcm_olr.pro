PRO igcm_olr, latoff, lonoff,frange=frange,nogrid=nogrid,psfile=psfile,makepng=makepng,$
              ppoints=ppoints,files=files,porb=porb,calcphase=calcphase,delobs=delobs,diur=diur,$
              pcurve=pcurve,phase=phase,movie=movie,rot=rot,incl=incl,printOLR=printOLR,$
              obl=obl,einfo=einfo
;
; General keywords:
;   frange: use to force frange[0,1] as min,max for contour plots
;   nogrid: if set, does not put map grid on contour plot  
;   psfile: if ppoint not set, plot olr, otherwise phase curve
;   makepng: used to make png of flux map in isotropic projection,
;            saved in local directory as "flux_map.png"
; Keywords relevant for single map of emitted flux:
;   latoff = latitude offset in deg for map center
;   lonoff = longitude offset in deg for map center 
; Keywords relevant for emitted flux phase curve:
; (phase curve will be calculated by either setting ppoints OR files keywords)  
;   ppoints: number of points to calculate in orbital phase curve (if
;            not reading separate fort.64 files)
;   files = strarr with paths for fort.64 files from last orbit
;   porb: Porb/Prot (same as definition in igcm), used to calculate
;         sslon for non-synchronously rotating models. If set with
;         "files" keyword, will read day and sslon from the end of
;         each fort.64 file.
;   calcphase: if porb set but want to override read from fort.64
;              (e.g., if older file with bad formatting)
;   delobs: 180+delobs is longitude facing observer at beginning and
;           end of phase curve (delobs defaults to zero)
;   diur: if this simulation used a diurnally averaged forcing
;         pattern, then the sslon is always reported as 0; use this
;         keyword to override reported value, based on porb and
;         reported day value
;   pcurve: output (start and ends at orbital phase = 0,1 [transit],
;                   ppoints+1 points)
;   phase: output, array with phase values corresponding to pcurve
;   movie: if set, will save set of files for frames into local
;          subdirectory: frames/pc_movie_'+*+'.png
;          Can then make movie by:  convert -delay 10 pc_*.png movie.gif  
;   rot: set to rotate view of planet in movie (e.g. in case of
;        equatorial viewing of obliquity ne 0 planet)
;   incl: latitude of subobserver point (in degrees, 0 = equator), can
;         be any value within +/- obliquity of planet
;   printOLR: added by Erin
;   obl: obliquity, in degrees; setting this also means that instead
;        of a phase curve, THIS IS NOW BEING USED to calculate eclipse
;        depth as a function of time, assuming the viewing orientation
;        is always set to be during secondary eclipse, see next keyword
;   einfo: returns an array of [4,nfiles], with:
;          [0,*] is fmax-fmin (for the observed hemisphere)
;          [1,*] is the latitude of fmax
;          [2,*] is the latitude of fmin
;          [3,*] is the projected distance between fmax and fmin, in
;                normalized units of planet radii
;          ALSO NOTE that the following keywords will now output:
;            pcurve: the eclipse depth as a function of time
;            rot: an array of [nfiles], the observed angular
;              tilt [in degrees, clockwise from straight up] for the
;              brightest part of the planet (i.e. the "rot" used for
;              mapping adjusted by direction of brightest hemisphere)
;            incl: an [nfiles] array, the subobserver latitude (in degrees, and =sslat)
;          NOTE: overwrites any input values given to keywords incl & rot
;
; NOTE: for einfo, rot assumes that files starts with sslat~0, and goes
; positive, before turning over and becoming negative (consistent with
; finalorb output from code, but not necessarily general)
  
if not keyword_set(psfile) then begin
   if not keyword_set(makepng) then begin
                                ;Screen output
;   device,true_color=24,decomposed=0;,retain=2
      !p.font=-1
      window,0
   endif else window,xs=1000,ys=450
endif else begin
                                ;Postscript output
   if not keyword_set(ppoints) then begin
      if (size(psfile))[1] eq 2 then filename='plot.eps' else filename=psfile
      psopen,filename,/enc,x=25,y=20,/color
      !p.font=0
   endif
endelse
!x.style=1
!p.charsize=1.5

if not keyword_set(files) then begin  ; only one flux map to read in

   if not keyword_set(incl) then incl=0. ; in degrees, with 0=equatorial view
   if not keyword_set(rot) then rot=0.

   OpenR, 17, 'fort.64'
   ReadF, 17, nlat,nlon,nlev	  
   xy=fltarr(3,nlat,nlon)   
   ReadF, 17, xy
   Close, 17

   lon=reform(xy[0,*,*]) & lat=reform(xy[1,*,*]) & olr=reform(xy[2,*,*])

   if not keyword_set(makepng) then begin
      flon=fltarr(nlat,nlon+1)
      flon[*,0:nlon-1]=lon
      flon[*,nlon]=reform(lon[*,0]+360.)
      flat=fltarr(nlat,nlon+1)
      flat[*,0:nlon-1]=lat
      flat[*,nlon]=reform(lat[*,0])
      qtp=fltarr(nlat,nlon+1)
      qtp[*,0:nlon-1]=olr
      qtp[*,nlon]=reform(olr[*,0])
      loadct,5,/silent
   endif else begin
      flon=fltarr(nlat+2,nlon+1)
      flon[1:nlat,0:nlon-1]=lon
      flon[0,0:nlon-1]=lon[0,*]
      flon[nlat+1,0:nlon-1]=lon[0,*]
      flon[1:nlat,nlon]=reform(lon[*,0]+360.)
      flon[0,nlon]=lon[0,0]+360.
      flon[nlat+1,nlon]=lon[0,0]+360.
      flat=fltarr(nlat+2,nlon+1)
      flat[1:nlat,0:nlon-1]=lat
      flat[1:nlat,nlon]=reform(lat[*,0])
      flat[0,*]=90.
      flat[nlat+1,*]=-90.
      qtp=fltarr(nlat+2,nlon+1)
      qtp[1:nlat,0:nlon-1]=olr
      qtp[1:nlat,nlon]=reform(olr[*,0])
      qtp[0,*]=min(olr)
      qtp[nlat+1,*]=min(olr)
      loadct,3,/silent
   endelse

   qmin=min(qtp,max=qmax)
   if keyword_set(frange) then begin
      qmin=frange[0]
      qmax=frange[1]
   endif

   if qmin ne qmax then begin
      if not keyword_set(makepng) then $
         MAP_SET, /MILLER_CYLINDRICAL,latoff,lonoff ,/ISOTROPIC $
      else map_set,0,0,/isotropic,/noborder,position=[0,0,0.9,1]
      cbottom=0
      nlevels=45.
      step=(qmax-qmin)/nlevels
      mylevels=indgen(nlevels)*step+qmin
      ccolors=findgen(nlevels)*(255.-cbottom)/(nlevels-1)+cbottom
   
      contour, qtp,flon,flat, /overplot,/cell_fill,/closed,c_colors=ccolors,levels=mylevels
      if not keyword_set(nogrid) then begin
         if not keyword_set(makepng) then map_grid, /label,charsize=1 $
         else map_grid,/label,lonalign=1,latalign=1
      endif
      if not keyword_set(makepng) then $
         colorbar,position=[0.1,0.07,0.90,0.1],range=[qmin,qmax],format='(e8.2)',$
                  charsize=1.5,bottom=cbottom,ncolors=255-cbottom,divisions=5
   endif

   print, ' Min and Max OLR Fluxes (W/m^2): ',min(olr),max(olr)

   if (keyword_set(psfile)) and (not keyword_set(ppoints)) then begin
      psclose
      spawn,'gs -r300 -dEPSCrop -dTextAlphaBits=4 -sDEVICE=png16m -sOutputFile=plot.png -dBATCH -dNOPAUSE plot.eps'
      spawn,'convert plot.png eps3:plot.eps'
   endif

   if keyword_set(makepng) then write_png,'flux_map.png',tvrd(0,0,1800,900,/true)

   ; make a phase curve
   if keyword_set(ppoints) then begin
      if keyword_set(frange) then begin
         qmin=frange[0]
         qmax=frange[1]
         cbottom=0
         nlevels=45.
         step=(qmax-qmin)/nlevels
         mylevels=indgen(nlevels)*step+qmin
         ccolors=findgen(nlevels)*(255.-cbottom)/(nlevels-1)+cbottom
      endif
      phase=findgen(ppoints)/ppoints ; 0 is transit, 0.5 is secondary eclipse
      sobslon=180.-360.*phase        ; sub-observer point moves W as planet spins E
      if keyword_set(porb) then begin
         sobslon+=360.*phase*(1.-porb)
         sobslon = ((sobslon mod 360)+360) mod 360
      endif
      if keyword_set(movie) then begin
         for i=0,ppoints-1 do begin
            window,2,xsize=800,ysize=800
            map_set,incl,sobslon[i],rot,/orthographic,/isotropic
            if keyword_set(frange) then $
               contour,qtp,flon,flat,/overplot,/cell_fill,/closed,c_colors=ccolors,levels=mylevels $
            else $
               contour,qtp,flon,flat,/overplot,/cell_fill,/closed,nlevels=45,$
                       min_value=qmin,max_value=qmax
            if not keyword_set(nogrid) then map_grid, /label
            if i lt 10 then mn=strcompress('00'+string(long(i)),/remov) else $
               mn=strcompress('0'+string(long(i)),/remov)
            filename=strcompress('frames/pc_movie_'+mn+'.png',/remove_all)
            write_png,filename,tvrd(True=1)
         endfor
      endif

      pcurve=fltarr(ppoints)
      intpoints=nlon*5
      if incl eq 0 then begin  ; viewing along the equator, easy geometry
         for ip=0,ppoints-1 do begin
            amin=sobslon[ip]-90.
            if amin lt 0 then amin+=360.
            llon=findgen(intpoints)/(intpoints-1)*180.+amin
            ihigh=where(llon gt 360.)
            if ihigh[0] ne -1 then llon[ihigh]-=360.
            for ilat=0,nlat-1 do begin
               lolr=interpol(reform(olr[ilat,*]),reform(lon[ilat,*]),llon)
               pcurve[ip] += total(lolr*cos((llon-sobslon[ip])*!pi/180.)*(cos(lat[ilat,0]*!pi/180.))^2)
            endfor
         endfor
         pcurve=pcurve/total((cos(lat[*,0]*!pi/180.))^2)/total(cos((llon-sobslon[ip-1])*!pi/180.))
      endif else begin          ; viewing at a non-zero sub-observer point, difficult geometry!
         ; create an array of (x,y) points in projection plane to sample
         pxx=fltarr(intpoints,intpoints)
         pyy=pxx
         toss=(findgen(intpoints)-(intpoints-1)/2.)*2./(intpoints-1.)
         for i=0,intpoints-1 do begin
            pxx[*,i]=toss  ; pxx,pyy is x,y-coord of point [xi,yi], [0,0] is x,y=-1
            pyy[i,*]=toss
         endfor
         rho=sqrt(pxx^2+pyy^2)
         circ=where(rho le 1.)         ; reject points outside of circle (assuming radius=1)
         if circ[0] ne -1 then begin
            pxx=pxx[circ]
            pyy=pyy[circ]
            rho=rho[circ]
         endif
         cc=asin(rho/1.)  ; cc in radians
         ; use inverse formulas for orthographic projection
         plat=asin( cos(cc)*sin(incl*!pi/180.) $
                    + pyy*sin(cc)*cos(incl*!pi/180.)/rho )*180./!pi ; in degrees
         ; interpolate in lat,lon to get indices for olr array
         ; (use flon, flat, qtp, to get 0->360 coverage in lon)
         latind=interpol(findgen(nlat),reform(flat[*,0]),plat)
         toss=where(plat gt max(flat)) ; correct for latitude values too high
         if toss[0] ne -1 then latind[toss]=0                ; (b/c flat[0,0] is max latitude)
         toss=where(plat lt min(flat)) ; correct for latitude values too low
         if toss[0] ne -1 then latind[toss]=nlat-1           ; (b/c flat[nlat-1,0] is min latitude)
         for ip=0,ppoints-1 do begin
            plon=sobslon[ip] + atan( pxx*sin(cc), rho*cos(cc)*cos(incl*!pi/180.) $
                                     - pyy*sin(cc)*sin(incl*!pi/180.) )*180./!pi ; in degrees
            toss=where(plon lt 0.)
            if toss[0] ne -1 then plon[toss]+=360.
            toss=where(plon gt 360.)
            if toss[0] ne -1 then plon[toss]-=360.
            lonind=interpol(findgen(nlon+1),reform(flon[0,*]),plon)
            ; interpolate olr array at lat, lon points
            lolr=interpolate(qtp,latind,lonind,cubic=-0.5)
            ; total fluxes, divide by number of points sampled
            pcurve[ip]=total(lolr)/n_elements(lolr)
         endfor
      endelse
         
      loadct,0,/silent
      if not keyword_set(psfile) then begin
         window,2
      endif else begin
         if (size(psfile))[1] eq 2 then filename='plot.eps' else filename=psfile
         psopen,filename,/enc,x=20,y=15,/color
         !p.font=0
      endelse
      phase=[phase,phase[0]+1.]
      pcurve=[pcurve,pcurve[0]]
      plot,phase,pcurve,xtitle='Orbital phase',ytitle='Infrared flux [W/m^2]',/ynoz,psym=-4
      pmax=max(pcurve,maxind)
      xyouts,0.1,min(pcurve),strcompress('Sub-observer longitude of max flux: '+$
                                         string(180.-360.*phase[maxind])+' degrees'),charsize=1.5
      if keyword_set(psfile) then psclose
   endif

endif else begin  ; set of individual files to read in, snapshots throughout orbit

   nfiles=n_elements(files)
   pcurve=fltarr(nfiles)

   if keyword_set(porb) then begin
      day=dblarr(nfiles)
      phase=fltarr(nfiles)
   endif else begin
      phase=findgen(nfiles-1)/nfiles ; 0 is transit, 0.5 is secondary eclipse
      phase=[phase,1.]
      ang=180.-findgen(nfiles-1)*360./(nfiles-1) ;sub-observer point moves W as planet spins E
      ang=([ang,ang[0]]+360.) mod 360.
      if keyword_set(delobs) then begin
         ang=(ang+delobs+360.) mod 360.
         phase=phase-delobs/360.
      endif
   endelse

   OpenR, 17, files[nfiles-1]
   ReadF, 17, nlat,nlon,nlev	  
   xy=fltarr(3,nlat,nlon)   
   ReadF, 17, xy
   if keyword_set(porb) then begin
      toss=''
      readf,17,toss
      readf,17,toss
      starti=strpos(toss,'DAY')+4   ; where the value of DAY begins
      endi=strpos(toss,',',starti)
      dayf=double(strmid(toss,starti,endi-starti))
   endif
   Close, 17

   lon=reform(xy[0,*,*]) & lat=reform(xy[1,*,*])

   intpoints=nlon*5
   ; create an array of (x,y) points in projection plane to sample
   pxx=fltarr(intpoints,intpoints)
   pyy=pxx
   toss=(findgen(intpoints)-(intpoints-1)/2.)*2./(intpoints-1.)
   for i=0,intpoints-1 do begin
      pxx[*,i]=toss             ; pxx,pyy is x,y-coord of point [xi,yi], [0,0] is x,y=-1
      pyy[i,*]=toss
   endfor
   rho=sqrt(pxx^2+pyy^2)
   circ=where(rho le 1.)        ; reject points outside of circle (assuming radius=1)
   if circ[0] ne -1 then begin
      pxx=pxx[circ]
      pyy=pyy[circ]
      rho=rho[circ]
   endif
   cc=asin(rho/1.)              ; cc in radians

   flon=fltarr(nlat,nlon+1)
   flon[*,0:nlon-1]=lon
   flon[*,nlon]=reform(lon[*,0]+360.)
   flat=fltarr(nlat,nlon+1)
   flat[*,0:nlon-1]=lat
   flat[*,nlon]=reform(lat[*,0])

   if not keyword_set(obl) then begin
      if not keyword_set(incl) then incl=0
      if not keyword_set(rot) then rot=0

      if incl ne 0 then begin
         ; use inverse formulas for orthographic projection
         plat=asin( cos(cc)*sin(incl*!pi/180.) $
                    + pyy*sin(cc)*cos(incl*!pi/180.)/rho )*180./!pi ; in degrees
         ; interpolate in lat,lon to get indices for olr array
         ; (use flon, flat, qtp, to get 0->360 coverage in lon)
         latind=interpol(findgen(nlat),reform(flat[*,0]),plat)
         toss=where(plat gt max(flat)) ; correct for latitude values too high
         if toss[0] ne -1 then latind[toss]=0                ; (b/c flat[0,0] is max latitude)
         toss=where(plat lt min(flat)) ; correct for latitude values too low
         if toss[0] ne -1 then latind[toss]=nlat-1           ; (b/c flat[nlat-1,0] is min latitude)
      endif
   endif else begin
      einfo=fltarr(4,nfiles)
      ; subobserver latitude (=sslat for secondary eclipse)
      incl=fltarr(nfiles)
      ; rotation of N pole, from perspective of observer
      rot=fltarr(nfiles)
      ; (vary with time, so will calculate for each file)
   endelse

   for ip=0,nfiles-1 do begin  ; loop through set of files

      OpenR, 17, files[ip]
      ReadF, 17, nlat,nlon,nlev	  
      xy=fltarr(3,nlat,nlon)   
      ReadF, 17, xy
      if keyword_set(porb) then begin
         toss=''
         readf,17,toss
         readf,17,toss
         starti=strpos(toss,'DAY')+4 ; where the value of DAY begins
         endi=strpos(toss,',',starti)
         day[ip]=double(strmid(toss,starti,endi-starti))
         starti=strpos(toss,':') ; last print before sslon, sslat
         sslon=float(strmid(toss,starti+1,8))
         sslat=float(strmid(toss,starti+1+8,8))
         ; check in case of old file format:
;         print,'day, sslon, sslat:',day[ip],sslon,sslat
         if keyword_set(diur) then begin
            sslon=360.*day[ip]*(1./porb-1)
;            print,'effective sslon:',sslon
         endif
      endif
      Close, 17

      olr=reform(xy[2,*,*])

      qtp=fltarr(nlat,nlon+1)
      qtp[*,0:nlon-1]=olr
      qtp[*,nlon]=reform(olr[*,0])
      qmin=min(qtp,max=qmax)
      print,'Min/max flux [W/m^2]',qmin,qmax

      if keyword_set(frange) then begin
         qmin=frange[0]
         qmax=frange[1]
      endif

      if keyword_set(porb) then begin
         if not keyword_set(calcphase) then begin
            phase[ip]=(day[ip]-day[0])/porb
            sobslon=(sslon+180.-360.*(day[ip]-day[0])/porb) mod 360.
         endif else begin
            phase[ip]=(ip*(dayf-day[0])/nfiles)/porb
            sobslon = (sslon + 180. - 360.*ip*(dayf-day[0])/nfiles/porb) mod 360.
         endelse
         while (sobslon lt 0.) do sobslon+=360.
      endif else $
         sobslon = ang[ip]

      if keyword_set(obl) then begin
         incl[ip]=sslat  ; in degrees
         rot[ip]=acos(cos(obl*!pi/180.)/cos(sslat*!pi/180.))*180./!pi   ; in degrees
         if (phase[ip] gt 0.25) and (phase[ip] lt 0.75) then rot[ip]*=-1
      endif

      if keyword_set(movie) then begin
         loadct,5,/silent
         window,2,xsize=800,ysize=800
         if not keyword_set(obl) then $
            map_set,incl,sobslon,rot,/orthographic,/isotropic $
         else $
            map_set,incl[ip],sobslon,rot[ip],/orthographic,/isotropic
         contour,qtp,flon,flat,/overplot,/cell_fill,/closed,nlevels=45,$
                 min_value=qmin,max_value=qmax
         if not keyword_set(nogrid) then map_grid, /label
         if ip lt 10 then mn=strcompress('00'+string(long(ip)),/remov) else $
            mn=strcompress('0'+string(long(ip)),/remov)
         filename=strcompress('frames/pc_movie_'+mn+'.png',/remove_all)
         write_png,filename,tvrd(True=1)
         loadct,0,/silent
      endif

      if not keyword_set(obl) then begin
         if incl eq 0 then begin
            amin=sobslon-90.
            if amin lt 0 then amin+=360.
            llon=findgen(intpoints)/(intpoints-1)*180.+amin
            ihigh=where(llon gt 360.)
            if ihigh[0] ne -1 then llon[ihigh]-=360.
            for ilat=0,nlat-1 do begin
               lolr=interpol(reform(qtp[ilat,*]),reform(flon[ilat,*]),llon)
               pcurve[ip] += total(lolr*cos((llon-sobslon)*!pi/180.)*(cos(lat[ilat,0]*!pi/180.))^2)
            endfor
         endif else begin
            ; use inverse formulas for orthographic projection
            plon=sobslon + atan( pxx*sin(cc), rho*cos(cc)*cos(incl*!pi/180.) $
                                 - pyy*sin(cc)*sin(incl*!pi/180.) )*180./!pi ; in degrees
            toss=where(plon lt 0.)
            if toss[0] ne -1 then plon[toss]+=360.
            toss=where(plon gt 360.)
            if toss[0] ne -1 then plon[toss]-=360.
            ; interpolate in lat,lon to get indices for olr array
            ; (use flon, flat, qtp, to get 0->360 coverage in lon)
            lonind=interpol(findgen(nlon+1),reform(flon[0,*]),plon)
            ; interpolate olr array at lat, lon points
            lolr=interpolate(qtp,latind,lonind,cubic=-0.5)
            ; total fluxes, divide by number of points sampled
            pcurve[ip]=total(lolr)/n_elements(lolr)
         endelse
      endif else begin
         plat=asin( cos(cc)*sin(incl[ip]*!pi/180.) $
                    + pyy*sin(cc)*cos(incl[ip]*!pi/180.)/rho )*180./!pi ; in degrees
         latind=interpol(findgen(nlat),reform(flat[*,0]),plat)
         toss=where(plat gt max(flat)) ; correct for latitude values too high
         if toss[0] ne -1 then latind[toss]=0                ; (b/c flat[0,0] is max latitude)
         toss=where(plat lt min(flat)) ; correct for latitude values too low
         if toss[0] ne -1 then latind[toss]=nlat-1           ; (b/c flat[nlat-1,0] is min latitude)
         plon=sobslon + atan( pxx*sin(cc), rho*cos(cc)*cos(incl[ip]*!pi/180.) $
                                 - pyy*sin(cc)*sin(incl[ip]*!pi/180.) )*180./!pi ; in degrees
         toss=where(plon lt 0.)
         if toss[0] ne -1 then plon[toss]+=360.
         toss=where(plon gt 360.)
         if toss[0] ne -1 then plon[toss]-=360.
         lonind=interpol(findgen(nlon+1),reform(flon[0,*]),plon)
         lolr=interpolate(qtp,latind,lonind,cubic=-0.5)
         pcurve[ip]=total(lolr)/n_elements(lolr)

         ; find Fmax, Fmin for observed hemisphere
         fmin=min(lolr,minind)
         fmax=max(lolr,maxind)
         ; adjust rot by 180 if southern hemisphere brighter
         if plat[maxind] lt 0. then rot[ip]+=180.
         einfo[0,ip]=fmax-fmin
         einfo[1,ip]=plat[maxind]
         einfo[2,ip]=plat[minind]
         ; find the projected distance between fmax and fmin
         if not keyword_set(movie) then begin
            toss=map_proj_init(2,sphere_radius=1.,center_latitude=incl[ip],$
                               center_long=sobslon,rotation=rot[ip])
            fxy=map_proj_forward([sobslon,sobslon],[plat[maxind],plat[minind]],$
                                map_structure=toss)
         endif else $
            fxy=map_proj_forward([sobslon,sobslon],[plat[maxind],plat[minind]])
         einfo[3,ip]=sqrt((fxy[0,0]-fxy[0,1])^2+(fxy[1,0]-fxy[1,1])^2)

      endelse

   endfor  ; end of loop over files

   ; (normalized/corrected by int_-90^90 lat^2 and int_-90^90 lon)
   if not keyword_set(obl) then $
      if incl eq 0 then $
         pcurve=pcurve/total((cos(lat[*,0]*!pi/180.))^2)/total(cos((llon-sobslon)*!pi/180.))
   
   if not keyword_set(psfile) then begin
      window,2
   endif else begin
      if (size(psfile))[1] eq 2 then filename='plot.eps' else filename=psfile
      psopen,filename,/enc,x=20,y=15,/color
      !p.font=0
   endelse

   plot,phase,pcurve,xtitle='Orbital phase',ytitle='Infrared flux [W/m^2]',/ynoz
   pmax=max(pcurve,maxind)
   if keyword_set(psfile) then begin
      psclose
      !p.font=-1
      spawn,'gs -r300 -dEPSCrop -dTextAlphaBits=4 -sDEVICE=png16m -sOutputFile=plot.png -dBATCH -dNOPAUSE plot.eps'
      spawn,'convert plot.png eps3:plot.eps'
   endif

endelse

   if keyword_set(printOLR) then begin
   n=N_ELEMENTS(OLR[*,0])
   OPENW,25,'OLR.txt'
   for j=0,n-1 DO BEGIN
      PRINTF,25,OLR[j,0]
   endfor
   CLOSE,25
   endif
!x.style=0
!p.charsize=1

END

