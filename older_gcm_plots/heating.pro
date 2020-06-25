pro heating,nlat,nlon,nlev,psfile=psfile,xrange=xrange,oned=oned,zeros=zeros,$
            dots=dots,kps=kps,mbar=mbar,missinglat=missinglat,oldheader=oldheader,$
            sslon=sslon,sslat=sslat,watts=watts,plevel=plevel,mheat=mheat,ape=ape,$
            radea=radea,grav=grav,gascon=gascon

; plots SW and LW heating/cooling rates throughout the atmosphere
; reads fort.63
;
; required:
; nlat,nlon,nlev: the number of latitude, longitude, and level points
;          
; optional keywords:
; psfile: will save the plot to an eps file, with "psfile" as the
;         filename, unless set to 1, in which case defaults to "plot.eps"
; xrange: used for user-defined xrange on plots, xmin,xmax=xrange[0:1]
; oned: if set, will only plot global totals (no lines for different locations)
; zeros: if set, will plot points along xmin axis for where rates are =0
; dots: if set, will plot points where one level fluctates to postive
;       values (from neighboring negative values), or vice versa
; kps: if set, will use (K/s) as units for heating (default is W/m^2)
; watts: if set, will use (W) as units for heating, warning: requires
;        cp, radea, and grav to be correct ... also could be misleading
;        to plot this way
; plevel: instead of vertical profiles, plots the heating rate at the
;         input pressure level (older part of code, should
;         double-check before use)
; mbar: for all fort.63 produced after Jun 2011, this keyword is
;       unnecessary, for earlier files the pressure units were printed
;       in mbar
; missinglat: for use before the correction to igcm in which the first
;             entry in fort.63 was overwritten and the northernmost
;             latitude is missing, this cheats by copying information
;             from the southernmost latitude
; oldheader: use when reading in old files (before they included a
;            print of the day and sslon at beginning of file)
; sslon: substellar longitude, should be automatically read in from
;        file, but for older files wasn't printed to the header
;        (i.e., use this when keyword "oldheader" is used)
; sslat: substellar latitude, as above
; radea: radius of the planet, in (m), only required if keywords watts
;        or ape set
; grav: gravitational acceleration (m/s^2), only required if keywords
;     watts or ape set
; ape: (older code, double-check before use) plots APE (available
;      potential energy ... and associated quantities?)
; mheat: (older code, double-check before use) returns [2,nlev] array
;        with the total magnetic heating rate at each level (from
;        magheating.pro), both the "true" value [0,*] and what was
;        used with the tdrag_min limit [1,*]
;
;
; format of fort.63 is (as of Fall 2014):
;DAY:   2.00, SUBSTELLAR LON,LAT:   0.00   0.00
; 
; LATITUDE, LONGITUDE:  -81.6505907502981       0.000000000000000E+000
;           UPWARD FLUX (Wm-2)     DOWNWARD FLUX (Wm-2)     NET FLUX (Wm-2)           SW NET FLUX (Wm-2)
;         0.006      0.21033E+06           0.00000E+00           0.21033E+06           0.15074E+06
;       .... (nlev+1 entries total, last being surface) ........
;    101.645      0.33613E+07           0.33578E+07           0.35000E+04           0.00000E+00
;      PRESSURE (bar)      HEATING RATES: SW (K/DAY) LW (K/DAY)
;          0.00  0.34641E+05 -0.34234E+05
;       .... (nlev+1 entries total, last being surface) ........
;     101.645  0.00000E+00  0.00000E+00
; 
; LATITUDE, LONGITUDE:  -81.6505907502981        11.2500000000000     
;       etc., for nlat, nlon
;
; order is: max neg lat, lon=0 to lon=max; max pos lat, lon=0,max;
; etc.
;
;  

if keyword_set(watts) then begin
   if not keyword_set(gascon) then gascon=3523.
   cp=gascon/0.286
   if not keyword_set(radea) then radea=1.e8
   if not keyword_set(grav) then grav=8.
   print,'Careful, is radea, c_p, and grav correct?', radea,cp,grav
endif

entry = create_struct( 'lat', 0.0, $
                       'lon', 0.0, $
                       'p0',0.0, $
                       'fpres', fltarr(nlev), $
                       'sflux', fltarr(nlev), $  ; SW flux
                       'flux', fltarr(nlev), $   ; actual IR flux
                       'rflux', fltarr(nlev), $    ; cirrad flux
                       'pres', fltarr(nlev), $
                       'sw', fltarr(nlev), $  ;sw heating
                       'lw', fltarr(nlev))  ;lw heating
nent=nlat*nlon
rates=replicate(entry,nent)

read1=strarr(2)
data1=fltarr(5,nlev+1)
read4=''
data2=fltarr(3,nlev+1)
read5=''

openr,17,'fort.63'

if not keyword_set(oldheader) then begin
   toss=strarr(2)
   readf,17,toss
   print,toss[0]
   sslon=strmid(toss[0],35,4)
   sslat=strmid(toss[0],42,4)
endif

;(the first set of nlon entries are the southernmost latitude, if
; missing lat, copy these to create northernmost)
if keyword_set(missinglat) then begin
   ; read first latitude separately, before moving onto the other entries
   for i=0,nlon-1 do begin
      readf,17,read1
      rates[i+nlon].lat=float(strmid(read1[0],23,6))
      rates[i+nlon].lon=float(strmid(read1[0],47,6))
      readf,17,data1
      rates[i+nlon].fpres=reform(data1[0,0:nlev-1])
      rates[i+nlon].rflux=reform(data1[1,0:nlev-1]-data1[2,0:nlev-1])
      rates[i+nlon].flux=reform(data1[3,0:nlev-1])
      rates[i+nlon].sflux=reform(data1[4,0:nlev-1])
      readf,17,read4
      readf,17,data2
      readf,17,read5
      rates[i+nlon].pres=reform(data2[0,0:nlev-1]) ; in bar (or mbar, if keyword set)
      rates[i+nlon].p0=data2[0,nlev]
      if keyword_set(kps) then begin
         rates[i+nlon].sw=reform(data2[1,0:nlev-1])/24./3600. ; in K/s
         rates[i+nlon].lw=reform(data2[2,0:nlev-1])/24./3600. ; in K/s
      endif else begin
         if keyword_set(watts) then begin
            rates[i+nlon].sw=reform(data2[1,0:nlev-1])/24./3600.*cp ; in J/s/kg
            rates[i+nlon].lw=reform(data2[2,0:nlev-1])/24./3600.*cp ; in J/s/kg
         endif else begin
            rates[i+nlon].sw=rates[i+nlon].sflux ;W/m^2
            rates[i+nlon].lw=-rates[i+nlon].flux ;W/m^2
         endelse
      endelse
      ; fix to deal with rates too small the exponent gets incorrectly dropped (making them big again):
      test=where(data2[1,*] lt min(abs(data2[2,0:nlev-1])))
      if test[0] ne -1 then rates[i+nlon].sw[test[0]:nlev-1]=replicate(0.,nlev-test[0])
      rates[i].lat=-rates[i+nlon].lat
      rates[i].lon=rates[i+nlon].lon
      rates[i].fpres=rates[i+nlon].fpres
      rates[i].rflux=rates[i+nlon].rflux
      rates[i].flux=rates[i+nlon].flux
      rates[i].sflux=rates[i+nlon].sflux
      rates[i].pres=rates[i+nlon].pres
      rates[i].p0=rates[i+nlon].p0
      rates[i].sw=rates[i+nlon].sw
      rates[i].lw=rates[i+nlon].lw
   endfor
   istart=2*nlon
endif else begin
   ; can just read through all entries
   istart=0
endelse
for i=istart,nent-1 do begin
   readf,17,read1
   rates[i].lat=float(strmid(read1[0],23,6))
   rates[i].lon=float(strmid(read1[0],47,6))
   readf,17,data1
   rates[i].fpres=reform(data1[0,0:nlev-1])
   rates[i].rflux=reform(data1[1,0:nlev-1]-data1[2,0:nlev-1])
   rates[i].flux=reform(data1[3,0:nlev-1])
   rates[i].sflux=reform(data1[4,0:nlev-1])
   readf,17,read4
   readf,17,data2
   readf,17,read5
   rates[i].pres=reform(data2[0,0:nlev-1])  ; in bar (or mbar, if keyword set)
   rates[i].p0=data2[0,nlev]
   if keyword_set(kps) then begin
      rates[i].sw=reform(data2[1,0:nlev-1])/24./3600. ; in K/s
      rates[i].lw=reform(data2[2,0:nlev-1])/24./3600. ; in K/s
   endif else begin
      if keyword_set(watts) then begin
         rates[i].sw=reform(data2[1,0:nlev-1])/24./3600.*cp ; in J/s/kg
         rates[i].lw=reform(data2[2,0:nlev-1])/24./3600.*cp ; in J/s/kg
      endif else begin
         rates[i].sw=rates[i].sflux ;W/m^2
         rates[i].lw=-rates[i].flux ;W/m^2
      endelse
   endelse
   ; fix to deal with rates too small the exponent gets incorrectly dropped (making them big again):
   test=where(data2[1,0:nlev-1] lt min(abs(data2[2,0:nlev-1])))
   if test[0] ne -1 then rates[i].sw[test[0]:nlev-1]=replicate(0.,nlev-test[0])
endfor
close,17

if keyword_set(psfile) then begin
   if (size(psfile))[1] eq 2 then filename='plot.eps' else filename=psfile
   !x.thick=5
   !y.thick=5
   !p.thick=5
   if not (keyword_set(oned) or keyword_set(plevel)) then psopen,filename,/enc,/color,xsize=25,ysize=35 $
     else psopen,filename,/enc,/color,xsize=25,ysize=20
   !p.font=0
   cthick=6
endif else begin
;   device,true_color=24,decomposed=0;,retain=2
;  !p.font=-1
   cthick=3
endelse
!p.charsize=1.5

if not (keyword_set(plevel) or keyword_set(ape)) then begin

   if keyword_set(mbar) then begin
      rates.pres/=1.e3          ;convert to bar
      rates.fpres/=1.e3         ;convert to bar
   endif

   if not keyword_set(kps) then begin
      pmin=min(rates.fpres)
      if pmin eq 0 then begin
         ;correction for when not enough sig figs are printed and P=0 appears
         oomest=nlev*alog10(rates[0].pres[1:nlev-1]/rates[0].pres[0:nlev-2])
         oomest=median(oomest)
         p0est=rates[0].pres[nlev-1]/10.^(-1.*oomest/2./nlev)
         vert=p0est*get_sigma(oomest,nlev)
         vert=[vert[0]/2.,0.5*(vert[1:nlev-1]+vert[0:nlev-2])]
         pmin=min(vert,max=pmax)
      endif else begin
         pmax=max(rates.fpres)
         vert=rates[0].fpres
      endelse
   endif else begin
      pmin=min(rates.pres)
      if pmin eq 0 then begin
         ;correction for when not enough sig figs are printed and P=0 appears
         oomest=nlev*alog10(rates[0].pres[1:nlev-1]/rates[0].pres[0:nlev-2])
         oomest=median(oomest)
         p0est=rates[0].pres[nlev-1]/10.^(-1.*oomest/2./nlev)
         vert=p0est*get_sigma(oomest,nlev)
         pmin=min(vert,max=pmax)
      endif else begin
         pmax=max(rates.pres)
         vert=rates[0].pres
      endelse
   endelse

   loadct,0,/silent
   ttrate=fltarr(3,nlev)        ;lw, sw, total

   if not keyword_set(oned) then begin
      ratelw=fltarr(6,nlev)
      ratesw=fltarr(6,nlev)
      trate=fltarr(6,nlev)
      for i=0,5 do begin
      ;slice into equal area portions
         locs=where((cos(rates.lat*!pi/180.)*cos((rates.lon-sslon)*!pi/180.) le (1.-i/3.)) $
                    and (cos(rates.lat*!pi/180.)*cos((rates.lon-sslon)*!pi/180.) gt (2./3-i/3.)))
;   dA=(2.*!pi/nlon)*(!pi/nlat);*radea^2
;   area=4.*!pi/6.;*radea^2
;   area=total(dA*cos(rates[locs].lat*!pi/180.))
         area=total(cos(rates[locs].lat*!pi/180.))
         for j=0,nlev-1 do begin
            if not keyword_set(watts) then begin
               ratelw[i,j]=total(rates[locs].lw[j]*cos(rates[locs].lat*!pi/180.))/area ; int(rates dA)/int(dA)
               ratesw[i,j]=total(rates[locs].sw[j]*cos(rates[locs].lat*!pi/180.))/area
               trate[i,j]=total((rates[locs].lw[j]+rates[locs].sw[j])*cos(rates[locs].lat*!pi/180.))/area
            endif else begin
               ratelw[i,j]=total(rates[locs].lw[j]*cos(rates[locs].lat*!pi/180.)) $
                           *radea^2*(2.*!pi/nlon)*(!pi/nlat)/grav*dP[j]*1.e5 ; W
               ratesw[i,j]=total(rates[locs].sw[j]*cos(rates[locs].lat*!pi/180.)) $
                           *radea^2*(2.*!pi/nlon)*(!pi/nlat)/grav*dP[j]*1.e5 ; W
               trate[i,j]=total((rates[locs].lw[j]+rates[locs].sw[j])*cos(rates[locs].lat*!pi/180.)) $
                          *radea^2*(2.*!pi/nlon)*(!pi/nlat)/grav*dP[j]*1.e5 ; W
            endelse
         endfor
      endfor
   endif
   
   totarea=total(cos(rates.lat*!pi/180.))
   if not keyword_set(watts) then begin
      for j=0,nlev-1 do begin
         ttrate[0,j]=total(rates.lw[j]*cos(rates.lat*!pi/180.))/totarea
         ttrate[1,j]=total(rates.sw[j]*cos(rates.lat*!pi/180.))/totarea
         ttrate[2,j]=total((rates.lw[j]+rates.sw[j])*cos(rates.lat*!pi/180.))/totarea
      endfor
      if not keyword_set(kps) then begin
         xtitc='Average cooling flux [W/m^2]'
         xtith='Average heating flux [W/m^2]'
      endif else begin
         xtitc='Average cooling rate [K/s]'
         xtith='Average heating rate [K/s]'
      endelse
   endif else begin
      ttrate[0,0]=total(rates.lw[0]*cos(rates.lat*!pi/180.))*radea^2 $
                  *(2.*!pi/nlon)*(!pi/nlat)/grav*1.e5*vert[0]  ; W
      ttrate[1,0]=total(rates.sw[0]*cos(rates.lat*!pi/180.))*radea^2 $
                  *(2.*!pi/nlon)*(!pi/nlat)/grav*1.e5*vert[0]  ; W
      ttrate[2,0]=total((rates.lw[0]+rates.sw[0])*cos(rates.lat*!pi/180.))*radea^2 $
                  *(2.*!pi/nlon)*(!pi/nlat)/grav*1.e5*vert[0]  ; W
      for j=1,nlev-1 do begin
         ttrate[0,j]=ttrate[0,j-1]+total(rates.lw[j]*cos(rates.lat*!pi/180.)) $
                     *radea^2*(2.*!pi/nlon)*(!pi/nlat)/grav*1.e5*(vert[j]-vert[j-1]) ; W
         ttrate[1,j]=ttrate[1,j-1]+total(rates.sw[j]*cos(rates.lat*!pi/180.)) $
                     *radea^2*(2.*!pi/nlon)*(!pi/nlat)/grav*1.e5*(vert[j]-vert[j-1]) ; W
         ttrate[2,j]=ttrate[2,j-1]+total((rates.lw[j]+rates.sw[j])*cos(rates.lat*!pi/180.)) $
                     *radea^2*(2.*!pi/nlon)*(!pi/nlat)/grav*1.e5*(vert[j]-vert[j-1]) ; W
      endfor
      xtitc='Cumulative cooling rate, from P=0 [W]'
      xtith='Cumulative heating rate, from P=0 [W]'
   endelse

   if not keyword_set(oned) then begin
      xmin=min([reform(ratelw),reform(ratesw),reform(trate),reform(ttrate)],max=xmax)
      noz=where(ratelw ne 0)
      xmin=min(abs(ratelw[noz]))
   endif else begin
      xmax=max(abs(ttrate))
      noz=where(ttrate ne 0)
      xmin=min(abs(ttrate[noz]))
   endelse
   xmax=10.^(ceil(alog10(xmax)))
   xmin=10.^(floor(alog10(xmin)))
   if keyword_set(xrange) then begin
      xmin=xrange[0]
      xmax=xrange[1]
   endif

   gxmax=max(ttrate)
   noz=where(ttrate ne 0)
   gxmin=min(abs(ttrate[noz]))

;ccolor=findgen(6)*255./5.
   ccolor=[65,50,100,240,200,110]
   temparr=fltarr(nlev)

   if not keyword_set(oned) then begin
      plot,[0,0],[0,0],xr=[xmax,xmin],/xlog,yr=[pmax,pmin],/ylog,xstyle=1,/ystyle,$
           xtit=xtitc,position=[0.15,0.1,0.55,0.5],ytit='Pressure [bar]'
      xyouts,xmax/2.,pmax/3.,'Global average'
   endif else begin
      plot,[0,0],[0,0],xr=[xmax,xmin],/xlog,yr=[pmax,pmin],/ylog,xstyle=1,/ystyle,$
           xtit=xtitc,position=[0.15,0.1,0.55,0.9],ytit='Pressure [bar]'
      if keyword_set(dots) then al_legend,['SW rate','LW rate','total rate'],psym=[6,4,2],/bottom
   endelse
   if keyword_set(mheat) then begin
      loadct,5,/silent
      oplot,-ttrate[0,*],vert,linestyle=2,color=50
      oplot,-ttrate[1,*],vert,linestyle=1,color=50
      oplot,-ttrate[2,*],vert,linestyle=0,color=50
      oplot,-1.*(ttrate[2,*]+mheat[1,*]),vert,thick=8
      loadct,0,/silent
   endif else begin
      oplot,-ttrate[0,*],vert,linestyle=2
      oplot,-ttrate[1,*],vert,linestyle=1
      oplot,-ttrate[2,*],vert,linestyle=0
   endelse
   if keyword_set(dots) then begin
      lon=lonely(reform(ttrate[0,*]))
      oplot,-lon*reform(ttrate[0,*]),vert,psym=4
      lon=lonely(reform(ttrate[2,*]))
      oplot,-lon*reform(ttrate[2,*]),vert,psym=2
      if keyword_set(zeros) then begin
         zero=where(ttrate[0,*] eq 0) ;, complement=nzero)
         if zero[0] ne -1 then begin
            temparr[zero]=xmin
            oplot,temparr,vert,psym=4
            temparr[zero]=0.
         endif
         zero=where(ttrate[1,*] eq 0) ;, complement=nzero)
         if zero[0] ne -1 then begin
            temparr[zero]=xmin
            oplot,temparr,vert,psym=6
            temparr[zero]=0.
         endif
         zero=where(ttrate[2,*] eq 0) ;, complement=nzero)
         if zero[0] ne -1 then begin
            temparr[zero]=xmin
            oplot,temparr,vert,psym=2
            temparr[zero]=0.
         endif
      endif
   endif

   if not keyword_set(oned) then begin
      plot,[0,0],[0,0],xr=[xmin,xmax],/xlog,yr=[pmax,pmin],/ylog,ytickname=replicate(' ',10),$
           xtit=xtith,position=[0.55,0.1,0.95,0.5],/noerase,/xstyle,/ystyle
   endif else begin
      plot,[0,0],[0,0],xr=[xmin,xmax],/xlog,yr=[pmax,pmin],/ylog,ytickname=replicate(' ',10),$
           xtit=xtith,position=[0.55,0.1,0.95,0.9],/noerase,/xstyle,/ystyle
   endelse
   if keyword_set(mheat) then begin
      loadct,5,/silent
      oplot,ttrate[0,*],vert,linestyle=2,color=50
      oplot,ttrate[1,*],vert,linestyle=1,color=50
      oplot,ttrate[2,*],vert,linestyle=0,color=50
      oplot,mheat[0,*],vert,linestyle=3,color=100  ; "true" value
      oplot,mheat[1,*],vert,color=100              ; value used
      oplot,ttrate[2,*]+mheat[1,*],vert,thick=8
      al_legend,['radiative','optical','infrared','magnetic','total'],$
                linestyle=[0,1,2,0,0],color=[50,50,50,100,0],/bottom,/right,charsize=1.25
      loadct,0,/silent
   endif else begin
      oplot,ttrate[0,*],vert,linestyle=2
      oplot,ttrate[1,*],vert,linestyle=1
      oplot,ttrate[2,*],vert,linestyle=0
      if not keyword_set(watts) then $
         al_legend,['SW rate','LW rate','total rate'],linestyle=[1,2,0],/bottom,/right,charsize=1.25 $
      else al_legend,['SW rate','LW rate','total rate'],linestyle=[1,2,0],/right,charsize=1.25
   endelse
   if keyword_set(dots) then begin
      lon=lonely(reform(ttrate[0,*]))
      oplot,lon*reform(ttrate[0,*]),vert,psym=4
      lon=lonely(reform(ttrate[2,*]))
      oplot,lon*reform(ttrate[2,*]),vert,psym=2
      if keyword_set(zeros) then begin
         zero=where(ttrate[0,*] eq 0) ;, complement=nzero)
         if zero[0] ne -1 then begin
            temparr[zero]=xmin
            oplot,temparr,vert,psym=4
            temparr[zero]=0.
         endif
         zero=where(ttrate[1,*] eq 0) ;, complement=nzero)
         if zero[0] ne -1 then begin
            temparr[zero]=xmin
            oplot,temparr,vert,psym=6
            temparr[zero]=0.
         endif
         zero=where(ttrate[2,*] eq 0) ;, complement=nzero)
         if zero[0] ne -1 then begin
            temparr[zero]=xmin
            oplot,temparr,vert,psym=2
            temparr[zero]=0.
         endif
      endif
   endif

   if not keyword_set(oned) then begin
      plot,[0,0],[0,0],xr=[xmax,xmin],/xlog,yr=[pmax,pmin],/ylog,xstyle=1,/noerase,$
           ytit='Pressure [bar]',position=[0.15,0.5,0.55,0.9],xtickname=replicate(' ',10),/ystyle
      for i=5,0,-1 do begin
         if i eq 5 then loadct,5,/silent
         if i eq 4 then loadct,4,/silent
         if i eq 1 then loadct,5,/silent
         oplot,-ratelw[i,*],vert,color=ccolor[i],linestyle=2
         oplot,-trate[i,*],vert,color=ccolor[i]
         if keyword_set(dots) then begin
            lon=lonely(reform(ratelw[i,*]))
            oplot,-reform(ratelw[i,*])*lon,vert,color=ccolor[i],psym=4
            lon=lonely(reform(trate[i,*]))
            oplot,-reform(trate[i,*])*lon,vert,color=ccolor[i],psym=2
            if keyword_set(zeros) then begin
               zero=where(ratesw[i,*] eq 0) ;, complement=nzero)
               if zero[0] ne -1 then begin	
                  temparr[zero]=xmin
                  oplot,temparr,vert,psym=6,color=ccolor[i]
                  temparr[zero]=0.
               endif
               zero=where(ratelw[i,*] eq 0) ;, complement=nzero)
               if zero[0] ne -1 then begin
                  temparr[zero]=xmin
                  oplot,temparr,vert,psym=4,color=ccolor[i]
                  temparr[zero]=0.
               endif
               zero=where(trate[i,*] eq 0) ;, complement=nzero)
               if zero[0] ne -1 then begin
                  temparr[zero]=xmin
                  oplot,temparr,vert,psym=2,color=ccolor[i]
                  temparr[zero]=0.
               endif
            endif
         endif
      endfor
      tcolors=[8,72,88,184,216,232]
      loadct,38,/silent
      if keyword_set(psfile) then mu0='!9'+string("155B)+'!X!D*!N' $ ;"
      else mu0='!4'+string("154B)+'!X!D!3'+string("52B)+'!X!N' ;
      al_legend,[mu0+' = 1 to 2/3',mu0+' = 2/3 to 1/3',mu0+' = 1/3 to 0',$
                 mu0+' = 0 to -1/3',mu0+' = -1/3 to -2/3',mu0+' = -2/3 to -1'],$
                /bottom,outline_color=0.,charsize=1.25,textcolors=tcolors
      loadct,0,/silent
      axis,xaxis=1,/xlog,xr=[xmax,xmin],/xstyle ;,xtitle=xtitc

      plot,[0,0],[0,0],xr=[xmin,xmax],/xlog,yr=[pmax,pmin],/ylog,ytickname=replicate(' ',10),xstyle=1,$
           position=[0.55,0.5,0.95,0.9],/noerase,xtickname=replicate(' ',10),/ystyle
;   loadct,25,/silent
      for i=5,0,-1 do begin
         if i eq 5 then loadct,5,/silent
         if i eq 4 then loadct,4,/silent
         if i eq 1 then loadct,5,/silent
         oplot,ratesw[i,*],vert,color=ccolor[i],linestyle=1
         oplot,ratelw[i,*],vert,color=ccolor[i],linestyle=2
         oplot,trate[i,*],vert,color=ccolor[i]
         if keyword_set(dots) then begin
            lon=lonely(reform(ratelw[i,*]))
            oplot,reform(ratelw[i,*])*lon,vert,color=ccolor[i],psym=4
            lon=lonely(reform(trate[i,*]))
            oplot,reform(trate[i,*])*lon,vert,color=ccolor[i],psym=2
            if keyword_set(zeros) then begin
               zero=where(ratesw[i,*] eq 0) ;, complement=nzero)
               if zero[0] ne -1 then begin
                  temparr[zero]=xmin
                  oplot,temparr,vert,psym=6,color=ccolor[i]
                  temparr[zero]=0.
               endif
               zero=where(ratelw[i,*] eq 0) ;, complement=nzero)
               if zero[0] ne -1 then begin
                  temparr[zero]=xmin
                  oplot,temparr,vert,psym=4,color=ccolor[i]
                  temparr[zero]=0.
               endif
               zero=where(trate[i,*] eq 0) ;, complement=nzero)
               if zero[0] ne -1 then begin
                  temparr[zero]=xmin
                  oplot,temparr,vert,psym=2,color=ccolor[i]
                  temparr[zero]=0.
               endif
            endif
         endif
      endfor
      loadct,0,/silent
      if keyword_set(dots) then al_legend,['SW rate','LW rate','total rate'],psym=[6,4,2],/bottom,/right
      axis,xaxis=1,/xlog,xr=[xmin,xmax],/xstyle ;,xtitle=xtith
   endif

   pos=where(ttrate[2,*] gt 0.,complement=neg)
   if not keyword_set(watts) then begin
      print,'Average SW: ',total(ttrate[1,*]*vert)/total(vert)
      print,'Average LW: ',total(ttrate[0,*]*vert)/total(vert)
      print,'Average heating: ',total(ttrate[2,pos]*vert[pos])/total(vert[pos])
      print,'Average cooling: ',total(ttrate[2,neg]*vert[neg])/total(vert[neg])
      print,''
      print,'Average total rate: ',total(ttrate[2,*]*vert)/total(vert)
   endif else begin
      print,'Total SW: ',total(ttrate[1,*])
      print,'Total LW: ',total(ttrate[0,*])
      print,'Total heating: ',total(ttrate[2,pos])
      print,'Total cooling: ',total(ttrate[2,neg])
      print,''
      print,'Total rate: ',total(ttrate[2,*])
   endelse

   loadct,0,/silent

endif else begin

   slat=(rates.lat)[sort(rates.lat)]
   slat=slat[uniq(slat)]
   slon=(rates.lon)[sort(rates.lon)]
   slon=slon[uniq(slon)]

   if keyword_set(ape) then begin
      OpenR, 17, 'fort.26'
      ReadF, 17, nlat,nlon,nlev	  
      xy=fltarr(6,nlat,nlon,nlev)   
      ReadF, 17, xy
      Close, 17

      dm=fltarr(nlat,nlon,nlev)
      for ilev=0,nlev-1 do begin
         for ilat=0,nlat-1 do begin
            dm[ilat,*,ilev]=(radea^2/grav)*dP[ilev]*1.e5*cos(xy[1,ilat,0,0]*!pi/180.)*(!pi/nlat)*(2.*!pi/nlon) ;kg
         endfor
      endfor
      tm=1./(total(dm/xy[5,*,*,*])/total(dm))
      print,'Tm=',tm
      ke=xy[3,*,*,*]^2 *dm  ; J
      heat=fltarr(nlat,nlon,nlev)

      for i=0,nent-1 do begin
         lati=where(slat eq rates[i].lat)
         loni=where(slon eq rates[i].lon)
         ; note, the latitudes here go from southernmost to northernmost,
         ; the opposite of the pattern in xy (from fort.26), so flip them
         if keyword_set(watts) then begin
            heat[nlat-1-lati,loni,*]=(rates[i].sw+rates[i].lw)*cos(rates[i].lat*!pi/180.) $
                              *radea^2*(2.*!pi/nlon)*(!pi/nlat)/grav*dP*1.e5  ; W
         endif else begin
            print,'Error: Set keyword watts'
            stop
         endelse
      endfor         
      print,'Net radiative heating:',total(heat),' W'
      print,'Total radiative APE gen:',total((1.-tm/reform(xy[5,*,*,*]))*heat),' W'

   endif

   if keyword_set(plevel) then begin
      heat=fltarr(nlat-1,nlon+1)
      swh=fltarr(nlat-1,nlon+1)

      for i=0,nent-1 do begin
         lati=where(rates[i].lat eq slat)
         loni=where(rates[i].lon eq slon)
         if keyword_set(watts) then begin
            heat[lati,loni]=(rates[i].sw[plevel]+rates[i].lw[plevel])*cos(rates[i].lat*!pi/180.) $
                            *radea^2*(2.*!pi/nlon)*(!pi/nlat)/grav*dP[plevel]*1.e5 ; W                         
            swh[lati,loni]=rates[i].sw[plevel]*cos(rates[i].lat*!pi/180.) $
                           *radea^2*(2.*!pi/nlon)*(!pi/nlat)/grav*dP[plevel]*1.e5 ; W                         
         endif else begin
            heat[lati,loni]=rates[i].sw[plevel]+rates[i].lw[plevel]
            swh[lati,loni]=rates[i].sw[plevel]
         endelse
      endfor
      heat[*,nlon]=heat[*,0]
      swh[*,nlon]=swh[*,0]

      if keyword_set(kps) then ctitle='Heating rate [K/s]' else $
         if keyword_set(watts) then ctitle='Heating rate [W]' else $
            ctitle='Heating rate [W/m^2]'

      lat=fltarr(nlat-1,nlon+1)
      for ilat=0,nlat-2 do lat[ilat,*]=slat[ilat]
      lon=fltarr(nlat-1,nlon+1)
      for ilon=0,nlon-1 do lon[*,ilon]=slon[ilon]
      lon[*,nlon]=slon[0]

      qmin=min(heat,max=qmax)
      MAP_SET, /Miller_CYLINDRICAL,0,0,/ISOTROPIC ;, TITLE=ctitle
      print,ctitle
      loadct, 3,/silent
      cbottom=0.
      nlevels=45.
      step=(qmax-qmin)/nlevels
      mylevels=indgen(nlevels)*step+qmin
      ccolors=findgen(nlevels)*(255.-cbottom)/(nlevels-1)+cbottom
      contour,heat,lon,lat, /overplot,/cell_fill,c_colors=ccolors,levels=mylevels,$
              c_thick=cthick
      toss=where(swh ne 0.)
      if toss[0] ne -1 then begin
         ncontours=5
         myclevels=max(swh)*findgen(ncontours+1)/ncontours
         loadct,1,/silent
         contour, swh,lon,lat, /overplot,c_thick=cthick,levels=myclevels,/follow,$
                  c_charsize=1.25,color=175
         loadct,3,/silent
      endif
      if qmin ne qmax then begin
         colorbar,position=[0.1,0.07,0.90,0.1],range=[qmin,qmax],format='(e9.2)',$
                  charsize=2,bottom=cbottom,ncolors=255-cbottom,divisions=4
      endif
      loadct,0,/silent
;   igcm,3,plevel,vfrac=0.2
;   contour, swh,lon,lat, /overplot,c_thick=3,color=0,levels=mylevels
      print,qmin,qmax
      print,myclevels
   endif

endelse

if keyword_set(psfile) then begin
  !x.thick=1
  !y.thick=1
  !p.thick=1
  psclose
  if keyword_set(plevel) then begin
     spawn,'gs -r300 -dEPSCrop -dTextAlphaBits=4 -sDEVICE=png16m -sOutputFile=plot.png -dBATCH -dNOPAUSE plot.eps'
     spawn,'convert plot.png eps3:plot.eps'
  endif
endif

;stop

end
