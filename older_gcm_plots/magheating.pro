pro magheating,p0,oom,bfield,tdrag_min,levplot,tmap=tmap,ncontours=ncontours,$
               psfile=psfile,coverplot=coverplot,vert=vert,ga=ga,radea=radea,mh=mh,$
               ape=ape,kediss=kediss,wms=wms,metals=metals

gascon=3523.  ; J/kg/K
if not keyword_set(radea) then radea = 1.e8  ; m
if not keyword_set(ga) then ga=8.  ; m/s
print,'           CHECK: gascon, radea, ga:',gascon,radea,ga
print,'  DOUBLE-CHECK: bfield is',bfield
;tdrag_min = 0.005 * (2.*!pi)/ww

; radea in m
; ga in m/s
; p0 in bar
;
; for horizontal plot, mh returns a [2,nlat,nlon] array with the true
; [0] and used [1] mag heating rates.  (In the case of tdrag_min=0
; these would be the same.)
; for a vertical plot, mh is [2,nlev]

; Read File with 3D fields
  OpenR, 17, 'fort.26'
  ReadF, 17, nlat,nlon,nlev	  
;	  print,nlat,nlon,nlev 
	  xy=fltarr(6,nlat,nlon,nlev)   
  ReadF, 17, xy
  Close, 17

;get surface presssures
  openr,18,'fort.50'
  readf,18,nlat,nlon
  ab=fltarr(3,nlat,nlon)
  readf,18,ab
  close,18

if not keyword_set(psfile) then begin
  device,true_color=24,decomposed=0,retain=2
  cthick=3
endif else begin
   if (size(psfile))[1] eq 2 then filename='plot.eps' else filename=psfile
   if keyword_set(zave) then xsz=40 else xsz=25
   psopen,filename,/enc,/color,xsize=xsz,ysize=20;,bits_per_pixel=24
  !p.font=0
  !x.thick=6
  !y.thick=6
  !p.thick=6
  !p.charsize=2
  cthick=6
endelse

sigma=get_sigma(oom,nlev)
pres=p0*sigma
dP=fltarr(nlev)
for ilev=1,nlev-2 do dP[ilev]=0.5*(pres[ilev+1]-pres[ilev-1]) ; to match swdp
dP[0]=0.5*(pres[1])
dP[nlev-1]=0.5*(p0-pres[nlev-2])

if not (keyword_set(vert) or keyword_set(ape)) then begin

   latoff=0
   lonoff=0

   ; Pick vertical level (lower layer is nlev-1)
   nl=levplot                   ; from function command line
    
   lon=reform(xy[0,*,*,nl]) & lat=reform(xy[1,*,*,nl]) 
   u=reform(xy[3,*,*,nl]) & v=reform(xy[4,*,*,nl]) & temp=reform(xy[5,*,*,nl])
   sp=reform(ab[2,*,*])
   sp=(sp+1.)*p0

   snl=sigma[nl]
   dmass=dP[nl]*1.e5/ga*radea^2  ; kg

   mh=fltarr(2,nlat,nlon)
   mh[0,*,*]=u^2                 ;(m/s)^2
   mh[1,*,*]=u^2                 ;(m/s)^2
   for ilat=0,nlat-1 do begin
      ang_f=abs(cos(!pi/2.-xy[1,ilat,0,0]*!pi/180.))
      rho=sp[ilat,*]*snl/gascon/temp[ilat,*]*1.e2 ; g/cm^3
      tdrag=tmag(rho,temp[ilat,*],bfield,metals=metals)/ang_f         ; s
      if not keyword_set(wms) then $
         mh[0,ilat,*]=mh[0,ilat,*]/tdrag $
                      *dmass*cos(lat[ilat,0]*!pi/180.)*(!pi/nlat)*(2.*!pi/nlon) $; W
      else $
         mh[0,ilat,*]=mh[0,ilat,*]/tdrag*dmass/radea^2  ;W/m^2
      for ilon=0,nlon-1 do tdrag[ilon]=max([tdrag[ilon],tdrag_min])
      if not keyword_set(wms) then $
         mh[1,ilat,*]=mh[1,ilat,*]/tdrag $
                      *dmass*cos(lat[ilat,0]*!pi/180.)*(!pi/nlat)*(2.*!pi/nlon) $; W
      else $
         mh[1,ilat,*]=mh[1,ilat,*]/tdrag*dmass/radea^2  ;W/m^2
   endfor
   qmin=min(mh[1,*,*],max=qmax)
   qtp=fltarr(nlat,nlon+1)
   qtp[*,0:nlon-1]=reform(mh[1,*,*])
   qtp[*,nlon]=reform(mh[1,*,0])
   
   flon=fltarr(nlat,nlon+1)
   flon[*,0:nlon-1]=lon
   flon[*,nlon]=reform(lon[*,0]+360.)
   flat=fltarr(nlat,nlon+1)
   flat[*,0:nlon-1]=lat
   flat[*,nlon]=reform(lat[*,0])

   if not keyword_set(wms) then units=' W' else units=' W/m^2'

   if not (keyword_set(tmap) or keyword_set(coverplot)) then begin
      MAP_SET, /Miller_CYLINDRICAL,latoff,lonoff,/ISOTROPIC
      loadct, 1,/silent
      cbottom=0.
      nlevels=45.
      step=(qmax-qmin)/nlevels
      mylevels=indgen(nlevels)*step+qmin
      ccolors=findgen(nlevels)*(255.-cbottom)/(nlevels-1)+cbottom

      contour, qtp,flon,flat, /overplot,/cell_fill,/closed,c_colors=ccolors,levels=mylevels
      if qmin ne qmax then begin
         colorbar,position=[0.1,0.07,0.90,0.1],range=[qmin,qmax],format='(e8.1)',$
                  charsize=2,bottom=cbottom,ncolors=255-cbottom
      endif
      loadct,0,/silent
     
      test=where(mh[0,*,*] ne mh[1,*,*])
      if test[0] ne -1 then begin
         oqtp=fltarr(nlat,nlon+1)
         oqtp[*,0:nlon-1]=reform(mh[0,*,*])
         oqtp[*,nlon]=reform(mh[0,*,0])
         oqtp=oqtp-qtp
         nlevels=4
         qmax=max(oqtp,min=qmin)
         mylevels=[0.25,0.5,0.75,1.]*qmax
         contour,oqtp,flon,flat,/over,levels=mylevels,c_thick=3
         print,'Max extra heating:',qmax,units
      endif
   endif else begin
      !p.thick=3
      if keyword_set(tmap) then igcm,3,nl,vfrac=0.2
      !p.thick=6
      if not keyword_set(ncontours) then ncontours=5
      mylevels=qmax*findgen(ncontours+1)/ncontours
      loadct,1,/silent
      contour, qtp,flon,flat, /overplot,c_thick=cthick,color=125,levels=mylevels,thick=8 ;,/follow,c_charsize=1.5
      loadct,0,/silent
      !p.thick=1
      print,'Min, max mag heating:',qmin,qmax,units
      print,'contours at % of max:',mylevels[1:ncontours-1]/qmax
      print,'Total mag heating at this level:',total(qtp)
      test=where(mh[0,*,*] ne mh[1,*,*])
      if test[0] ne -1 then begin
         oqtp=fltarr(nlat,nlon+1)
         oqtp[*,0:nlon-1]=reform(mh[0,*,*])
         oqtp[*,nlon]=reform(mh[0,*,0])
         oqtp=oqtp-qtp
         oqmax=max(oqtp,min=oqmin)
         mylevels=mylevels*oqmax/qmax
         loadct,1,/silent
         contour,oqtp,flon,flat,/over,levels=mylevels,c_thick=cthick,color=200
         loadct,0,/silent
         print,'Max extra heating:',oqmax,units
         print,'Total extra heating:',total(mh[0,*,*]-mh[1,*,*]),units
      endif
   endelse

endif else begin

   lon=reform(xy[0,*,*,*]) & lat=reform(xy[1,*,*,*]) 
   u=reform(xy[3,*,*,*]) & v=reform(xy[4,*,*,*]) & temp=reform(xy[5,*,*,*])
   sp=reform(ab[2,*,*])
   sp=(sp+1.)*p0

   if not keyword_set(ape) then begin

      mh=fltarr(2,nlev)
      
      ke=u^2              ;(m/s)^2
      for ilev=0,nlev-1 do begin
         for ilat=0,nlat-1 do begin
            ang_f=abs(cos(!pi/2.-xy[1,ilat,0,0]*!pi/180.))
            rho=reform(sp[ilat,*]*sigma[ilev]/gascon/temp[ilat,*,ilev]*1.e2) ; g/cm^3
            tdrag=tmag(rho,reform(temp[ilat,*,ilev]),bfield,metals=metals)/ang_f      ; s
            mh[0,ilev]+=total(reform(ke[ilat,*,ilev])/tdrag) $           ; m^2/s^3  (or W/kg)
                        *cos(lat[ilat,0,0]*!pi/180.)*(!pi/nlat)*(2.*!pi/nlon) $
                        *(radea^2/ga)*dP[ilev]*1.e5 ; W
            for ilon=0,nlon-1 do tdrag[ilon]=max([tdrag[ilon],tdrag_min])
            mh[1,ilev]+=total(reform(ke[ilat,*,ilev])/tdrag) $ ; m^2/s^3  (or W/kg)
                        *cos(lat[ilat,0,0]*!pi/180.)*(!pi/nlat)*(2.*!pi/nlon) $
                        *(radea^2/ga)*dP[ilev]*1.e5 ; W
         endfor
      endfor
      plot,mh[1,*],sigma*p0,/ylog,yr=[p0,min(sigma)*p0],/xlog,xtit='Magnetic heating [W]',ytit='Pressure [bar]'
      oplot,mh[0,*],sigma*p0,linestyle=2

      print,'Total magnetic heating: ',total(mh[1,*]),' W'

   endif else begin

      kediss=fltarr(nlat,nlon,nlev)
      apegen=fltarr(nlat,nlon,nlev)
      dm=fltarr(nlat,nlon,nlev)
      for ilev=0,nlev-1 do begin
         for ilat=0,nlat-1 do begin
            dm[ilat,*,ilev]=(radea^2/ga)*dP[ilev]*1.e5*cos(lat[ilat,0,0]*!pi/180.)*(!pi/nlat)*(2.*!pi/nlon) ;kg
         endfor
      endfor
      tm=1./(total(dm/temp)/total(dm))
      print,'Tm=',tm

      ke=u^2              ;(m/s)^2
      for ilev=0,nlev-1 do begin
         for ilat=0,nlat-1 do begin
            ang_f=abs(cos(!pi/2.-xy[1,ilat,0,0]*!pi/180.))
            rho=reform(sp[ilat,*]*sigma[ilev]/gascon/temp[ilat,*,ilev]*1.e2) ; g/cm^3
            tdrag=tmag(rho,reform(temp[ilat,*,ilev]),bfield,metals=metals)/ang_f       ; s
            for ilon=0,nlon-1 do tdrag[ilon]=max([tdrag[ilon],tdrag_min])
            kediss[ilat,*,ilev]=-1.*ke[ilat,*,ilev]/tdrag*dm[ilat,*,ilev]  ; W
            apegen[ilat,*,ilev]=(1.-tm/temp[ilat,*,ilev])*(ke[ilat,*,ilev]/tdrag)*dm[ilat,*,ilev]  ; W
         endfor
      endfor
      print,'Magnetic APE gen:',total(apegen),' W'
      print,'Magnetic KE diss:',total(kediss),' W'
      
   endelse

endelse

if keyword_set(psfile) then begin
   !p.font=-1
   !x.thick=1
   !y.thick=1
   !p.thick=1
   !p.charsize=1
   psclose
   if not keyword_set(vert) then begin
      spawn,'gs -r300 -dEPSCrop -dTextAlphaBits=4 -sDEVICE=png16m -sOutputFile=plot.png -dBATCH -dNOPAUSE plot.eps'
      spawn,'convert plot.png eps3:plot.eps'
   endif
endif

end
