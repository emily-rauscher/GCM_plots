pro vertvel,ga,gascon,radea,oom,p0,psfile=psfile,psurf=psurf,ww=ww,noplot=noplot,om=om,file=file

; Calculates vertical velocities throughout the atmosphere, from
; horizontal velocities in fort.26 and using the continuity equation,
; with surface pressure information from fort.50.
;
; Inputs:
;    ga: surface gravity, in m/s^2
;    gascon: specific gas constant, in SI
;    radea: planet radius, in m
;    oom: vertical extent of atmosphere, in orders of magnitude in pressure
;    p0: surface pressure, in bar
;
; Optional outputs:
;    om: (nlat,nlon,nlev) array with vertical velocity in bar/s
;    ww: (nlat,nlon,nlev) array with vertical velocity in m/s
;    file: if keyword set, will create fort.26.uvw
;          which will be an updated version of fort.26, containing
;          vertical velocity (in m/s) as an additional column,
;          after the u and v entries
;
; Optional keywords:
;    psfile: if set, will create plot.eps (or can set = string to use
;            as filename)
;    psurf: set to read in 2D surface pressure information from
;           fort.50
;    noplot: if set, no plot will be produced


if not keyword_set(noplot) then begin
   if not keyword_set(psfile) then begin
      device,true_color=24,decomposed=0,retain=2
      !p.font=-1
      !p.charsize=1.5
   endif else begin
      if (size(psfile))[1] eq 2 then filename='plot.eps' else filename=psfile
      psopen,filename,/enc,/color,xsize=35,ysize=25 ;,bits_per_pixel=24
      !p.font=0
      !p.charsize=1.5
      !p.thick=4
      !x.thick=4
      !y.thick=4
   endelse
endif

  fstring=''
  OpenR, 17, 'fort.26'
  ReadF, 17, nlat,nlon,nlev	  
	  xy=fltarr(6,nlat,nlon,nlev)   
  ReadF, 17, xy
  readf,17,fstring
  readf,17,fstring
  Close, 17
  lon=reform(xy[0,*,*,*]) & lat=reform(xy[1,*,*,*]) 
  u=reform(xy[3,*,*,*]) & v=reform(xy[4,*,*,*]) 
  temp=reform(xy[5,*,*,*])

if keyword_set(psurf) then begin
   if not keyword_set(p0) then begin
      print,'MUST SET P0 IF USING PSURF!'
      stop
   endif
  openr, 18, 'fort.50'
  readf, 18, nlat,nlon
  ab=fltarr(3,nlat,nlon)
  readf, 18, ab
  close, 18
  sp=reform(ab[2,*,*])
  sp = (sp+1.)*p0
endif

sigma=get_sigma(oom,nlev)
pres=p0*sigma

;integrate continuity for omega (=0 at p=0)
om=fltarr(nlat,nlon,nlev)
dudx=fltarr(nlat,nlon,2)  ; 1=this level, 0=previous level
dvdy=fltarr(nlat,nlon,2)
for ilev=0,nlev-1 do begin
   ;need to worry about poles, and 360->0
   for ilat=0,nlat-1 do begin
      dudx[ilat,*,1]=deriv(cos(lat[ilat,0,0]*!pi/180.)*reform(lon[0,*,0]*!pi/180),reform(u[ilat,*,ilev]))/radea
      dudx[ilat,0,1]=(u[ilat,1,ilev]-u[ilat,nlon-1,ilev]) $
                     /(cos(lat[ilat,0,0]*!pi/180.)*radea*(2.*!pi+!pi/180.*(lon[0,1,0]-lon[0,nlon-1,0])))
      dudx[ilat,nlon-1,1]=(u[ilat,0,ilev]-u[ilat,nlon-2,ilev]) $
                          /(cos(lat[ilat,0,0]*!pi/180.)*radea*(2.*!pi+!pi/180.*(lon[0,0,0]-lon[0,nlon-2,0])))
   endfor
   for ilon=0,nlon/2-1 do begin
      dvdy[*,ilon,1]=deriv(reform(lat[*,0,0]*!pi/180.),reform(v[*,ilon,ilev]))/radea
      dvdy[0,ilon,1]=(-v[0,ilon+nlon/2,ilev]-v[1,ilon,ilev])/(!pi-(lat[0,0,0]+lat[1,0,0])*!pi/180.)/radea
      dvdy[nlat-1,ilon,1]=(v[nlat-1,ilon+nlon/2,ilev]+v[nlat-2,ilon,ilev]) $
                             /(!pi+(lat[nlat-1,0,0]+lat[nlat-2,0,0])*!pi/180.)/radea
   endfor
   for ilon=nlon/2,nlon-1 do begin
      dvdy[*,ilon,1]=deriv(reform(lat[*,0,0]*!pi/180.),reform(v[*,ilon,ilev]))/radea
      dvdy[0,ilon,1]=(-v[0,ilon-nlon/2,ilev]-v[1,ilon,ilev])/(!pi-(lat[0,0,0]+lat[1,0,0])*!pi/180.)/radea
      dvdy[nlat-1,ilon,1]=(v[nlat-1,ilon-nlon/2,ilev]+v[nlat-2,ilon,ilev]) $
                             /(!pi+(lat[nlat-1,0,0]+lat[nlat-2,0,0])*!pi/180.)/radea
   endfor
;   if ilev eq 0 then
;   om[*,*,ilev]=0.-pres[ilev]*(dudx[*,*,1]+dvdy[*,*,1]) else $
   if not keyword_set(psurf) then begin
      if ilev eq 0 then om[*,*,ilev]=0.-0.5*(dudx[*,*,1]+dvdy[*,*,1])*pres[ilev] else $  ;bar/s
         om[*,*,ilev]=om[*,*,ilev-1]-0.5*(dudx[*,*,1]+dudx[*,*,0]+dvdy[*,*,1]+dvdy[*,*,0])*(pres[ilev]-pres[ilev-1])
   endif else begin
      if ilev eq 0 then om[*,*,ilev]=0.-0.5*(dudx[*,*,1]+dvdy[*,*,1])*sigma[ilev]*sp[*,*] else $  ;bar/s
         om[*,*,ilev]=om[*,*,ilev-1]-0.5*(dudx[*,*,1]+dudx[*,*,0]+dvdy[*,*,1]+dvdy[*,*,0]) $
                      *(sigma[ilev]-sigma[ilev-1])*sp[*,*]
   endelse
   dudx[*,*,0]=dudx[*,*,1]
   dvdy[*,*,0]=dvdy[*,*,1]
endfor

;vertical velocity
ww=fltarr(nlat,nlon,nlev)
dw=fltarr(nlev)
nw=fltarr(nlev)
gw=fltarr(nlev)

for ilev=0,nlev-1 do begin
   if not keyword_set(psurf) then begin
      ww[*,*,ilev]=-1.*om[*,*,ilev]*gascon*temp[*,*,ilev]/pres[ilev]/ga ;m/s
   endif else begin
      ww[*,*,ilev]=-1.*om[*,*,ilev]*gascon*temp[*,*,ilev]/sigma[ilev]/sp[*,*]/ga ;m/s
   endelse
   dw[ilev]=sqrt(total((ww[*,0:nlon/4.,ilev])^2*cos(lat[*,0:nlon/4.,ilev]*!pi/180.) $
             +(ww[*,nlon*3./4.:nlon-1,ilev])^2*cos(lat[*,nlon*3./4.:nlon-1,ilev]*!pi/180.)) $
              /total(cos(lat[*,0:nlon/4.,ilev]*!pi/180.)+cos(lat[*,nlon*3./4.:nlon-1,ilev]*!pi/180.)))
   nw[ilev]=sqrt(total((ww[*,nlon/4.:nlon*3./4.,ilev])^2*cos(lat[*,nlon/4.:nlon*3./4.,ilev]*!pi/180.)) $
                    /total(cos(lat[*,nlon/4.:nlon*3./4.,ilev]*!pi/180.)))
   gw[ilev]=sqrt(total((ww[*,*,ilev])^2*cos(lat[*,*,ilev]*!pi/180.))/total(cos(lat[*,*,ilev]*!pi/180.)))
endfor

if not keyword_set(noplot) then begin

   if keyword_set(psurf) then begin
      yrange=[max(sp),min(sigma)*min(sp)] 
      if min(sigma)*min(sp) gt p0/10.^(oom) then yrange[1]=p0/10.^(oom)
   endif else begin
      yrange=[p0,min(pres)]
   endelse

   !p.multi=[0,2,2,1,0]
   xmax=max([dw,nw,gw],min=xmin)
   plot,gw,pres,/ylog,yr=yrange,/ystyle,xr=[xmin,xmax],xtit='RMS vertical velocity [m/s]',ytit='Pressure [bar]',$
        xmargin=[12,3]
   oplot,dw,pres,linestyle=1
   oplot,nw,pres,linestyle=2
   al_legend,['Global average','Dayside average','Nightside average'],linestyle=[0,1,2],charsize=1.25,/bottom,/right
   
   xmax=max(ww[nlat/2.,*,*],min=xmin)
   plot,[0,0],[p0,p0],/ylog,yr=yrange,xr=[xmin,xmax],xtit='vertical velocity [m/s]',$
        title='Equatorial longitudes',ystyle=8,ytickname=replicate(' ',8),xmargin=[3,12]
   axis,yaxis=1,/ylog,yr=yrange,/ystyle,/save,ytit='Pressure [bar]'
   loadct,25,/silent
   ccolors=findgen(nlon)*255/(nlon-1)
   if not keyword_set(psurf) then begin
      for i=0,nlon-1,4 do oplot,ww[nlat/2.,i,*],pres,color=ccolors[i]
   endif else begin
      for i=0,nlon-1,4 do oplot,ww[nlat/2.,i,*],sp[nlat/2.,i]*sigma,color=ccolors[i]
   endelse
   colorbar,range=[0.,360.],charsize=1,/vertical,position=[0.54,0.6,0.55,0.9],divisions=4,/right
   loadct,0,/silent
;xyouts,xmin*0.5,p0*0.5,'Longitudes along the equator',charsize=1.25

   xmax=max(ww[*,0,*],min=xmin)
   plot,[0,0],[p0,p0],/ylog,yr=yrange,/ystyle,xr=[xmin,xmax],xtit='vertical velocity [m/s]',ytit='Pressure [bar]',$
        title='Substellar latitudes',xmargin=[12,3]
   loadct,25,/silent
   ccolors=findgen(nlat)*255/(nlat-1)
   if not keyword_set(psurf) then begin
      for i=nlat-1,0,-1 do oplot,ww[i,0.,*],pres,color=ccolors[i]
   endif else begin
      for i=nlat-1,0,-1 do oplot,ww[i,0.,*],sp[i,0.]*sigma,color=ccolors[i]
   endelse
   colorbar,range=[-90.,90.],charsize=1,/vertical,/reverse,position=[0.45,0.1,0.46,0.4]
   loadct,0,/silent
;xyouts,xmin*0.5,p0*0.5,'Latitudes through the substellar point',charsize=1.25

   xmax=max(ww[*,nlon/2.-1,*],min=xmin)
   plot,[0,0],[p0,p0],/ylog,yr=yrange,xr=[xmin,xmax],xtit='vertical velocity [m/s]',$
        title='Antistellar latitudes',ystyle=8,ytickname=replicate(' ',8),xmargin=[3,12]
   axis,yaxis=1,/ylog,yr=yrange,/ystyle,/save,ytit='Pressure [bar]'
   loadct,25,/silent
   ccolors=findgen(nlat)*255/(nlat-1)
   if not keyword_set(psurf) then begin
      for i=nlat-1,0,-1 do oplot,ww[i,nlon/2.-1,*],pres,color=ccolors[i]
   endif else begin
      for i=nlat-1,0,-1 do oplot,ww[i,nlon/2.-1,*],sp[i,nlon/2.-1]*sigma,color=ccolors[i]
   endelse
   colorbar,range=[-90.,90.],charsize=1.,/vertical,/reverse,position=[0.54,0.1,0.55,0.4],/right
   loadct,0,/silent
;xyouts,xmin*0.9,p0*0.5,'Latitudes through the antistellar point',charsize=1.25

   if keyword_set(psfile) then psclose

   loadct,0,/silent
   !p.multi=0
   !p.charsize=1
   !p.thick=1
   !x.thick=1
   !y.thick=1

endif

if keyword_set(file) then begin

   ab=fltarr(7,nlat,nlon,nlev)
   ab[0,*,*,*]=xy[0,*,*,*]
   ab[1,*,*,*]=xy[1,*,*,*]
   ab[2,*,*,*]=xy[2,*,*,*]
   ab[3,*,*,*]=xy[3,*,*,*]
   ab[4,*,*,*]=xy[4,*,*,*]
   ab[5,*,*,*]=ww
   ab[6,*,*,*]=xy[5,*,*,*]
   openw,19,'fort.26.uvw'
   printf,19,nlat,nlon,nlev
   printf,19,ab
   printf,19,fstring
   close,19
   
endif

end


