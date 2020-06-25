pro rcb,alpha,lev,oom=oom,p0=p0,psurf=psurf,psfile=psfile,trend=trend,xlog=xlog,nskip=nskip,prange=prange,trange=trange

if not keyword_set(psfile) then begin
   ;Screen output
   device,true_color=24,decomposed=0,retain=2
   !p.font=-1
   window,0
endif else begin
   ;Postscript output
   if (size(psfile))[1] eq 2 then filename='plot.eps' else filename=psfile
   psopen,filename,/enc,x=20,y=15,/color
   !p.font=0
   !x.thick=5
   !y.thick=5
endelse
!x.style=1
!p.charsize=1.5

;grab fort.26 info
openr, 17, 'fort.26'
readf,17,nlat,nlon,nlev
xy=fltarr(6,nlat,nlon,nlev)  ;0: lon, 1: lat, 2: level, 3: u, 4: v, 5: T
readf,17,xy
close,17

   akap=2./7

;calculate sigma
if keyword_set(oom) then begin
  sigma=get_sigma(oom,nlev)
endif else begin
  sigma=findgen(nlev)/nlev+0.5/nlev
endelse

if keyword_set(psurf) then begin
   if not keyword_set(p0) then begin
      print,'MUST SET P0 IF USING PSURF!'
      abort
   endif
  openr, 18, 'fort.50'
  readf, 18, nlat,nlon
  ab=fltarr(3,nlat,nlon)
  readf, 18, ab
  close, 18
  sp=reform(ab[2,*,*])
  sp = (sp+1.)*p0
endif

if not keyword_set(trend) then begin

   prcb=fltarr(nlat,nlon+1)
   trcb=prcb
   lamb=(1.+alpha)/4.
   prcb[*,0:nlon-1]=sigma[lev]*sp*(akap/(lamb-akap))^(1./(1.+alpha))
   trcb[*,0:nlon-1]=reform(xy[5,*,*,lev])*(lamb/(lamb-akap))^0.25

   prcb[*,nlon]=prcb[*,0]
   trcb[*,nlon]=trcb[*,0]

   flon=fltarr(nlat,nlon+1)
   flon[*,0:nlon-1]=reform(xy[0,*,*,0])
   flon[*,nlon]=reform(xy[0,*,0,0])+360.
   flat=fltarr(nlat,nlon+1)
   flat[*,0:nlon-1]=reform(xy[1,*,*,0])
   flat[*,nlon]=reform(xy[1,*,0,0])

   pmin=min(prcb,max=pmax)
   tmin=min(trcb,max=tmax)
   print,'pmin, pmax:',pmin,pmax
   print,'tmin, tmax:',tmin,tmax

   if pmin ne pmax then begin
      window,0
      ctload,10,/silent,/reverse
      MAP_SET, /MILLER_CYLINDRICAL,0,0,/ISOTROPIC 
      contour, prcb, flon,flat,/overplot ,/cell_fill,/closed ,NLEVELS=45,$
               min_value=pmin, max_value=pmax
      map_grid, /label,charsize=1  ,color=0
      colorbar,  position=[0.1,0.07,0.90,0.10], range=[pmin,pmax], format='(f5.1)', $ ;, format='(e8.2)', $
                 charsize=1.5,divisions=5,color=255
      loadct,0,/silent
      window,2
      loadct,4,/silent
      MAP_SET, /MILLER_CYLINDRICAL,0,0,/ISOTROPIC 
      contour, trcb, flon,flat,/overplot ,/cell_fill,/closed ,NLEVELS=45,$
               min_value=tmin, max_value=tmax
      map_grid, /label,charsize=1  ,color=0
      colorbar,  position=[0.1,0.07,0.90,0.10], range=[tmin,tmax], format='(f6.0)', $
                 charsize=1.5,divisions=5,color=255
      loadct,0,/silent
   endif

endif else begin

   if not keyword_set(trange) then begin
      xmin=min(xy[5,*,*,*],max=xmax)
      xmax=xmin*(1./sigma[0])^akap
   endif else begin
      xmin=trange[0]
      xmax=trange[1]
   endelse
   if not keyword_set(prange) then begin
      pmax=max(sp,min=pmin) 
      pmin*=sigma[0]
   endif else begin
      pmax=prange[0]
      pmin=prange[1]
   endelse
   plot,[1,1],[1,1],xtit='Potential temperature [K]',xr=[xmin,xmax],ytit='Pressure [bar]',yr=[pmax,pmin],/ylog,/ystyle,xlog=xlog
   loadct,13,/silent  
;   loadct,25,/silent  
;  red (day side, eqtr): 255, 255
;  yellow (d.s., poles) : 200, 200
;  light blue/green (n.s, poles): 85, 120?
;  dark blue (night side, eqtr): 50, 50
   poles=where(abs(reform(xy[1,*,0,0]) ) gt 30.,complement=eqtr)
   ns=where(((reform(xy[0,0,*,0]) ge 90.) and (reform(xy[0,0,*,0]) le 270.)),complement=ds)
   if not keyword_set(nskip) then nskip=1
   for ilon=0,nlon-1,nskip do begin
      toss2=where(ds eq ilon)
      for ilat=0,nlat-1,nskip do begin
         toss=where(poles eq ilat)
         if toss2[0] eq -1 then begin ; ns
            if toss[0] eq -1 then pcolor=50 $  ; ns, eqtr
                             else pcolor=120    ; ns, poles
         endif else begin
            if toss[0] eq -1 then pcolor=255 $  ; ds, eqtr
                             else pcolor=200    ; ds, poles
         endelse
         oplot,xy[5,ilat,ilon,*]*(p0/(sp[ilat,ilon]*sigma))^akap,sp[ilat,ilon]*sigma,color=pcolor
;         oplot,xy[5,ilat,ilon,*],sp[ilat,ilon]*sigma,color=pcolor
      endfor
   endfor
   al_legend,['Day side, equator','Day side, poles','Night side, equator','Night side, poles'],color=[255,200,50,120],linestyle=[0,0,0,0],$
             /right,/bottom,charsize=1,thick=3
   loadct,0,/silent
stop

endelse

if keyword_set(psfile) then begin
   !p.font=-1
   !x.thick=1
   !y.thick=1
   psclose
endif

end

;prcb=fltarr(nlat,nlon+1)+1.
;if keyword_set(p0) then prcb*=p0
;for ilat=0,nlat-1 do begin
;   for ilon=0,nlon-1 do begin
;      if keyword_set(psurf) then begin
;;        ALOG gives better deriv answer:
;         dtdp=deriv(alog(sigma*sp[ilat,ilon]),alog(reform(xy[5,ilat,ilon,*])))
;      endif else begin
;         dtdp=deriv(alog(sigma),alog(reform(xy[5,ilat,ilon,*])))
;      endelse
;      test=where(dtdp gt 0.286)
;      if test[0] ne -1 then begin
;         print,test
;         if test[0] eq 0 then begin
;            if n_elements(test) gt 1 then begin
;               if test[1] eq 1 then begin
;                  prcb[ilat,ilon]=sigma[test[0]]
;                  if keyword_set(psurf) then prcb[ilat,ilon]*=sp[ilat,ilon] else $
;                     if keyword_set(p0) then prcb[ilat,ilon]*=p0
;               endif else begin
;                  prcb[ilat,ilon]=sigma[test[1]]
;                  if keyword_set(psurf) then prcb[ilat,ilon]*=sp[ilat,ilon] else $
;                     if keyword_set(p0) then prcb[ilat,ilon]*=p0
;               endelse
;            endif
;         endif else begin
;            prcb[ilat,ilon]=sigma[test[0]]
;            if keyword_set(psurf) then prcb[ilat,ilon]*=sp[ilat,ilon] else $
;               if keyword_set(p0) then prcb[ilat,ilon]*=p0
;         endelse
;      endif
;   endfor
;endfor

;aprcb=0.
;for i=0,nlat-1 do aprcb=aprcb+total(prcb[i,0:nlon-1]*cos(xy[1,i,0,0]*!pi/180.))
;norm=total(cos(xy[1,*,*,0]*!pi/180.))
;print,norm
;aprcb/=norm
;print,'Area-weighted average RCB pressure:',aprcb
