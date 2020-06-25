pro plottmag,p0,oom,bfield,tmagmin=tmagmin,psfile=psfile,xrange=xrange,yrange=yrange

  gascon=3523.

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

if keyword_set(psfile) then begin
   if (size(psfile))[1] eq 2 then filename='plot.eps' else filename=psfile
   if keyword_set(zave) then xsz=40 else xsz=25
   psopen,filename,/enc,/color,xsize=xsz,ysize=20 ;,bits_per_pixel=24
   !p.font=0
   !x.thick=6
   !y.thick=6
   !p.thick=3
   !p.charsize=2
endif

   lon=reform(xy[0,*,*,*]) & lat=reform(xy[1,*,*,*]) 
   u=reform(xy[3,*,*,*]) & v=reform(xy[4,*,*,*]) & temp=reform(xy[5,*,*,*])
   sp=reform(ab[2,*,*])
   sp=(sp+1.)*p0

   sigma=get_sigma(oom,nlev)

   rho=fltarr(nlat,nlon,nlev)
   tau=fltarr(nlat,nlon,nlev)

   for ilev=0,nlev-1 do begin
      for ilat=0,nlat-1 do begin
         rho[ilat,*,ilev]=reform(sp[ilat,*]*sigma[ilev]/gascon/temp[ilat,*,ilev]*1.e2)                  ; g/cm^3
         tdrag=tmag(reform(rho[ilat,*,ilev]),reform(temp[ilat,*,ilev]),bfield)/abs(sin(lat[ilat,0,0]))  ; s
         tau[ilat,*,ilev]=tdrag
      endfor
   endfor

   if keyword_set(yrange) then begin
      tmin=yrange[0]
      tmax=yrange[1]
   endif else tmin=min(temp,max=tmax)
   if keyword_set(xrange) then begin
      rmin=xrange[0]
      rmax=xrange[1]
   endif else rmin=min(rho,max=rmax)
   ltau=alog10(tau)
   pmin=min(ltau,max=pmax,/nan)
   print,pmin,pmax

   plot,[0,0],[0,0],xr=[rmin,rmax],/xlog,yr=[tmin,tmax],xtit='Density [cgs]',ytit='Temperature [K]',$
        /xstyle,/ystyle
   ctload,13,/silent,/reverse
   for ilev=0,nlev-1 do begin
      for ilat=0,nlat-1 do begin
         for ilon=0,nlon-1 do begin
            ;(to max out color for Inf values)
            if (ltau[ilat,ilon,ilev] gt pmax) then pcolor=255 else pcolor=255.*(ltau[ilat,ilon,ilev]-pmin)/(pmax-pmin)
            oplot,[rho[ilat,ilon,ilev],rho[ilat,ilon,ilev]],[temp[ilat,ilon,ilev],temp[ilat,ilon,ilev]],$
                  color=pcolor,psym=1
         endfor
      endfor
   endfor
   colorbar,range=[pmin,pmax],format='(f5.2)',charsize=1.5,color=255,position=[0.25,0.85,0.85,0.9]

   loadct,0,/silent

   if keyword_set(tmagmin) then begin
      out=where(tau lt tmagmin)
      if out[0] ne -1 then begin
         for i=0,n_elements(out)-1 do oplot,[rho[out[i]],rho[out[i]]],[temp[out[i]],temp[out[i]]],psym=4
         ctload,13,/silent,/reverse
         for i=0,n_elements(out)-1 do oplot,[rho[out[i]],rho[out[i]]],[temp[out[i]],temp[out[i]]],$
                  color=255.*(ltau[out[i]]-pmin)/(pmax-pmin),psym=3,symsize=6
         loadct,0,/silent
      endif
   endif

if keyword_set(psfile) then begin
   !p.font=-1
   !x.thick=1
   !y.thick=1
   !p.thick=1
   !p.charsize=1
   psclose
;   spawn,'gs -r300 -dEPSCrop -dTextAlphaBits=4 -sDEVICE=png16m -sOutputFile=plot.png -dBATCH -dNOPAUSE plot.eps'
;   spawn,'convert plot.png eps3:plot.eps'
endif
   
end
