pro lcool,prcb,kth,alpha,p0,oom

;grab fort.26 info
openr, 17, 'fort.26'
readf,17,nlat,nlon,nlev
xy=fltarr(6,nlat,nlon,nlev)  ;0: lon, 1: lat, 2: level, 3: u, 4: v, 5: T
readf,17,xy
close,17

;grab psurf
openr, 18, 'fort.50'
readf, 18, nlat,nlon
ab=fltarr(3,nlat,nlon)
readf, 18, ab
close, 18
sp=reform(ab[2,*,*])
sp = (sp+1.)*p0

sigma=get_sigma(oom,nlev)

Lconst=16.*5.67e-8*6.67e-11*0.268/3.
Larr=fltarr(nlat,nlon)  ; W/kg

for ilat=0,nlat-1 do begin
   coslat=cos(xy[1,ilat,0,0]*!pi/180.)
   for ilon=0,nlon-1 do begin
      if prcb[ilat,ilon] eq sp[ilat,ilon] then begin
         Tbot=xy[5,ilat,ilon,pind]
         Pbot=sp[ilat,ilon]*sigma(nlev-1)
         Larr[ilat,ilon]=Lconst*coslat*(2.*!pi/nlon)*(!pi/nlat) $
                         *(Tbot+(Tbot/Pbot)*(prcb[ilat,ilon]-Pbot)*0.286)^4 $
                         /(prcb[ilat,ilon]*1.e5) $            ; convert to Pa
                         /(kth*(prcb[ilat,ilon]/1.)^alpha*1.e-3) ; convert to m^3/kg
      endif else begin
         pind=where(sp[ilat,ilon]*sigma eq prcb[ilat,ilon])
         if ((pind[0] eq -1) or (n_elements(pind) gt 1)) then stop
         Larr[ilat,ilon]=Lconst*coslat*(2.*!pi/nlon)*(!pi/nlat) $
                         *xy[5,ilat,ilon,pind]^4 $
                         /(prcb[ilat,ilon]*1.e5) $            ; convert to Pa
                         /(kth*(prcb[ilat,ilon]/1.)^alpha*1.e-3) ; convert to m^3/kg
      endelse
   endfor
endfor

print,'Total L/M [W/kg]: ',total(Larr)

end
