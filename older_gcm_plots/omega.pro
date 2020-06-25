pro omega

; SEE VERTVEL.PRO FOR A MORE DETAILED PROGRAM, THIS ONE IS NO LONGER USED

OOM = 5.34363
p0=220.612
RADEA = 9.44e7
GA = 9.42
GASCON = 4593.

  OpenR, 17, 'fort.26'
  ReadF, 17, nlat,nlon,nlev	  
	  xy=fltarr(6,nlat,nlon,nlev)   
  ReadF, 17, xy
  Close, 17
  lon=reform(xy[0,*,*,*]) & lat=reform(xy[1,*,*,*]) 
  u=reform(xy[3,*,*,*]) & v=reform(xy[4,*,*,*]) 
  temp=reform(xy[5,*,*,*])

sigma=get_sigma(oom,nlev)

;integrate continuity for omega (=0 at p=0)
; careful: increasing the lat index is going south, not north
om=fltarr(nlat,nlon,nlev)
for ilev=0,nlev-1 do begin
   if ilev gt 0 then dP=p0*(sigma[ilev]-sigma[ilev-1]) else dP=p0*sigma[0]
   for ilon=0,nlon-1 do begin
      ilo=ilon-1
      ilonn=ilon+1
      if ilon eq 0 then ilo=nlon-1
      if ilon eq nlon-1 then ilonn=0
      ;first do poles (should be ilat -> ilat,ilon+nlon/2), with flow now negative)
      ;careful with +nlon/2
      if ilon+nlon/2 gt nlon-1 then iilon=ilon-nlon/2 else iilon=ilon+nlon/2
      if ilev gt 0 then begin
         om[0,ilon,ilev]=(-0.5/radea)$
                          *[(u[0,ilonn,ilev]-u[0,ilo,ilev])/(cos(lat[0,ilon,ilev]*!pi/180.)*(2.*!pi/nlon))$
                            +(-1.*v[0,iilon,ilev]-v[1,ilon,ilev])/(!pi/nlat)] *dP + om[0,ilon,ilev-1]
         om[nlat-1,ilon,ilev]=(-0.5/radea)$
                          *[(u[nlat-1,ilonn,ilev]-u[nlat-1,ilo,ilev])$
                            /(cos(lat[nlat-1,ilon,ilev]*!pi/180.)*(2.*!pi/nlon))$
                            +(v[nlat-2,ilon,ilev]+v[nlat-1,iilon,ilev])/(!pi/nlat)] *dP + om[nlat-1,ilon,ilev-1]
      endif else begin
         om[0,ilon,ilev]=(-0.5/radea) $
                          *[(u[0,ilonn,ilev]-u[0,ilo,ilev])/(cos(lat[0,ilon,ilev]*!pi/180.)*(2.*!pi/nlon))$
                            +(-1.*v[0,iilon,ilev]-v[1,ilon,ilev])/(!pi/nlat)] *dP 
         om[nlat-1,ilon,ilev]=(-0.5/radea) $
                          *[(u[nlat-1,ilonn,ilev]-u[nlat-1,ilo,ilev])$
                            /(cos(lat[nlat-1,ilon,ilev]*!pi/180.)*(2.*!pi/nlon))$
                            +(v[nlat-2,ilon,ilev]+v[nlat-1,iilon,ilev])/(!pi/nlat)] *dP 
      endelse
      for ilat=1,nlat-2 do begin
         if ilev gt 0 then om[ilat,ilon,ilev]=(-0.5/radea)$
                          *[(u[ilat,ilonn,ilev]-u[ilat,ilo,ilev])/(cos(lat[ilat,ilon,ilev]*!pi/180.)*(2.*!pi/nlon))$
                            +(v[ilat-1,ilon,ilev]-v[ilat+1,ilon,ilev])/(!pi/nlat)] *dP + om[ilat,ilon,ilev-1] $
         else om[ilat,ilon,ilev]=(-0.5/radea) $
                          *[(u[ilat,ilonn,ilev]-u[ilat,ilo,ilev])/(cos(lat[ilat,ilon,ilev]*!pi/180.)*(2.*!pi/nlon))$
                            +(v[ilat-1,ilon,ilev]-v[ilat+1,ilon,ilev])/(!pi/nlat)] *dP 
      endfor
   endfor
endfor
om*=1.e5  ;to convert units to: (kg/m/s^2)/s

dzdt=om
for i=0,nlev-1 do begin
   dzdt[*,*,i]=-1.*dzdt[*,*,i]*(gascon/ga)*temp[*,*,i]/(p0*sigma[i]*1.e5) ;should be in m/s
endfor

stop

end
