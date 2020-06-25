pro torque,tdrag,dp,ga,radea

; calculates north/south, east/west torque at bottom layer
;
; tdrag: drag timescale (in s) at bottom level
; dp: delta P of lowest level in bar, = 0.5*(psurf-pres[nlev-2])
; ga: grav acc in m/s^2
; radea: planet radius in m

latoff=0
lonoff=0

; Read File with 3D fields

  OpenR, 17, 'fort.26'
  ReadF, 17, nlat,nlon,nlev	  
;	  print,nlat,nlon,nlev 
	  xy=fltarr(6,nlat,nlon,nlev)   
  ReadF, 17, xy
  Close, 17

   nl=nlev-1
    
   lon=reform(xy[0,*,*,nl]) & lat=reform(xy[1,*,*,nl])
   u=reform(xy[3,*,*,nl]) & v=reform(xy[4,*,*,nl])
   temp=reform(xy[5,*,*,nl])

   dm=fltarr(nlat,nlon)
   for ilat=0,nlat-1 do dm[ilat,*]=cos(lat[ilat,*]*!pi/180.)
   dm*=((radea^2/ga)*dp*1.e5*(!pi/nlat)*(2.*!pi/nlon))

   nst=(radea/tdrag)*dm*v
   ewt=(radea/tdrag)*dm*u

   print,'Rate of N momentum given to interior: ',total(nst),' N m'
   print,'Rate of E momentum given to interior: ',total(ewt),' N m'

;stop

end
