pro get_temp,level=level,long=long

; returns temp at all latitudes for given sigma level and given longitude

  OpenR, 17, 'fort.26'
    ReadF, 17, nlat, nlon, nlev
    xy=fltarr(6,nlat,nlon,nlev)
    READF, 17, xy
  Close, 17
  lon=reform(xy[0,*,long,level]) & lat=reform(xy[1,*,long,level])
  u=reform(xy[3,*,*,level]) & v=reform(xy[4,*,long,level])
  temp=reform(xy[5,*,long,level])
  n= N_ELEMENTS(temp)
  print,temp
  OPENW, 18,'PhotoTemp.txt'
  for j=0,n-1 DO BEGIN
    PRINTF, 18, temp[j]
  endfor
  CLOSE,18
  
; really only care about the temperature one, left the other ones in
; in case needed in the future

end
