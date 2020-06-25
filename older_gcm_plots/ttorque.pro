pro ttorque,p0,lonss=lonss
  
; purpose: to calculate the torque applied by the star on the atmosphere
; uses surfp.pro to read from fort.50 and plot the surface pressure field

  surfp,p0=p0,sp=sp,lon=lon,lat=lat
  colat=(90.-lat)*!pi/180.
  lon*=!pi/180.
  nlat=(size(sp))[1]
  nlon=(size(sp))[2]
  dlat=median(colat[1:nlat-1,0]-colat[0:nlat-2,0])
  dlon=mean(lon[0,1:nlon-1]-lon[0,0:nlon-2])
  
  y22=spher_harm(colat,lon,2,2)

  p22=total(conj(y22)*sp*sin(colat))*dlat*dlon ; units: bar

  if not keyword_set(llons) then lonss=0.
  ta=imaginary(p22*exp(2.*complex(0,1)*lonss)) ; units: bar

  
end
