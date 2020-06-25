pro test,lat,lon

  nlat=(size(lat))[1]
  nlon=(size(lat))[2]
  obliq=90.*!pi/180.
  
  for itime=0,100 do begin

     sslat=obliq*sin(2.*!pi*itime/100.)
     print,sslat

     H=acos(-tan(lat)*tan(sslat))
     if sslat gt 0 then begin
        pday=where(lat gt !pi/2.-sslat)
        pnight=where(lat lt -!pi/2.+sslat)
        H[pday]=!pi
        H[pnight]=0.
     endif else begin
        pday=where(lat lt -!pi/2.-sslat)
        pnight=where(lat gt !pi/2.+sslat)
        H[pday]=!pi
        H[pnight]=0.
     endelse

     ftoa=(sin(lat)*sin(sslat)*H+cos(lat)*cos(sslat)*sin(H))
     print,total(ftoa*cos(lat))*2.*!pi^2/nlat/nlon

     contour,ftoa,lon,lat,/cell_fill,/follow,levels=findgen(100)/62.,$
             xr=[0,2.*!pi],yr=[-1,1]*!pi/2,/xstyle,/ystyle

     wait,0.5
     
  endfor
  
end

