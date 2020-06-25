pro make_radeq

;want 3D set of T,pressure,altitude 
;        as fort.26 equivalent: lon,lat,lev,u,v,temp
;        and vertical.txt:      sigma,press,altitude
;(assume no dynamics, so u=v=0 and no variation in p_surf)
;

;SET THESE BEFORE RUNNING:
oom=7.
psurf=100.
R=4593.
g=9.42
TGR=1200.

;steals lat,lon,lev info from fort.26 in same directory
openr,17,'fort.26'
readf,17,nlat,nlon,nlev
print,nlat,nlon,nlev
xy=fltarr(6,nlat,nlon,nlev)
readf,17,xy
close,17

;overwrite and set u=v=0
xy[3,*,*,*]=replicate(0.,nlat,nlon,nlev)
xy[4,*,*,*]=replicate(0.,nlat,nlon,nlev)

;Calculate T at each point:
;set values from fort.7
resttt=[670.,678.,687.,698.,711.,725.,742.,760.,782.,806.,$
833.,864.,899.,938.,983.,1036.,1094.,1159.,1240.,1345.,$
1428.,1488.,1531.,1546.,1524.,1453.,1401.,1366.,1341.,1321.,$
1307.,1295.,1287.,1280.,1275.,1271.,1269.,1269.,1269.,1270.,$
1270.,1272.,1274.,1277.,1283.]
redtep=[1085.,1088.,1094.,1101.,1110.,1120.,1131.,1142.,1152.,1159.,$
1162.,1157.,1140.,1109.,1059.,972.,801.,606.,480.,316.,$
194.,121.,97.,147.,283.,480.,623.,720.,788.,839.,$
875.,899.,915.,926.,933.,936.,938.,938.,938.,938.,$
937.,936.,934.,931.,926.]
if (n_elements(resttt) ne nlev) or (n_elements(redtep) ne nlev) then begin
   print,'ERROR: n_elements(resttt or redtep) ne nlev'
   stop
endif


;set T profiles
for i=0,nlev-1 do begin
   for j=0,nlon-1 do begin
      for k=0,nlat-1 do begin
         alon=xy[0,k,j,i]
         alat=xy[1,k,j,i]
         if (alon ge 90.) and (alon le 270.) then begin
            xy[5,k,j,i]=resttt[i]
         endif else begin
            xy[5,k,j,i]=( (resttt[i])^4 + ((resttt[i]+redtep[i])^4 - (resttt[i])^4) $
                                          *cos(alat*!pi/180.)*cos(alon*!pi/180.) )^0.25
         endelse
      endfor
   endfor
endfor

sigma=get_sigma(oom,nlev)
z=fltarr(nlat,nlon,nlev)

;Calculate z at each point:
;set altitude of first level (up from base=p0=[sigma=1], where T=TGR
;and z=0; no variation in psurf to account for here)
z[*,*,nlev-1]=(R/g)*0.5*(xy[5,*,*,nlev-1]+TGR)*alog(1./sigma[nlev-1])

;integrate hydrostatic to solve for higher levels
for i=nlev-2,0,-1 do begin
   z[*,*,i]=z[*,*,i+1]+(R/g)*0.5*(xy[5,*,*,i]+xy[5,*,*,i+1])*alog(sigma[i+1]/sigma[i])
endfor

;print results
openw,18,'fort.26.radeq'
openw,19,'vertical.txt.radeq'
printf,18,nlat,nlon,nlev
printf,19,nlat,nlon,nlev

for i=0,nlev-1 do begin
   for j=0,nlon-1 do begin
      for k=0,nlat-1 do begin
         printf,18,xy[*,k,j,i]
         printf,19,sigma[i],sigma[i]*psurf,z[k,j,i]
      endfor
   endfor
endfor

printf,19,'R,g,TGR,p0,OOM:',R,g,TGR,psurf,oom

close,18
close,19

;stop

end
