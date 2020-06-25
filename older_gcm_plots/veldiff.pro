pro veldiff,alpha,theta,lat,lon,r=r,kappa=kappa,nom=nom,p0=p0

;reads fort.26, makes simplified velocity profile (horizontal
;component vs pressure for one 2D location)

;can input multiple alpha,theta

;alpha is angle from normal (0: ray goes straight up, 90: horizontal)
;theta is angle of propogation (0:E, 90:N, 180:W, 270:S)
;lat is latitude of profile
;lon in longitude of profile
;ALL ANGLES INPUT IN DEGREES

;r is specific gas constant (R)
;kappa is R/c_p=1+1/gamma

if not keyword_set(r) then r=4593.          ;value for C&S05
if not keyword_set(kappa) then kappa=0.321  ;value for C&S05
if not keyword_set(nom) then nom=5.34363
if not keyword_set(p0) then p0=220.16

openr,17,'fort.26'
  readf,17,nlat,nlon,nlev
  xy=fltarr(6,nlat,nlon,nlev)
  readf,17,xy
close,17

if nom ne 1 then begin
   vert=fltarr(nlev)
   stp=-1.*nom/nlev
   vert[nlev-1]=10.^(stp/2.)
   for i=nlev-2,0,-1 do vert[i]=vert[i+1]*10.^(stp)
   vert*=p0
endif else begin
   vert=findgen(nlev)+1
   vert*=p0
endelse

ntheta=n_elements(theta)
nalpha=n_elements(alpha)
winds=fltarr(ntheta,nalpha,nlev)

;find indices for lat,lon location
print,'location chosen:',lat,lon
toss=min(abs(xy[0,0,*,0]-lon),loni)
toss=min(abs(xy[1,*,0,0]-lat),lati)
print,'location used:',xy[1,lati,0,0],xy[0,0,loni,0]   

;get wind and temp vertical profiles
u=reform(xy[3,lati,loni,*]) & v=reform(xy[4,lati,loni,*]) & temp=reform(xy[5,lati,loni,*])
;calculate local sound speeds
lss=sqrt(r*temp/(1.-kappa))

for j=0,ntheta-1 do begin

   ;use theta to get appropriate u,v components
   wind=u*cos(theta[j]*!pi/180.)+v*sin(theta[j]*!pi/180.)
   print,'theta:',theta[j]
   
   for k=0,nalpha-1 do begin

      ;use alpha to get horizontal component
      wind*=sin(alpha[k]*!pi/180.)
      print,'   alpha:',alpha[k]

      ;put in units of local sound speed
      wind/=lss

      winds[j,k,*]=wind

   endfor
endfor

if nom ne 1 then $
     plot,winds[0,0,*],vert,yr=[max(vert),min(vert)],/ystyle,xr=[min(winds),max(winds)],$
       ytit='pressure [bar]',xtit='projected wind speed [v/c!Ds!N]',/ylog $
else plot,winds[0,0,*],vert,yr=[max(vert),min(vert)],/ystyle,xr=[min(winds),max(winds)],$
       ytit='pressure [bar]',xtit='projected wind speed [v/c!Ds!N]'

i=0
for j=0,ntheta-1 do begin
   for k=0,nalpha-1 do begin
      oplot,winds[j,k,*],vert,linestyle=i
      i+=1
   endfor
endfor

end
