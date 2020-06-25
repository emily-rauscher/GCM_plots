pro terminator_profiles,oom=oom,p0=p0,sigma=sigma,tprof=tprof

;reads from fort.26 in same directory
;can optionally output sigma and tprof
if not keyword_set(p0) then p0=1.

;east (lon=nlon/4.) terminator in purple
;west (lon=nlon*3./4) terminator in yellow

device,true_color=24,decomposed=0,retain=2
loadct,5

;grab fort.26 info
openr, 17, 'fort.26'
readf,17,nlat,nlon,nlev
xy=fltarr(6,nlat,nlon,nlev)  ;0: lon, 1: lat, 2: level, 3: u, 4: v, 5: T
readf,17,xy
close,17

tprof=fltarr(2,nlat,nlev)  ;0:east, 1:west

for i=0,nlat-1 do begin
   tprof[0,i,*]=xy[5,i,nlon/4.,*]
   tprof[1,i,*]=xy[5,i,nlon*3./4.,*]
endfor

;calculate sigma
if keyword_set(oom) then begin
  stp=-1.*oom/nlev
  sigma=fltarr(nlev)
  sigma[nlev-1]=10.^(stp/2.)
  for i=nlev-2,0,-1 do sigma[i]=sigma[i+1]*10.^(stp)
endif else begin
  sigma=findgen(nlev)/nlev+0.5/nlev
endelse

press=sigma*p0

xn=min(tprof)
xx=max(tprof)
xtit='Temperature [K]'
ytit='Pressure [bar]'

if keyword_set(oom) then $
   plot,tprof[0,0,*],press,xr=[xn,xx],yr=[max(press),min(press)],/ylog,xtitle=xtit,ytitle=ytit,/ystyle $
        else plot,tprof[0,0,*],press,xr=[xn,xx],yr=[1.,0.]*press,xtitle=xtit,ytitle=ytit

for i=0,nlat-1 do oplot,tprof[1,i,*],press,color=200
for i=0,nlat-1 do oplot,tprof[0,i,*],press,color=75

loadct,0

end
