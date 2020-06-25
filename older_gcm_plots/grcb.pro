pro grcb,psfile=psfile,prcb=prcb,avept=avept

; calculates 3D RCB from Guillot profiles
; and compares cooling rate for 3D vs 1D (average) RCBs

oom=4

tirr=1750.
tint=500.
ga=1.5e3   ;cgs
gascon=3523.*1.e4  ;cgs
akap=0.286
p0=1.e3
gnlev=50

alpha=1.
kth=1.e-2
kv=2.e-3

;use fort.26 file for lat,lon values (and nlat,nlon,nlev)
openr, 17, 'fort.26'
readf,17,nlat,nlon,nlev
xy=fltarr(6,nlat,nlon,nlev)  ;0: lon, 1: lat, 2: level, 3: u, 4: v, 5: T
readf,17,xy
close,17

sigma=get_sigma(oom,gnlev)
pres=p0*sigma
ymax=max(pres,min=ymin)
bvf=fltarr(nlat,nlon+1,gnlev)
cosa=reform(cos(xy[0,*,*,0]*!pi/180.)*cos(xy[1,*,*,0]*!pi/180.))
pcosa=fltarr(nlat,nlon+1)
pcosa[*,0:nlon-1]=cosa
pcosa[*,nlon]=cosa[*,0]

tnight=guillot(pres,1.,tint=tint,tirr=tirr,kth=kth,kv=kv,g=ga,alpha=alpha,/night)

;for ilon=0,nlon-1 do begin
for ilon=0,nlon do begin
   for ilat=0,nlat-1 do begin
      if pcosa[ilat,ilon] le 0 then temp=tnight else begin
         temp=guillot(pres,pcosa[ilat,ilon],tint=tint,tirr=tirr,kth=kth,kv=kv,g=ga,alpha=alpha)
      endelse
      dtdp=deriv(alog(pres),alog(temp))
      bvf[ilat,ilon,*]=(ga^2/gascon/temp)*(akap-dtdp)
   endfor
endfor

xmin=min(abs(bvf),max=xmax)
;plot,[0,0],[0,0],yr=[ymax,ymin],/ystyle,/ylog,xr=[xmin,xmax],$
;     xtit='N!E2!N [s!E-2!N]',ytit='Pressure [bar]',/xlog
;loadct,5,/silent
;for ilon=0,nlon-1 do begin
;   for ilat=0,nlat-1 do begin
;      oplot,bvf[ilat,ilon,*],sigma*p0
;      oplot,-bvf[ilat,ilon,*],sigma*p0,color=100
;   endfor
;endfor
;loadct,0,/silent

prcb=fltarr(nlat,nlon+1)+p0
trcb=fltarr(nlat,nlon+1)
toss=where(bvf le 0)
if toss[0] ne -1 then begin
   convind=array_indices(bvf,toss)
   for i=0,n_elements(toss)-1 do begin
      pconv=sigma[convind[2,i]]*p0
      if pconv lt prcb[convind[0,i],convind[1,i]] then begin
         prcb[convind[0,i],convind[1,i]]=pconv
         if convind[1,i] eq nlon then lonind=0 else lonind=convind[1,i]
         cosa=cos(xy[0,convind[0,i],lonind,0]*!pi/180.)*cos(xy[1,convind[0,i],lonind,0]*!pi/180.)
         if cosa le 0. then begin
            cosa=1.
            night=1
         endif else night=0
         temp=guillot(pconv,cosa,tint=tint,tirr=tirr,kth=kth,kv=kv,g=ga,alpha=alpha,nightside=night)
         trcb[convind[0,i],convind[1,i]]=temp
;         trcb[convind[0,i],convind[1,i]]=xy[5,convind[0,i],lonind,convind[2,i]]
      endif
   endfor
endif
toss=where(trcb eq 0.)
if toss[0] ne -1 then trcb[toss]=max(trcb)

if keyword_set(psfile) then begin
   if (size(psfile))[1] eq 2 then filename='plot.eps' else filename=psfile
   psopen,filename,/color,xsize=25,ysize=20,/enc;,bits_per_pixel=24
   !p.font=0
   !p.charsize=2
   !p.thick=4
   !x.thick=4
   !y.thick=4
endif else begin
   !p.font=-1
   !p.charsize=1.5
   device,true_color=24,decomposed=0,retain=2
   window,0
endelse

;pprcb=fltarr(nlat,nlon+1)
;pprcb[*,0:nlon-1]=prcb
;pprcb[*,nlon]=prcb[*,0]
pprcb=prcb
aprcb=alog10(pprcb)
flon=fltarr(nlat,nlon+1)
flon[*,0:nlon-1]=reform(xy[0,*,*,0])
flon[*,nlon]=reform(xy[0,*,0,0])+360.
flat=fltarr(nlat,nlon+1)
flat[*,0:nlon-1]=reform(xy[1,*,*,0])
flat[*,nlon]=reform(xy[1,*,0,0])
pmin=min(prcb,max=pmax)
if pmin ne pmax then begin
   ctload,10,/silent,/reverse
   MAP_SET, /MILLER_CYLINDRICAL,0,0,/ISOTROPIC,color=255;,title='P_rcb'
   cbottom=0
   nlevels=60
   step=(alog10(pmax)-alog10(pmin))/(nlevels-1)
   mylevels=indgen(nlevels)*step+alog10(pmin)
   ccolors=findgen(nlevels)*(255.-cbottom)/(nlevels-1)+cbottom

   contour, aprcb, flon,flat,/overplot ,/cell_fill,/closed,$
            levels=mylevels
               ;NLEVELS=45,min_value=alog10(pmin), max_value=alog10(pmax)
   ndiv=4
   divstep=(alog10(pmax)-alog10(pmin))/ndiv
   plevels=pmin*10.^(divstep/4.)*10.^(findgen(ndiv*2)/4.)
   pdiv=pmin*10.^(findgen(ndiv+1)*divstep)
;   map_grid, /label,charsize=1  ,color=0
   contour, pprcb, flon, flat,/over,/follow,levels=plevels,c_charsize=1.5,color=255
   colorbar,  position=[0.1,0.07,0.90,0.10], range=alog10([pmin,pmax]),$
              charsize=2,divisions=ndiv,color=255,ticknames=string(pdiv);,format='(f4.0)'
   loadct,0,/silent
endif

if keyword_set(psfile) then begin
   psclose
   !p.font=-1
   !p.charsize=1
   !p.thick=1
   !x.thick=1
   !y.thick=1
   spawn,'gs -r300 -dEPSCrop -dTextAlphaBits=4 -sDEVICE=png16m -sOutputFile=plot.png -dBATCH -dNOPAUSE plot.eps'
   spawn,'convert plot.png eps3:plot.eps'
endif

if keyword_set(psfile) then begin
   if (size(psfile))[1] eq 2 then filename='plot2.eps' else filename=psfile
   psopen,filename,/color,xsize=30,ysize=20,/enc;,bits_per_pixel=24
   !p.font=0
   !p.charsize=2
   !p.thick=4
   !x.thick=4
   !y.thick=4
endif else begin
   !p.font=-1
   !p.charsize=2
   device,true_color=24,decomposed=0,retain=2
   window,1
endelse
plot,tnight,pres,/ylog,yr=[max(prcb)*1.1,min(pres)],/ystyle,xr=[min(tnight),max(trcb,/nan)],/xstyle,$
     xtit='Temperature [K]',ytit='Pressure [bar]';,/xlog
nlines=n_elements(uniq(prcb(sort(prcb))))
mu=findgen(nlines)/(nlines-1)
for i=0,nlines-1 do begin
   temp=guillot(pres,mu[i],tint=tint,tirr=tirr,kth=kth,kv=kv,g=ga,alpha=alpha)
   oplot,temp,pres
endfor
oplot,trcb,prcb,psym=4
if keyword_set(avept) then begin
   temp=guillot(pres,1.,f=0.25,tint=tint,tirr=tirr,kth=kth,kv=kv,g=ga,alpha=alpha)
   loadct,5,/silent
   oplot,temp,pres,color=100,thick=8
   dtdp=deriv(alog(pres),alog(temp))
   avebvf=(ga^2/gascon/temp)*(akap-dtdp)
   zero=where(avebvf le 0.)
   if zero[0] ne -1 then begin
      tave=temp[zero[0]]
      pave=pres[zero[0]]
      oplot,[tave,tave],[pave,pave],psym=6,color=100,symsize=1.5,thick=6
      print,tave,pave
   endif
   loadct,0,/silent
endif
if keyword_set(psfile) then begin
   psclose
   !p.font=-1
   !p.charsize=1
   !p.thick=1
   !x.thick=1
   !y.thick=1
endif

Lconst=16.*5.67e-8*6.67e-11*0.268/3.
LM=0.
for ilat=0,nlat-1 do begin
   coslat=cos(xy[1,ilat,0,0]*!pi/180.)
   for ilon=0,nlon-1 do begin
      LM=LM + Lconst*coslat*(2.*!pi/nlon)*(!pi/nlat) $
         *trcb[ilat,ilon]^4/(prcb[ilat,ilon]*1.e5) $
         /(kth*(prcb[ilat,ilon]/1.)^alpha*1.e-3)
   endfor
endfor
print,'Total L/M [W/kg]: ',LM

if keyword_set(avept) then begin
   print,'uniform RCB would give [W/kg]:',Lconst*4.*!pi*tave^4/(kth*pave^alpha*1.e-3)/(pave*1.e5)
endif

end
