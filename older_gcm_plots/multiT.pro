pro multiT

;device,true_color=24, decomposed=0,retain=2

openr, 17, 'fort.26'
readf,17,nlat,nlon,nlev
xy=fltarr(6,nlat,nlon,nlev)  ;0: lon, 1: lat, 2: level, 3: u, 4: v, 5: T
readf,17,xy
close,17
print,nlat,nlon

oom=5.34363
stp=-1.*oom/nlev
sigma=fltarr(nlev)
sigma[nlev-1]=10.^(stp/2.)
for i=nlev-2,0,-1 do sigma[i]=sigma[i+1]*10.^(stp)
pres=sigma*220.16

for k=0,nlat-1 do begin

    name=strcompress('tprofs_lat'+string(xy[1,k,0,0]),/remov)
    psopen,name,/enc,/color,/landscape
    !p.font=0

    xmin=min(xy[5,k,*,*])
    xmax=max(xy[5,k,*,*])
    plot,xy[5,k,0,*],pres,yr=[max(pres),min(pres)],xr=[xmin,xmax],$
         /ystyle,/xstyle,/ylog,xthick=7,ythick=7,thick=7
    xyouts,xmin+100.,0.1,strcompress('lat='+string(xy[1,k,0,0]),/remov)

    ibeg=1
    iend=nlon-1
    istep=5
    nlines=(iend-ibeg+1.)/istep+1
    print,nlines

    yout=10.^(findgen(nlines)*2./(nlines-1)-0.5)
    xyouts,xmin+100.,yout[0],strcompress('lon='+string(xy[0,0,0,0]),/remov)
    
    loadct,13
    j=1
    for i=ibeg,iend,istep do begin
        cind=(i-ibeg+1)*255./(iend-ibeg+1)
        oplot,xy[5,k,i,*],pres,color=cind,thick=7
        xyouts,xmin+100.,yout[j],strcompress('lon='+string(xy[0,0,i,0]),/remov),color=cind
        j++
    endfor

    loadct,0
    !p.font=-1
    psclose

endfor

;0: lon, 1: lat, 2: level, 3: u, 4: v, 5: T
;psopen,'tprofs_lon0',/enc,/color,/landscape
;!p.font=0
;xmin=min(xy[5,*,0,*])
;xmax=max(xy[5,*,0,*])
;plot,xy[5,0,0,*],pres,yr=[max(pres),min(pres)],xr=[xmin,xmax],$
;     /ystyle,/xstyle,/ylog,xthick=7,ythick=7,thick=7
;ibeg=1
;iend=nlat-1
;istep=5
;nlines=(iend-ibeg+1.)/istep+1
;nlines=ceil(nlines)
;print,nlines
;yout=10.^(findgen(nlines)*2./(nlines-1)-0.5)
;xyouts,xmin+100.,yout[0],strcompress('lat='+string(xy[1,0,0,0]),/remov)
;loadct,13
;j=1
;for i=ibeg,iend,istep do begin
;    cind=(i-ibeg+1)*255./(iend-ibeg+1)
;    oplot,xy[5,i,0,*],pres,color=cind,thick=7
;    xyouts,xmin+100.,yout[j],strcompress('lat='+string(xy[1,i,0,0]),/remov),color=cind
;    j++
;endfor
;loadct,0
;!p.font=-1
;psclose

end
