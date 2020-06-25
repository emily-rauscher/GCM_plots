pro tadv,psfile=psfile

  p0=1.e2
  oom=5
  gascon=3523.
  radea=1.e8
;  ga=8.
;  ww=2.1e-5

print,'Check: are parameters correct?',p0,oom,gascon,radea

openr, 17, 'fort.26'
readf,17,nlat,nlon,nlev
xy=fltarr(6,nlat,nlon,nlev)  ;0: lon, 1: lat, 2: level, 3: u, 4: v, 5: T
readf,17,xy
close,17

sigma=get_sigma(oom,nlev)

openr, 18, 'fort.50'
readf, 18, nlat,nlon
ab=fltarr(3,nlat,nlon)
readf, 18, ab
close, 18
sp=reform(ab[2,*,*])
sp = (sp+1.)*p0*1.e5  ;SI

rhocgs=fltarr(nlat,nlon,nlev)
for ilev=0,nlev-1 do rhocgs[*,*,ilev]=reform(sp[*,*]*sigma[ilev]/gascon/xy[5,*,*,ilev])*1.e-3  ;cgs

rho=reform(rhocgs,nlat*nlev*nlon)
temp=reform(xy[5,*,*,*],nlat*nlev*nlon)

zadv=abs(radea/reform(xy[3,*,*,*],nlat*nlev*nlon))
lz=alog10(zadv)
zmax=max(lz,min=zmin)
;zmax=zmin+4.
madv=abs(radea/reform(xy[4,*,*,*],nlat*nlev*nlon))
lm=alog10(madv)
mmax=max(lm,min=mmin)
;mmax=mmin+4.
;lmax=max([lz,lm],min=lmin)

if keyword_set(psfile) then begin
   psopen,'plot.eps',/enc,/color,xsize=25,ysize=35
   !p.font=0
   !x.thick=4
   !y.thick=4
   !p.charsize=2
endif else begin
   device,true_color=24,decomposed=0,retain=2
endelse

;plot,[0,0],[0,0],xtit='Density [cgs]',ytit='Temperature [K]',psym=3,xr=minmax(rho),yr=minmax(temp),/xlog
plot,rho,temp,/xlog,/ynoz,xtit='Density [cgs]',ytit='Temperature [K]',psym=3,position=[0.15,0.1,0.95,0.545]
loadct,13,/silent
for i=0,n_elements(rho)-1 do begin
   if lz[i] ge zmax then lcolor=255 else lcolor=255-(lz[i]-zmin)*255/(zmax-zmin)
   oplot,[rho[i],rho[i]],[temp[i],temp[i]],color=lcolor,psym=3
endfor
colorbar,range=[zmin,zmax],format='(f5.2)',position=[0.2,0.495,0.82,0.525],charsize=1.5,/reverse

plot,rho,temp,/xlog,/ynoz,ytit='Temperature [K]',psym=3,position=[0.15,0.545,0.95,0.99],$
     /noerase,xtickn=replicate(' ',6)
for i=0,n_elements(rho)-1 do begin
   if lm[i] ge mmax then lcolor=255 else lcolor=255-(lm[i]-mmin)*255/(mmax-mmin)
   oplot,[rho[i],rho[i]],[temp[i],temp[i]],color=lcolor,psym=3
endfor
colorbar,range=[mmin,mmax],format='(f5.2)',position=[0.2,0.94,0.82,0.97],charsize=1.5,/reverse

loadct,0,/silent

if keyword_set(psfile) then begin
   !p.font=-1
   !x.thick=1
   !y.thick=1
   !p.charsize=1
   psclose
endif

;stop

end
