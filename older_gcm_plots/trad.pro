pro trad,psfile=psfile

  p0=1.e2  ;bar
  oom=5
  gascon=3523.  ;SI
  akap=0.286
  radea=1.e8  ;m
  ga=8.  ;SI
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
tsec=fltarr(nlat,nlon,nlev)
for ilev=0,nlev-1 do begin
   rhocgs[*,*,ilev]=reform(sp[*,*]*sigma[ilev]/gascon/xy[5,*,*,ilev])*1.e-3 ;cgs
   tsec[*,*,ilev]=(gascon/akap)*sp[*,*]*sigma[ilev]/ga/5.67e-8/(xy[5,*,*,ilev])^3
endfor

rho=reform(rhocgs,nlat*nlev*nlon)
temp=reform(xy[5,*,*,*],nlat*nlev*nlon)
taur=reform(tsec,nlat*nlev*nlon)
ltau=alog10(taur)
lmin=min(ltau,max=lmax)

if keyword_set(psfile) then begin
   psopen,'plot.eps',/enc,/color,/landscape
   !p.font=0
   !x.thick=4
   !y.thick=4
   !p.charsize=2
endif else begin
   device,true_color=24,decomposed=0,retain=2
endelse

plot,rho,temp,/xlog,psym=3,/ynoz,xtit='Density [cgs]',ytit='Temperature [K]',charsize=2,xthick=4,ythick=4
loadct,13,/silent
for i=0,n_elements(rho)-1 do oplot,[rho[i],rho[i]],[temp[i],temp[i]],color=255-(ltau[i]-lmin)*255./(lmax-lmin),$
                                   symsize=2
colorbar,range=[lmin,lmax],format='(f5.2)',position=[0.2,0.88,0.84,0.92],charsize=1.5,/reverse
loadct,0,/silent

if keyword_set(psfile) then begin
   !p.font=-1
   !x.thick=1
   !y.thick=1
   !p.charsize=1
   psclose
endif


end
