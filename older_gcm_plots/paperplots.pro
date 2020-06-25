pro paperplots

;(for 58w)
;run getfiles script first

;for T,wind maps at L3,15,27:
; (may need to play with fraction keyword in partvelvec,
;  and/or xsize,ysize for ps file size)
;igcm,3,2,psfile='L3.eps'
;igcm,3,14,psfile='L15.eps'
;igcm,3,26,psfile='L27.eps'
;igcm,3,22,psfile='L23.eps'
;igcm,3,10,psfile='L11.eps'
;igcm,3,12,psfile='L13.eps'
;igcm,3,13,psfile='L14.eps'

;for [u] profile:
;igcm_uz,nom=5.34363,psfile='uz.eps',p0=220.16

;for T-sigma profiles:
Tprofiles,oom=5.34363,tprof=tprof,sigma=sigma
tnight=[511.,524.,539.,555.,570.,585.,598.,618.,653.,690.,$
        730.,769.,817.,890.,963.,1049.,1125.,1190.,1258.,1299.,$
        1345.,1387.,1437.,1477.,1497.,1780.,1782.,1787.,1791.,1793.,$
        1796.,1801.,1808.]
tday=[replicate(1000.,12),994.,956.,918.,880.,842.,804.,766.,728.,$
      690.,652.,614.,576.,538.,replicate(0.,8)]
tday+=tnight
xmin=min([reform(tprof,6*33),tnight,tday],max=xmax)
pres=sigma*220.16
psopen,'Tprof.eps',/enc,xsize=25,ysize=20
  !p.font=0
  !p.charsize=2
  plot,tprof[0,*],pres,xr=[xmin,xmax],yr=[max(pres),min(pres)],/ylog,/ystyle,$
       xthick=8,ythick=8,xtitle='temperature [K]',ytitle='pressure [bar]',thick=8,/xstyle
  for i=1,5 do oplot,tprof[i,*],pres,linestyle=i,thick=8
  oplot,tday,pres,thick=14
  oplot,tnight,pres,thick=14
psclose

;for N^2 profiles:
;N2,sigma,tprof,4593.,9.42,0.321,psfile='N2.eps'

;for KE diag:
;KE_diag,220.16,33,1450,oom=5.34363,psfile='KE_diag.eps'
;KE_diag,220.16,33,800,oom=5.34363,psfile='KE_diag.eps'
;KE_diag,220.16,33,800,oom=5.34363,dayrange=[0.,100.],psfile='KE_diag_start.eps'

;for jump:
;DV_map,6,psfile='L7div.eps',/moreinfo,vfrac=0.3
;DV_map,7,psfile='L8div.eps',/moreinfo,vfrac=0.3
;DV_map,8,psfile='L9div.eps',/moreinfo,vfrac=0.3
;DV_map,9,psfile='L10div.eps',/moreinfo,vfrac=0.3
;DV_map,10,psfile='L11div.eps',/moreinfo,vfrac=0.0003
;DV_map,12,psfile='L13div.eps',/moreinfo,vfrac=0.0003
;DV_map,13,psfile='L14div.eps',/moreinfo,vfrac=0.0003
;psopen,'jump.eps',/enc,/color,xsize=25,ysize=20
;  !P.FONT=0
;  !X.STYLE=1
;  !P.CHARSIZE=1.5
;  openr,20,'fort.26'
;  readf,20,nlat,nlon,nlev
;  xy=fltarr(6,nlat,nlon,nlev)
;  readf,20,xy
;  close,20
;  lev=13
;  lon=reform(xy[0,*,*,lev]) & lat=reform(xy[1,*,*,lev])
;  u=reform(xy[3,*,*,lev]) & v=reform(xy[4,*,*,lev]) & temp=reform(xy[5,*,*,lev])
;  flon=lon*360./lon[0,nlon-1] & flat=lat*90./lat[0,nlat-1]
;  MAP_SET, /Miller_CYLINDRICAL,0,0,/ISOTROPIC
;  loadct,4
;  contour,temp,flon,flat,/over,/cell_fi,/close,nlevels=55
;  colorbar,position=[0.1,0.07,0.9,0.1],range=[min(temp),max(temp)],format='(i5)',charsize=2
;  partvelvec,u,v,flon,flat,/over,fraction=0.9,color=0
;  map_grid,/label,color=255
;  openr,21,'fort.51'
;  readf,21,nlat,nlon,nlev
;  ab=fltarr(5,nlat,nlon,nlev)
;  readf,21,ab
;  close,21
;  div=reform(ab[3,*,*,lev])
;  contour,div,flon,flat,/over,color=255,levels=[1.e-4,1.5e-4,2.e-4,2.5e-4],thick=5
;  loadct,2
;  contour,div,flon,flat,/over,color=240.,levels=[-2.5e-4,-2.e-4,-1.5e-4,-1.e-4],thick=5
;  loadct,0
;psclose

;for variability (must do this manually):
;$cp -f ~/research/Tests/runs/fortR0058o3.26.day776 fort.26
;igcm,3,14,prange=[930.,1336.],psfile='L15day776.eps'
;$cp -f ~/research/Tests/runs/fortR0058o3.26.day781 fort.26
;igcm,3,14,prange=[930.,1336.],psfile='L15day781.eps'
;$cp -f ~/research/Tests/runs/fortR0058o3.26.day791 fort.26
;igcm,3,14,prange=[930.,1336.],psfile='L15day791.eps'
;$cp -f ~/research/Tests/runs/fortR0058o3.26.day796 fort.26
;igcm,3,14,prange=[930.,1336.],psfile='L15day796.eps'
;$cp -f ~/research/Tests/runs/fortR0058o3.26.day800 fort.26
;igcm,3,14,prange=[930.,1336.],psfile='L15day800.eps'

;for erates:
;erates, 220.16,2.06e-5,33,800

end
