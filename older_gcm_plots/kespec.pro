pro kespec,lev=lev,oom=oom,p0=p0,spec=spec,psfile=psfile

; plot KE spectrum (can be total or zonal wave no., and by level or
; integrated).  will plot integrated, unless lev set (1=top)
; reads from fort.29

if not keyword_set(oom) then oom=5
if not keyword_set(p0) then p0=100.

if not keyword_set(psfile) then begin
; Screen output
  device,true_color=24,decomposed=0,retain=2
  !p.font=-1
  !p.charsize=1.5
endif else begin
   if (size(psfile))[1] eq 2 then filename='plot.eps' else filename=psfile
   psopen,filename,/enc,/color,xsize=25,ysize=25;,bits_per_pixel=24
  !p.font=0
  !p.charsize=2
endelse

; Read File with 3D fields
  OpenR, 17, 'fort.29'
  ReadF, 17, twn,zwn,nlev
  nlev=nlev-1
;	  print,nlat,nlon,nlev 
  if (twn lt zwn) then begin
     print,'zonal wave no. gt total wave no. (?)'
     abort
  endif
    xy=fltarr(nlev+2,twn+zwn)
  ReadF, 17, xy
  Close, 17

sigma=get_sigma(oom,nlev)
dsigma=fltarr(nlev)
for ilev=0,nlev-1 do begin
   if ilev gt 0 and ilev lt nlev-1 then dsigma[ilev]=sigma[ilev+1]-sigma[ilev-1]
   if ilev eq 0 then dsigma[ilev]=sigma[0]
   if ilev eq nlev-1 then dsigma[ilev]=1.-sigma[nlev-1]
endfor

twaveno=reform(xy[0,0:twn-1])
;print,twaveno
zwaveno=reform(xy[0,twn:zwn+twn-1])
;print,zwaveno

if not keyword_set(lev) then begin
   itke=reform(xy[nlev+1,0:twn-1])
   izke=reform(xy[nlev+1,twn:zwn+twn-1])
   plot,twaveno,itke,xtitle='Total wavenumber',ytitle='Kinetic energy [J/m^2]',xr=[1,zwn],/xlog,charsize=1.5,/ylog,/ystyle,/xstyle
   for i=1,nlev do oplot,twaveno,xy[i,0:twn-1]*dsigma,linestyle=2
   spec=itke
endif else begin
   tke=reform(xy[lev,0:twn-1])
   plot,twaveno,tke,xtitle='Total wavenumber',ytitle='Kinetic energy',xr=[1,zwn],/xlog,charsize=1.5,/ylog,/ystyle,/xstyle
   xyouts,1.3,max(tke)*.2,strcompress('L'+string(lev)),charsize=1.5
   spec=tke
endelse

if keyword_set(psfile) then psclose

;stop

end
