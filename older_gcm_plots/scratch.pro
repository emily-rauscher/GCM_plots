pro scratch,psfile=psfile

if keyword_set(psfile) then begin
   psopen,'plot.eps',/enc,xsize=25,ysize=20
   !p.font=0
   !p.charsize=2
   !p.thick=4
   !x.thick=5
   !y.thick=5
endif
;!p.position=[0.2,0.15,0.95,0.95]
;!P.MULTI=[0,2,2]

rho=10.^(-5.175*findgen(901)/900.-2.825)
temp=findgen(901.)*3.+400.
leta=fltarr(901,901)

for i=0,900 do begin
   toss=tmag(replicate(rho[i],901),temp,1.,eta=etai)
   leta[i,*]=alog10(etai)
;   toss=tmag(replicate(rho[i],901),temp,1.,eta=etai,/simp)
endfor

lmin=min(leta,max=lmax,/nan)
print,lmin,lmax
nlevels=40
step=(lmax-lmin)/nlevels
llevels=indgen(nlevels)*step+lmin
;llevels=indgen(nlevels)*0.5+9.
llevels=[9,9.5,10.,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,$
         15.5,16,17,18,19,20,22,24,26,28,30]
;         15.5,16,16.5,17,$
;         17.5,18,18.5,19,19.5,20,22,24,26,28,30]
;22.5,25,27.5,30];21,22,23,24,25,26,27,28,29,30]

;loadct,1
contour,leta,rho,temp,/xlog,/xstyle,/ystyle,levels=llevels,/follow,$
        xtit='Density [cgs]',ytit='Temperature [K]',c_charsize=1.5

if keyword_set(psfile) then begin
;!P.MULTI=0
   !p.charsize=1
   !p.font=-1
   !p.thick=1
   !x.thick=1
   !y.thick=1
   psclose
endif

end
