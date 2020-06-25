pro freedman,met,kv=kv

; to read in freedman opacities ... and plot or output or something
;
; met = string(0.0, 0.3, or -0.3), metallicity and identifier of which
;        file to read
; kv = visible opacity (assumed constant, ignore whether Ross or Pla)
;       used to "solve" for IR component

file=strcompress('~/research/data/Freedman_opacities_MH'+met+'.txt',/remove_all)

data=rd_tfile(file,5,-1,/convert)
; data[0,*] is temperature in K
; data[1,*] is pressure in 0.1Pa (cgs)
; data[2,*] is density in g/cm^3
; data[3,*] is Rosseland mean opacity in cm^2/g
; data[4,*] is Planck mean opacity in cm^2/g

if keyword_set(kv) then begin
   test=where(data[3,*] lt kv,complement=otest)
   if test[0] ne -1 then data[3,test]=1./(1./data[3,test]-1/kv)
   if otest[0] ne -1 then data[3,otest]=0.
   data[4,*]=data[4,*]-kv
   test=where(data[4,*] lt 0.)
   if test[0] ne -1 then data[4,test]=0.
   print,'---------------------Opacities "corrected" to be only IR'
endif

test=where(data[3:4,*] ne 0.)
mind=min((data[3:4,*])[test])

plot,data[3,*],data[1,*]/1.e6,psym=3,/xlog,/ylog,yr=[max(data[1,*])/1.e6,min(data[1,*])/1.e6],$
     xr=[mind,max(data[3:4,*])],$
     xtit='mean opacities [cm^2/g]',ytit='Pressure [bar]',charsize=1.5,tit=met
oplot,data[4,*],data[1,*]/1.e6,psym=3

;hot=where(data[0,*] ge 500.)
;oplot,data[3,hot],data[1,hot]/1.e6,psym=4

;pres=300.*get_sigma(6,30)
;tg=guillot(pres,1.,/night,tint=2000)

toss=sort(data[1,*])
presi=uniq(data[1,toss])
pres=(data[1,toss])[presi]
pres=pres(sort(pres))
templ=[replicate(400.,4),500.,650.,800.,1000.,1400.,1800.,2300.,2300.]*0.9

fcount=0
for i=0,n_elements(pres)-1 do begin
   toss=where((data[1,*] eq pres[i]) and (data[0,*] gt templ[i]))
   oplot,data[3,toss],data[1,toss]/1.e6,psym=1
   oplot,data[4,toss],data[1,toss]/1.e6,psym=7
   fcount=fcount+n_elements(toss)
endfor

al_legend,['Rosseland opacities','Planck opacities','(at T outside of interest)'],psym=[1,7,3],/bottom

fdata=fltarr(2,fcount) ;Rosseland
bdata=fltarr(2,fcount) ;Planck
dcount=0
for i=0,n_elements(pres)-1 do begin
   toss=where((data[1,*] eq pres[i]) and (data[0,*] gt templ[i]))
   count=n_elements(toss)
   fdata[0,dcount:dcount+count-1]=data[1,toss]
   fdata[1,dcount:dcount+count-1]=data[3,toss]
   bdata[0,dcount:dcount+count-1]=data[1,toss]
   bdata[1,dcount:dcount+count-1]=data[4,toss]
   dcount=dcount+count
endfor

;remove zeros
fdind=where(fdata[1,*] ne 0.)
bdind=where(bdata[1,*] ne 0.)

;fit powerlaw to fdata: fdata[1,*]=fit[0]*(fdata[0,*])^fit[1]+fit[2]
fit=comfit(reform(fdata[0,fdind]/1.e6),reform(fdata[1,fdind]),[1.e-1,1.,0.],yfit=toplot,/geo)
print,fit
oplot,toplot,fdata[0,fdind]/1.e6
bfit=comfit([reform(fdata[0,fdind]/1.e6),reform(bdata[0,bdind]/1.e6)],[reform(fdata[1,fdind]),reform(bdata[1,bdind])],[1.e-1,1.,0.],yfit=bplot,/geo)
print,bfit
oplot,bplot,[reform(fdata[0,fdind]/1.e6),reform(bdata[0,bdind]/1.e6)]

xyouts,mind*3.,40.,'kappa = k (P/[1 bar])^alpha:'
k=fit[0]
alpha=fit[1]
oplot,k*(bdata[0,bdind]/1.e6)^alpha,bdata[0,bdind]/1.e6,thick=3
oplot,(bfit[0])*(bdata[0,bdind]/1.e6)^(bfit[1]),bdata[0,bdind]/1.e6,thick=3
xyouts,mind*3.,120.,strcompress('k:'+string(k)+', '+string(bfit[0])+' [cm^2/g]')
xyouts,mind*3.,360.,strcompress('alpha:'+string(alpha)+', '+string(bfit[1]))


;stop

end
