pro N2,p0,oom,ga,gascon,akap,psfile=psfile,allplots=allplots,tprof=tprof,tlimit=tlimit,nlimit=nlimit,prcb=prcb,pot=pot

;N^2=(g/theta)*d(theta)/d(z), theta=T(P_s/P)^(R/cp), R/cp=kappa
;d(p)/d(z)=-rho g, rho=P/RT
;
;p0 in bar
;ga in m/s^2
;gascon in J/kg/K
;akap unitless
;
; use tprof to feed in [n,nlev] array to overplot N2 values
; tlimit is [nlev] array giving most stable possible profile, used to
;    determine where to set boundary in N2


;grab fort.26 info
openr, 17, 'fort.26'
readf,17,nlat,nlon,nlev
xy=fltarr(6,nlat,nlon,nlev)  ;0: lon, 1: lat, 2: level, 3: u, 4: v, 5: T
readf,17,xy
close,17

;grab psurf
openr, 18, 'fort.50'
readf, 18, nlat,nlon
ab=fltarr(3,nlat,nlon)
readf, 18, ab
close, 18
sp=reform(ab[2,*,*])
sp = (sp+1.)*p0


if keyword_set(psfile) then begin
   if not keyword_set(allplots) then begin
      if (size(psfile))[1] eq 2 then filename='plot.eps' else filename=psfile
      psopen,filename,/color,xsize=25,ysize=20 ,/enc;,bits_per_pixel=24
   endif else begin
      if (size(psfile))[1] eq 2 then filename='plot.ps' else filename=psfile
      psopen,filename,/color,xsize=25,ysize=20 ;,/enc,bits_per_pixel=24
   endelse
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
!x.style=1

sigma=get_sigma(oom,nlev)
ymax=p0
ymin=min(sigma)*p0
bvf=fltarr(nlat,nlon,nlev)

for ilon=0,nlon-1 do begin
   for ilat=0,nlat-1 do begin
      pres=sp[ilat,ilon]*sigma    ;*1.e5        ; to convert to SI (Pa) (cancels out so not necessary)
      dtdp=deriv(alog(pres),alog(reform(xy[5,ilat,ilon,*])))
      bvf[ilat,ilon,*]=(ga^2/gascon/xy[5,ilat,ilon,*])*(akap-dtdp)
   endfor
endfor

if keyword_set(tprof) then begin
   sz=size(tprof)
   if (sz[0] eq 2) and (sz[2] eq nlev) then begin
      obvf=tprof
      for i=0,sz[1]-1 do begin
         dtdp=deriv(alog(p0*sigma),alog(reform(tprof[i,*])))
         obvf[i,*]=(ga^2/gascon/tprof[i,*])*(akap-dtdp)         
      endfor
   endif else begin
      print,'tprof unacceptable size'
      stop
   endelse
endif
if keyword_set(tlimit) then begin
   bvflim=fltarr(nlev)
   dtdp=deriv(alog(p0*sigma),alog(tlimit))
   bvflim=(ga^2/gascon/tlimit)*(akap-dtdp)
endif else begin
   if not keyword_set(nlimit) then begin
      negbvf=where(bvf le 0.)
      if negbvf[0] ne -1 then begin
         mnegbvf=median(bvf[negbvf])
         print,'median negative value:',mnegbvf
         xnegbvf=max(abs(bvf[negbvf]))
         print,'maximum negative value:',-xnegbvf
         noise=where(abs(bvf) le xnegbvf)

         if keyword_set(allplots) then begin 
            if not keyword_set(psfile) then window,2 else erase
            loadct,5,/silent
            nbins=50.
            hist=myhistogram(bvf[noise],bins=bins,nbins=nbins,xtitle='N^2')
            gfit=gaussfit(bins,hist,coeffs,nterms=3)
            oplot,bins,gfit,psym=10,color=50
            mrange=where(hist ge max(hist)*0.75)
            zfit=gaussfit(bins[0:max(mrange)],hist[0:max(mrange)],zcoeffs,nterms=3)
            oplot,bins,zfit,psym=10,color=150
            gauss=zcoeffs[0]*exp(-(bins^2)/2./(zcoeffs[2]+abs(zcoeffs[1]))^2)
            oplot,bins,gauss,color=150 ;,thick=3
            onesig=zcoeffs[2]+abs(zcoeffs[1])
            print,'1-sigma:',onesig
            oop=where(gauss[0:nbins-5] gt hist[0:nbins-5])
            mxnoise=max(bins[oop])
            print,'highest value before excess:',mxnoise
            loadct,0,/silent
         endif
      endif
   endif
endelse

xmin=min(abs(bvf),max=xmax)

if not keyword_set(psfile) then window,0,xsize=1000,ysize=700 else erase
plot,[0,0],[0,0],yr=[ymax,ymin],/ystyle,/ylog,xr=[xmin,xmax],$
     xtit='N!E2!N [s!E-2!N]',ytit='Pressure [bar]',/xlog
loadct,5,/silent
for ilon=0,nlon-1 do begin
   for ilat=0,nlat-1 do begin
      oplot,bvf[ilat,ilon,*],sigma*sp[ilat,ilon]
;      oplot,bvf[ilat,ilon,0:nlev-2],sigma[0:nlev-2]*sp[ilat,ilon]
;      oplot,[bvf[ilat,ilon,nlev-1],bvf[ilat,ilon,nlev-1]],$
;            [sigma[nlev-1]*sp[ilat,ilon],sigma[nlev-1]*sp[ilat,ilon]],psym=4
   endfor
endfor
for ilon=0,nlon-1 do begin
   for ilat=0,nlat-1 do begin
      oplot,-bvf[ilat,ilon,*],sigma*sp[ilat,ilon],color=100
;      oplot,-bvf[ilat,ilon,0:nlev-2],sigma[0:nlev-2]*sp[ilat,ilon],color=100
;      oplot,[-bvf[ilat,ilon,nlev-1],-bvf[ilat,ilon,nlev-1]],$
;            [sigma[nlev-1]*sp[ilat,ilon],sigma[nlev-1]*sp[ilat,ilon]],psym=4,color=100
   endfor
endfor
if keyword_set(tprof) then begin
   for i=0,sz[1]-1 do oplot,obvf[i,*],sigma*p0,color=50,thick=7
   for i=0,sz[1]-1 do oplot,-obvf[i,*],sigma*p0,color=175,thick=7
endif
if keyword_set(tlimit) then begin
   oplot,bvflim,sigma*p0,color=50,thick=5
   al_legend,['Positive N^2','Negative N^2','Stability limit'],color=[0,100,50],linestyle=[0,0,0]
endif
if keyword_set(nlimit) then begin
   oplot,[nlimit,nlimit],[ymax,ymin],color=50,thick=5
   al_legend,['Positive N^2','Negative N^2','Stability limit'],color=[0,100,50],linestyle=[0,0,0]
endif

loadct,0,/silent

if keyword_set(pot) then begin
   pot=reform(xy[5,*,*,*])
   for ilev=0,nlev-1 do pot[*,*,ilev]*=(1./sigma[ilev])^akap
   window,2
   xmin=min(pot,max=xmax)
   plot,[1,1],[1,1],xtit='Potential temperature [K]',ytit='Pressure [bar]',/ylog,yr=[p0,min(sigma)*p0],$
        xr=[xmin,xmax],/xlog
   for ilat=0,nlat-1 do begin
      for ilon=0,nlon-1 do begin
         oplot,pot[ilat,ilon,*],p0*sigma
      endfor
   endfor
endif

if keyword_set(allplots) then begin

   if (keyword_set(tlimit) or keyword_set(nlimit)) then begin

      if keyword_set(nlimit) then blimit=nlimit else blimit=bvflim
      prcb=sp
      for ilat=0,nlat-1 do begin
         for ilon=0,nlon-1 do begin
            toss=where((bvf[ilat,ilon,*] le blimit)); or (bvf[ilat,ilon,*] le 0))
            if ((toss[0] ne -1) and (max(toss) eq nlev-1)) then begin
               check=1
               convind=nlev
               itoss=n_elements(toss)-1
               while check eq 1 do begin
                  if toss[itoss] eq convind-1 then begin
                     convind-=1
                     itoss-=1
                  endif else check=0
               endwhile
               prcb[ilat,ilon]*=sigma[convind]
            endif
         endfor
      endfor
      if not keyword_set(psfile) then window,2 else erase
      pprcb=fltarr(nlat,nlon+1)
      pprcb[*,0:nlon-1]=prcb
      pprcb[*,nlon]=prcb[*,0]
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
         MAP_SET, /MILLER_CYLINDRICAL,0,0,/ISOTROPIC;,title='Lowest P convective',color=255
         cbottom=0
         nlevels=60
         step=(alog10(pmax)-alog10(pmin))/(nlevels-1)
         mylevels=indgen(nlevels)*step+alog10(pmin)
         ccolors=findgen(nlevels)*(255.-cbottom)/(nlevels-1)+cbottom
         contour, aprcb, flon,flat,/overplot ,/cell_fill,/closed,$
                  levels=mylevels
         ndiv=5
         divstep=(alog10(pmax)-alog10(pmin))/ndiv
         plevels=pmin*10.^(divstep/4.)*10.^(findgen(ndiv*2)/4.)
         pdiv=pmin*10.^(findgen(ndiv+1)*divstep)
         contour, pprcb, flon, flat,/over,/follow,levels=plevels,c_charsize=1.5
         colorbar,  position=[0.1,0.07,0.90,0.10], range=alog10([pmin,pmax]),$
                    charsize=1.5,divisions=5,color=255,ticknames=string(pdiv) ;,format='(e8.2)'
         loadct,0,/silent
      endif
;stop

   endif else begin

      oplot,[-mnegbvf,-mnegbvf],[ymax,ymin],linestyle=2
      oplot,[xnegbvf,xnegbvf],[ymax,ymin],linestyle=1
      loadct,5,/silent
      if keyword_set(psfile) then lthick=6 else lthick=3
      oplot,[mxnoise,mxnoise],[ymax,ymin],color=50,linestyle=2,thick=lthick
      oplot,[onesig,onesig],[ymax,ymin],color=50,linestyle=1
      if not keyword_set(psfile) then begin
         lcolors=[255,100,255,255,50,50]
         lthick=[1,1,1,1,1,3]
      endif else begin
         lcolors=[0,100,0,0,50,50]
         lthick=[3,3,3,3,3,6]
      endelse
      al_legend,['Positive N^2','Negative N^2',$
                 'Maximum negative value of N^2','Median negative value of N^2',$
                 'One sigma deviation of N^2 from zero','Value where positive N^2 above gaussian fit'],$
                linestyle=[0,0,1,2,1,2],thick=lthick,color=lcolors,charsize=1
      loadct,0,/silent

      prcb=sp
      toss=where(bvf le mxnoise)
      convind=array_indices(bvf,toss)
      for i=0,n_elements(toss)-1 do begin
         pconv=sigma[convind[2,i]]*sp[convind[0,i],convind[1,i]]
         if pconv lt prcb[convind[0,i],convind[1,i]] then prcb[convind[0,i],convind[1,i]]=pconv
      endfor

      if not keyword_set(psfile) then window,3 else erase
      pprcb=fltarr(nlat,nlon+1)
      pprcb[*,0:nlon-1]=prcb
      pprcb[*,nlon]=prcb[*,0]
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
         MAP_SET, /MILLER_CYLINDRICAL,0,0,/ISOTROPIC,title='Lowest P convective',color=255
         cbottom=0
         nlevels=60
         step=(alog10(pmax)-alog10(pmin))/(nlevels-1)
         mylevels=indgen(nlevels)*step+alog10(pmin)
         ccolors=findgen(nlevels)*(255.-cbottom)/(nlevels-1)+cbottom
         
         contour, aprcb, flon,flat,/overplot ,/cell_fill,/closed,$
                  levels=mylevels
               ;NLEVELS=45,min_value=alog10(pmin), max_value=alog10(pmax)
         ndiv=5
         divstep=(alog10(pmax)-alog10(pmin))/ndiv
         plevels=pmin*10.^(divstep/4.)*10.^(findgen(ndiv*2)/4.)
         pdiv=pmin*10.^(findgen(ndiv+1)*divstep)
         contour, pprcb, flon, flat,/over,/follow,levels=plevels,c_charsize=1.5
;      map_grid, /label,charsize=1  ,color=0
         colorbar,  position=[0.1,0.07,0.90,0.10], range=alog10([pmin,pmax]),$
                    charsize=1.5,divisions=5,color=255,ticknames=string(pdiv) ;,format='(e8.2)'
         loadct,0,/silent
      endif

      oprcb=sp
      sprcb=sp
      sprcbmed=sp
      sprcbmax=sp
      mprcb=sp
      cs=0
      for ilat=0,nlat-1 do begin
         if ilat lt nlat/2 then latc=ccolors[ilat] else latc=ccolors[nlat-1-ilat]
         lat=where(convind[0,*] eq ilat)
         for ilon=0,nlon-1 do begin
            spot=where(convind[1,lat] eq ilon)
            if spot[0] ne -1 then begin
               nspot=n_elements(spot)
               ispot=0
               deepest=max((convind[2,lat])[spot])
               pfconv=(convind[2,lat])(spot[ispot])
               if mean(bvf[ilat,ilon,pfconv:nlev-2]) lt mxnoise then $
                  mprcb[ilat,ilon]=sigma[pfconv]*sp[ilat,ilon]
               spotcheck=0
               if pfconv le nlev-4 then begin
                                ; don't use final point for smoothing, as it's always high
                  sbvf=smooth(reform(bvf[ilat,ilon,pfconv:nlev-2]),3,/edge_tr)
                  toss=where(sbvf le mxnoise)
                  if toss[0] ne -1 then begin
                     ntoss=n_elements(toss)
                     itoss=0
                     while (spotcheck eq 0) and (itoss lt ntoss) do begin
                        if (ntoss-itoss eq nlev-1-pfconv) then begin
                           sprcb[ilat,ilon]=sigma[pfconv+toss[itoss]]*sp[ilat,ilon]
                           spotcheck=1
                        endif else itoss+=1
                     endwhile
                  endif
                  spotcheck=0
                  toss=where(sbvf le mnegbvf)
                  if toss[0] ne -1 then begin
                     ntoss=n_elements(toss)
                     itoss=0
                     while (spotcheck eq 0) and (itoss lt ntoss) do begin
                        if (ntoss-itoss eq nlev-1-pfconv) then begin
                           sprcbmed[ilat,ilon]=sigma[pfconv+toss[itoss]]*sp[ilat,ilon]
                           spotcheck=1
                        endif else itoss+=1
                     endwhile
                  endif
                  spotcheck=0
                  toss=where(sbvf le xnegbvf)
                  if toss[0] ne -1 then begin
                     ntoss=n_elements(toss)
                     itoss=0
                     while (spotcheck eq 0) and (itoss lt ntoss) do begin
                        if (ntoss-itoss eq nlev-1-pfconv) then begin
                           sprcbmax[ilat,ilon]=sigma[pfconv+toss[itoss]]*sp[ilat,ilon]
                           spotcheck=1
                        endif else itoss+=1
                     endwhile
                  endif
               endif
               spotcheck=0
               while (spotcheck eq 0) and (ispot lt nspot) do begin
                  pfconv=(convind[2,lat])(spot[ispot])
                  if (pfconv+nspot-ispot eq nlev) or ((pfconv+nspot-ispot eq nlev-1) and (deepest eq nlev-2)) then begin
                     oprcb[ilat,ilon]=sigma[pfconv]*sp[ilat,ilon]
                     spotcheck=1
                  endif else begin
                     ispot+=1
;                  print,(convind[2,lat])[spot]
                  endelse
               endwhile
            endif else cs+=1
         endfor
      endfor
      print,'Percent of columns totally stable:',cs/(nlat*nlon)

      if not keyword_set(psfile) then window,6 else erase
      pprcb=fltarr(nlat,nlon+1)
      pprcb[*,0:nlon-1]=oprcb
      pprcb[*,nlon]=oprcb[*,0]
      aprcb=alog10(pprcb)
      pmin=min(oprcb,max=pmax)
      if pmin ne pmax then begin
         ctload,10,/silent,/reverse
         MAP_SET, /MILLER_CYLINDRICAL,45,0,/ISOTROPIC,title='P_RCB when convective to bottom-1',color=255
         cbottom=0
         nlevels=100
         step=(alog10(pmax)-alog10(pmin))/(nlevels-1)
         mylevels=indgen(nlevels)*step+alog10(pmin)
         ccolors=findgen(nlevels)*(255.-cbottom)/(nlevels-1)+cbottom

         contour, aprcb, flon,flat,/overplot ,/cell_fill,/closed,$
                  levels=mylevels
                                ;min_value=alog10(pmin), max_value=alog10(pmax),NLEVELS=45
         ndiv=5
         divstep=(alog10(pmax)-alog10(pmin))/ndiv
         plevels=pmin*10.^(divstep/4.)*10.^(findgen(ndiv*2)/4.)
         pdiv=pmin*10.^(findgen(ndiv+1)*divstep)
;      contour, pprcb, flon, flat,/over,/follow,levels=plevels,c_charsize=1.5
         map_grid, /label,charsize=1  ,color=0
         colorbar,  position=[0.1,0.07,0.90,0.10], range=alog10([pmin,pmax]), $
                    charsize=1.5,divisions=5,color=255,ticknames=string(pdiv) ;,format='(e8.2)'
         loadct,0,/silent
      endif

      if not keyword_set(psfile) then window,7 else erase
      pprcb=fltarr(nlat,nlon+1)
      pprcb[*,0:nlon-1]=sprcb
      pprcb[*,nlon]=sprcb[*,0]
      aprcb=alog10(pprcb)
      pmin=min(oprcb,max=pmax)
      if pmin ne pmax then begin
         ctload,10,/silent,/reverse
         MAP_SET, /MILLER_CYLINDRICAL,45,0,/ISOTROPIC,title='P_RCB when vertically smoothed',color=255
         cbottom=0
         nlevels=100
         step=(alog10(pmax)-alog10(pmin))/(nlevels-1)
         mylevels=indgen(nlevels)*step+alog10(pmin)
         ccolors=findgen(nlevels)*(255.-cbottom)/(nlevels-1)+cbottom

         contour, aprcb, flon,flat,/overplot ,/cell_fill,/closed,$
                  levels=mylevels
               ;min_value=alog10(pmin), max_value=alog10(pmax),NLEVELS=45
         ndiv=5
         divstep=(alog10(pmax)-alog10(pmin))/ndiv
         plevels=pmin*10.^(divstep/4.)*10.^(findgen(ndiv*2)/4.)
         pdiv=pmin*10.^(findgen(ndiv+1)*divstep)
;      contour, pprcb, flon, flat,/over,/follow,levels=plevels,c_charsize=1.5
         map_grid, /label,charsize=1  ,color=0
         colorbar,  position=[0.1,0.07,0.90,0.10], range=alog10([pmin,pmax]), $
                    charsize=1.5,divisions=5,color=255,ticknames=string(pdiv) ;,format='(e8.2)'
         loadct,0,/silent
      endif
      if not keyword_set(psfile) then window,8 else erase
      pprcb=fltarr(nlat,nlon+1)
      pprcb[*,0:nlon-1]=mprcb
      pprcb[*,nlon]=mprcb[*,0]
      aprcb=alog10(pprcb)
      pmin=min(oprcb,max=pmax)
      if pmin ne pmax then begin
         ctload,10,/silent,/reverse
         MAP_SET, /MILLER_CYLINDRICAL,45,0,/ISOTROPIC,title='P_RCB when vertically averaged',color=255
         cbottom=0
         nlevels=100
         step=(alog10(pmax)-alog10(pmin))/(nlevels-1)
         mylevels=indgen(nlevels)*step+alog10(pmin)
         ccolors=findgen(nlevels)*(255.-cbottom)/(nlevels-1)+cbottom

         contour, aprcb, flon,flat,/overplot ,/cell_fill,/closed,$
                  levels=mylevels
               ;min_value=alog10(pmin), max_value=alog10(pmax),NLEVELS=45
         ndiv=5
         divstep=(alog10(pmax)-alog10(pmin))/ndiv
         plevels=pmin*10.^(divstep/4.)*10.^(findgen(ndiv*2)/4.)
         pdiv=pmin*10.^(findgen(ndiv+1)*divstep)
;      contour, pprcb, flon, flat,/over,/follow,levels=plevels,c_charsize=1.5
         map_grid, /label,charsize=1  ,color=0
         colorbar,  position=[0.1,0.07,0.90,0.10], range=alog10([pmin,pmax]), $
                    charsize=1.5,divisions=5,color=255,ticknames=string(pdiv) ;,format='(e8.2)'
         loadct,0,/silent
      endif

;average profiles (or separated?)
      if not keyword_set(psfile) then window,9 else erase
      !P.MULTI=[0,3,1,1,0]
;loadct,6,/silent
;lcolors=[64,128,192]
      i60=where(abs(xy[1,*,0,0]) ge 60.)
      i30=where((abs(xy[1,*,0,0]) lt 60.) and (abs(xy[1,*,0,0]) ge 30.))
      i00=where((abs(xy[1,*,0,0]) lt 30.) and (abs(xy[1,*,0,0]) ge 0.))
      plot,[0,0],[0,0],yr=[ymax,min(prcb)],/ystyle,/ylog,xr=[-xnegbvf,xnegbvf],$
           xtit='N!E2!N [s!E-2!N]',ytit='Pressure [bar]'
      loadct,25,/silent
      for ilat=0,n_elements(i60)-1 do begin
         for ilon=0,nlon-1 do begin
            oplot,bvf[i60[ilat],ilon,*],sigma*sp[i60[ilat],ilon],color=255./(nlon-1)*ilon
         endfor
      endfor
      loadct,0,/silent
      oplot,[mxnoise,mxnoise],[ymax,ymin],linestyle=2
      xyouts,-xnegbvf*0.8,ymax/2,'Polar',charsize=1
      plot,[0,0],[0,0],yr=[ymax,min(prcb)],/ystyle,/ylog,xr=[-xnegbvf,xnegbvf],$
           xtit='N!E2!N [s!E-2!N]',ytit='Pressure [bar]'
      loadct,25,/silent
      for ilat=0,n_elements(i30)-1 do begin
         for ilon=0,nlon-1 do begin
            oplot,bvf[i30[ilat],ilon,*],sigma*sp[i30[ilat],ilon],color=255./(nlon-1)*ilon
         endfor
      endfor
      loadct,0,/silent
      oplot,[mxnoise,mxnoise],[ymax,ymin],linestyle=2
      xyouts,-xnegbvf*0.8,ymax/2,'Mid-lats',charsize=1
      plot,[0,0],[0,0],yr=[ymax,min(prcb)],/ystyle,/ylog,xr=[-xnegbvf,xnegbvf],$
           xtit='N!E2!N [s!E-2!N]',ytit='Pressure [bar]'
      loadct,25,/silent
      for ilat=0,n_elements(i00)-1 do begin
         for ilon=0,nlon-1 do begin
            oplot,bvf[i00[ilat],ilon,*],sigma*sp[i00[ilat],ilon],color=255./(nlon-1)*ilon
         endfor
      endfor
      loadct,0,/silent
      oplot,[mxnoise,mxnoise],[ymax,ymin],linestyle=2
      xyouts,-xnegbvf*0.8,ymax/2,'Equatorial',charsize=1
      !P.MULTI=0
   endelse

endif
;stop

if keyword_set(psfile) then begin
   psclose
   !p.font=-1
   !p.charsize=1
   !p.thick=1
   !x.thick=1
   !y.thick=1

   if not keyword_set(allplots) then begin
      spawn,'gs -r300 -dEPSCrop -dTextAlphaBits=4 -sDEVICE=png16m -sOutputFile=plot.png -dBATCH -dNOPAUSE plot.eps'
      spawn,'convert plot.png eps3:plot.eps'
   endif

endif

end
