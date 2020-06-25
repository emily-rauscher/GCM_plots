function zguillot, press,lat,tint=tint,kth=kth,alpha=alpha,tirr=tirr,g=g,kv=kv,$
                   psfile=psfile,complex=complex,noplot=noplot
; input press: pressure in bars
;       lat: latitude in degrees
;
;  THIS CALCULATES AVERAGE TEMPERATURES, BUT REALLY FLUX IS THE BETTER
;  QUANTITY TO AVERAGE AND THE AVERAGE TEMPERATURE DOESN'T HAVE MUCH MEANING

if not keyword_set(tint) then tint=500. ; in K
if not keyword_set(g) then g=8.e2 ; in cgs
if not keyword_set(tirr) then tirr=2077.5 ; in K, =2077.5 for fig 2, =2013 for solc=9.31e5
if not keyword_set(kth) then kth=1.e-2 ; in cgs
if not keyword_set(kv) then kv=4.e-3 ; in cgs
if not keyword_set(alpha) then alpha=0

nlev=n_elements(press)
nlat=n_elements(lat)
nlon=nlat  ; only need longitudes from -90 to 90
temp=fltarr(nlat,nlev)

mlat=max(lat)*!pi/180.
lon=findgen(nlon)/(nlon-1)*2.*mlat-mlat  ; just short of +/- pi/2
clon=cos(lon)

if not keyword_set(complex) then begin

   tempg=fltarr(nlon,nlev)
   tnight=guillot(press,1.,tint=tint,kth=kth,alpha=alpha,tirr=tirr,g=g,kv=kv,/night)

   for ilat=0,nlat-1 do begin
      clat=cos(lat[ilat]*!pi/180.)
      for ilon=0,nlon-1 do begin
         tg=guillot(press,clat*clon[ilon],tint=tint,kth=kth,alpha=alpha,tirr=tirr,g=g,kv=kv)
         tempg[ilon,*]=tg
      endfor
      for ilev=0,nlev-1 do begin
         tday=(1./!pi)*int_tabulated(lon,tempg[*,ilev])
         temp[ilat,ilev]=0.5*(tnight[ilev]+tday)
      endfor
   endfor

endif else begin

; tau= kth/g * P, convert from bar to dyne/cm2 with 1.e6
; tau= kth/g *P * (P/P_ref)^alpha
   tau=(kth/g)*1.e6*press^(alpha+1) ;unitless

   gamma=kv/kth
   gamma=gamma/(press^alpha)
; (kth = k_o * (press/1 bar)^alpha, so gamma unitless

   tnight=((3./4)*tint^4*(2./3 + tau/(alpha+1)))^0.25

   lintgrl=fltarr(nlat,nlon,nlev)
   for ilat=0,nlat-1 do begin
      clat=cos(!pi/180.*lat[ilat])   
      for ilon=0,nlon-1 do begin
         func=[0,(press^alpha)*exp(-gamma*tau/clat/clon[ilon])] ;unitless
         funcx=[0,press*1.e6]
         intgrl=fltarr(nlev)
         for i=0,nlev-1 do begin
     ;integration wrong if func decreases too quickly
            test=int_tabulated(funcx[0:i+1],func[0:i+1])
            if ((i gt 0) and (test lt intgrl[i-1])) then $
               intgrl[i]=intgrl[i-1]+func[i+1]*(funcx[i+1]-funcx[i]) $
            else intgrl[i]=test
         endfor
         lintgrl[ilat,ilon,*]=intgrl
      endfor
   endfor

   for ilat=0,nlat-1 do begin
      clat=cos(!pi/180.*lat[ilat])
      for ilev=0,nlev-1 do begin
         bfunc=(replicate(tnight[ilev]^4,nlon) + 0.5*tirr^4*clat*clon $
                + (gamma[ilev]/4.)*tirr^4*exp(-gamma[ilev]*tau[ilev]/clat/clon) $
                + 0.75*tirr^4*(kth/g)*clat*clon*lintgrl[ilat,*,ilev])^0.25
         bintgrl=int_tabulated(lon,bfunc)
         tday=(1./!pi)*bintgrl
         temp[ilat,ilev]=0.5*(tnight[ilev]+tday)
      endfor
   endfor
endelse

if not keyword_set(noplot) then begin
   if not keyword_set(psfile) then begin
      device,true_color=24,decomposed=0,retain=2
      !p.font=-1
      !p.charsize=1.5
   endif else begin
      if (size(psfile))[1] eq 2 then filename='plot.eps' else filename=psfile
      psopen,filename,/enc,/color,xsize=25,ysize=25 ;,bits_per_pixel=24
      !p.font=0
      !p.charsize=2
   endelse

   tmax=max(temp,min=tmin)
   nlevels=45
   step=(tmax-tmin)/nlevels
   tlevels=indgen(nlevels)*step+tmin

   loadct,2,/silent
   contour, temp,lat,press, /cell_fill,levels=tlevels,/ylog,yr=[max(press),min(press)],/xstyle,/ystyle,$
            xtit='Latitude [degrees]',ytit='Pressure [bar]',ymargin=[4,4]
   contour, temp,lat,press, /over,/follow,nlevels=nlev
   colorbar,position=[0.23,0.965,0.9,0.98],range=[tmin,tmax],format='(i5)',charsize=2
   print,tmin,tmax
   loadct,0,/silent

   if keyword_set(psfile) then begin
      psclose
      spawn,'gs -r300 -dEPSCrop -dTextAlphaBits=4 -sDEVICE=png16m -sOutputFile=plot.png -dBATCH -dNOPAUSE plot.eps'
      spawn,'convert plot.png eps3:plot.eps'
   endif
endif

return,temp

end
