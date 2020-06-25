pro Tprofiles,oom=oom,tprof=tprof,sigma=sigma,xmin = xmin,xmax = xmax,unstable = unstable,$
              p0=p0, psurf=psurf,pot=pot,psfile=psfile,ave=ave,noplot=noplot,eqtr=eqtr,$
              xlog=xlog,otprof=otprof,oxyouts=oxyouts,allprof=allprof,lats=lats,lskip=lskip

;reads from fort.26 in same directory
;
; if keyword psurf set, then reads fort.50 and adjusts local pressures accordingly
;can set tprof to export tprofiles, and sigma to export sigma
;can set unstable to mark where profiles becomes convectively unstable
; use otprof to feed in [n,nlev] array to overplot (in solid, thick black lines)
; use eqtr to see profiles along equator or lats to see profiles along lon=0
; use lskip to plot every lskip profile instead of all

if not keyword_set(lskip) then lskip=1
  
;grab fort.26 info
openr, 17, 'fort.26'
readf,17,nlat,nlon,nlev
xy=fltarr(6,nlat,nlon,nlev)  ;0: lon, 1: lat, 2: level, 3: u, 4: v, 5: T
readf,17,xy
close,17

tprof=fltarr(6,nlev)
;0: substellar, 1: antistellar, 2: east eqtr terminator, 
;3: west eqtr terminator, 4: north pole, 5: southpole
;
;substellar is ave of  [T,nlat/2.-1:nlat/2.,0        ,L]
;antistellar is ave of [T,nlat/2.-1:nlat/2.,nlon/2.  ,L]
;east is ave of        [T,nlat/2.-1:nlat/2.,nlon/4.  ,L]
;west is ave of        [T,nlat/2.-1:nlat/2.,nlon*3./4,L]
;north is ave of       [T,0     ,* ,L]
;south is ave of       [T,nlat-1,* ,L]

;calculate location T profiles
if not keyword_set(ave) then begin
   loci=intarr(2,6)
   loci[*,0]=[nlat/2.,0]
   loci[*,1]=[nlat/2.,nlon/2.]
   loci[*,2]=[nlat/2.,nlon/4.]
   loci[*,3]=[nlat/2.,nlon*3./4.]
   loci[*,4]=[0,nlon/4.]
   loci[*,5]=[nlat-1,nlon/4.]
   for i=0,nlev-1 do begin
      tprof[0,i]=xy[5,loci[0,0],loci[1,0],i]
      tprof[1,i]=xy[5,loci[0,1],loci[1,1],i]
      tprof[2,i]=xy[5,loci[0,2],loci[1,2],i]
      tprof[3,i]=xy[5,loci[0,3],loci[1,3],i]
      tprof[4,i]=xy[5,loci[0,4],loci[1,4],i]
      tprof[5,i]=xy[5,loci[0,5],loci[1,5],i]
   endfor
endif else begin
   for i=0,nlev-1 do begin
      tprof[0,i]=mean(xy[5,nlat/2.-1:nlat/2.,0        ,i])
      tprof[1,i]=mean(xy[5,nlat/2.-1:nlat/2.,nlon/2.  ,i])
      tprof[2,i]=mean(xy[5,nlat/2.-1:nlat/2.,nlon/4.  ,i])
      tprof[3,i]=mean(xy[5,nlat/2.-1:nlat/2.,nlon*3./4,i])
      tprof[4,i]=mean(xy[5,0     ,* ,i])
      tprof[5,i]=mean(xy[5,nlat-1,* ,i])
   endfor
endelse

;calculate sigma
if keyword_set(oom) then begin
  stp=-1.*oom/nlev
  sigma=fltarr(nlev)
  sigma[nlev-1]=10.^(stp/2.)
  for i=nlev-2,0,-1 do sigma[i]=sigma[i+1]*10.^(stp)
endif else begin
  sigma=findgen(nlev)/nlev+0.5/nlev
endelse

if keyword_set(psurf) then begin
   if not keyword_set(p0) then begin
      print,'MUST SET P0 IF USING PSURF!'
      stop
   endif
  openr, 18, 'fort.50'
  readf, 18, nlat,nlon
  ab=fltarr(3,nlat,nlon)
  readf, 18, ab
  close, 18
  sp=reform(ab[2,*,*])
  sp = (sp+1.)*p0
endif

if not keyword_set(xmin) then begin
   xmin=min(tprof)
   if keyword_set(eqtr) then begin
      txmin=min(xy[5,nlat/2.,*,*])
      if txmin lt xmin then xmin=txmin
   endif
   if keyword_set(otprof) then begin
      txmin=min(otprof)
      if txmin lt xmin then xmin=txmin
   endif
endif
if not keyword_set(xmax) then begin
   xmax=max(tprof)
   if keyword_set(eqtr) then begin
      txmax=max(xy[5,nlat/2.,*,*])
      if txmax gt xmax then xmax=txmax
   endif
   if keyword_set(otprof) then begin
      txmax=max(otprof)
      if txmax gt xmax then xmax=txmax
   endif
endif

if keyword_set(oom) then yrange = [1., 10.^(-oom)] else yrange = [1., 0]
if keyword_set(p0) then begin
   vert=sigma*p0
   verttit='Pressure [bar]'
endif else begin
   vert = sigma
   verttit = 'sigma'
endelse
if keyword_set(p0) then yrange*=p0
if keyword_set(psurf) then yrange=[max(sp),min(sp)*min(sigma)]

if keyword_set(psfile) then begin
   if (size(psfile))[1] eq 2 then filename='plot.eps' else filename=psfile
   psopen,filename,/enc,/color,xsize=25,ysize=20;,bits_per_pixel=24
   !p.font=0
   !p.charsize=2
   !p.thick=4
   !y.thick=4
   !x.thick=4
   xmrg=[10,3]
   cthick=12
endif else begin
  device,true_color=24,decomposed=0,retain=2
  xmrg=[12,3]
  cthick=8
endelse   

if not keyword_set(noplot) then begin
   if not keyword_set(psurf) then begin
      if keyword_set(oom) then $
         plot,tprof[0,*],vert,xr=[xmin,xmax],yrange = yrange,/ylog,/ystyle, $
              xtitle='Temperature [K]',ytitle=verttit,xmargin=xmrg,xlog=xlog,/xstyle
      oplot,tprof[1,*],vert,linestyle=4
      oplot,tprof[2,*],vert,linestyle=5
      oplot,tprof[3,*],vert,linestyle=3
      oplot,tprof[4,*],vert,linestyle=1
      oplot,tprof[5,*],vert,linestyle=2
   endif else begin
      if keyword_set(oom) then $
         plot,tprof[0,*],sigma*sp[loci[0,0],loci[1,0]],xr=[xmin,xmax],yrange = yrange,/ylog,/ystyle, $
              xtitle='Temperature [K]',ytitle=verttit,xmargin=xmrg,xlog=xlog,/xstyle
      oplot,tprof[1,*],sigma*sp[loci[0,1],loci[1,1]],linestyle=4
      oplot,tprof[2,*],sigma*sp[loci[0,2],loci[1,2]],linestyle=5
      oplot,tprof[3,*],sigma*sp[loci[0,3],loci[1,3]],linestyle=3
      oplot,tprof[4,*],sigma*sp[loci[0,4],loci[1,4]],linestyle=1
      oplot,tprof[5,*],sigma*sp[loci[0,5],loci[1,5]],linestyle=2
   endelse      
   
   if keyword_set(eqtr) then begin
      loadct,25,/silent
      ccolors=indgen(nlon)/nlon*255.
      if not keyword_set(psurf) then begin
         for i=0,nlon-1,lskip do oplot,xy[5,nlat/2.,i,*],vert,color=ccolors[i]
      endif else begin
         for i=0,nlon-1,lskip do oplot,xy[5,nlat/2.,i,*],sigma*sp[nlat/2.,i],color=ccolors[i]
      endelse
      if keyword_set(unstable) then begin
            akap=0.286
            print,'Using AKAP=',akap
            for i=0,nlon-1,lskip do begin
               if not keyword_set(psurf) then test=deriv(alog(sigma),alog(reform(xy[5,nlat/2.,i,*]))) $
               else test=deriv(alog(sigma*sp[nlat/2.,i]),alog(reform(xy[5,nlat/2.,i,*])))
               convi=where(test ge akap)
               if convi[0] ne -1 then begin
                  ;find continuous convective zone(s) and plot
                  nconv=n_elements(convi)
                  ic=0
                  eloop:
                  cstart=convi[ic]
                  while (ic lt nconv-1) do begin
                     if (convi[ic+1] eq convi[ic]+1) then ic+=1 else break
                  endwhile
                  cend=convi[ic]
                  oplot,xy[5,nlat/2.,i,cstart:cend],vert[cstart:cend],color=ccolors[i],thick=cthick
                  if cend ne max(convi) then begin
                     ic+=1
                     goto, eloop
                  endif
               endif
            endfor
      endif
      colorbar,position=[0.25,0.175,0.45,0.2],range=[0,360],divisions=4,/top,charsize=1.25,$
               title='Longitude (along equator)'
      loadct,0,/silent
   endif else begin
      if keyword_set(lats) then begin
         loadct,25,/silent
         ccolors=indgen(nlat)/nlat*255.
         if not keyword_set(psurf) then begin
            for i=0,nlat-1,lskip do oplot,xy[5,i,0,*],vert,color=ccolors[nlat-1-i]
         endif else begin
            for i=0,nlat-1,lskip do oplot,xy[5,i,0,*],sigma*sp[nlat/2.,i],$
                                    color=ccolors[nlat-1-i]
         endelse
         if keyword_set(unstable) then begin
            akap=0.286
            print,'Using AKAP=',akap
            for i=0,nlat-1,lskip do begin
               if not keyword_set(psurf) then test=deriv(alog(sigma),alog(reform(xy[5,i,0,*]))) $
               else test=deriv(alog(sigma*sp[i,0]),alog(reform(xy[5,i,0,*])))
               convi=where(test ge akap)
               if convi[0] ne -1 then begin
                  ;find continuous convective zone(s) and plot
                  nconv=n_elements(convi)
                  ic=0
                  cloop:
                  cstart=convi[ic]
                  while (ic lt nconv-1) do begin
                     if (convi[ic+1] eq convi[ic]+1) then ic+=1 else break
                  endwhile
                  cend=convi[ic]
                  oplot,xy[5,i,0,cstart:cend],vert[cstart:cend],color=ccolors[nlat-1-i],thick=cthick
                  if cend ne max(convi) then begin
                     ic+=1
                     goto, cloop
                  endif
               endif
            endfor
         endif
         colorbar,position=[0.25,0.175,0.45,0.2],range=[-90,90.],divisions=4,/top,$
                  charsize=1.25,title='Latitude (at longitude=0)'
         loadct,0,/silent
      endif
   endelse

   if keyword_set(allprof) then begin
      if not keyword_set(psurf) then begin
         for ilat=0,nlat-1 do begin
            for i=0,nlon-1 do begin
               oplot,xy[5,ilat,i,*],vert
            endfor
         endfor
      endif else begin
         for ilat=0,nlat-1 do begin
            for i=0,nlon-1 do begin
               oplot,xy[5,ilat,i,*],sigma*sp[nlat/2.,i]
            endfor
         endfor
      endelse
   endif

   if keyword_set(otprof) then begin
      sz=size(otprof)
;      loadct,5,/silent
      if (sz[0] eq 2) and (sz[2] eq nlev) then begin
         for i=0,sz[1]-1 do oplot,otprof[i,*],vert,thick=8;,color=100
      endif else begin
         print,'otprof unacceptable size'
         stop
      endelse
      loadct,0,/silent
   endif

   if keyword_set(oxyouts) then begin
      xyouts,xmax-(xmax-xmin)*0.25,yrange[1]*3.,oxyouts
   endif
   
;   if keyword_set(unstable) then begin
;      for i = 0, 5 do oplot, convec[i, *], vert, psym = 2
      if keyword_set(pot) then begin
         pot=fltarr(6,nlev)
         for i=0,5 do pot[i,*]=tprof[i,*]*(1./vert)^(akap)
         pmin=min(pot,max=pmax)
         window,2
         plot,pot[0,*],vert,xr=[pmin,pmax],yrange=yrange,/ylog,xtitle='Potential temperature [K]',ytit=verttit,$
              charsize=1.5,/xlog,/ystyle
         for i=1,5 do oplot,pot[i,*],vert,linestyle=i
         if keyword_set(eqtr) then begin
            loadct,25,/silent
            for i=0,nlon-1 do oplot,xy[5,nlat/2.,i,*]*(1./vert)^(akap),vert,color=ccolors[i]         
            loadct,0,/silent
         endif
      endif
;   endif
   if keyword_set(psfile) then begin
      !p.font=-1
      !p.charsize=2
      !p.thick=1
      !y.thick=1
      !x.thick=1
      psclose
   endif   
endif

;stop

end
