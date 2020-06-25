pro KE_diag, p0, nl, ndays, oom=oom, dayrange=dayrange,psfile=psfile,tdrag=tdrag,cpt=cpt,trim=trim

; plots the kinetic energy (winds) of each pressure level as a
;    function of time
; reads from fort.52
; requires: psopen.pro, psclose.pro, colorbar.pro
;
; p0: surface pressure in bar
; nl: number of vertical levels
; ndays: number of days in simulation
; oom: vertical pressure extent in orders of magnitude, if not set
;      assumes linear pressure levels
; dayrange=[begday,endday] to plot
; psfile: outputs an eps file instead of plotting to the window
; tdrag: see description below, maybe not fully implemented
; cpt: overplots the total enthalpy as a function of time
; trim: doesn't plot days 0:2, during which there are often no winds
;
; tdrag (in seconds) an array of tfrc[nl], if set will
; plot KE diss as a function of time in run ... not sure this has been tested

if keyword_set(trim) then dayrange=[3,ndays]

if not keyword_set(psfile) then begin
; Screen output
  device,true_color=24,decomposed=0,retain=2
;  set_plot, 'x'
;  window, 0
  !p.font=-1
  !p.charsize=1.5
  cbarc=255
  pxmarg=[12,12]
endif else begin
   if (size(psfile))[1] eq 2 then filename='plot.eps' else filename=psfile
   psopen,filename,/enc,/color,xsize=25,ysize=20;,bits_per_pixel=24
  !p.font=0
  !p.charsize=1.5
  cbarc=0
  pxmarg=[12,8]
endelse
  !x.style=1


  loadct,4,/silent

  openr,17,'fort.52'
    xy=fltarr(2,nl,ndays+1)   ;0: kinetic, 1: kinetic + cpT
  readf,17,xy
  close,17

  xy*=(p0*1.e5)

  days=findgen(ndays+1)
  if keyword_set(dayrange) then days=days[dayrange[0]:dayrange[1]]
  ke=reform(xy[0,*,days])
  ke=transpose(ke)
  te=reform(xy[1,*,days])
  te=transpose(te)

  if keyword_set(tdrag) then begin
     if n_elements(tdrag) ne nl then begin
        print,'tdrag must contain NL elements'
        stop
     endif
     if keyword_set(dayrange) then begin
        numdays=dayrange[1]-dayrange[0]+1
     endif else begin
        numdays=ndays+1
     endelse
     ked=fltarr(numdays)
     noz=where(tdrag ne 0.)
     nnoz=n_elements(noz)
     for i=0,nnoz-1 do begin
        ked=ked+ke[*,noz[i]]/tdrag[noz[i]]
     endfor
     nonz=where(ked ne 0)
     minked=min(ked[nonz],max=maxked)

     ysty=8
     print,'KE diss on last day:',ked[numdays-1],' W'
  endif else ysty=1

  if keyword_set(cpt) then begin
     ysty=9
     ymrgn=[4,4]
  endif else ymrgn=[4,2]

;  print,minmax(ke),minmax(te)
  nonz=where(ke ne 0)
  minlke=min(alog10((ke)[nonz]),max=maxlke)

  if keyword_set(oom) then begin
     levs=fltarr(nl)
     stp=-1.*oom/nl
     levs[nl-1]=10.^(stp/2.)
     for i=nl-2,0,-1 do levs[i]=levs[i+1]*10.^(stp)
     levs*=p0
     contour,alog10(ke),days,levs,/cell_fill,nlevels=55,yr=[max(levs),min(levs)],ymargin=ymrgn,$
             min_value=minlke,max_value=maxlke,ystyle=ysty,/ylog,/xstyle,xmargin=pxmarg,$
             xtitle='Planet Day',ytitle='Pressure [bar]' ;,title='log(KE)'
;     contour,alog10(te-ke),days,levs,/overplot,nlevels=10,color=0
     if keyword_set(tdrag) then begin
        loadct,0,/silent
        axis,yaxis=1,/ylog,yr=[minked,maxked],ytit='KE drag loss [W]',/save
        oplot,days,ked,thick=4
     endif
     if keyword_set(cpt) then begin
        colorbar,range=[minlke,maxlke],position=[0.15,0.925,0.875,0.95],format='(f4.1)',color=cbarc,/top,charsize=1.25
        loadct,0,/silent
        norme=total(xy[1,*,0]-xy[0,*,0])  ;c_p T in atmosphere for initial conditions
        plote=total(te-ke,2)
        mine=min(plote,max=maxe)
        axis,yaxis=1,yr=[mine,maxe],ytit='Total c_p T [J]',/save,/ystyle,ylog=0
        oplot,days,plote,thick=4
     endif else begin
;        colorbar,range=[minlke,maxlke],position=[0.935,0.15,0.96,0.85],$
        colorbar,range=[minlke,maxlke],position=[0.91,0.15,0.94,0.85],$
                 format='(f4.1)',color=cbarc,/vertical,/right,charsize=1.25
     endelse
  endif else begin
     levs=findgen(nl)+1
     levs*=p0
     contour,alog10(ke),days,levs,/cell_fill,nlevels=55,yr=[nl-1,0],/ystyle,$
             min_value=minlke,max_value=maxlke,$
             xtitle='Day',ytitle='pressure [bar]';,title='log(KE)'
;     contour,alog10(te),days,levs,/overplot,nlevels=10
  endelse

  loadct,0,/silent


if keyword_set(psfile) then begin
   psclose
   spawn,'gs -r300 -dEPSCrop -dTextAlphaBits=4 -sDEVICE=png16m -sOutputFile=plot.png -dBATCH -dNOPAUSE plot.eps'
   spawn,'convert plot.png eps3:plot.eps'
endif
!x.style=0

;stop

end
