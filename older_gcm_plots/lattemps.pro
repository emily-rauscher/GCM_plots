pro lattemps,nl,lats=lats,ts=ts

  ; reads from fort.26 at level=nl to get average T profile from pole to pole,
  ;   including min/max variations from average value
  ; outputs:
  ;   lats: [nlat] array holding unique latitude values from fort.26
  ;   ts: [3,nlat] array holding average [0] and min [1], max [2] values
  ;       of temperature, as function of latitude
  
   OpenR, 17, 'fort.26'
   ReadF, 17, nlat,nlon,nlev	  
   xy=fltarr(6,nlat,nlon,nlev)   
   ReadF, 17, xy
   Close, 17

   temp=reform(xy[5,*,*,nl])
   lats=reform(xy[1,*,0,nl])

   ts=fltarr(3,nlat)
   for ilat=0,nlat-1 do begin
      ts[0,ilat]=mean(temp[ilat,*])
      ts[1,ilat]=min(temp[ilat,*],max=toss)
      ts[2,ilat]=toss
   endfor

   xmin=min(ts,max=xmax)
   plot,ts[0,*],lats,xr=[xmin,xmax],xtitle='Temperature [K]',ytit='Latitude [degrees]'
   oplot,ts[1,*],lats,linestyle=2
   oplot,ts[2,*],lats,linestyle=2
   
end
