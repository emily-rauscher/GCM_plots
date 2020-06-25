pro level_temp_ave,runnum,days=days

  ;we'll assume L15
nl=15

  files=strcompress('../Tests/runs/fortR00'+string(runnum)+'.26.*',/remove)
  ff=findfile(files)  ;presumably these are 091,096,101,106,110

temps=fltarr(nl,5)

for i=0,4 do begin

  openr,17,ff[i]
  readf,17,nlat,nlon,nlev
  xy=fltarr(6,nlat,nlon,nlev)
  readf,17,xy
  close,17

  for j=0,nl-1 do begin
    temps[j,i]=mean(xy[5,*,*,j])
  endfor

endfor

if keyword_set(days) then xx=days else xx=[91,96,101,106,110]
ymax=min(temps)
ymin=max(temps)  ;plot from high to low temp, b/c follows height

plot,xx,temps[0,*],yr=[ymin,ymax],xtitle='day',ytitle='T (K)',psym=-4
for i=1,nl-1 do oplot,xx,temps[i,*],psym=-4

end
