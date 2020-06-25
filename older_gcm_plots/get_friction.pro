function get_friction,k,sb,oom,nl

s=get_sigma(oom,30)
;print,s

fric=fltarr(nl)
sf=fltarr(nl)
for i=0,nl-1 do begin
    sf0=(s[i]-sb)/(1.-sb)
    sf[i]=MAX([0,sf0])
;    print, i,s[i],sb,sf0,sf[i]
endfor

;print, sf
return,sf

end 

