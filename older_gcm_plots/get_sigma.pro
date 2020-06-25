function get_sigma,oom,nl

sigma=fltarr(nl)
stp=-1.*oom/nl
sigma[nl-1]=10.^(stp/2.)
for i=nl-2,0,-1 do sigma[i]=sigma[i+1]*10.^(stp)

return,sigma

end
