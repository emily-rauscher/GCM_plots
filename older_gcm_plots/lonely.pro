function lonely,array

nel=n_elements(array)

pos=where(array gt 0.)
neg=where(array lt 0.)

parr=intarr(nel)
narr=intarr(nel)

parr[pos]=1.
narr[neg]=1.

parr=[0,parr,0]
narr=[0,narr,0]

parr=(parr*shift(parr,1)+parr*shift(parr,-1)) ;should >0 where not lonely
narr=(narr*shift(narr,1)+narr*shift(narr,-1)) ;should >0 where not lonely

farr=parr+narr
farr=farr[1:nel]

lon=where((farr eq 0) and (array ne 0))

farr=intarr(nel)
if lon[0] ne -1 then farr[lon]=1.

return,farr

end
