pro plottemp,pkv,palpha=palpha,pkth=pkth

pres=100.*get_sigma(7,45)

ga=2.19e3  ;cgs
ptirr=1722.
if not keyword_set(pkth) then pkth=0.279
if not keyword_set(palpha) then palpha=0.196

;kv=1.e-6, 1.e-4, 1.e-2

;tint=0 (1.e-30), 500.

;f = 1 (substellar), 0.5 (dayside average), 0.25 (global average)
; /night

prange=[max(pres),min(pres)]

count=0

for i=0,n_elements(pkv)-1 do begin
   for j=0,n_elements(palpha)-1 do begin
      for k=0,n_elements(pkth)-1 do begin
;         window,count
         print,count

         temp=guillot(pres,1.,tint=1.e-30,tirr=ptirr,kth=pkth[k],alpha=palpha[j],g=ga,kv=pkv[i],f=1.)
         plot,temp,pres,/ylog,yr=prange,xr=[500.,3000.],xstyle=8
         temp=guillot(pres,1.,tint=1.e-30,tirr=ptirr,kth=pkth[k],alpha=palpha[j],g=ga,kv=pkv[i],f=0.5)
         oplot,temp,pres,linestyle=1
         temp=guillot(pres,1.,tint=1.e-30,tirr=ptirr,kth=pkth[k],alpha=palpha[j],g=ga,kv=pkv[i],f=0.25)
         oplot,temp,pres,linestyle=2
         temp=guillot(pres,1.,tint=1.e-30,tirr=ptirr,kth=pkth[k],alpha=palpha[j],g=ga,kv=pkv[i],/night)
         oplot,temp,pres,linestyle=3
         temp=guillot(pres,1.,tint=500.,tirr=ptirr,kth=pkth[k],alpha=palpha[j],g=ga,kv=pkv[i],f=1.)
         oplot,temp,pres
         temp=guillot(pres,1.,tint=500.,tirr=ptirr,kth=pkth[k],alpha=palpha[j],g=ga,kv=pkv[i],f=0.5)
         oplot,temp,pres,linestyle=1
         temp=guillot(pres,1.,tint=500.,tirr=ptirr,kth=pkth[k],alpha=palpha[j],g=ga,kv=pkv[i],f=0.25)
         oplot,temp,pres,linestyle=2
         temp=guillot(pres,1.,tint=500.,tirr=ptirr,kth=pkth[k],alpha=palpha[j],g=ga,kv=pkv[i],/night)
         oplot,temp,pres,linestyle=3

         gamma=(pkv[i]/pkth[k])*pres^(-palpha[j])
         axis,xaxis=1,/xlog,xr=[min(gamma),max(gamma)],/save,xtick_get=v
         oplot,gamma,pres,thick=3
         print,'gamma:',min(gamma),' to',max(gamma)

         print,'kv:',pkv[i]
         print,'alpha:',palpha[j]
         print,'kth:',pkth[k]

         count+=1
      endfor
   endfor
endfor


end
