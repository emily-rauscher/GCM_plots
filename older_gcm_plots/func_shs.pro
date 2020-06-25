pro func_shs,X,A,F,pder

; X: longitudes or latitudes, IN RADIANS
; A: coefficients for fitting function
; F: fluxes (calculated)
; pder: partial derivatives

  F = A[2] * ( cos(X-A[0]) )^A[3] + A[1]

  if n_params() ge 4 then $
     pder = [ [A[2]*A[3]*(cos(X-A[0]))^(A[3]-1)*sin(X-A[0])], $
              [replicate(1.,n_elements(X))], $
              [(cos(X-A[0]))^A[3]], $
              [alog(cos(X-A[0]))*A[2]*(cos(X-A[0]))^A[3]] ]

end
