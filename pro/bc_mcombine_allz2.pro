;+
; NAME:
;   bc_mcombine
;
; PURPOSE:
;    Define a composite spectrum as a (reddened) linear combination 
;    of instantaneous burst models.  Reddening of the composite  
;    spectrum is controlled by the first parameter of the array a.  The 
;    reddening is performed using the Charlot & Fall law.   
;
; CALLING SEQUENCE:
;    y = bc_mcombine(x, a, mlib=mlib)
;
; INPUTS:
;    x:  array containing wavelength [npix]
;    a:  array of coefficients defining the model composite [N + 1]
;          a(0):      Tau_V   (must be > 0)
;          a(1...N):  linear coefficients for each model (must be >=0)
;    mlib: a 2-d array containing the model fluxes [npix, N]
;  
; OUTPUTS:
;    y:  returned composite template
; 
; REVISION HISTORY: 
;    2007-Sept-18 Written by Christy Tremonti (tremonti@as.arizona.edu)
;
;-
;------------------------------------------------------------------------------

function bc_mcombine_allZ2, x, a, mlib=mlib;, savedata=savedata

light_factor = 100.
vel_factor = 100.

y = mlib # (a[2:*] * light_factor)

; Redshift
xred = x * (a[0]*vel_factor / 3e5 + 1)

; Redden using the Charlot & Fall law 
; F_obs = F_int * exp(-Tau_V * (lambda / 5500 A)^-0.7)

klam = (xred / 5500.0)^(-0.7)
e_tau_lam = exp(-1. * klam * a[1])

; Create a linear combination of the templates
y = y * e_tau_lam

y = interpol(y,xred,x)

return, y

end
