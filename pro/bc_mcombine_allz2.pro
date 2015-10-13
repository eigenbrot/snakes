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

function bc_mcombine_allZ2, x, a, mlib=mlib, savedata=savedata

light_factor = 100.

y = mlib # (a[2:*] * light_factor)

; Redshift
xred = x * (a[0]*100. / 3e5 + 1)

; Redden using the Charlot & Fall law 
; F_obs = F_int * exp(-Tau_V * (lambda / 5500 A)^-0.7)

klam = (xred / 5500.0)^(-0.7)
e_tau_lam = exp(-1. * klam * a[1])

; Create a linear combination of the templates
y = y * e_tau_lam

y = interpol(y,xred,x)

if n_elements(savedata) eq 0 then savedata = 0

if tag_exist(savedata, 'flux',/quiet) then begin
   redidx = where(x ge 5250)
   blueidx = where(x lt 5250)
   hklow = 3920
   hkhigh = 4000
   hkidx = where(x gt hklow and x lt hkhigh)
   chisq = total((y - savedata.flux)^2/savedata.err^2)/(n_elements(y) - n_elements(a) - 1)
   
   redchi = total((y[redidx] - savedata.flux[redidx])^2/savedata.err[redidx]^2)
   bluechi = total((y[blueidx] - savedata.flux[blueidx])^2/savedata.err[blueidx]^2)
   hkchi = total((y[hkidx] - savedata.flux[hkidx])^2/savedata.err[hkidx]^2)

   SNR = mean(savedata.flux[savedata.lightidx]/savedata.err[savedata.lightidx])
   MMWA = total(savedata.agearr*a[1:*]*light_factor/savedata.norm) $
          / total(a[1:*]*light_factor/savedata.norm)

   redd = exp(-a[0]*(savedata.wave[savedata.lightidx]/5500)^(-0.7))
   light_weight = mean(savedata.flux[savedata.lightidx,*] * $
                       rebin(redd,n_elements(savedata.lightidx),$
                             n_elements(savedata.agearr)), $
                       dimension=1) * a[1:*] * light_factor
   MLWA = total(light_weight * savedata.agearr) / total(light_weight)

   MMWZ = total(savedata.Z * a[1:*] * light_factor / savedata.norm) $
          / total(a[1:*] * light_factor / savedata.norm)

   MLWZ = total(light_weight * savedata.Z) / total(light_weight)

   printf, savedata.lun, -1, a[1:*]/savedata.norm, MMWA, MLWA, $
           MMWZ, MLWZ, a[0], SNR, chisq, redchi, bluechi, hkchi, format=savedata.fmt
endif

return, y

end
