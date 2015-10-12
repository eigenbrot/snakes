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

function bc_mcombine_sky, x, a, mlib=mlib, savestep=savestep

; Create a linear combination of the templates

y = mlib[*,1:*] # (a[3:*] * 1000.)

; Redden using the Charlot & Fall law 
; F_obs = F_int * exp(-Tau_V * (lambda / 5500 A)^-0.7)

xred = x*(a[0]*100./3e5 + 1)

klam = (xred / 5500.0)^(-0.7)
e_tau_lam = exp(-a[1] * klam)
y = y * e_tau_lam

y = interpol(y,xred,x)
;y = y + (mlib[*,0] * a[2] * 1000.)

if n_elements(savestep) eq 0 then savestep = 0

if tag_exist(savestep, 'flux',/quiet) then begin
   yfit = bc_mcombine(savestep.wave, a, mlib=savestep.mlib, savestep=0)

   blueidx = savestep.blueidx
   redidx = savestep.redidx
   hkidx = savestep.hkidx

   chisq = total((yfit - savestep.flux)^2/savestep.err^2)/(n_elements(yfit) + n_elements(a) - 1)
   
   redchi = total((yfit[redidx] - savestep.flux[redidx])^2/savestep.err[redidx]^2)/$
            (n_elements(redidx) + n_elements(a) - 1)
   bluechi = total((yfit[blueidx] - savestep.flux[blueidx])^2/savestep.err[blueidx]^2)/$
             (n_elements(blueidx) + n_elements(a) - 1)
   hkchi = total((yfit[hkidx] - savestep.flux[hkidx])^2/savestep.err[hkidx]^2)/$
           (n_elements(hkidx) + n_elements(a) - 1)

   SNR = -99
   MMWA = total(savestep.agearr*a[1:*]*1./savestep.norm) $
          / total(a[1:*]*1./savestep.norm)

   redd = exp(-a[0]*(savestep.wave[savestep.lightidx]/5500)^(-0.7))
   light_weight = mean(savestep.mlib[savestep.lightidx,*] * $
                       rebin(redd,n_elements(savestep.lightidx),$
                             n_elements(savestep.agearr)), $
                       dimension=1) * a[1:*]
   MLWA = total(light_weight * savestep.agearr) / total(light_weight)

   printf, savestep.lun, 99, a[1:*]/savestep.norm, MMWA, MLWA, a[0],$
           SNR, chisq, redchi, bluechi, hkchi, -99, format=savestep.fmt
endif

return, y

end
