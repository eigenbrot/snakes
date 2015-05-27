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

function bc_tau_combine, x, a, mlib=mlib, ages=ages, savestep=savestep

; Create a linear combination of the templates
tau_V = a[0]
psi0 = a[1]
tau_sf = a[2]
t_form = a[3]

logt = alog10([99,ages,t_form])
tdiff = logt[1:*] - logt[0:-2]
borders = 10.^(logt[1:*] - tdiff/2.)
borders[0] = 1e-99
borders[-1] = t_form

while (borders[1:*] - borders[0:-2])[-1] lt 0 do begin
   useidx = where(borders[1:*] - borders[0:-2] ge 0)
   useidx = [useidx,useidx[-1]+1]
   borders = borders[useidx]
   borders[-1] = t_form
   mlib = mlib[*,useidx[0:-2]]
endwhile

mass = tau_sf * (exp(borders[1:*]/tau_sf) - exp(borders[0:-2]/tau_sf))

weighted_age = tau_sf * alog(0.5 * (exp(borders[1:*]/tau_sf) + $
                                    exp(borders[0:-2]/tau_sf)))

psi = exp(weighted_age/tau_sf)
mass *= psi0
psi *= psi0

galaxy = mlib # mass

; Redden using the Charlot & Fall law 
; F_obs = F_int * exp(-Tau_V * (lambda / 5500 A)^-0.7)

klam = (x / 5500.0)^(-0.7)
e_tau_lam = exp(-tau_V * klam)
galaxy *= e_tau_lam

if n_elements(savestep) eq 0 then savestep = 0

if tag_exist(savestep, 'flux',/quiet) then begin
   chisq = total((y - savestep.flux)^2/savestep.err^2)
   
   SNR = -99
   MMWA = total(savestep.agearr*a[1:*]) $
          / total(a[1:*])

   redd = exp(-a[0]*(savestep.wave[savestep.lightidx]/5500)^(-0.7))
   light_weight = mean(savestep.flux[savestep.lightidx,*] * $
                       rebin(redd,n_elements(savestep.lightidx),$
                             n_elements(savestep.agearr)), $
                       dimension=1) * a[1:*]
   MLWA = total(light_weight * savestep.agearr) / total(light_weight)

   printf, savestep.lun, 2, mass, MMWA, MLWA, a[0],$
           SNR, chisq, -99, -99, format=savestep.fmt
endif

return, galaxy

end
