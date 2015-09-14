;+
; NAME:
;   bc_continuum
;
; PURPOSE:
;   Fit data with a continuum model which is a linear combination of 
;   Bruzual & Charlot 2003 models resampled to the correct velocity
;   dispersion.  Reddening is treated as a free parameter using the
;   Charlot and Fall Law: F_obs = F_int * exp(-Tau_V * (wave/5500)^-0.7)
;   The program returns Tau_V and the light fraction of each model at
;   5500 A.  (Note that A_V = 1.086 * Tau_V.)
;
; CALLING SEQUENCE:
;    s = bc_continuum(model, restwl, flux, err, vdisp, emmaskw=, yfit=)
;
; INPUTS:
;    model: a structure containing the Bruzual & Charlot models
;           re-sampled to a logarithmic wavelength scale (needed for
;           easy convolution to correct velocity dispersion.)  The structure
;           format is {wave: [npix], flux: [npix, nage], age: [nage]}
;    restwl: rest wavelength of data
;    flux: array of data fluxes
;    err: array of flux errors      
;    vdisp: velocity dispersion of galaxy in km/s
;
; OPTIONAL INPUTS:
;    emmaskw: width in km/s of emission line mask.  Default is 400 km/s.
;  
; OUTPUTS:
;    s: structure containing tags, "tau_v", "tau_v_err", "light_frac", 
;       "light_frac_err", and "model_age".  The reddening is stored in
;       "tau_v".  The light fraction contributed by each model to the
;       final fit is measured at 5500 A.  The tag "model_age" contains
;       the ages of the  Bruzual & Charlot templates used.
;    y: best fit continuum model
;
; NOTES: the "light_frac" is not normalized to 1.  The coefficients add
;        up to equal the value of the unreddened data at 5500 A.  To 
;        normalize the light_frac simply do: 
;        norm = s.light_frac/total(s.light_frac)
;
;        Fitting is done by MPFIT
;        http://cow.physics.wisc.edu/%7Ecraigm/idl/idl.html
;
; REVISION HISTORY: 
;    2007-Sept-18 Written by Christy Tremonti (tremonti@as.arizona.edu)
;
;-
;------------------------------------------------------------------------------


function bc_continuum_tau, model, restwl, flux, err, vdisp, emmaskw=emmaskw, $
                           yfit = yfit, plotlabel = plotlabel, $
                           savestep=savestep, lun=lun, lightidx=lightidx, $
                           fmt=fmt

if n_elements(savestep) eq 0 then savestep = 0

; width of emission line masks in km/s
if not keyword_set(emmaskw) then emmaskw = 400.0

nmodels = n_elements(model.age)

;-----------------------------------------------------------------------------
; Constrain model fit coefs to be positive, let TauV be pos or neg
; (fitcoefs = [TauV, psi0, tau_sf, t_form])

parinfo = replicate({value:0.D, fixed:0, limited:[0,0], tied:'', $
                    limits:[0.0,0], step:0}, 4)

parinfo[0].limited = [0,1]
parinfo[0].limits = [-5.0,20.0]
parinfo[1].value = 10.0
parinfo[1].limited = [1,1]
parinfo[1].limits = [0,1e9]
parinfo[2].value = 0.2
parinfo[2].limited = [0,0]
parinfo[2].limits = [-12.0, 20.0]
parinfo[3].value = 10
parinfo[3].limited = [0,0]
parinfo[3].limits = [7, 9]

;-----------------------------------------------------------------------------
; Mask out bad data regions 

npix = n_elements(restwl)
quality = fltarr(npix) + 1

outside_model = where(restwl le min(model.wave) or restwl ge max(model.wave))
if outside_model[0] ne -1 then quality[outside_model] = 0

bad = where(finite(flux) ne 1 or finite(err) ne 1 or err eq 0)
if bad[0] ne -1 then quality[bad] = 0

;------------------------------------------------------------------------------
; Mask out emission lines

;     OII       OII     H8       NeIII     Hg        Hd      Hb      OIII 
em= [3726.03, 3728.82, 3889.05, 3869.06, 4101.73, 4340.46, 4861.33, 4959.91, $
;    OIII     He I     OI        NII      Ha       NII      SII      SII 
    5006.84, 5875.67, 6300.30, 6548.04, 6562.82, 6583.41, 6716.44, 6730.81]
; bad sky lines
sk = [5569., 5882.6]

dz = emmaskw / 3e5 ; clipping interval
dzsk = 1000. / 3e5

for ii = 0, n_elements(em) - 1 do begin 
  maskout = where(restwl gt em[ii]*(1-dz) and restwl lt em[ii]*(1+dz))
  if maskout[0] ne -1 then quality[maskout] = 0
endfor

for ii = 0, n_elements(sk) - 1 do begin 
  maskout = where(restwl gt sk[ii]*(1-dzsk) and restwl lt sk[ii]*(1+dzsk))
  if maskout[0] ne -1 then quality[maskout] = 0
endfor

ok = where(quality eq 1)

;-----------------------------------------------------------------------------
; Convolve models to velocity dispersion of data and interpolate to
; match data

bc03_pix = 70.0 ; size of 1 model pixel in km/s 
bc03_vdisp = 75.0 ; approximate velocity dispersion of BC03 models
 
;Deconvolve template instrumental resolution
if vdisp lt bc03_vdisp then vdisp_add = 0 $
else vdisp_add = sqrt(vdisp^2 - bc03_vdisp^2)  
sigma_pix = vdisp_add / bc03_pix
  
custom_lib = dblarr(npix, nmodels)
for ii = 0, nmodels - 1 do begin
  cflux = gconv(model.flux[*,ii], sigma_pix) ; convolve with gaussian
  custom_lib[*,ii] = interpol(cflux, model.wave, restwl)
endfor
if outside_model[0] ne -1 then custom_lib[outside_model, *] = 0.0

;-------------------------------------------------------------------------------

if savestep eq 1 then $
   savestep = {flux: flux[ok],$
               err: err[ok], $
               agearr: model.age/1e9, $
               norm: model.norm, $
               lightidx: lightidx, $
               lun: lun, $
               wave: restwl[ok], $
               fmt: fmt}

fitcoefs = mpfitfun('bc_tau_combine', restwl[ok], flux[ok], err[ok], $
                    parinfo = parinfo, $
                    functargs = {mlib: custom_lib[ok,*], ages: model.age/1e9, $
                                 savestep: savestep}, $
                    perror=perror, niter=niter, status=status, $
                    maxiter = 500, /NAN)

print, 'CONTINUUM FIT ITERATIONS: ', strtrim(niter, 2)
print, 'CONTINUUM_FIT EXIT STATUS: ', strtrim(status, 2)
print, fitcoefs
; fit to full spectrum including masked pixels
yfit = bc_tau_combine(restwl, fitcoefs, mlib=custom_lib, ages=model.age/1e9)

chisq = total((yfit - flux)^2/err^2)
redchi = total((yfit - flux)^2/err^2)/(n_elements(flux) + n_elements(fitcoefs) - 1)

psi0 = fitcoefs[1]
tau_sf = 1./fitcoefs[2]
t_form = fitcoefs[3]

logt = alog10([99,model.age/1e9,t_form])
tdiff = logt[1:*] - logt[0:-2]
borders = 10.^(logt[1:*] - tdiff/2.)
borders[0] = 1e-99
borders[-1] = t_form

while (borders[1:*] - borders[0:-2])[-1] lt 0 do begin
   useidx = where(borders[1:*] - borders[0:-2] ge 0)
   useidx = [useidx,useidx[-1]+1]
   borders = borders[useidx]
   borders[-1] = t_form
   custom_lib = custom_lib[*,useidx[0:-2]]
   nmodels -= 1
endwhile

mass = psi0 * tau_sf * (exp(borders[1:*]/tau_sf) - $
                                    exp(borders[0:-2]/tau_sf))
weighted_age = tau_sf * alog(0.5 * (exp(borders[1:*]/tau_sf) + $
                                    exp(borders[0:-2]/tau_sf)))
while n_elements(mass) lt 10 do begin
   mass = [mass, -99]
endwhile

while n_elements(weighted_age) lt 10 do begin
   weighted_age = [weighted_age, -99]
endwhile

; structure containing fit coefs
coefs = {tauv: fitcoefs[0], tauv_err: perror[0], $
         psi0: psi0, psi0_err: perror[1], $
         tau_sf: tau_sf, tau_sf_err: perror[2], $
         t_form: t_form, t_form_err: perror[3], $
         model_age: model.age, chisq: chisq, redchi: redchi, $
         mass: mass, weighted_age: weighted_age}

;---------------------------------------------------------------------------
; Plot spectrum, best fit, and individual stellar components
defplotcolors
ymax = max(yfit) * 1.1
xmin = min(restwl) * 0.98
xmax = max(restwl) * 1.02
;xtitle='Wavelength (Angstroms)'
plot, restwl, flux, xtickformat='(A1)', /nodata,$
      ytitle = 'Flux', yrange = [1, ymax], xrange = [xmin,xmax], $
      position = [0.15,0.3,0.95,0.99], charsize=1.0, charthick=1.0, /xs, /ys, /t3d

errsmooth = 5
oband, restwl[0:*:errsmooth], (flux-err)[0:*:errsmooth], (flux+err)[0:*:errsmooth], color=!gray

; Show masked regions in green
galfit = flux  + 'NaN'
galfit[ok] = flux[ok]
masked = flux
masked[ok] = 'NaN'
thick=1.0
oplot, restwl, masked, color=!green, thick=thick, /t3d
oplot, restwl, galfit, thick=thick, color=!black, /t3d

oplot, restwl, yfit, color = !red, thick=thick, /t3d

if status ge 5 then $
    xyouts, 0.20, 0.9, "FAILED", charsize = 2, color = !red, /norm, /t3d

if keyword_set(plotlabel) then $
   xyouts, 0.2, 0.93, plotlabel, color = !black, /norm, /t3d
 
lidx = where(restwl ge 5450 and restwl le 5550)
for i=0, nmodels - 1 do begin
   print, '^^^', mass[i]
   yi = mass[i] * custom_lib[*,i] * $
         exp(-coefs.tauv*(restwl/5500.0)^(-0.7))
   oplot, restwl, yi, color = !blue, thick=thick, linestyle=2, /t3d
   ;; xyouts, 0.2, 0.88 - 0.02*i, string('f_',i,' = ',$
   ;;                          mean(yi[lidx]),$
   ;;                          format='(A2,I02,A3,E10.3)'),$
   ;;         charsize=0.6, /norm, /t3d
   
endfor

plot, restwl, (galfit - yfit)/err, xtitle='Wavelength (Angstroms)', ytitle='Residuals/error', $
      position=[0.15,0.15,0.95,0.3], xrange=[xmin,xmax], yrange=[-5,5], $
      yminor=1, yticks=4, charsize=1.0, charthick=1.0, thick=thick, $
      /xs, /ys, /noerase, /t3d
;, ytickv=[-200,0,200]
xyouts, 0.2, 0.9, 'Tau V = ' + string(coefs.tauv, format = '(F10.2)'), $
        /norm, /t3d
xyouts, 0.2, 0.87, 'psi0 = ' + string(coefs.psi0, format = '(F10.2)'), $
        /norm, /t3d
xyouts, 0.2, 0.84, 'tau_sf = ' + string(coefs.tau_sf, format = '(F10.2)'), $
        /norm, /t3d
xyouts, 0.2, 0.81, 't_form = ' + string(coefs.t_form, format = '(F10.2)'), $
        /norm, /t3d

;; for i = 0, n_elements(coefs.light_frac) - 1 do begin
;;    xyouts, 0.8, 0.55 - i*0.02, string(model.age[i]/1e9, ': ', coefs.light_frac[i], $
;;                                       format='(F6.3,A2,F10.3)'), $
;;            charsize=0.6, alignment=0.0, /norm, /t3d
;; endfor

print, fitcoefs, format = '(4F6.1)'

return, coefs

end
