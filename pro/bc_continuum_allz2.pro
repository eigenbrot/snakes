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


function bc_continuum_allZ2, model, restwl, flux, err, vdisp, emmaskw=emmaskw, $
                            yfit = yfit, plotlabel = plotlabel, bluefit = bluefit, $
                            savestep=savestep, lun=lun, lightidx=lightidx, fmt=fmt

print, vdisp

light_factor = 100.

if n_elements(savestep) eq 0 then savestep = 0

; width of emission line masks in km/s
if not keyword_set(emmaskw) then emmaskw = 400.0

dims = size(model.flux, /dimensions)
npix = n_elements(restwl)
nmodels = dims[1]

;-----------------------------------------------------------------------------
; Constrain model fit coefs to be positive, let TauV be pos or neg
; (fitcoefs = [TauV, model_coefs[*]])

parinfo = replicate({value:1.D, fixed:0, limited:[0,0], tied:'', $
                    limits:[0.0,0], step:0}, nmodels + 1)

parinfo[0].limited = [1,1]
parinfo[0].limits = [0,20]
parinfo[1:*].limited = [1,0]
parinfo[1:*].limits = [0,0]

;-----------------------------------------------------------------------------
; Mask out bad data regions 

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
if keyword_set(bluefit) then begin
   fitidx = where(restwl[ok] lt 5250)
   fitflux = flux[ok[fitidx]]
   fiterr = err[ok[fitidx]]
   fitwave = restwl[ok[fitidx]]
   fitlib = custom_lib[ok[fitidx],*]
endif else begin
   fitflux = flux[ok]
   fiterr = err[ok]
   fitwave = restwl[ok]
   fitlib = custom_lib[ok,*]
endelse

if savestep eq 1 then begin
   savedata = {flux: fitflux,$
               err: fiterr, $
               agearr: model.age/1e9, $
               Z: model.Z, $
               norm: model.norm, $
               lightidx: lightidx, $
               lun: lun, $
               wave: fitwave, $
               fmt: fmt}
endif else begin
   savedata = 0
endelse

fitcoefs = mpfitfun('bc_mcombine_allZ2', fitwave, fitflux, fiterr, $
                    parinfo = parinfo, $
                    functargs = {mlib: fitlib, savedata: savedata}, $
                    perror=perror, niter=niter, status=status, $
                    errmsg=errmsg, maxiter = 50000, xtol=1d-10, ftol=1d-10, /NAN)

print, 'CONTINUUM FIT ITERATIONS: ', strtrim(niter, 2)
print, 'CONTINUUM_FIT EXIT STATUS: ', strtrim(status, 2)
print, 'CONTINUUM_FIT ERRMSG: ', errmsg

; fit to full spectrum including masked pixels
redidx = where(restwl ge 5250)
blueidx = where(restwl lt 5250)
hklow = 3920
hkhigh = 4000
hkidx = where(restwl gt hklow and restwl lt hkhigh)

yfit = bc_mcombine_allZ2(restwl, fitcoefs, mlib=custom_lib)

bluefree = (n_elements(blueidx) - n_elements(fitcoefs) - 1)

chisq = total((yfit - flux)^2/err^2)/(n_elements(flux) - n_elements(fitcoefs) - 1)
redchi = total((yfit[redidx] - flux[redidx])^2/err[redidx]^2)/$
         (n_elements(redidx) - n_elements(fitcoefs) - 1)
bluechi = total((yfit[blueidx] - flux[blueidx])^2/err[blueidx]^2)/$
          (n_elements(blueidx) - n_elements(fitcoefs) - 1)
hkchi = total((yfit[hkidx] - flux[hkidx])^2/err[hkidx]^2)/$
        (n_elements(hkidx) - n_elements(fitcoefs) - 1)
;redchi = total((yfit - flux)^2/err^2)/(n_elements(flux) + n_elements(fitcoefs) - 1)


SNR = mean(flux[lightidx]/err[lightidx])

MMWA = total(model.age/1e9*fitcoefs[1:*]*light_factor/model.norm) $
          / total(fitcoefs[1:*]*light_factor/model.norm)

redd = exp(-fitcoefs[0]*(restwl[lightidx]/5500)^(-0.7))
light_weight = mean(model.flux[lightidx,*] * rebin(redd,n_elements(lightidx),$
                                                   n_elements(model.age)), $
                    dimension=1) * fitcoefs[1:*] * light_factor
MLWA = total(light_weight * model.age/1e9) / total(light_weight)

MMWZ = total(model.Z*fitcoefs[1:*]*light_factor/model.norm) $
          / total(fitcoefs[1:*]*light_factor/model.norm)
MLWZ = total(light_weight * model.Z) / total(light_weight)

; structure containing fit coefs
coefs = {tauv: fitcoefs[0], tauv_err: perror[0], $
         light_frac: fitcoefs[1:*]*light_factor, light_frac_err: perror[1:*]*light_factor, $
         model_age: model.age, chisq: chisq, $
         redchi: redchi, bluechi: bluechi, hkchi: hkchi, $
         bluefree: bluefree, MMWA: MMWA, MLWA: MLWA, MMWZ: MMWZ, MLWZ: MLWZ, SNR: SNR}

;---------------------------------------------------------------------------
; Plot spectrum, best fit, and individual stellar components

;defplotcolors
smoothkern = 5

blueymax = max(yfit[blueidx]) / 0.8
ymax = max(yfit) * 1.1
ymax = max([ymax,blueymax])
xmin = min(restwl) * 0.98
xmax = max(restwl) * 1.02
;xtitle='Wavelength (Angstroms)'
plot, restwl, flux, xtickformat='(A1)', /nodata,$
      ytitle = 'Flux', yrange = [1, ymax], xrange = [xmin,xmax], $
      position = [0.15,0.3,0.95,0.99], charsize=1.0, charthick=1.0, /xs, /ys, /t3d

vline, 5250., color=!gray, linestyle=2
vline, hklow, color=!gray, linestyle=2
vline, hkhigh, color=!gray, linestyle=2

for i=0, nmodels - 1 do begin
   yi = coefs.light_frac[i] * custom_lib[*,i] * $
         exp(-fitcoefs[0]*(restwl/5500.0)^(-0.7))
    oplot, restwl, smooth(yi,smoothkern), color = !lblue, thick=thick, linestyle=0, /t3d
    ;; xyouts, 0.2, 0.88 - 0.02*i, string('f_',i,' = ',$
    ;;                         mean(yi[lightidx]),$
    ;;                         format='(A2,I02,A3,E10.3)'),$
    ;;        charsize=0.6, /norm, /t3d
endfor

errsmooth = 5
oband, restwl[0:*:errsmooth], (flux-err)[0:*:errsmooth], (flux+err)[0:*:errsmooth], color=!gray

; Show masked regions in green
galfit = flux  + 'NaN'
galfit[ok] = flux[ok]
masked = flux
masked[ok] = 'NaN'
thick=1.0
oplot, restwl, smooth(galfit,smoothkern,/NAN), thick=thick, color=!black, /t3d
oplot, restwl, smooth(masked,smoothkern,/NAN), color=!green, thick=thick, /t3d

oplot, restwl, smooth(yfit,smoothkern), color = !red, thick=thick, /t3d

if status ge 5 then $
    xyouts, 0.20, 0.9, "FAILED", charsize = 2, color = !red, /norm, /t3d
 
plot, restwl, smooth((galfit - yfit)/err,smoothkern,/NAN), xtitle='Wavelength (Angstroms)', $
      ytitle='Residuals/error', $
      position=[0.15,0.15,0.95,0.3], xrange=[xmin,xmax], yrange=[-5,5], $
      yminor=1, yticks=4, charsize=1.0, charthick=1.0, thick=thick, $
      /xs, /ys, /noerase, /t3d
;, ytickv=[-200,0,200]
if keyword_set(plotlabel) then $
   xyouts, 0.2, 0.955, plotlabel, color = !black, /norm, /t3d

xyouts, 0.2, 0.90, 'Tau V = ' + string(fitcoefs[0], format = '(F5.2)'), $
        /norm, /t3d
xyouts, 0.2, 0.88, 'V_disp = ' + string(vdisp, format = '(F8.2)'), $
        /norm, /t3d
xyouts, 0.2, 0.86, 'MLWA = ' + string(MLWA, format = '(F8.2)'), $
        /norm, /t3d
xyouts, 0.2, 0.84, 'MMWA = ' + string(MMWA, format = '(F8.2)'), $
        /norm, /t3d
xyouts, 0.2, 0.82, 'SNR = ' + string(SNR, format = '(F8.2)'), $
        /norm, /t3d

xyouts, 0.4, 0.90, 'Chisq = ' + string(chisq, format = '(F8.2)'), $
        /norm, /t3d
xyouts, 0.4, 0.88, 'redChi = ' + string(redchi, format = '(F8.2)'), $
        /norm, /t3d
xyouts, 0.4, 0.86, 'bluChi = ' + string(bluechi, format = '(F8.2)'), $
        /norm, /t3d
xyouts, 0.4, 0.84, 'HK_Chi = ' + string(hkchi, format = '(F8.2)'), $
        /norm, /t3d
xyouts, 0.4, 0.82, 'MLWZ = ' + string(MLWZ, format = '(F8.2)') + '(Z/Z_sol)', $
        /norm, /t3d
xyouts, 0.4, 0.8, 'MMWZ = ' + string(MMWZ, format = '(F8.2)'), $
        /norm, /t3d

;; for i = 0, n_elements(coefs.light_frac) - 1 do begin
;;    xyouts, 0.8, 0.55 - i*0.02, string(model.age[i]/1e9, ': ', coefs.light_frac[i],$
;;                                       format='(F6.3,A2,F10.3)'), $
;;            charsize=0.6, alignment=0.0, /norm, /t3d
;; endfor

print, fitcoefs, format = '(15F6.1)'

return, coefs

end
