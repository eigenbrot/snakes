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


function bc_continuum_alldust, model, restwl, flux, err, vdidx, emmaskw=emmaskw, $
                               yfit = yfit, plotlabel = plotlabel, fitregion = fitregion, $
                               savestep=savestep, lun=lun, lightidx=lightidx, fmt=fmt,$
                               velstart=velstart, chivec=chivec

light_factor = 100.
vel_factor = 100.
if n_elements(savestep) eq 0 then savestep = 0
if n_elements(velstart) eq 0 then velstart = 528.

; width of emission line masks in km/s
if not keyword_set(emmaskw) then emmaskw = 1000.0

dims = size(model.flux, /dimensions)
npix = n_elements(restwl)
nmodels = dims[2]

;-----------------------------------------------------------------------------
; Constrain model fit coefs to be positive, let TauV be pos or neg
; (fitcoefs = [TauV, model_coefs[*]])

parinfo = replicate({value:1.D, fixed:0, limited:[0,0], tied:'', $
                    limits:[0.0,0]}, nmodels * 2 + 1)

parinfo[0].limited = [1,1]
parinfo[0].limits = [(velstart - 100.)/vel_factor, (velstart + 100.)/vel_factor]
parinfo[0].value = velstart/vel_factor
parinfo[1:nmodels].limited = [1,0]
parinfo[1:nmodels].limits = [0,20]
parinfo[nmodels+1:*].limited = [1,0]
parinfo[nmodels+1:*].limits = [0,0]

;-----------------------------------------------------------------------------
; Mask out bad data regions 

quality = fltarr(npix) + 1

outside_model = where(restwl le min(model.wave) or restwl ge max(model.wave))
if outside_model[0] ne -1 then quality[outside_model] = 0

bad = where(finite(flux) ne 1 or finite(err) ne 1 or err eq 0)
if bad[0] ne -1 then quality[bad] = 0

;------------------------------------------------------------------------------
; Mask out emission lines

sk =    [6300.,        5890., 5683.8, 5577.,      5461., 5199.,      4983., 4827.32, 4665.69, 4420.23, 4358., 4165.68, 4047.0]
sknam = ['[OI] (atm)', 'NaD', 'NaI',  'OI (atm)', 'HgI', 'NI (atm)', 'NaI', 'HgI',   'NaI',   'NaI',   'HgI', 'NaI',   'HgI']

sk2 = [6300., 5890., 5577.]

;; em=     [3727.3,  4959.,    5006.8,   6563.8, 6716.0]
;; emnam = ['[OII]', '[OIII]', '[OIII]', 'Ha',   'S2']

em = [6563.8, 4861.,    4959., 5006.8, 6716.0, 6583.41, 6548.04]
emnam = ['Ha', 'Hb', '[OIII]','[OIII]',   'S2',   'NII',  'NII']

;; em = [6563.8, 6716.0, 6583.41, 6548.04]
;; emnam = ['Ha', 'S2', 'NII', 'NII']

abs =    [3933.7, 3968.5, 4304.4,   5175.3, 5894.0, 4861., 4341., 4102.]
absnam = ['H',    'K',    'G band', 'Mg',   'Na',   'HB',  'HG',  'HD']
;sk = [5569., 5882.6]
HPS = 5914.
HPS_wid = 230.

;shift emission to velocity guess. The mask width allows for slop
em *= (velstart/3e5 + 1)

dz = emmaskw / 3e5 ; clipping interval
dzsk = 1500. / 3e5

for ii = 0, n_elements(em) - 1 do begin 
  maskout = where(restwl gt em[ii]*(1-dz) and restwl lt em[ii]*(1+dz))
  if maskout[0] ne -1 then quality[maskout] = 0
endfor

for ii = 0, n_elements(sk2) - 1 do begin 
  maskout = where(restwl gt sk2[ii]*(1-dzsk) and restwl lt sk2[ii]*(1+dzsk))
  if maskout[0] ne -1 then quality[maskout] = 0
endfor

ok = where(quality eq 1)

;-----------------------------------------------------------------------------
; Interpolate to match data
  
custom_lib = dblarr(npix, nmodels)
for ii = 0, nmodels - 1 do begin
   print, size(model.flux, /dimensions)
   custom_lib[*,ii] = interpol(model.flux[vdidx,*,ii], model.wave, restwl)
endfor
if outside_model[0] ne -1 then custom_lib[outside_model, *] = 0.0

;-------------------------------------------------------------------------------
if keyword_set(fitregion) then begin
   fitidx = where(restwl[ok] ge fitregion[0] and restwl[ok] le fitregion[1])
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

if savestep eq 1 then $
   savestep = {flux: fitflux,$
               err: fiterr, $
               agearr: model.age/1e9, $
               norm: model.norm, $
               lightidx: lightidx, $
               lun: lun, $
               wave: fitwave, $
               fmt: fmt}

fitcoefs = mpfitfun('bc_mcombine_alldust', fitwave, fitflux, fiterr, $
                    parinfo = parinfo, $
                    functargs = {mlib: fitlib, savestep: savestep}, $
                    perror=perror, niter=niter, status=status, $
                    maxiter = 50000, /NAN)

print, 'CONTINUUM FIT ITERATIONS: ', strtrim(niter, 2)
print, 'CONTINUUM_FIT EXIT STATUS: ', strtrim(status, 2)

; structure containing fit coefs
coefs = {vsys: fitcoefs[0]*vel_factor, vsys_error: perror[0]*vel_factor,$
         tauv: fitcoefs[1:nmodels], tauv_err: perror[1:nmodels], $
         light_frac: fitcoefs[nmodels+1:*]*light_factor, light_frac_err: perror[nmodels+1:*]*light_factor, $
         model_age: reform(model.age[vdidx,*]), chisq: 0.0D, $
         redchi: 0.0D, bluechi: 0.0D, hkchi: 0.0D, $
         MMWA: 0.0D, MLWA: 0.0D, MMWT: 0.0D, MLWT: 0.0D, SNR: 0.0D}

; fit to full spectrum including masked pixels
redidx = where(restwl ge 5400)
blueidx = where(restwl lt 5400)
hklow = 3920
hkhigh = 4000
hkidx = where(restwl gt hklow and restwl lt hkhigh)

yfit = bc_mcombine_alldust(restwl, fitcoefs, mlib=custom_lib)

coefs.chisq = total((yfit - flux)^2/err^2)/(n_elements(flux) - n_elements(fitcoefs) - 1)
coefs.redchi = total((yfit[redidx] - flux[redidx])^2/err[redidx]^2)/$
               (n_elements(redidx) - n_elements(fitcoefs) - 1)
coefs.bluechi = total((yfit[blueidx] - flux[blueidx])^2/err[blueidx]^2)/$
                (n_elements(blueidx) - n_elements(fitcoefs) - 1)
coefs.hkchi = total((yfit[hkidx] - flux[hkidx])^2/err[hkidx]^2)/$
              (n_elements(hkidx) - n_elements(fitcoefs) - 1)
;redchi = total((yfit - flux)^2/err^2)/(n_elements(flux) + n_elements(fitcoefs) - 1)


coefs.SNR = mean(flux[lightidx]/err[lightidx])

coefs.MMWA = total(reform(model.age[vdidx,*])/1e9*coefs.light_frac/reform(model.norm[vdidx,*])) $
          / total(coefs.light_frac/reform(model.norm[vdidx,*]))

redd = exp(-coefs.tauv(restwl[lightidx]/5500)^(-0.7))
light_weight = mean(custom_lib[lightidx,*]*rebin(redd,n_elements(lightidx),$
                                                 n_elements(reform(model.age[vdidx,*]))),$
                    dimension=1) * coefs.light_frac

coefs.MLWA = total(light_weight * reform(model.age[vdidx,*])/1e9) / total(light_weight)


;; coefs.MMWZ = total(model.Z[vdidx,*]*coefs.light_frac/model.norm[vdidx,*]) $
;;           / total(coefs.light_frac/model.norm[vdidx,*])
;; coefs.MLWZ = total(light_weight * model.Z[vdidx,*]) / total(light_weight)

coefs.MMWT = total(coefs.tauv*coefs.light_frac/reform(model.norm[vdidx,*])) $
             / total(coefs.light_frac/reform(model.norm[vdidx,*]))
coefs.MLWT = total(light_weight * coefs.tauv) / total(light_weight)


;---------------------------------------------------------------------------
; Plot spectrum, best fit, and individual stellar components

defplotcolors
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
         exp(-coefs.tauv[i]*(restwl/5500.0)^(-0.7))
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

if status eq 5 then $
    xyouts, 0.20, 0.9, "MAXITER", charsize = 2, color = !red, /norm, /t3d
if status ge 6 then $
    xyouts, 0.20, 0.9, "FAILED", charsize = 2, color = !red, /norm, /t3d
 
plot, restwl, smooth((galfit - yfit)/err,smoothkern,/NAN), xtitle='Wavelength (Angstroms)', $
      ytitle='Residuals/error', $
      position=[0.15,0.15,0.95,0.3], xrange=[xmin,xmax], yrange=[-5,5], $
      yminor=1, yticks=4, charsize=1.0, charthick=1.0, thick=thick, $
      /xs, /ys, /noerase, /t3d
;, ytickv=[-200,0,200]
if keyword_set(plotlabel) then $
   xyouts, 0.2, 0.955, plotlabel, color = !black, /norm, /t3d

;; xyouts, 0.2, 0.90, 'Tau V = ' + string(fitcoefs[0], format = '(F5.2)'), $
;;         /norm, /t3d
xyouts, 0.2, 0.88, 'V_disp = ' + string(vdidx, format = '(F8.2)'), $
        /norm, /t3d
xyouts, 0.2, 0.86, 'MLWA = ' + string(coefs.MLWA, format = '(F8.2)'), $
        /norm, /t3d
xyouts, 0.2, 0.84, 'MMWA = ' + string(coefs.MMWA, format = '(F8.2)'), $
        /norm, /t3d
xyouts, 0.2, 0.82, 'SNR = ' + string(coefs.SNR, format = '(F8.2)'), $
        /norm, /t3d
xyouts, 0.2, 0.80, 'VSYS = ' + string(coefs.vsys, format = '(F8.2)'), $
        /norm, /t3d

xyouts, 0.4, 0.90, 'Chisq = ' + string(coefs.chisq, format = '(F8.2)'), $
        /norm, /t3d
xyouts, 0.4, 0.88, 'redChi = ' + string(coefs.redchi, format = '(F8.2)'), $
        /norm, /t3d
xyouts, 0.4, 0.86, 'bluChi = ' + string(coefs.bluechi, format = '(F8.2)'), $
        /norm, /t3d
xyouts, 0.4, 0.84, 'HK_Chi = ' + string(coefs.hkchi, format = '(F8.2)'), $
        /norm, /t3d
xyouts, 0.4, 0.82, 'MMWT = ' + string(coefs.MMWT, format = '(F8.2)'), $
        /norm, /t3d
xyouts, 0.4, 0.80, 'MLWT = ' + string(coefs.MLWT, format = '(F8.2)'), $
        /norm, /t3d

xyouts, 0.85, 0.56, '*', /norm, /t3d
for i = 0, n_elements(coefs.light_frac) - 1 do begin
   xyouts, 0.8, 0.55 - i*0.02, string(model.age[vdidx,i]/1e9, ': ', light_weight[i], ', ', coefs.tauv[i],$
                                      format='(F6.3,A2,F10.3,A2,F10.3)'), $
           charsize=0.6, alignment=0.0, /norm, /t3d
endfor

chivec = (flux - yfit)/err

print, fitcoefs, format = '(15F6.1)'

return, coefs

end
