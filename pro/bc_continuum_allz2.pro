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

function bc_continuum_allZ2, model, restwl, flux, err, vdidx, emmaskw=emmaskw, $
                             yfit = yfit, fitregion = fitregion, velstart = velstart, $
                             savestep=savestep, lun=lun, lightidx=lightidx, fmt=fmt, $
                             chivec=chivec, fixedV=fixedV

light_factor = 100.
vel_factor = 100.

if n_elements(savestep) eq 0 then savestep = 0
if n_elements(velstart) eq 0 then velstart = 528.

; width of emission line masks in km/s
if n_elements(emmaskw) eq 0 then begin
   emmaskw = 1000.0
   print, 'EMMASKW NOT SET'
   stop
endif

dims = size(model.flux, /dimensions)
npix = n_elements(restwl)
nmodels = dims[2]

;-----------------------------------------------------------------------------
; Constrain model fit coefs to be positive, let TauV be pos or neg
; (fitcoefs = [TauV, model_coefs[*]])

parinfo = replicate({value:1.D, fixed:0, limited:[0,0], tied:'', $
                    limits:[0.0,0], step:0, relstep:0}, nmodels + 2)

if keyword_set(fixedV) then begin 
   parinfo[0].value = fixedV/vel_factor
   parinfo[0].fixed = 1
   fixedVBool = 1
endif else begin
   parinfo[0].limited = [1,1]
   parinfo[0].limits = [(velstart - 200.)/vel_factor, (velstart + 200.)/vel_factor]
   parinfo[0].fixed = 0
   parinfo[0].value = velstart/vel_factor
   fixedVBool = 0
endelse

parinfo[1].limited = [1,1]
parinfo[1].limits = [0.,20.]
parinfo[2:*].limited = [1,0]
parinfo[2:*].limits = [0.,0.]

;Put in some age/metallicity priors
;default: don't fit
;; parinfo[2:*].fixed = 1
;; parinfo[2:*].value = 0

;Let's use only the middle range of metallicity
;; ids = where(model.Z[0,*] gt 0.01 and model.Z[0,*] lt 0.03)
;; parinfo[ids+2].fixed = 0
;; parinfo[ids+2].value = 1

;And take out the 0.4 solar, leaving us with only two metallicities
;; id4 = where(model.Z[0,*] gt 0.1 and model.Z[0,*] lt 0.5)
;; parinfo[id4+2].fixed = 1
;; parinfo[id4+2].value = 0

; Only fit solar metallicity
;; ids = where(model.Z[0,*] ge 0.9 and model.Z[0,*] le 1.1)
;; parinfo[ids+2].fixed = 0
;; parinfo[ids+2].value = 1.D

;; ;
;; ; t > 8 Gyr: Z_sol = 0.02 - 0.4
;; id8 = where(model.age[0,*] gt 8d9 and model.Z[0,*] ge 0.02 and model.Z[0,*] le 0.41)
;; parinfo[id8+2].fixed = 0
;; parinfo[id8+2].value = 1.D

;; ; 3 < t < 8: Z_sol = 0.02 - 0.4
;; id3 = where(model.age[0,*] le 8d9 and model.age[0,*] gt 3d9 and model.Z[0,*] ge 0.02 and model.Z[0,*] le 0.41)
;; parinfo[id3+2].fixed = 0
;; parinfo[id3+2].value = 1.D

;; ; 1 < t < 3 Gyr; Z_sol = 0.2 - 1
;; id1 = where(model.age[0,*] le 3d9 and model.age[0,*] gt 1d9 and model.Z[0,*] ge 0.2 and model.Z[0,*] le 1.1)
;; parinfo[id1+2].fixed = 0
;; parinfo[id1+2].value = 1.D

;; ; t < 1 Gyr: Z_sol = 0.4 - 2.5
;; id0 = where(model.age[0,*] le 1d9 and model.Z[0,*] ge 0.4 and model.Z[0,*] le 2.51)
;; parinfo[id0+2].fixed = 0
;; parinfo[id0+2].value = 1.D

;-----------------------------------------------------------------------------
; Mask out bad data regions 

quality = fltarr(npix) + 1

outside_model = where(restwl le min(model.wave) or restwl ge max(model.wave))
if outside_model[0] ne -1 then quality[outside_model] = 0

bad = where(finite(flux) ne 1 or finite(err) ne 1 or err eq 0)
if bad[0] ne -1 then quality[bad] = 0


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
if keyword_set(fixedV) then begin 
   em *= (fixedV/3e5 + 1)
endif else begin
   em *= (velstart/3e5 + 1)
endelse

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

if savestep eq 1 then begin
   savedata = {flux: fitflux,$
               err: fiterr, $
               agearr: model.age[vdidx,*]/1e9, $
               Z: model.Z[vdidx,*], $
               norm: model.norm[vdidx,*], $
               lightidx: lightidx, $
               lun: lun, $
               wave: fitwave, $
               fmt: fmt}
endif else begin
   savedata = 0
endelse

fitcoefs = mpfitfun('bc_mcombine_allz2', fitwave, fitflux, fiterr, $
                    parinfo = parinfo, $
                    functargs = {mlib: fitlib, savedata: savedata}, $
                    perror=perror, niter=niter, status=status, $
                    errmsg=errmsg, maxiter = 50000, xtol=1d-10, ftol=1d-10, /NAN)

print, 'CONTINUUM FIT ITERATIONS: ', strtrim(niter, 2)
print, 'CONTINUUM_FIT EXIT STATUS: ', strtrim(status, 2)
print, 'CONTINUUM_FIT ERRMSG: ', errmsg

; structure containing fit coefs
coefs = {vsys: fitcoefs[0]*vel_factor, vsys_error: perror[0]*vel_factor, $
         tauv: fitcoefs[1], tauv_err: perror[1], $
         velstart: velstart, fixedvbool: fixedVBool, emmaskw: emmaskw, $
         light_frac: fitcoefs[2:*]*light_factor, $
         light_frac_err: perror[2:*]*light_factor, $
         model_age: model.age[vdidx,*], chisq: 0.0D, $
         redchi: 0.0D, bluechi: 0.0D, hkchi: 0.0D, $
         bluefree: 0L, MMWA: 0.0D, MLWA: 0.0D, $
         MMWZ: 0.0D, MLWZ: 0.0D, SNR: 0.0D}

; fit to full spectrum including masked pixels
redidx = where(restwl ge 5400)
blueidx = where(restwl lt 5400)
hklow = 3920
hkhigh = 4000
hkidx = where(restwl gt hklow and restwl lt hkhigh)

yfit = bc_mcombine_allz2(restwl, fitcoefs, mlib=custom_lib)

coefs.bluefree = (n_elements(blueidx) - n_elements(fitcoefs) - 1)

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


coefs.MMWZ = total(model.Z[vdidx,*]*coefs.light_frac/model.norm[vdidx,*]) $
          / total(coefs.light_frac/model.norm[vdidx,*])
coefs.MLWZ = total(light_weight * model.Z[vdidx,*]) / total(light_weight)

galfit = flux  + 'NaN'
galfit[ok] = flux[ok]

chivec = (flux - yfit)/err

print, fitcoefs, format = '(60F6.1)'
print, coefs.vsys
print, coefs.tauv
print, coefs.hkchi
print, coefs.MLWZ
print, coefs.MMWZ
print, '$$$$$$$$$$$$$$$$$$$$$', coefs.MLWA

return, coefs

end
