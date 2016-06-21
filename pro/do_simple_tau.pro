
pro do_simple_tau, datafile, errorfile, output, location=location, $
                   model=model, plot=plot, $
                   wavemin=wavemin, wavemax=wavemax, lightmin=lightmin, $
                   lightmax=lightmax, multimodel=multimodel, savefiber=savefiber

; read in models
if keyword_set(multimodel) then begin
   readcol, model, metals, models, format='F,A'
   m = mrdfits(models[0], 1)
endif else begin
   if not n_elements(model) then model=$
      '/d/monk/eigenbrot/WIYN/14B-0456/anal/models/bc03_solarZ_ChabIMF.fits'
   m = mrdfits(model, 1)
   metal = -99.9
endelse

m.flux *= rebin(transpose(m.norm), size(m.flux, /dimensions))

; read in data
data = MRDFITS(datafile,0,header)
error = MRDFITS(errorfile,0)

numfibers = N_ELEMENTS(data[0,*])
wavesize = N_ELEMENTS(data[*,0])

cdelt = float(sxpar(header,'CDELT1'))
crval = float(sxpar(header,'CRVAL1'))
crpix = float(sxpar(header,'CRPIX1'))
print, 'CDELT1 = ',cdelt
print, 'CRVAL1 = ',crval
print, 'CRPIX1 = ',crpix
wave = (FINDGEN(wavesize) - crpix) * cdelt + crval
;vdisp = 377. ; measured velocity dispersion

if keyword_set(location) then begin
   readcol, location, apnum, fiber_radii, ra, dec, rkpc, zkpc
   sizeidx = [0.937,1.406,1.875,2.344,2.812]
endif

vdisp = [493., 589., 691., 796., 966.]/2.355
size_borders = [19, 43, 62, 87, 109] ; The last one is needed to prevent indexing errors
size_switch = 0

if n_elements(wavemin) eq 0 then $
   wavemin = min(wave)
if n_elements(wavemax) eq 0 then $
   wavemax = max(wave)

idx = where(wave ge wavemin and wave le wavemax)
wave = wave[idx]

if n_elements(lightmin) eq 0 then $
   lightmin = 5450
if n_elements(lightmax) eq 0 then $
   lightmax = 5550

lightidx = where(wave ge lightmin and wave le lightmax)

agearr = m.age/1e9
numages = N_ELEMENTS(agearr)
colarr = STRARR(numages)


FOR k=0, numages - 1 DO BEGIN
   colarr[k] = string(agearr[k],' Gyr',format='(F6.3,A4)')
ENDFOR

t3d, /reset;, translate=[-1,-1,0], rotate=[0,0,180]
fmt = '(I11,'+string(numages+2)+'E13.3,F7.2,F12.3,2E12.3,F10.3,E12.3,2F13.3)'
openw, lun, output, /get_lun
printf, lun, '# Generated on ',systime()
printf, lun, '# Data file: ',datafile
printf, lun, '# Error file: ',errorfile
printf, lun, '# Model file: ',model,format='(A14,A90)'
printf, lun, '# Fiber Num',colarr,'MMWA [Gyr]','MLWA [Gyr]',$
        'Tau_V','S/N','Chisq','redChi','Z/Z_sol','psi_0','tau_sf','t_form',$
        format='(A-11,'+string(numages+2)+'A13,A7,3A12,A10,A12,2A13)'
printf, lun, '#'

if n_elements(savefiber) ne 0 then begin
   savename = (strsplit(output,'.',/extract))[0] + '_steps.dat'
   openw, savelun, savename, /get_lun
   printf, savelun, '# Generated on ',systime()
   printf, savelun, '# Data file: ',datafile
   printf, savelun, '# Error file: ',errorfile
   printf, savelun, '# Model file: ',model,format='(A14,A90)'
   printf, savelun, '# Fiber Num',colarr,'MMWA [Gyr]','MLWA [Gyr]',$
           'Tau_V','S/N','Chisq','redChi','Z/Z_sol',$
           format='(A-11,'+string(numages+2)+'A13,A7,3A12,2A10)'
   printf, savelun, '#'
endif else begin
   savefiber = -1
   savelun = 0
endelse

if keyword_set(plot) then begin
   plotfile = (strsplit(output,'.',/extract))[0] + '.ps'
   dfpsplot, plotfile, /color, /times, /landscape
endif

fitsfile = (strsplit(output,'.',/extract))[0] + '.fits'
outputarray = {TAUV: 0.0D, TAUV_ERR: 0.0D, $
               psi0: 0.0D, psi0_err: 0.0D, $
               tau_sf: 0.0D, tau_sf_err: 0.0D, $
               t_form: 0.0D, t_form_err: 0.0D, $
               MODEL_AGE: dblarr(10),$
               CHISQ: 0.0D, REDCHI: 0.0D, $
               mass: dblarr(10), weighted_age: dblarr(10)}
outputarray = replicate(outputarray, numfibers)

L_sun = 3.826e33 ;ergs s^-1
dist_mpc = 10.062
flux_factor = 1d17 ;to avoid small number precision errors
tau = 2*!DPI

for i = 0, numfibers - 1 DO BEGIN
   
   print, 'Grabbing aperture '+string(i+1,format='(I3)')
   flux = data[idx,i]*flux_factor
   err = error[idx,i]*flux_factor
   
   if keyword_set(multimodel) then begin
      m = mrdfits(models[i], 1)
      m.flux *= rebin(transpose(m.norm), size(m.flux, /dimensions))
      metal = metals[i]
      print, 'Using mode '+models[i]
   endif

   vd = 200.

   if keyword_set(location) then begin
      lidx = where(sizeidx eq fiber_radii[i])
      vd = vdisp[lidx]
      plotlabel = string('Aperture',i+1,'r=',rkpc[i],'z=',zkpc[i],$
                         format='(A8,I4,A3,F5.2,A3,F5.2)')
      print, plotlabel
   endif else begin
      print, i, size_borders
      if i eq size_borders[0] then begin
         size_switch += 1
         size_borders = size_borders[1:*]
         vd = vdisp[size_switch]
         plotlabel = string('Aperture',i+1,format='(A8,I4)')
      endif
   endelse

   if i eq savefiber then begin
      savestep = 1
   endif else begin
      savestep = 0
   endelse
   
   print, '#########', vd

; fit continuum
   coef = bc_continuum_tau(m, wave, flux, err, vd, $
                           plotlabel=plotlabel,$
                           yfit=continuum, $
                           savestep=savestep, lun=savelun, $
                           lightidx=lightidx, fmt=fmt)

;; ; measure absorption line indices
;;    icoef = absline_index(wave, flux, err)
;;    mcoef = absline_index(wave, continuum, tag='_model') ; measure off model

;; ; fit emission lines
;;  lcoef = emlinefit(wave, flux, continuum, err, inst_res = 65.0, yfit = emfit)
 
;; ; Create output structure & save
;;    s = create_struct(coef, icoef, mcoef, lcoef)
;mwrfits, s, 'NGC_test.fits', /create

   outputarray[i] = coef
   
   goodidx = where(coef.mass ne -99)

   ;SNR = sqrt(total((flux[lightidx]/err[lightidx])^2)/n_elements(lightidx))
   SNR = mean(flux[lightidx]/err[lightidx])

   MMWA = total(coef.weighted_age[goodidx]*coef.mass[goodidx]) $
          / total(coef.mass[goodidx])

   redd = exp(-coef.tauv*(wave[lightidx]/5500)^(-0.7))

   light_weight = mean(m.flux[lightidx,goodidx,0] * $
                       rebin(redd,n_elements(lightidx),$
                             n_elements(goodidx)), $
                       dimension=1) * coef.mass[goodidx]
   MLWA = total(light_weight * coef.weighted_age[goodidx]) / total(light_weight)

   printf, lun, i+1, coef.mass, MMWA, MLWA, coef.tauv,$
           SNR, coef.chisq, coef.redchi, metal, coef.psi0, $
           coef.tau_sf, coef.t_form, format=fmt

ENDFOR

if keyword_set(plot) then dfpsclose

free_lun, lun
mwrfits, outputarray, fitsfile, /create

end
