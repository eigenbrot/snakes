
pro do_simple_allZ2, datafile, errorfile, output, location=location, $
                     model=model, fitregion=fitregion, velstart=velstart, $
                     wavemin=wavemin, wavemax=wavemax, lightmin=lightmin, $
                     lightmax=lightmax, multimodel=multimodel, savestep=savestep
;defplotcolors
; read in models
if not n_elements(model) then model=$
   '/d/monk/eigenbrot/WIYN/14B-0456/anal/models/allZ2/allz2_test.fits'

m = mrdfits(model, 1)

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

agearr = m.age[0,*]/1e9
numages = N_ELEMENTS(agearr)
colarr = STRARR(numages)

FOR k=0, numages - 1 DO BEGIN
   colarr[k] = string(agearr[k],floor(k/10),' Gyr',format='(F7.3,I2,A4)')
ENDFOR

t3d, /reset;, translate=[-1,-1,0], rotate=[0,0,180]
fmt = '(I11,'+string(numages + 5)+'F13.3,2F11.3,4E25.17)'
openw, lun, output, /get_lun
printf, lun, '# Generated on ',systime()
printf, lun, '# Data file: ',datafile
printf, lun, '# Error file: ',errorfile
printf, lun, '# Model file: ',model,format='(A14,A90)'
printf, lun, '# Fiber Num',colarr,'MMWA [Gyr]','MLWA [Gyr]',$
        'MMWZ [Z_sol]','MLWZ [Z_sol]','V_sys','Tau_V','S/N','Chisq',$
        'redChi','blueChi','HKChi',$
        format='(A-11,'+string(numages + 5)+'A13,2A11,4A25)'
printf, lun, '#'

fitsfile = (strsplit(output,'.',/extract))[0] + '.coef.fits'
outputarray = {VSYS: 0.0D, VSYS_ERROR: 0.0D,TAUV: 0.0D, TAUV_ERR: 0.0D, $
               LIGHT_FRAC: dblarr(numages),$
               LIGHT_FRAC_ERR: dblarr(numages), $
               MODEL_AGE: fltarr(numages), $
               CHISQ: 0.0D, REDCHI: 0.0D, BLUECHI: 0.0D, HKCHI: 0.0D, $
               BLUEFREE: 0L, MMWA: 0.0D, MLWA: 0.0D, MMWZ: 0.0D, MLWZ: 0.0D, SNR: 0.0D}
outputarray = replicate(outputarray, numfibers)

chifile = (strsplit(output,'.',/extract))[0] + '.chi.fits'
chiarray = fltarr(numfibers, n_elements(wave))

yfitfile = (strsplit(output,'.',/extract))[0] + '.fit.fits'
yfitarray = fltarr(n_elements(wave), numfibers)

L_sun = 3.826e33 ;ergs s^-1
dist_mpc = 10.062
flux_factor = 1d17 ;to avoid small number precision errors
tau = 2*!DPI

for i = 0, numfibers - 1 DO BEGIN
   
   print, 'Grabbing fiber '+string(i+1,format='(I3)')
   flux = data[idx,i]*flux_factor
   err = error[idx,i]*flux_factor
   
   if keyword_set(location) then begin
      vdidx = where(sizeidx eq fiber_radii[i])
   endif else begin
      print, i, size_borders
      if i eq size_borders[0] then begin
         size_switch += 1
         size_borders = size_borders[1:*]
      endif
      vdidx = size_switch
   endelse

   if keyword_set(savestep) then begin
      savename = 'steps/' + (strsplit(output,'.',/extract))[0] + '_' + string(i,format='(I02)') + '_steps.dat'
      openw, savelun, savename, /get_lun
      printf, savelun, '# Generated on ',systime()
      printf, savelun, '# Data file: ',datafile
      printf, savelun, '# Error file: ',errorfile
      printf, savelun, '# Model file: ',model,format='(A14,A90)'
      printf, savelun, '# Fiber Num',colarr,'MMWA [Gyr]','MLWA [Gyr]',$
              'MMWZ [Z_sol]','MLWZ [Z_sol]','V_sys','Tau_V','S/N','Chisq',$
              'redChi','blueChi','HKChi',$
              format='(A-11,'+string(numages + 5)+'A13,2A11,4A25)'
      printf, lun, '#'
   endif else begin
      savestep = 0
      savelun = 0
   endelse

; fit continuum
   coef = bc_continuum_allZ2(m, wave, flux, err, vdidx, $
                             fitregion=fitregion, $
                             yfit=yfit, velstart=velstart, $
                             savestep=savestep, lun=savelun, $
                             lightidx=lightidx, fmt=fmt, $
                             chivec=chivec)
   
   if savestep then close, savelun
;; ; measure absorption line indices
;;    icoef = absline_index(wave, flux, err)
;;    mcoef = absline_index(wave, continuum, tag='_model') ; measure off model

;; ; fit emission lines
;;  lcoef = emlinefit(wave, flux, continuum, err, inst_res = 65.0, yfit = emfit)
 
;; ; Create output structure & save
;;    s = create_struct(coef, icoef, mcoef, lcoef)
;mwrfits, s, 'NGC_test.fits', /create

   ;; print, size(outputarray[i].vsys, /type), size(coef.vsys, /type)
   ;; print, size(outputarray[i].tauv, /type), size(coef.tauv, /type)
   ;; print, size(outputarray[i].light_frac, /type), size(coef.light_frac, /type)
   ;; print, size(outputarray[i].light_frac, /dimensions), size(coef.light_frac, /dimensions)
   ;; print, size(outputarray[i].light_frac_err, /dimensions), size(coef.light_frac_err, /dimensions)
   ;; print, size(outputarray[i].model_age, /type), size(coef.model_age, /type)
   ;; print, size(outputarray[i].model_age, /dimensions), size(coef.model_age, /dimensions)
   ;; print, size(outputarray[i].chisq, /type), size(coef.chisq, /type)
   ;; print, size(outputarray[i].MLWA, /type), size(coef.MLWA, /type)
   ;; print, size(outputarray[i].MMWA, /type), size(coef.MMWA, /type)
   ;; print, size(outputarray[i].SNR, /type), size(coef.SNR, /type)
   ;; print, size(outputarray[i].redchi, /type), size(coef.redchi, /type)
   ;; print, size(outputarray[i].bluechi, /type), size(coef.bluechi, /type)

   print, coef.bluefree
   outputarray[i] = coef
   chiarray[i,*] = chivec
   yfitarray[*,i] = yfit/flux_factor

   if i eq 0 then begin
      printf, lun, '# Blue_free: ', coef.bluefree
      printf, lun, '#'
   endif

   printf, lun, i+1, coef.light_frac/m.norm, coef.MMWA, $
           coef.MLWA, coef.MMWZ, coef.MLWZ, coef.vsys, coef.tauv, coef.SNR, $
           coef.chisq, coef.redchi, coef.bluechi, coef.hkchi, format=fmt
   print, i

ENDFOR


free_lun, lun
mwrfits, outputarray, fitsfile, /create
mwrfits, transpose(chiarray), chifile, /create
mwrfits, yfitarray, yfitfile, /create
print, m.norm

end
