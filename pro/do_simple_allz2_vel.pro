
pro do_simple_allZ2_vel, datafile, errorfile, output, location=location, $
                         fitfile=fitfile, fitregion=fitregion, $
                         wavemin=wavemin, wavemax=wavemax, lightmin=lightmin, $
                         lightmax=lightmax, multimodel=multimodel, savestep=savestep

;;;Here are the aperture subset assignments we will use. This is only
;;;temporary. 1.19.16
;
P1subset = [6,8,10,26,32] - 1
P2subset = [5,8,17,22,23,36] - 1
P3subset = [5,8,13,15,20,25] - 1
P4subset = [6,38,43] - 1
P5subset = [14,16] - 1
P6subset = [8,17] - 1

Pnum = fix((stregex(datafile, '_P([1-9])_', /subexpr, /extract))[1])

case Pnum of
   1: subidx = P1subset
   2: subidx = P2subset
   3: subidx = P3subset
   4: subidx = P4subset
   5: subidx = P5subset
   6: subidx = P6subset
   else: stop
endcase   

;defplotcolors
; read in models

fits = mrdfits(fitfile, 0)

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
wave = (DINDGEN(wavesize) - (crpix - 1)) * cdelt + crval
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

t3d, /reset;, translate=[-1,-1,0], rotate=[0,0,180]
fmt = '(I11,F13.3,F11.3,4E25.17)'
openw, lun, output, /get_lun
printf, lun, '# Generated on ',systime()
printf, lun, '# Data file: ',datafile
printf, lun, '# Error file: ',errorfile
printf, lun, '# Previous fit: ',fitfile,format='(A14,A90)'
printf, lun, '# Fiber Num','V_sys','S/N','Chisq',$
        'redChi','blueChi','HKChi',$
        format='(A-11,A13,A11,4A25)'
printf, lun, '#'

fitsfile = (strsplit(output,'.',/extract))[0] + '.coef.vel.fits'
outputarray = {VSYS: 0.0D, VSYS_ERROR: 0.0D,$
               CHISQ: 0.0D, REDCHI: 0.0D, BLUECHI: 0.0D, HKCHI: 0.0D, $
               BLUEFREE: 0L, SNR: 0.0D}
outputarray = replicate(outputarray, numfibers)

chifile = (strsplit(output,'.',/extract))[0] + '.chi.vel.fits'
chiarray = fltarr(numfibers, n_elements(wave))

yfitfile = (strsplit(output,'.',/extract))[0] + '.fit.vel.fits'
yfitarray = fltarr(n_elements(wave), numfibers)

L_sun = 3.826e33 ;ergs s^-1
dist_mpc = 10.062
flux_factor = 1d17 ;to avoid small number precision errors
tau = 2*!DPI

for i = 0, numfibers - 1 DO BEGIN
;foreach i, subidx DO BEGIN
  
   print, 'Grabbing fiber '+string(i+1,format='(I3)')
   flux = data[idx,i]*flux_factor
   err = error[idx,i]*flux_factor
   prevfit = fits[*,i]*flux_factor
   ;; plot, wave, flux
   ;; oplot, wave, prevfit
   ;; print, flux - prevfit
   ;; stop

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

   ;; if keyword_set(savestep) then begin
   ;;    savename = 'steps/' + (strsplit(output,'.',/extract))[0] + '_' + string(i,format='(I02)') + '_steps.dat'
   ;;    openw, savelun, savename, /get_lun
   ;;    printf, savelun, '# Generated on ',systime()
   ;;    printf, savelun, '# Data file: ',datafile
   ;;    printf, savelun, '# Error file: ',errorfile
   ;;    printf, savelun, '# Model file: ',model,format='(A14,A90)'
   ;;    printf, savelun, '# Fiber Num',colarr,'MMWA [Gyr]','MLWA [Gyr]',$
   ;;            'MMWZ [Z_sol]','MLWZ [Z_sol]','V_sys','Tau_V','S/N','Chisq',$
   ;;            'redChi','blueChi','HKChi',$
   ;;            format='(A-11,'+string(numages + 5)+'A13,2A11,4A25)'
   ;;    printf, lun, '#'
   ;; endif else begin
   ;;    savestep = 0
   ;;    savelun = 0
   ;; endelse

; fit continuum
   savestep=0
   coef = bc_continuum_allZ2_vel(prevfit, wave, flux, err, vdidx, $
                                 fitregion=fitregion, $
                                 yfit=yfit, $
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

   printf, lun, i+1, coef.vsys, coef.SNR, $
           coef.chisq, coef.redchi, coef.bluechi, coef.hkchi, format=fmt
   print, i

ENDFOR
;ENDFOREACH


free_lun, lun
mwrfits, outputarray, fitsfile, /create
mwrfits, transpose(chiarray), chifile, /create
mwrfits, yfitarray, yfitfile, /create

end
