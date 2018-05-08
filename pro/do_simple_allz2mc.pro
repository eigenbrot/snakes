
pro do_simple_allz2mc, datafile, errorfile, output, NMC=NMC, location=location, $
                       velocity=velocity, emmaskw=emmaskw, mcerrfile=mcerrfile, $
                       maskBalm=maskBalm, parinfofile=parinfofile, $
                       model=model, fitregion=fitregion, velstart=velstart, $
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
   else: print, "Couldn't identify pointing"
endcase   

;defplotcolors
; read in models
if not n_elements(model) then model=$
   '/d/monk/eigenbrot/WIYN/14B-0456/anal/models/allZ2/allz2_test.fits'

if not n_elements(NMC) then NMC=30

m = mrdfits(model, 1)

; read in data
data = MRDFITS(datafile,0,header)
error = MRDFITS(errorfile,0)

if n_elements(mcerrfile) gt 0 then begin
   mcerr = MRDFITS(mcerrfile,0)
endif else begin
   mcerr = error
endelse 

numfibers = N_ELEMENTS(data[0,*])
wavesize = N_ELEMENTS(data[*,0])

cdelt = float(sxpar(header,'CDELT1'))
crval = float(sxpar(header,'CRVAL1'))
crpix = float(sxpar(header,'CRPIX1'))
print, 'CDELT1 = ',cdelt
print, 'CRVAL1 = ',crval
print, 'CRPIX1 = ',crpix
wave = (DINDGEN(wavesize) - crpix) * cdelt + crval
;vdisp = 377. ; measured velocity dispersion

if keyword_set(location) then begin
   readcol, location, apnum, fiber_radii, ra, dec, rkpc, zkpc
   sizeidx = [0.937,1.406,1.875,2.344,2.812]
endif

if n_elements(velstart) eq 0 then velstart = 528.

if keyword_set(velocity) then begin
   readcol, velocity, apnum, vobs, dvobs, vstar, known_V, vhans, dvhans
   velstart = 12345.
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

fitsfile = strjoin((strsplit(output,'.',/extract))[0:-2],'.') + '.coef.fits'
outputarray = {VSYS: 0.0D, VSYS_ERROR: 0.0D,TAUV: 0.0D, TAUV_ERR: 0.0D, $
               VELSTART: velstart, FIXEDVBOOL: 0, emmaskw: 0.0, $
               LIGHT_FRAC: dblarr(numages),$
               LIGHT_FRAC_ERR: dblarr(numages), $
               LIGHT_WEIGHT: dblarr(numages), LIGHT_WEIGHT_ERR: dblarr(numages), $
               MODEL_AGE: fltarr(numages), $
               CHISQ: 0.0D, REDCHI: 0.0D, BLUECHI: 0.0D, HKCHI: 0.0D, $
               BLUEFREE: 0L, TOTFREE: 0L, REDFREE: 0L, HKFREE: 0L,$
               MMWA: 0.0D, MLWA: 0.0D, dMLWA: 0.0D, MMWZ: 0.0D, $
               MLWZ: 0.0D, dMLWZ: 0.0D, SNR: 0.0D}
outputarray = replicate(outputarray, numfibers)

chifile = (strsplit(output,'.',/extract))[0] + '.chi.fits'
chiarray = fltarr(numfibers, n_elements(wave))

yfitfile = (strsplit(output,'.',/extract))[0] + '.fit.fits'
yfitarray = fltarr(n_elements(wave), numfibers)

L_sun = 3.826e33 ;ergs s^-1
dist_mpc = 10.062
flux_factor = 1d17 ;to avoid small number precision errors
tau = 2*!DPI

;smooth models to in the same way we will smooth the data
mdims = size(m.flux, /dimensions)
for j = 0, mdims[0] - 1 do begin
   for k = 0, mdims[2] - 1 do begin
      print, size(m.flux[j,*,k], /dimensions)
      m.flux[j,*,k] = convol(reform(m.flux[j,*,k]),[0.7,1,0.7],/edge_mirror)
   endfor
endfor

if file_test('MCdir',/directory) eq 0 then file_mkdir, 'MCdir'

for i = 0, numfibers - 1 DO BEGIN
;foreach i, subidx DO BEGIN
  
   print, 'Grabbing fiber '+string(i+1,format='(I3)')
   flux = data[idx,i]*flux_factor
   err = error[idx,i]*flux_factor
   mce = mcerr[idx,i]*flux_factor
   
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

   if keyword_set(velocity) then begin
      fixedV = known_V[i]
   endif else begin
      fixedV = 0
   endelse

   if keyword_set(parinfofile) then begin
      parinfo = mrdfits(parinfofile, i+1)
   endif else begin
      parinfo = 0
   endelse
   

   if keyword_set(savestep) then begin
      savename = 'steps/' + savestep + '_' + string(i+1,format='(I02)') + '_steps.dat'
      openw, savelun, savename, /get_lun
      printf, savelun, '# Generated on ',systime()
      printf, savelun, '# Data file: ',datafile
      printf, savelun, '# Error file: ',errorfile
      printf, savelun, '# Model file: ',model,format='(A14,A90)'
      printf, savelun, '# Fiber Num',colarr,'MMWA [Gyr]','MLWA [Gyr]',$
              'MMWZ [Z_sol]','MLWZ [Z_sol]','V_sys','Tau_V','S/N','Chisq',$
              'redChi','blueChi','HKChi',$
              format='(A-11,'+string(numages + 5)+'A13,2A11,4A25)'
      printf, savelun, '#'
   endif else begin
      savestep = 0
      savelun = 0
   endelse

   ;Do the Monte Carlo stuff
   MCfile = 'MCdir/' + strjoin((strsplit(output,'.',/extract))[0:-2],'.') + string('.MC',i+1,'.fits',format='(A3,I03,A5)')
   MCarray = {TAUV: 0.0D, $
              LIGHT_FRAC: dblarr(numages),$
              LIGHT_WEIGHT: dblarr(numages), $
              CHISQ: 0.0D, REDCHI: 0.0D, BLUECHI: 0.0D, HKCHI: 0.0D, $
              MMWA: 0.0D, MLWA: 0.0D, MMWZ: 0.0D, MLWZ: 0.0D}
   MCarray = replicate(MCarray, NMC)

   for NN = 0, NMC - 1 do begin
      df = randomn(seed, n_elements(flux))*mce
      tmpflux = convol(flux+df, [0.7,1,0.7],/edge_mirror)
      tcoef = bc_continuum_allZ2(m, wave, tmpflux, err, vdidx, $
                                 fitregion=fitregion, fixedV=fixedV, $
                                 yfit=yfit, velstart=velstart, $
                                 maskBalm=maskBalm, emmaskw=emmaskw, $
                                 savestep=savestep, lun=savelun, $
                                 lightidx=lightidx, fmt=fmt, $
                                 chivec=chivec, parinfo=parinfo)

      MCarray[NN].tauv = tcoef.tauv
      MCarray[NN,*].light_frac = tcoef.light_frac
      MCarray[NN,*].light_weight = tcoef.light_weight
      MCarray[NN].MLWA = tcoef.MLWA
      MCarray[NN].MMWA = tcoef.MMWA
      MCarray[NN].MLWZ = tcoef.MLWZ
      MCarray[NN].MMWZ = tcoef.MMWZ
      MCarray[NN].CHISQ = tcoef.chisq
      MCarray[NN].redchi = tcoef.redchi
      MCarray[NN].bluechi = tcoef.bluechi
      MCarray[NN].hkchi = tcoef.hkchi
      
   endfor
   
   outputarray[i].vsys = tcoef.vsys
   outputarray[i].vsys_error = tcoef.vsys_error
   outputarray[i].velstart = tcoef.velstart
   outputarray[i].fixedvbool = tcoef.fixedvbool
   outputarray[i].emmaskw = tcoef.emmaskw
   outputarray[i].model_age = tcoef.model_age
   outputarray[i].chisq = tcoef.chisq
   outputarray[i].redchi = tcoef.redchi
   outputarray[i].bluechi = tcoef.bluechi
   outputarray[i].hkchi = tcoef.hkchi
   outputarray[i].bluefree = tcoef.bluefree
   outputarray[i].totfree = tcoef.totfree
   outputarray[i].redfree = tcoef.redfree
   outputarray[i].hkfree = tcoef.hkfree
   outputarray[i].SNR = tcoef.snr
   outputarray[i].tauv = mean(MCarray.tauv)
   outputarray[i].tauv_err = stddev(MCarray.tauv)
   outputarray[i].light_frac = mean(MCarray.light_frac, dimension=2)
   outputarray[i].light_frac_err = stddev(MCarray.light_frac, dimension=2)
   outputarray[i].light_weight = mean(MCarray.light_weight, dimension=2)
   outputarray[i].light_weight_err = stddev(MCarray.light_weight, dimension=2)   
   outputarray[i].MLWA = mean(MCarray.MLWA)
   outputarray[i].dMLWA = stddev(MCarray.MLWA)
   outputarray[i].MLWZ = mean(MCarray.MLWZ)
   outputarray[i].dMLWZ = stddev(MCarray.MLWZ)
   outputarray[i].MMWA = mean(MCarray.MMWA)
   outputarray[i].MMWZ = mean(MCarray.MMWZ)

   mwrfits, MCarray, MCfile, /create

   if savestep then free_lun, savelun

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

   chiarray[i,*] = chivec
   yfitarray[*,i] = yfit/flux_factor

   if i eq 0 then begin
      printf, lun, '# Blue_free: ', tcoef.bluefree
      printf, lun, '#'
   endif

ENDFOR
;ENDFOREACH


free_lun, lun
mwrfits, outputarray, fitsfile, /create
mwrfits, transpose(chiarray), chifile, /create
mwrfits, yfitarray, yfitfile, /create

end
