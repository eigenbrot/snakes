
pro do_simple_sky, datafile, errorfile, output, location=location, $
                   model=model, plot=plot, bluefit=bluefit, skysub=skysub,$
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
OGwave = (FINDGEN(wavesize) - crpix) * cdelt + crval
;vdisp = 377. ; measured velocity dispersion

if keyword_set(location) then begin
   readcol, location, apnum, fiber_radii, ra, dec, rkpc, zkpc
   sizeidx = [0.937,1.406,1.875,2.344,2.812]
endif

vdisp = [493., 589., 691., 796., 966.]/2.355
size_borders = [19, 43, 62, 87, 109] ; The last one is needed to prevent indexing errors
size_switch = 0

if n_elements(wavemin) eq 0 then $
   wavemin = min(OGwave)
if n_elements(wavemax) eq 0 then $
   wavemax = max(OGwave)

idx = where(OGwave ge wavemin and OGwave le wavemax)
wave = OGwave[idx]

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
fmt = '(I11,'+string(numages+2)+'E13.3,F7.2,F12.3,4F12.3,3F10.3)'
openw, lun, output, /get_lun
printf, lun, '# Generated on ',systime()
printf, lun, '# Data file: ',datafile
printf, lun, '# Error file: ',errorfile
printf, lun, '# Model file: ',model,format='(A14,A90)'
printf, lun, '# Fiber Num',colarr,'MMWA [Gyr]','MLWA [Gyr]',$
        'Tau_V','S/N','Chisq','redChi','blueChi','HKChi','Z/Z_sol','Skyfrac',$
        format='(A-11,'+string(numages+2)+'A13,A7,5A12,3A10)'

if n_elements(savefiber) ne 0 then begin
   savename = (strsplit(output,'.',/extract))[0] + '_steps.dat'
   openw, savelun, savename, /get_lun
   printf, savelun, '# Generated on ',systime()
   printf, savelun, '# Data file: ',datafile
   printf, savelun, '# Error file: ',errorfile
   printf, savelun, '# Model file: ',model,format='(A14,A90)'
   printf, savelun, '# Fiber Num',colarr,'MMWA [Gyr]','MLWA [Gyr]',$
           'Tau_V','S/N','Chisq','redChi','blueChi','HKChi','Z/Z_sol',$
           format='(A-11,'+string(numages+2)+'A13,A7,5A12,2A10)'   
   printf, savelun, '#'

   startfiber = savefiber - 1
   endfiber = savefiber - 1

endif else begin
   savefiber = -1
   savelun = 0
   startfiber = 0
   endfiber = numfibers - 1
endelse

if keyword_set(plot) then begin
   plotfile = (strsplit(output,'.',/extract))[0] + '.ps'
   dfpsplot, plotfile, /color, /times, /landscape
endif

fitsfile = (strsplit(output,'.',/extract))[0] + '.fits'
outputarray = {TAUV: 0.0D, TAUV_ERR: 0.0D, LIGHT_FRAC: dblarr(10),$
               LIGHT_FRAC_ERR: dblarr(10), SKYFRAC: 0.0D, MODEL_AGE: dblarr(10),$
               CHISQ: 0.0D, REDCHI: 0.0D, BLUECHI: 0.0D, HKCHI: 0.0D, $
               BLUEFREE: 0L, MMWA: 0.0D, MLWA: 0.0D, SNR: 0.0D}
outputarray = replicate(outputarray, numfibers)

chifile = (strsplit(output,'.',/extract))[0] + '.chi.fits'
chiarray = fltarr(numfibers, n_elements(wave))

L_sun = 3.826e33 ;ergs s^-1
dist_mpc = 10.062
flux_factor = 1d17 ;to avoid small number precision errors
tau = 2*!DPI

for i = startfiber, endfiber DO BEGIN
   
   print, 'Grabbing fiber '+string(i+1,format='(I3)')
   flux = data[idx,i]*flux_factor
   err = error[idx,i]*flux_factor
   
   if keyword_set(location) then begin
      lidx = where(sizeidx eq fiber_radii[i])
      vd = vdisp[lidx]
      plotlabel = string('Aperture',i+1,'r=',rkpc[i],'z=',zkpc[i],$
                         format='(A8,I4,A3,F6.2,A3,F5.2)')
      print, plotlabel

   endif else begin
      print, i, size_borders
      if i eq size_borders[0] then begin
         size_switch += 1
         size_borders = size_borders[1:*]
      endif
      vd = vdisp[size_switch]
      plotlabel = string('Aperture',i+1,format='(A8,I4)')
   endelse

   if keyword_set(multimodel) then begin
      m = mrdfits(models[i], 1)
      metal = metals[i]
      plotlabel += '!c!c'
      plotlabel += string('Z/Z_sol = ',metal,format='(A10,F8.4)')
      print, 'Using model '+models[i]
   endif

   if i eq savefiber - 1 then begin
      savestep = 1
   endif else begin
      savestep = 0
   endelse

; fit continuum
   print, vd
   coef = bc_continuum_sky(m, wave, flux, err, vd, $
                           plotlabel=plotlabel, $
                           bluefit=bluefit, $
                           yfit=continuum, $
                           savestep=savestep, lun=savelun, $
                           lightidx=lightidx, fmt=fmt, $
                           chivec=chivec)

   if keyword_set(skysub) then begin
      sky = interpol(m.flux[*,0]*m.skynorm, m.wave, OGwave)
      sky *= coef.skyfrac
      print, '%%%%%%%%%',mean(data[*,i]), mean(sky)
      data[*,i] -= sky
   endif

;; ; measure absorption line indices
;;    icoef = absline_index(wave, flux, err)
;;    mcoef = absline_index(wave, continuum, tag='_model') ; measure off model

;; ; fit emission lines
;;  lcoef = emlinefit(wave, flux, continuum, err, inst_res = 65.0, yfit = emfit)
 
;; ; Create output structure & save
;;    s = create_struct(coef, icoef, mcoef, lcoef)
;mwrfits, s, 'NGC_test.fits', /create

   print, coef.bluefree
   outputarray[i] = coef
   chiarray[i,*] = chivec
   print, size(chiarray, /dimensions), n_elements(chivec)
;   print, chivec
   ;SNR = sqrt(total((flux[lightidx]/err[lightidx])^2)/n_elements(lightidx))
   ;; SNR = mean(flux[lightidx]/err[lightidx])

   ;; MMWA = total(agearr*coef.light_frac*1./m.norm) $
   ;;        / total(coef.light_frac*1./m.norm)

   ;; redd = exp(-coef.tauv*(wave[lightidx]/5500)^(-0.7))
   ;; light_weight = mean(m.flux[lightidx,*] * rebin(redd,n_elements(lightidx),$
   ;;                                                n_elements(agearr)), $
   ;;                     dimension=1) * coef.light_frac
   ;; MLWA = total(light_weight * agearr) / total(light_weight)

   if i eq startfiber then begin
      printf, lun, '# Blue_free: ', coef.bluefree
      printf, lun, '#'
   endif

   printf, lun, i+1, coef.light_frac/m.norm, coef.MMWA, coef.MLWA, coef.tauv,$
           coef.SNR, coef.chisq, coef.redchi, coef.bluechi, coef.hkchi, metal,coef.skyfrac,$
           format=fmt

ENDFOR

if keyword_set(plot) then dfpsclose
if keyword_set(skysub) then writefits, skysub, data, header

chiplot = (strsplit(output,'.',/extract))[0] + '.chi.ps'
dfpsplot, chiplot, /color, /times, /landscape
meanchi = smooth(mean(chiarray,dimension=1),5,/NAN)
stdchi = smooth(stddev(chiarray,dimension=1),5,/NAN)
plot, wave, meanchi, xtitle='Wavelength', ytitle='<Chi>', /nodata, $
      xrange = [min(wave), max(wave)], yrange=[-10,10], /t3d, /xs, /ys

oband, wave[0:*:5], (meanchi-stdchi)[0:*:5], $
       (meanchi+stdchi)[0:*:5], color=!gray, /noclip, /t3d

oplot, wave, meanchi, color=!black, /t3d

sk =    [6300.,        5890., 5683.8, 5577.,      5461., 5199.,      4983., 4827.32, 4665.69, 4420.23, 4358., 4165.68, 4047.0]
sknam = ['[OI] (atm)', 'NaD', 'NaI',  'OI (atm)', 'HgI', 'NI (atm)', 'NaI', 'HgI',   'NaI',   'NaI',   'HgI', 'NaI',   'HgI']
em = [6563.8,  6716.0]
emnam = ['Ha', 'S2']

abs =    [3933.7, 3968.5, 4304.4,   5175.3, 5894.0, 4861., 4341., 4102.]
absnam = ['H',    'K',    'G band', 'Mg',   'Na',   'HB',  'HG',  'HD']
for s=0, n_elements(sk) - 1 do begin
   ypos = abs(interpol(meanchi, wave, sk[s])*1.7)
   xyouts, sk[s], ypos, sknam[s], alignment=0.5, charsize=0.5, /data
endfor

for e=0, n_elements(em) - 1 do begin
   ypos = abs(interpol(meanchi, wave, em[e])*1.7)
   xyouts, em[e], ypos, emnam[e], alignment=0.5, charsize=0.5, /data, color=!blue
endfor

for a=0, n_elements(abs) - 1 do begin
   ypos = abs(interpol(meanchi, wave, abs[a])*1.7)
   xyouts, abs[a], ypos, absnam[a], alignment=0.5, charsize=0.5, /data, color=!red
endfor
dfpsclose

free_lun, lun
if savefiber ne -1 then free_lun, savelun
mwrfits, outputarray, fitsfile, /create
mwrfits, transpose(chiarray), chifile, /create
print, m.norm

end
