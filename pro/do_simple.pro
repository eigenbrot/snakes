
pro do_simple, datafile, errorfile, output, model=model, plot=plot, $
               wavemin=wavemin, wavemax=wavemax

; read in models
if not n_elements(model) then model=$
   '/d/monk/eigenbrot/WIYN/14B-0456/anal/models/bc03_solarZ_ChabIMF.fit'
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
vdisp = 377. ; measured velocity dispersion

if n_elements(wavemin) eq 0 then $
   wavemin = min(wave)

if n_elements(wavemax) eq 0 then $
   wavemax = max(wave)

idx = where(wave ge wavemin and wave le wavemax)
wave = wave[idx]

agearr = m.age/1e6
numages = N_ELEMENTS(agearr)
colarr = STRARR(numages)


FOR k=0, numages - 1 DO BEGIN
   colarr[k] = string(agearr[k],' Myr',format='(I6,A4)')
ENDFOR

t3d, /reset;, translate=[-1,-1,0], rotate=[0,0,180]
fmt = '(I11,'+string(numages)+'F10.7)'
openw, lun, output, /get_lun
printf, lun, '# Fiber Num',colarr,format='(A-11,'+string(numages)+'A10)'
printf, lun, '#'

if keyword_set(plot) then begin
   plotfile = (strsplit(output,'.',/extract))[0] + '.ps'
   dfpsplot, plotfile, /color, xsize=9, ysize=9/1.5, /times ;/landscape
endif

for i = 74, 74  DO BEGIN
   
   print, 'Grabbing fiber '+string(i,format='(I3)')
   flux = data[idx,i]*1e19
   err = error[idx,i]*1e19
   
; fit continuum
   coef = bc_continuum(m, wave, flux, err, vdisp, $
                       plotlabel=string('Fiber',i+1,format='(A5,I3)'),$
                       yfit=continuum)

;; ; measure absorption line indices
;;    icoef = absline_index(wave, flux, err)
;;    mcoef = absline_index(wave, continuum, tag='_model') ; measure off model

;; ; fit emission lines
;;  lcoef = emlinefit(wave, flux, continuum, err, inst_res = 65.0, yfit = emfit)

;; ; Create output structure & save
;;    s = create_struct(coef, icoef, mcoef, lcoef)
;mwrfits, s, 'NGC_test.fits', /create

;   coef.light_frac *= m.norm
   coef.light_frac /= total(coef.light_frac)

   printf, lun, i+1, coef.light_frac, format=fmt

ENDFOR

if keyword_set(plot) then dfpsclose

free_lun, lun
;help, s, /struct

end
