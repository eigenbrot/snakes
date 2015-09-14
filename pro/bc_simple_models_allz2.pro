

; Write out BC03 model grid used in SDSS code re-binned to log-lambda

pro bc_simple_models_allZ2, output, plot=plot, age=age

if keyword_set(plot) then begin
   plotfile = (strsplit(output,'.',/extract))[0] + '.ps'
   dfpsplot, plotfile, /color, /times, /landscape
endif

angstrom = '!6!sA!r!u!9 %!6!n'
lambda = 'l'
;lambda = '!4'+string("153B)+'!X'
; Ages in Gyr

Z_arr = [0.0001, 0.0004, 0.004, 0.008, 0.02, 0.05]/0.02

if n_elements(age) eq 0 then $
   age=[0.005, 0.025, 0.100, 0.286, 0.640, 0.904, 1.434, 2.500, 5.000, 10.000]

isedpath = '/usr/users/tremonti/models/Padova1994/chabrier/'

; Set up log-lambda wavelength grid
npix = round((alog10(8500) - alog10(3400)) / 1e-4)
logwl = findgen(npix) * 1e-4 + alog10(3400)

tmp = im_read_bc03(isedpath=isedpath, metallicity=0, age=age,$
                   minwave=3400, maxwave=8500, bc03_extras=ssp_extras)

; Define output structure
m = {wave: 10.0^logwl, flux: fltarr(npix, n_elements(age) * 6), $
     age: [age,age,age,age,age,age]*1e9, Z: fltarr(n_elements(age) * 6), $
     id: string('All metallicities, Chabrier IMF',format='(A4,I1,A14)'), $
     norm: fltarr(n_elements(age) * 6), m_remaining: fltarr(n_elements(age)*6)}

; normalization wavelengths
nwl = where(m.wave gt 5450 and m.wave lt 5550)

for zz = 0, 5 do begin
   bc03 = im_read_bc03(isedpath=isedpath, metallicity=zz, age=age,$
                       minwave=3400, maxwave=8500, bc03_extras=ssp_extras)

; Interpolate to log-lambda & normalize
   for ii = 0, n_elements(age) - 1 do begin
      linterp, bc03.wave, bc03.flux[*,ii], m.wave, spec
      
      print, min(spec), max(spec), mean(spec), age[ii]
      plot, m.wave, spec, xtitle='Wavelength ( '+angstrom+' )', $
            ytitle = string('F!D!X',lambda,'!N (L',sunsymbol(),$
                            '/M',sunsymbol(),' ',angstrom,')')

      xyouts, 0.5, 0.98, string(age[ii],' Gyr',format='(F7.3,A4)'), /norm
      hline, median(spec[nwl])
      vline, m.wave[nwl[0]]
      vline, m.wave[nwl[-1]]
      m.norm[ii+(n_elements(age) * zz)] = median(spec[nwl])
      spec = spec / median(spec[nwl]) ; normalize at 5500 A
      m.flux[*,ii+(n_elements(age) * zz)] = spec
      m.Z[ii + (n_elements(age) * zz)] = Z_arr[zz]
      junk = min(abs(bc03.age/1e9 - age[ii]), idx)
      m.m_remaining[ii+(n_elements(age) * zz)] = ssp_extras[idx].m_
   endfor
endfor

print, m.norm
mwrfits, m, output, /create
if keyword_set(plot) then dfpsclose

end
