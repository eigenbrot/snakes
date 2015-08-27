

; Write out BC03 model grid used in SDSS code re-binned to log-lambda

pro bc_simple_models_sky, output, skyfile, metallicity=metallicity, plot=plot, age=age

if keyword_set(plot) then begin
   plotfile = (strsplit(output,'.',/extract))[0] + '.ps'
   dfpsplot, plotfile, /color, /times, /landscape
endif

if n_elements(metallicity) eq 0 then metallicity = 4

angstrom = '!6!sA!r!u!9 %!6!n'
lambda = 'l'
;lambda = '!4'+string("153B)+'!X'
; Ages in Gyr

if n_elements(age) eq 0 then $
   age=[0.005, 0.025, 0.100, 0.286, 0.640, 0.904, 1.434, 2.500, 5.000, 10.000]

isedpath = '/usr/users/tremonti/models/Padova1994/chabrier/'

bc03 = im_read_bc03(isedpath=isedpath, metallicity=metallicity, age=age,$
                    minwave=3400, maxwave=8500, bc03_extras=ssp_extras)

; Set up log-lambda wavelength grid
npix = round((alog10(8500) - alog10(3400)) / 1e-4)
logwl = findgen(npix) * 1e-4 + alog10(3400)

; Define output structure
m = {wave: 10.0^logwl, flux: fltarr(npix, n_elements(age)+1), age: bc03.age, $
     id: string('Z = ',metallicity,', Chabrier IMF',format='(A4,I1,A14)'), $
     norm: fltarr(n_elements(age)), m_remaining: fltarr(n_elements(age))}

; normalization wavelengths
nwl = where(m.wave gt 5450 and m.wave lt 5550)

; Read in the sky spectrum
sdata = mrdfits(skyfile,0,shead)
for ss = 0, n_elements(shead) - 1 do begin
   if strpos(shead[ss], 'CDELT1') GE 0 then begin 
      cdelt = float((strsplit(shead[ss],'=',/extract))[1])
      print, 'FOUND CDELT1'
   endif
   if strpos(shead[ss], 'CRVAL1') GE 0 then begin
      crval = float((strsplit(shead[ss],'=',/extract))[1])
      print, 'FOUND CRVAL'
   endif
endfor

swave = findgen(n_elements(sdata))*cdelt + crval
linterp, swave, sdata, m.wave, sspec
plot, m.wave, sspec, xtitle='Wavelength ( '+angstrom+' )', $
      ytitle = string('F!D!X',lambda,'!N (L',sunsymbol(),$
                      '/M',sunsymbol(),' ',angstrom,')')
xyouts, 0.5, 0.98, 'SKY', /norm

hline, median(sspec[nwl])
vline, m.wave[nwl[0]]
vline, m.wave[nwl[-1]]
sspec = sspec / median(sspec[nwl])
m.flux[*,0] = sspec


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
   m.norm[ii] = median(spec[nwl])
   spec = spec / median(spec[nwl]) ; normalize at 5500 A
   m.flux[*,ii+1] = spec
   junk = min(abs(bc03.age/1e9 - age[ii]), idx)
   m.m_remaining[ii] = ssp_extras[idx].m_
endfor

print, m.norm
mwrfits, m, output, /create
if keyword_set(plot) then dfpsclose

end
