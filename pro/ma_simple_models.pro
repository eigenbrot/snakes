pro ma_simple_models, output, Z, plot=plot

if keyword_set(plot) then begin
   plotfile = (strsplit(output,'.',/extract))[0] + '.ps'
   dfpsplot, plotfile, /color, /times, /landscape
endif

angstrom = '!6!sA!r!u!9 %!6!n'
lambda = 'l'
;lambda = '!4'+string("153B)+'!X'

modelpaths = ['/d/monk/eigenbrot/WIYN/14B-0456/anal/DFK/MA11/MA11_models/raw/ssp_M11_MILES.chaz10m4.rhb',$
              '/d/monk/eigenbrot/WIYN/14B-0456/anal/DFK/MA11/MA11_models/raw/ssp_M11_MILES.chaz0001.rhb',$
              '/d/monk/eigenbrot/WIYN/14B-0456/anal/DFK/MA11/MA11_models/raw/ssp_M11_MILES.chaz001',$
              '/d/monk/eigenbrot/WIYN/14B-0456/anal/DFK/MA11/MA11_models/raw/ssp_M11_MILES_revisednearIRslope.chaz002',$
              '/d/monk/eigenbrot/WIYN/14B-0456/anal/DFK/MA11/MA11_models/raw/ssp_M11_MILES.chaz004']

Z_arr = [0.0001, 0.001, 0.01, 0.02, 0.04]/0.02

; Set up log-lambda wavelength grid
npix = round((alog10(8500) - alog10(3400)) / 1e-4)
logwl = findgen(npix) * 1e-4 + alog10(3400)

zid = where(Z*0.02/0.02 eq Z_arr)
print, 'Zid: ', zid
print, 'Reading ',modelpaths[zid]
readcol, string(modelpaths[zid],format='(A)'), Mage, MZ, Mwave, Mflux

Mage *= 1e9

;Put us in L_sol units
Mflux = Mflux / 3.826e33

uniA = Mage[uniq(Mage, sort(Mage))]
nssp = n_elements(uniA)

; Define output structure
m = {wave: 10.0^logwl, flux: fltarr(npix, nssp), $
     age: uniA, Z: Z, $
     id: string('Ma11 models, MILES, Chabrier IMF',format='(A30)'), $
     norm: fltarr(nssp)}

; normalization wavelengths
nwl = where(m.wave gt 5450 and m.wave lt 5550)

;Get singular wavelength vector
wid = where(Mage eq Mage[0])
Fwave = Mwave[wid]
nwave = n_elements(Fwave)

for ii = 0, nssp - 1 do begin
   aid = where(Mage eq uniA[ii])
   linterp, Fwave, Mflux[aid], m.wave, spec
   plot, m.wave, spec, xtitle='Wavelength ( '+angstrom+' )', $
         ytitle = string('F!D!X',lambda,'!N (L',sunsymbol(),$
                         '/M',sunsymbol(),' ',angstrom,')')
   xyouts, 0.5, 0.98, string(uniA[ii]/1e9,' Gyr',format='(F7.3,A4)'), /norm
   hline, median(spec[nwl])
   vline, m.wave[nwl[0]]
   vline, m.wave[nwl[-1]]
   m.norm[ii] = median(spec[nwl])
   spec = spec / median(spec[nwl]) ; normalize at 5500 A
   m.flux[*,ii] = spec
endfor

mwrfits, m, output, /create
if keyword_set(plot) then dfpsclose

end
