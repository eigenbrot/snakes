;

function mconv, y, sig
  d1 = n_elements(y)
  x = findgen(d1)
  yp = dblarr(d1)
  norm = sqrt(2.*!DPI)*sig

  for i = 0, d1 - 1 do begin
     div = ((x-i)/sig)^2.
     close = where(div lt 50.)
     kern = exp(-0.5*div[close])/norm[close]
     yp[i] = total(y[close]*kern)/total(kern)
  endfor
  
  return, yp
end

pro bc_dfk_models, dfkfile, output, dispdata, Z, plot=plot

if keyword_set(plot) then begin
   plotfile = (strsplit(output,'.',/extract))[0] + '.ps'
   dfpsplot, plotfile, /color, /times, /landscape
endif

angstrom = '!6!sA!r!u!9 %!6!n'
lambda = 'l'
;lambda = '!4'+string("153B)+'!X'

readcol, dfkfile, ages, weights, dfkgrp

isedpath = '/usr/users/tremonti/models/Padova1994/chabrier/'
numdfk = 4

Z_arr = [0.0001, 0.0004, 0.004, 0.008, 0.02, 0.05]/0.02
numZ = n_elements(Z)

; Set up log-lambda wavelength grid
npix = round((alog10(8500) - alog10(3400)) / 1e-4)
logwl = findgen(npix) * 1e-4 + alog10(3400)

; Define output structure
m = {wave: 10.0^logwl, flux: fltarr(5, npix, numdfk * numZ), $
     age: fltarr(5, numdfk * numZ), Z: fltarr(5, numdfk * numZ), $
     id: string('multi metallicity, Chabrier IMF',format='(A14)'), $
     norm: fltarr(5, numdfk * numZ)}

; normalization wavelengths
nwl = where(m.wave gt 5450 and m.wave lt 5550)

; Get dispersion information
v_data = MRDFITS(dispdata,0,v_header)
numdisp = n_elements(v_data[0,*])
v_wavesize = n_elements(v_data[*,0])
v_cdelt = float(sxpar(v_header,'CDELT1'))
v_crval = float(sxpar(v_header,'CRVAL1'))
v_crpix = float(sxpar(v_header,'CRPIX1'))
print, 'V_CDELT1 = ',v_cdelt
print, 'V_CRVAL1 = ',v_crval
print, 'V_CRPIX1 = ',v_crpix
v_wave = (FINDGEN(v_wavesize) - v_crpix) * v_cdelt + v_crval
vdisp_vec = dblarr(n_elements(m.wave),numdisp)
for dd = 0, numdisp - 1 do $
   vdisp_vec[*,dd] = interpol(v_data[*,dd],v_wave,m.wave)/2.355

bc03_pix = 70.0 ; size of 1 model pixel in km/s 
bc03_vdisp = 75.0 ; approximate velocity dispersion of BC03 models

for zz = 0, numZ - 1 do begin
   for ff = 0, 4 do begin
      for gg = 1, numdfk do begin
         zid = where(Z[zz]*0.02/0.02 eq Z_arr)
         print, 'Zid: ', zid

         bc03data = im_read_bc03(isedpath=isedpath, metallicity=zid, age=ages/1e9,$
                                 minwave=3400, maxwave=8500)
         
         nwave = (size(bc03data.flux, /dimensions))[0]
         nssp = (size(bc03data.flux, /dimensions))[1]
         print, 'found', nssp, ' SSPs'
         print, nwave, nssp
      
         idx = where(dfkgrp eq gg)
         tw = transpose(rebin(weights[idx], n_elements(idx), nwave))
         ta = transpose(rebin(ages[idx], n_elements(idx), nwave))
         
         print, size(bc03data.flux[*,idx], /dimensions)
         print, size(tw, /dimensions)
         
         tdata = total(bc03data.flux[*,idx]*tw, 2)/total(tw,2)
         dfkage = total(tw*ta)/total(tw)
         
         linterp, bc03data.wave, tdata, m.wave, tspec
         
         vdisp_add = sqrt(vdisp_vec[*,ff]^2 - bc03_vdisp^2)  
         sigma_pix = vdisp_add / bc03_pix
         spec = mconv(tspec,sigma_pix)
         
         plot, m.wave, spec, xtitle='Wavelength ( '+angstrom+' )', $
               ytitle = string('F!D!X',lambda,'!N (L',sunsymbol(),$
                               '/M',sunsymbol(),' ',angstrom,')')
         
         xyouts, 0.5, 0.98, string(dfkage/1e9,' Gyr, ', Z[zz], ' Z/Z_sol, ',$
                                   ff+2,"'' dfkgrp ",gg,$
                                   format='(F7.3,A6,F7.3,A10,I1,A10,I1)'), /norm
         
         hline, median(spec[nwl])
         vline, m.wave[nwl[0]]
         vline, m.wave[nwl[-1]]
         
         m.age[ff, gg-1 + zz*numdfk] = dfkage
         m.norm[ff, gg-1 + zz*numdfk] = median(spec[nwl])
         m.Z[ff,gg-1 + zz*numdfk] = Z[zz]
         m.flux[ff,*,gg-1 + zz*numdfk] = spec / median(spec[nwl])
         
      endfor
   endfor
endfor
mwrfits, m, output, /create
if keyword_set(plot) then dfpsclose

end
