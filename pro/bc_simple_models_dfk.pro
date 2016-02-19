; Write out BC03 model grid used in SDSS code re-binned to log-lambda

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

pro bc_simple_models_dfk, output, dispdata, plot=plot, constdisp=constdisp

if keyword_set(plot) then begin
   plotfile = (strsplit(output,'.',/extract))[0] + '.ps'
   dfpsplot, plotfile, /color, /times, /landscape
endif

angstrom = '!6!sA!r!u!9 %!6!n'
lambda = 'l'
;lambda = '!4'+string("153B)+'!X'
; Ages in Gyr

Z_arr = [1.0]
age=[0.0025, 0.0641, 2.5, 9.75]

isedpath = '/d/monk/eigenbrot/WIYN/14B-0456/anal/DFK/models/DFK_models.txt'

; Set up log-lambda wavelength grid
npix = round((alog10(8500) - alog10(3400)) / 1e-4)
logwl = findgen(npix) * 1e-4 + alog10(3400)

; Define output structure
m = {wave: 10.0^logwl, flux: fltarr(5, npix, n_elements(age) * n_elements(Z_arr)), $
     age: fltarr(5, n_elements(age) * n_elements(Z_arr)), Z: fltarr(5, n_elements(age) * n_elements(Z_arr)), $
     id: string('DFK basis set, solar metallicity'), $
     norm: fltarr(5, n_elements(age) * n_elements(Z_arr))}

; normalization wavelengths
nwl = where(m.wave gt 5450 and m.wave lt 5550)

; Get dispersion information
if n_elements(constdisp) eq 0 then begin
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
endif

bc03_pix = 70.0 ; size of 1 model pixel in km/s 
bc03_vdisp = 75.0 ; approximate velocity dispersion of BC03 models

for ff = 0, 4 do begin
   m.age[ff,*] = age*1e9

   for zz = 0, 0 do begin
      readcol, isedpath, dwave, Y, I1, I2, O
      da = [[Y],[I1],[I2],[O]]
      
    ; Interpolate to log-lambda & normalize
      for ii = 0, n_elements(age) - 1 do begin
         linterp, dwave, da[*,ii], m.wave, tspec
      
         ;Deconvolve template instrumental resolution
         print, 'Convolving...'
         if n_elements(constdisp) eq 0 then begin
            vdisp_add = sqrt(vdisp_vec[*,ff]^2 - bc03_vdisp^2)  
            sigma_pix = vdisp_add / bc03_pix
            spec = mconv(tspec,sigma_pix)
         endif else begin
            dispdiff = (constdisp[ff])^2 - bc03_vdisp^2
            if dispdiff gt 0 then begin
               vdisp_add = sqrt((constdisp[ff])^2 - bc03_vdisp^2)
               sigma_pix = vdisp_add / bc03_pix
               spec = gconv(tspec, sigma_pix)
            endif else begin
               print, 'Skipping ff=',ff,', disp of ',constdisp[ff]/2.355,' is too low'
               spec = tspec
            endelse
         endelse

         print, min(spec), max(spec), mean(spec), age[ii]
         plot, m.wave, spec, xtitle='Wavelength ( '+angstrom+' )', $
               ytitle = string('F!D!X',lambda,'!N (L',sunsymbol(),$
                               '/M',sunsymbol(),' ',angstrom,')')

         xyouts, 0.5, 0.98, string(age[ii],' Gyr, ', Z_arr[zz], ' Z/Z_sol, ',ff+2,"''",$
                                   format='(F7.3,A6,F7.3,A10,I1,A2)'), /norm
         hline, median(spec[nwl])
         vline, m.wave[nwl[0]]
         vline, m.wave[nwl[-1]]
         m.norm[ff, ii+(n_elements(age) * zz)] = median(spec[nwl])
         spec = spec / median(spec[nwl]) ; normalize at 5500 A
         m.flux[ff,*,ii+(n_elements(age) * zz)] = spec
         m.Z[ff, ii + (n_elements(age) * zz)] = Z_arr[zz]
         
      endfor
   endfor
endfor
print, m.norm
mwrfits, m, output, /create
if keyword_set(plot) then dfpsclose

end
