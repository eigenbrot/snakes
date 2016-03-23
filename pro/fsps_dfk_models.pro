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

pro fsps_dfk_models, output, dispdata, Z, plot=plot

if keyword_set(plot) then begin
   plotfile = (strsplit(output,'.',/extract))[0] + '.ps'
   dfpsplot, plotfile, /color, /times, /landscape
endif

angstrom = '!6!sA!r!u!9 %!6!n'
lambda = 'l'
;lambda = '!4'+string("153B)+'!X'

;readcol, dfkfile, ages, weights, dfkgrp
numdfk = 4
dfklims = [[0.9,5.2],$
           [5.5,404],$
           [456.5,5750],$
           [6000,13500]]/1e3


modelpaths = ['/d/monk/eigenbrot/WIYN/14B-0456/anal/DFK/FSPS/FSPS_models/raw/fsps_MILES_cha00190z.ssp',$
              '/d/monk/eigenbrot/WIYN/14B-0456/anal/DFK/FSPS/FSPS_models/raw/fsps_MILES_cha00077z.ssp']

Z_arr = [1,0.4]
numZ = n_elements(Z)

; Set up log-lambda wavelength grid
npix = round((alog10(8500) - alog10(3400)) / 1e-4)
logwl = findgen(npix) * 1e-4 + alog10(3400)

; Define output structure
m = {wave: 10.0^logwl, flux: fltarr(5, npix, numdfk * numZ), $
     age: fltarr(5, numdfk * numZ), Z: fltarr(5, numdfk * numZ), $
     id: string('FSPS models, MILES, Chabrier IMF',format='(A30)'), $
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

AA_pix = 0.9 ;size of pixel, in Angstroms. This doesn't matter
AA_vdisp = 2.54/2.355 ;resolution (sigma) in Angstroms
fsps_vdisp = AA_vdisp/m.wave*3e5 ;resolution in km/s
final_pix = 70.0 ;Spacing, in km/s of final model grid

for zz = 0, numZ - 1 do begin
   
   zid = where(Z[zz]*0.02/0.02 eq Z_arr)
   print, 'Zid: ', zid
   print, 'Reading ',modelpaths[zid]
   
   readcol, string(modelpaths[zid],format='(A)'), Mage, MZ, Mwave, Mflux
   
   ;Put us in L_sol units
   ;Mflux = Mflux / 3.826e33

   ;Compute weights
   uniA = Mage[uniq(Mage, sort(Mage))]
   nssp = n_elements(uniA)
   weights = findgen(nssp) * 0.0
   weights[0] = 0.5*(uniA[0] + uniA[1])
   weights[-1] = 0.5*(uniA[-1] - uniA[-2])
   for ww = 1, nssp - 2 do begin
      weights[ww] = 0.5*(uniA[ww+1] - uniA[ww-1])
   endfor

   ;Get singular wavelength vector
   wid = where(Mage eq Mage[0])
   Fwave = Mwave[wid]
   nwave = n_elements(Fwave)

   ;fill flux array
   flux = findgen(nwave,nssp)*0.0

   for aa = 0, nssp - 1 do begin
      aid = where(Mage eq uniA[aa])
      flux[*,aa] = Mflux[aid]
   endfor

   print, 'found', nssp, ' SSPs'
   print, nwave, nssp
         
   for ff = 0, 4 do begin
      for gg = 0, numdfk - 1 do begin
         gidx = where(uniA ge dfklims[0,gg] and uniA le dfklims[1,gg])
         if gidx[0] eq -1 then continue

         print, uniA[gidx]
         tw = transpose(rebin(weights[gidx], n_elements(gidx), nwave))
         ta = transpose(rebin(uniA[gidx], n_elements(gidx), nwave))
         
         print, size(flux[*,gidx], /dimensions)
         print, size(tw, /dimensions)
         
         tdata = total(flux[*,gidx]*tw, 2)/total(tw,2)
         dfkage = total(tw*ta)/total(tw)
         
         linterp, Fwave, tdata, m.wave, tspec
         
         print, n_elements(vdisp_vec[*,ff]), n_elements(fsps_vdisp)
         vdisp_add = sqrt(vdisp_vec[*,ff]^2 - fsps_vdisp^2)  
         sigma_pix = vdisp_add / final_pix
         spec = mconv(tspec,sigma_pix)
         
         plot, m.wave, spec, xtitle='Wavelength ( '+angstrom+' )', $
               ytitle = string('F!D!X',lambda,'!N (L',sunsymbol(),$
                               '/M',sunsymbol(),' ',angstrom,')')
         
         xyouts, 0.5, 0.98, string(dfkage,' Gyr, ', Z[zz], ' Z/Z_sol, ',$
                                   ff+2,"'' dfkgrp ",gg,$
                                   format='(F7.3,A6,F7.3,A10,I1,A10,I1)'), /norm
         
         hline, median(spec[nwl])
         vline, m.wave[nwl[0]]
         vline, m.wave[nwl[-1]]
         
         m.age[ff, gg + zz*numdfk] = dfkage*1e9
         m.norm[ff, gg + zz*numdfk] = median(spec[nwl])
         m.Z[ff,gg + zz*numdfk] = Z[zz]
         m.flux[ff,*,gg + zz*numdfk] = spec / median(spec[nwl])
         
      endfor
   endfor
endfor
mwrfits, m, output, /create
if keyword_set(plot) then dfpsclose

end
