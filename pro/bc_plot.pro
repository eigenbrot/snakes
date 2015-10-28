pro bc_plot, coeffile, fitfile, datafile, errorfile, model, location=location

wavemin=3800.
wavemax=6800.
flux_factor=1d17

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

idx = where(wave ge wavemin and wave le wavemax)
restwl = wave[idx]

if keyword_set(location) then begin
   readcol, location, apnum, fiber_radii, ra, dec, rkpc, zkpc
   sizeidx = [0.937,1.406,1.875,2.344,2.812]
endif

m = mrdfits(model,1)
nmodels = (size(m.flux,/dimensions))[1]

plotfile = (strsplit(fitfile,'.',/extract))[0] + '.ps'
dfpsplot, plotfile, /color, /times, /landscape

coef_arr = mrdfits(coeffile,1)
yfits = mrdfits(fitfile,0)

vdisp = [493., 589., 691., 796., 966.]/2.355
size_borders = [19, 43, 62, 87, 109] ; The last one is needed to prevent indexing errors
size_switch = 0

redidx = where(restwl ge 5400)
blueidx = where(restwl lt 5400)
hklow = 3920
hkhigh = 4000
hkidx = where(restwl gt hklow and restwl lt hkhigh)
npix = n_elements(restwl)

defplotcolors
smoothkern = 5
thick=1.0

for i = 0, numfibers - 1 DO BEGIN

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

   quality = fltarr(npix) + 1

   bad = where(finite(flux) ne 1 or finite(err) ne 1 or err eq 0)
   if bad[0] ne -1 then quality[bad] = 0

   outside_model = where(restwl le min(m.wave) or restwl ge max(m.wave))
   if outside_model[0] ne -1 then quality[outside_model] = 0

   sk =    [6300.,        5890., 5683.8, 5577.,      5461., 5199.,      4983., 4827.32, 4665.69, 4420.23, 4358., 4165.68, 4047.0]
   sknam = ['[OI] (atm)', 'NaD', 'NaI',  'OI (atm)', 'HgI', 'NI (atm)', 'NaI', 'HgI',   'NaI',   'NaI',   'HgI', 'NaI',   'HgI']
   
   sk2 = [6300., 5890., 5577.]
   
;; em=     [3727.3,  4959.,    5006.8,   6563.8, 6716.0]
;; emnam = ['[OII]', '[OIII]', '[OIII]', 'Ha',   'S2']

   em2 = [6563.8,  6716.0, 6583.41, 6548.04]

   em = [6563.8,  6716.0, 6583.41, 6548.04, 4959., 5006.8]
   emnam = ['Ha', 'S2', 'NII', 'NII', '[OIII]', '[OIII]']

   abs =    [3933.7, 3968.5, 4304.4,   4341., 5175.3, 5894.0, 4861.,  4102., 3820.4]
   absnam = ['H',    'K',    'G band', 'HG',  'Mg',   'Na',   'HB',   'HD',  'L band']
;sk = [5569., 5882.6]
   HPS = 5914.
   HPS_wid = 230.
   
   dz = 1500. / 3e5           ; clipping interval
   dzsk = 1600. / 3e5
   
   for ii = 0, n_elements(em2) - 1 do begin 
      maskout = where(restwl gt em2[ii]*(1-dz) and restwl lt em2[ii]*(1+dz))
      if maskout[0] ne -1 then quality[maskout] = 0
   endfor
   
   for ii = 0, n_elements(sk2) - 1 do begin 
      maskout = where(restwl gt sk2[ii]*(1-dzsk) and restwl lt sk2[ii]*(1+dzsk))
      if maskout[0] ne -1 then quality[maskout] = 0
   endfor
   
   ok = where(quality eq 1)

   ;-----------------------------------------------------------------------------
   ; Convolve models to velocity dispersion of data and interpolate to
   ; match data

   bc03_pix = 70.0              ; size of 1 model pixel in km/s 
   bc03_vdisp = 75.0            ; approximate velocity dispersion of BC03 models
   
   ;Deconvolve template instrumental resolution
   if vd lt bc03_vdisp then vdisp_add = 0 $
   else vdisp_add = sqrt(vd^2 - bc03_vdisp^2)  
   sigma_pix = vdisp_add / bc03_pix
  
   custom_lib = dblarr(npix, nmodels)
   for ii = 0, nmodels - 1 do begin
      cflux = gconv(m.flux[*,ii], sigma_pix) ; convolve with gaussian
      custom_lib[*,ii] = interpol(cflux, m.wave, restwl)
   endfor
   if outside_model[0] ne -1 then custom_lib[outside_model, *] = 0.0

   yfit = yfits[*,i]*flux_factor
   coefs = coef_arr[i]
   
   blueymax = max(yfit[blueidx]) / 0.8
   ymax = max(yfit) * 1.1
   ymax = max([ymax,blueymax,max(flux)*1.1])
   xmin = min(restwl) * 0.98
   xmax = max(restwl) * 1.02
;xtitle='Wavelength (Angstroms)'
   plot, restwl, alog10(flux), xtickformat='(A1)', /nodata,$
         ytitle = 'Log Flux + 17 [erg/s/cm^2/A]', yrange = [-0.5, 2.6], xrange = [xmin,xmax], $
         position = [0.15,0.3,0.95,0.99], charsize=1.0, charthick=1.0, /xs, /ys, /t3d
   
   vline, 5400., color=!gray, linestyle=2
   vline, hklow, color=!gray, linestyle=2
   vline, hkhigh, color=!gray, linestyle=2
   
   for j=0, n_elements(coefs.light_frac) - 1 do begin
      xred = restwl * (coefs.vsys/3e5 + 1)
      yi = coefs.light_frac[j] * custom_lib[*,j] * $
           exp(-coefs.tauv*(restwl/5500.0)^(-0.7))
      yi = interpol(yi,xred,restwl)
      oplot, restwl, alog10(smooth(yi,smoothkern)), $
             color = !lblue, thick=thick, linestyle=0, /t3d
      ;; xyouts, 0.2, 0.88 - 0.02*j, string('f_',j,' = ',$
      ;;                         mean(yi[lightidx]),$
      ;;                         format='(A2,I02,A3,E10.3)'),$
      ;;        charsize=0.6, /norm, /t3d
   endfor
   
   errsmooth = 5
   oband, restwl[0:*:errsmooth], alog10((smooth(flux-err,smoothkern))[0:*:errsmooth]), $
          alog10((smooth(flux+err,smoothkern))[0:*:errsmooth]), color=!gray, /t3d, /noclip
   
; Show masked regions in green
   galfit = flux  + 'NaN'
   galfit[ok] = flux[ok]
   masked = flux
   masked[ok] = 'NaN'
   oplot, restwl, alog10(smooth(galfit,smoothkern,/NAN)), thick=thick, color=!black, /t3d
   oplot, restwl, alog10(smooth(masked,smoothkern,/NAN)), color=!cyan, thick=thick*4, /t3d
   
   oplot, restwl, alog10(smooth(yfit,smoothkern)), color = !red, thick=thick, /t3d
      
   for s=0, n_elements(sk) - 1 do begin
      tidx = where(restwl ge sk[s] - 10. and restwl le sk[s] + 10.)
      ypos = alog10(max([max(flux[tidx], /NAN), max(yfit[tidx], /NAN)], /NAN)) + 0.1
;      ypos = alog10(interpol(flux, restwl, sk[s]))*1.03
      xyouts, sk[s], ypos, sknam[s], alignment=0.5, charsize=0.5, /data
   endfor

   for e=0, n_elements(em) - 1 do begin
      tidx = where(restwl ge em[e] - 10. and restwl le em[e] + 10.)
      ypos = alog10(max([max(flux[tidx], /NAN), max(yfit[tidx], /NAN)], /NAN)) + 0.1
;      ypos = alog10(interpol(flux, restwl, em[e]))*1.03
      xyouts, em[e]*(coefs.vsys/3e5 + 1), ypos, emnam[e], alignment=0.5, charsize=0.5, /data, color=!blue
   endfor
   
   prevy = 99.
   for a=0, n_elements(abs) - 1 do begin
      tidx = where(restwl ge abs[a] - 10. and restwl le abs[a] + 10.)
      ypos = alog10(min([min(abs(flux[tidx]), /NAN), min(abs(yfit[tidx]), /NAN)], /NAN)) - 0.1
      if a eq 3 and abs(ypos - prevy) le 0.05 then ypos -= 0.04
      prevy = ypos
;      ypos = alog10(interpol(flux, restwl, abs[a]))*0.9
      xyouts, abs[a]*(coefs.vsys/3e5 + 1), ypos, absnam[a], alignment=0.5, charsize=0.5, /data, color=!red
   endfor
   
   chivec = (galfit - yfit)/err
   plot, restwl, smooth(chivec,smoothkern,/NAN), xtitle='Wavelength (Angstroms)', $
         ytitle='Residuals/error', $
         position=[0.15,0.15,0.95,0.3], xrange=[xmin,xmax], yrange=[-6,6], $
         yminor=1, yticks=2, charsize=1.0, charthick=1.0, thick=thick, $
         /xs, /ys, /noerase, /t3d;, ytickv=[-200,0,200]
   if keyword_set(plotlabel) then $
      xyouts, 0.2, 0.96, plotlabel, color = !black, /norm, /t3d
   
   xyouts, 0.2, 0.94, 'SNR = ' + string(coefs.SNR, format = '(F8.2)'), $
           /norm, /t3d
   xyouts, 0.2, 0.92, 'V = ' + string(coefs.vsys, format = '(F8.2)') + ' km/s', $
           /norm, /t3d
   xyouts, 0.2, 0.90, 'V_disp = ' + string(vd, format = '(F8.2)'), $
           /norm, /t3d
   xyouts, 0.2, 0.88, 'Tau V = ' + string(coefs.tauv, format = '(F5.2)'), $
           /norm, /t3d

   xyouts, 0.4, 0.94, 'MLWA = ' + string(coefs.MLWA, format = '(F8.2)'), $
           /norm, /t3d
   xyouts, 0.4, 0.92, 'MMWA = ' + string(coefs.MMWA, format = '(F8.2)'), $
           /norm, /t3d
   xyouts, 0.4, 0.90, 'MLWZ = ' + string(coefs.MLWZ, format = '(F8.2)') + '(Z/Z_sol)', $
           /norm, /t3d
   xyouts, 0.4, 0.88, 'MMWZ = ' + string(coefs.MMWZ, format = '(F8.2)'), $
           /norm, /t3d

   xyouts, 0.6, 0.94, 'Chisq = ' + string(coefs.chisq, format = '(F8.2)'), $
           /norm, /t3d
   xyouts, 0.6, 0.92, 'redChi = ' + string(coefs.redchi, format = '(F8.2)'), $
           /norm, /t3d
   xyouts, 0.6, 0.90, 'bluChi = ' + string(coefs.bluechi, format = '(F8.2)'), $
           /norm, /t3d
   xyouts, 0.6, 0.88, 'HK_Chi = ' + string(coefs.hkchi, format = '(F8.2)'), $
           /norm, /t3d

   xyouts, 0.5, 0.05, systime(), alignment=0.5, /norm, /t3d
   

ENDFOR

dfpsclose

end
