function simple_make_parinfo, coeffile, output, fixweights=fixweights, fixtau=fixtau, zeroweights=zeroweights, zerotau=zerotau, fixvec=fixvec

light_factor = 100.
vel_factor = 100.

c = mrdfits(coeffile,1)
dims = size(c.light_frac, /dimensions)
nap = dims[1]
npars = dims[0] + 2

weights = c.light_frac/light_factor

parinfo = replicate({value:1.D, fixed:0, limited:[0,0], tied:'', $
                    limits:[0.0,0], step:0, relstep:0}, nap, npars)

;set all the proper values
parinfo[*,0].value = c.VSYS/vel_factor
parinfo[*,1].value = c.TAUV
parinfo[*,1].limited = [1,1]
parinfo[*,1].limits = [0., 20.]
for i = 0, nap - 1 do begin
   for j = 0, npars - 3 do begin 
      parinfo[i,j+2].value = weights[j,i]
   endfor
endfor
parinfo[*,2:*].limited = [1,0]
parinfo[*,2:*].limits = [0., 0.]

;Always fix V
parinfo[*,0].fixed = 1

if keyword_set(fixtau) then parinfo[*,1].fixed=1
if keyword_set(fixweights) then parinfo[*,2:*].fixed=1

if keyword_set(zerotau) then parinfo[*,1].value = 0.0
if keyword_set(zeroweights) then parinfo[*,2:*].value = 0.0

if keyword_set(fixvec) then begin
   if n_elements(fixvec) eq npars then begin
      for a = 0, nap - 1 do begin
         for p = 0, npars - 1 do begin
            parinfo[a,p].fixed = fixvec[p]
         endfor
      endfor
   endif else begin
      print, 'Oops, fixvec had length ',n_elements(fixvec),' and it probably should have been ', npars
      stop
   endelse
endif

mwrfits, parinfo[0,*], output, /create
for a = 1, nap - 1 do begin
   mwrfits, parinfo[a,*], output
endfor

return, parinfo

end
