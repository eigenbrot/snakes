pro brightSurf, image, num_ap, ap_vec, fluxes,$
              NOREAD=noread, DRAW=draw, OUTPUT=output

on_error, 2

;Read in the data image
IF (NOT KEYWORD_SET(NOREAD)) THEN data = MRDFITS(image) ELSE data = image

;Find the center of mass of the image. Hopefully this is the center of the
;pin-hole/fiber extended source
center = centroid(data)  ;[x_cent,y_cent]

dims = SIZE(data, /DIMENSIONS)
rlimit = MIN(dims)
Alimit = (!PI*rlimit^2)/num_ap

;Create an array that has the distance of each pixel from the center
distances = DBLARR(dims)
FOR j=0, dims[0] - 1 DO BEGIN
   FOR k=0, dims[1] - 1 DO BEGIN
      distances[j,k] = SQRT((center[0] - j)^2 + (center[1] - k)^2)
   ENDFOR
ENDFOR

fluxes = []
r1 = 0.0
r2 = SQRT(Alimit/!PI)
ap_vec = []

FOR i=0, num_ap - 1 DO BEGIN
   pixels = data[WHERE((distances GT r1) AND (distances LE r2), /NULL)]
   fluxes = [fluxes,TOTAL(pixels)]
   ap_vec = [ap_vec,SQRT(0.5*(r2^2 + r1^2))]
   IF N_ELEMENTS(pixels) LE Alimit/2 THEN BREAK
   r3 = SQRT(2*r2^2 - r1^2)
   r1 = r2
   r2 = r3
ENDFOR


IF KEYWORD_SET(DRAW) THEN PLOT, ap_vec, fluxes
IF KEYWORD_SET(OUTPUT) THEN FORPRINT, ap_vec, fluxes,  TEXTOUT=output, /NOCOMMENT

end
