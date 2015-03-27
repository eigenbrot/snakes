pro annulize, image, num_an, r_vec, fluxes,$
              NOREAD=noread, DRAW=draw, OUTPUT=output,$
              EXTEN=exten
;annulize takes in a FITS image and lays down a user-specified amount of
;annuli with equal area. It then creates a file containing total counts per
;annulus as a function of annulus radius.

;on_error, 2

;Read in the data image
IF NOT KEYWORD_SET(NOREAD) THEN BEGIN
   IF KEYWORD_SET(EXTEN) THEN data=MRDFITS(image,exten) ELSE data=MRDFITS(image)
ENDIF ELSE data = image

;Find the center of mass of the image. Hopefully this is the center of the
;pin-hole/fiber extended source
center = centroid(data)  ;[x_cent,y_cent]

print, center

;We want our annuli to all have the same area
dims = SIZE(data, /DIMENSIONS)
;rlimit = MAX(dims)/2
;Alimit = !PI*rlimit^2
;area = Alimit/num_an

;Create an array that has the distance of each pixel from the center
distances = DBLARR(dims)
FOR j=0, dims[0] - 1 DO BEGIN
   FOR k=0, dims[1] - 1 DO BEGIN
      distances[j,k] = SQRT((center[0] - j)^2 + (center[1] - k)^2)
   ENDFOR
ENDFOR

rlimit = MAX(distances)
print, rlimit
Alimit = !PI*rlimit^2
area = Alimit/num_an

area_array = INTARR(dims) + 1

fluxes = DBLARR(num_an)
r1 = 0.0
r2 = SQRT(area/!PI)
r_vec = DBLARR(num_an)
FOR i=0, num_an - 1 DO BEGIN
   idxs = WHERE((distances LE r2) AND (distances GT r1))
   flux = TOTAL(data[idxs])
   an_area = TOTAL(area_array[idxs])
   print, i,r1,r2,flux,an_area
   fluxes[i] = flux/an_area
   r_mid = SQRT(0.5*(r2^2 + r1^2)) ;middle radius based on area
   r_vec[i] = r_mid
   r1 = r2
   r2 = SQRT((area/!PI)+r1^2)
ENDFOR


IF KEYWORD_SET(DRAW) THEN PLOT, r_vec, fluxes
IF KEYWORD_SET(OUTPUT) THEN FORPRINT, r_vec, fluxes,  TEXTOUT=output, /NOCOMMENT

end
