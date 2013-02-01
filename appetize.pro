pro appetize, image, num_ap, ap_vec, fluxes,$
              NOREAD=noread, DRAW=draw, OUTPUT=output,$
              EXTEN=exten

on_error, 2

;Read in the data image
IF NOT KEYWORD_SET(NOREAD) THEN BEGIN
   IF KEYWORD_SET(EXTEN) THEN data=MRDFITS(image,exten) ELSE data=MRDFITS(image)
ENDIF ELSE data = image

;Find the center of mass of the image. Hopefully this is the center of the
;pin-hole/fiber extended source
center = centroid(data)  ;[x_cent,y_cent]

dims = SIZE(data, /DIMENSIONS)
;rlimit = MAX(dims)/2
;ap_step = rlimit/num_ap

;Create an array that has the distance of each pixel from the center
distances = DBLARR(dims)
FOR j=0, dims[0] - 1 DO BEGIN
   FOR k=0, dims[1] - 1 DO BEGIN
      distances[j,k] = SQRT((center[0] - j)^2 + (center[1] - k)^2)
   ENDFOR
ENDFOR

rlimit = MAX(distances)
ap_step = rlimit/num_ap
;print, rlimit, ap_step

area_array = DBLARR(dims) + 1

fluxes = DBLARR(num_ap)
ap_edge = 0.0
ap_vec = DBLARR(num_ap)
FOR i=0, num_ap - 1 DO BEGIN
   ap_edge += ap_step
   idxs = WHERE(distances LE ap_edge)
   flux = TOTAL(data[idxs])
   ap_area = TOTAL(area_array[idxs])
   fluxes[i] = flux
   ap_vec[i] = ap_edge
ENDFOR


IF KEYWORD_SET(DRAW) THEN PLOT, ap_vec, fluxes
IF KEYWORD_SET(OUTPUT) THEN FORPRINT, ap_vec, fluxes,  TEXTOUT=output, /NOCOMMENT

end
