function numPeaks, data, cutoff
;Given a vector of DATA and a CUTOFF value, numPeaks finds how many separate
;peaks there are in DATA that extend above the CUTOFF value. It also returns
;the location of the peaks.

count = 0
maxes = [0]
idx = WHERE(data GT cutoff)

;WHERE returns [-1] when it doesn't find anything
IF N_ELEMENTS(idx) GT 1 THEN BEGIN 
   ;If there are any data GT cutoff then there must be at least one peak
   count ++
   ;Ironically, if you set this to data[idx[0]] then if the first
   ;element is actually the max it will not be found.
   max_val = 0
ENDIF

WHILE N_ELEMENTS(idx) GT 1 DO BEGIN
   
   ;If we find a new max in this peak
   IF data[idx[0]] GT max_val THEN BEGIN
      max_val = data[idx[0]]
      maxes[-1] = idx[0]
   ENDIF

   ;A new peak occurs when the next data value greater than cutoff
   ;is NOT adjacent to the previous value
   IF idx[0] + 1 NE idx[1] THEN BEGIN
      count ++
      max_val = data[idx[1]]
      maxes = [maxes,idx[1]]
   ENDIF

   idx = idx[1:-1]
ENDWHILE

;So this function returns a list where the first element is the number of
;peaks and the rest of the elements are the locations of all the peaks.
RETURN, [count, maxes]

END
