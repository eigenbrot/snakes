FUNCTION mode, data, DEGREE=degree
; Calculates the mode of a data array. Taken from 
; http://www.dfanning.com/code_tips/mode.html
; with modifications to allow for non-integer data


;This is how we deal with non integer data types. Because the call to
;histogram will bin the data as integers we multiply by 1000 to effectively
;eliminate the loss of information this causes. If your data are truely
;significant out to the thousands place then I guess you're boned.
IF NOT(KEYWORD_SET(DEGREE)) THEN degree=1000

array = data*degree

distfreq = HISTOGRAM(array, MIN=MIN(array))

maxfreq = MAX(distfreq)

mode = WHERE(distfreq EQ maxfreq, count) + MIN(array)

IF count NE 1 THEN print, "watch out!"

RETURN, mode[0]/degree

END
