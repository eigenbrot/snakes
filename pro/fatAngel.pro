pro fatAngel, num_ap, findstr, EXTEN=exten,$
              OUTPUT=output, DEBUG=debug
;fatAngel takes a bunch of files that are hopefully related in a meaningful
;way and creates a list of ring radius and ring width as a function of input
;angle. The input angle is read from the FITS header, so it's your job
;not to fuck that up.


;First we take Manhattan
;on_error, 2


;Find all the data files we'll be working on. Ideally this is a bunch
;of files that correspond to different angles for the same laser. I guess you
;could do multiple lasers, too, but you wouldn't be able to then
;extract the data by laser. That's a job for a batch file.
file_list = FILE_SEARCH(findstr,COUNT=numfiles)

widths = DBLARR(numfiles)
radii = DBLARR(numfiles)
angles = FLTARR(numfiles)
aradii = DBLARR(numfiles)

FOR i=0, N_ELEMENTS(file_list) - 1 DO BEGIN

   print, file_list[i]
   ;We read in the image here because we want to do some 
   ;statistics on it later and don't want the data structure
   ;stuck inside annulize
   data = MRDFITS(file_list[i],exten,head,/SILENT)
;   sigma = STDDEV(data)

   ;Get the angle from the FITS header
   angle_str = head[WHERE(STRMATCH(head,'ANGLE*') EQ 1)]
   angles[i] = FLOAT(STRMID(angle_str,STRPOS(angle_str,'.') - 3,5))

   mode = MODE(data)
   print, mode
   data -= mode

   ;ANNULIZE!!

   IF KEYWORD_SET(DEBUG) THEN BEGIN
      window, 0
      annulize, data, num_ap, r_vec, fluxes,/NOREAD, /DRAW
   ENDIF ELSE annulize, data, num_ap,r_vec,fluxes,/NOREAD
   
   CDF = TOTAL(fluxes,/CUMULATIVE)
   IF KEYWORD_SET(DEBUG) THEN BEGIN
      window, 1
      plot, r_vec, CDF
   ENDIF

   CDF = CDF/MAX(CDF)

   rm_idx = (WHERE(CDF GE 0.50))[0]
   r1_idx = (WHERE(CDF GE 0.05))[0]
   r2_idx = (WHERE(CDF GE 0.95))[0]

   r1 = r_vec[r1_idx]
   r2 = r_vec[r2_idx]
   rm = r_vec[rm_idx]

   IF KEYWORD_SET(DEBUG) THEN BEGIN
      print, "    rm     r1    r2"
      print, rm, r1, r2
      print, TOTAL(fluxes[r1_idx:rm_idx]), TOTAL(fluxes[rm_idx:r2_idx])
      vline, [r1,rm,r2]
      wset, 0
      vline, [r1,rm,r2]
   ENDIF

   widths[i] = r2 - r1
   radii[i] = rm
   aradii[i] = SQRT(0.5*(r1^2 + r2^2))
ENDFOR

;Sort by angle so we can easily look at the output data
sort_idx = SORT(angles)

IF KEYWORD_SET(OUTPUT) THEN FORPRINT, angles[sort_idx], radii[sort_idx],$
 widths[sort_idx], aradii[sort_idx],$
 TEXTOUT=output, COMMENT="      angle          radius           width"

;Then we take Berlin
END

