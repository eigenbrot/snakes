pro fatAngel, num_ap, findstr, THRESHOLD=threshold,$
              EXTEN=exten, OUTPUT=output, DEBUG=debug
;fatAngel takes a bunch of files that are hopefully related in a meaningful
;way and creates a list of ring radius and ring width as a function of input
;angle. The input angle is read from the FITS header, so it's your job
;not to fuck that up.


;First we take Manhattan
on_error, 2


IF NOT KEYWORD_SET(THRESHOLD) THEN threshold = 4

;Find all the data files we'll be working on. Ideally this is a bunch
;of files that correspond to different angles for the same laser. I guess you
;could do multiple lasers, too, but you wouldn't be able to then
;extract the data by laser. That's a job for a batch file.
file_list = FILE_SEARCH(findstr,COUNT=numfiles)

widths = DBLARR(numfiles)
radii = DBLARR(numfiles)
angles = FLTARR(numfiles)

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

   ;ANNULIZE!!

   IF KEYWORD_SET(DEBUG) THEN BEGIN
      window, 0
      annulize, data, num_ap, r_vec, fluxes, /NOREAD, /DRAW
   ENDIF ELSE annulize, data, num_ap,r_vec,fluxes,/NOREAD

   ;To find the boundaries of the ring we find peaks in the 2nd
   ;derivative of the flux vector
   kernel = [-0.5,0,0.5] ;this kernel just takes a derivative
   ddf = CONVOL(CONVOL(fluxes,kernel),kernel)
   
   IF KEYWORD_SET(DEBUG) THEN BEGIN
      window, 1
      plot,r_vec,ddf
   ENDIF

   ;peakinfo = [#_of_peaks,peak1_location,peak2_location,...]
   peakinfo = numPeaks(ddf,MAX(ddf)/threshold)
   
   IF KEYWORD_SET(DEBUG) THEN print, peakinfo

   CASE peakinfo[0] OF
      ;2 peaks means we have a well defined ring
      2: BEGIN
         r1 = r_vec[peakinfo[1]]
         r2 = r_vec[peakinfo[2]]
      END
      ;1 peak means the inner radius has not been resolved yet
      1: BEGIN
         r1 = 0
         r2 = r_vec[peakinfo[1]]
      END
      3: BEGIN
         r1 = r_vec[peakinfo[1]]
         r2 = r_vec[peakinfo[2]]
         print, "Three"
      END
      ELSE: BEGIN
         print, "This is weird, there are either too many or none peaks"
         r1 = 0
         r2 = 0
      END
   ENDCASE

   IF KEYWORD_SET(DEBUG) THEN print, r1, r2

   widths[i] = r2 - r1
   radii[i] = SQRT(0.5*(r1^2 + r2^2)) ;Found by equal areas
ENDFOR

;Sort by angle so we can easily look at the output data
sort_idx = SORT(angles)

FORPRINT, angles[sort_idx], radii[sort_idx], widths[sort_idx],$
 TEXTOUT=output, COMMENT="      angle          radius           width"

;Then we take Berlin
END

