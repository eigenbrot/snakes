pro FReD, num_ap, Ddata, Fdata, OUTPUT=output, FL=focal_length, FR=focal_ratio

;on_error, 2

;annulize, Ddata, num_ap, Dr_vec, Dpdf
;annulize, Fdata, num_ap, Fr_vec, Fpdf

appetize, Ddata, num_ap, Dr_vec, Dflux
appetize, Fdata, num_ap, Fr_vec, Fflux

;IF MAX(ABS(Dr_vec - Fr_vec)) GT 0 THEN BEGIN
;   print, "**************"
;   print, "WARNING! the r_vectors for these images are not the same!"
;   print, "Something weird is going on here. Ending program..."
;   print, "**************"
;ENDIF

;r_vec = Dr_vec

;Dflux = TOTAL(Dpdf,/CUMULATIVE)
;Fflux = TOTAL(Fpdf,/CUMULATIVE)


IF NOT(KEYWORD_SET(FL)) THEN focal_length = 50 ;mm
IF NOT(KEYWORD_SET(FR)) THEN focal_ratio = 4.2

;Turn pixels into an actual length. The SBIG STL-1001E has 24micron square
;pixels
Dr_vec *= 0.024 ;->mm
Fr_vec *= 0.024

;Normalize them all so we are talking about EE
Dflux = Dflux/MAX(Dflux)
Fflux = Fflux/MAX(Fflux)


;Now we need to use the difference between the ideal beam and the direct beam
;to create a correction that we can apply to the both the direct and fiber
;beams
Fcorrection = DBLARR(N_ELEMENTS(Fflux))
Dr_corr = DBLARR(N_ELEMENTS(Dflux))

i=0
FOR j=0, N_ELEMENTS(Fcorrection) - 1 DO BEGIN

   ;we're assuming here that r_d will be zero for the direct beam
   ;(it had damn well better be!) and so for the direct beam correction
   ;dr = r_i
   Dr_i = SQRT(Dflux[j]*(focal_length/(2*focal_ratio))^2)
   Fr_i = SQRT(Fflux[j]*(focal_length/(2*focal_ratio))^2)

   Dr_corr[j] = Dr_i

   WHILE Dflux[i] LT Fflux[j] DO i += 1

;As of 1.12.2011 I think this is wrong. It should be i+1 and i, not i and i-1
;use the python version

   print, "WARNING, THIS PROGRAM IS PROBABLY MAKING A SUBTLE MISTAKE"
   print, "LOOK AT THE COMMENTS ON LINE 56 OF THE CODE"
   m = (Dflux[i] - Dflux[i-1])/(Dr_vec[i] - Dr_vec[i-1])
   r_d = (Fflux[j]-Dflux[i-1])/m + Dr_vec[i-1]

   dr2 = Fr_i^2 - r_d^2
   Fcorrection[j] = dr2
   IF (SQRT(ABS(Fr_vec[j]^2 - ABS(Fcorrection[j]))) LT $
       SQRT(ABS(Fr_vec[j-1]^2 - ABS(Fcorrection[j-1])))) $
       THEN Fcorrection[j] = Fcorrection[j-1]
ENDFOR

Fr_corr = SQRT(ABS(Fr_vec^2 - ABS(Fcorrection)))

;Now get the effective f-ratio for various values of EE

ratios = [.50,.80,.90,.95] ; EE values of [50,80,90,95]
Dratios = DBLARR(N_ELEMENTS(ratios))
Fratios = DBLARR(N_ELEMENTS(ratios))

FOR j=0, N_ELEMENTS(ratios) - 1 DO BEGIN
   ridx = WHERE(ABS(Dflux - ratios[j]) EQ MIN(ABS(Dflux - ratios[j])))
   d = 2*Dr_vec[ridx]
   fRatio = focal_length/d
   Dratios[j] = fRatio
ENDFOR

FOR j=0, N_ELEMENTS(ratios) - 1 DO BEGIN
   ridx = WHERE(ABS(Fflux - ratios[j]) EQ MIN(ABS(Fflux - ratios[j])))   
   d = 2*Fr_vec[ridx]
   fRatio = focal_length/d
   Fratios[j] = fRatio
ENDFOR

;We also need a complete vector of f numbers for various plots. This
;is easy because it only depends on the radius, i.e. r_vec
Df_vec = focal_length/(2*Dr_vec)
Ff_vec = focal_length/(2*Fr_vec)
Df_corr = focal_length/(2*Dr_corr)
Ff_corr = focal_length/(2*Fr_corr)

;And we also need the EFFECTIVE f-number
Dr_eff = Dr_vec/SQRT(Dflux)
Fr_eff = Fr_vec/SQRT(Fflux)
Dr_eff_corr = Dr_corr/SQRT(Dflux)
Fr_eff_corr = Fr_corr/SQRT(Fflux)

Df_eff = Df_vec*SQRT(Dflux)
Ff_eff = Ff_vec*SQRT(Fflux)
Df_eff_corr = Df_corr*SQRT(Dflux)
Ff_eff_corr = Ff_corr*SQRT(Fflux)

;Cut off the long tail so the important part of the data is emphazied on the plot
endpt = (WHERE(Dflux EQ MAX(Dflux)))[0]
endindex = CEIL(endpt)

IF KEYWORD_SET(OUTPUT) THEN FORPRINT, Dr_vec[0:endindex],$
                                      Fr_vec[0:endindex],$
                                      Df_vec[0:endindex],$
                                      Ff_vec[0:endindex],$
                                      Dr_corr[0:endindex],$
                                      Fr_corr[0:endindex],$
                                      Df_corr[0:endindex],$
                                      Ff_corr[0:endindex],$
                                      Dflux[0:endindex],$
                                      Fflux[0:endindex],$
                                      Df_eff[0:endindex],$
                                      Ff_eff[0:endindex],$
                                      Df_eff_corr[0:endindex],$
                                      Ff_eff_corr[0:endindex],$
                                      TEXTOUT=output,$
                                      FORMAT="(14E11.3)",$
COMMENT="#     Dr_vec     Fr_vec     Df_vec       Ff_vec    Dr_corr    Fr_corr   Df_corr   Ff_corr      D_EE        F_EE     Df_eff    Ff_eff   Df_eff_corr   Ff_eff_corr"

IF KEYWORD_SET(OUTPUT) THEN FORPRINT, ratios, Dratios, Fratios,$
                                      TEXTOUT="ratio_"+output,$
                       COMMENT="           EE      Direct       Fiber"

end
