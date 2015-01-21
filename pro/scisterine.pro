pro scisterine, inputdata, obstimes, starstoplot,$
                SUM=sum, THRESHOLD=threshold, SILENT=silent,$
                OUTPUT=output, NUMFREQ=numfreq, PMIN=pmin,$
                PMAX=pmax, FREQRATIO=freqratio
;+
;PURPOSE:
;     Scisterine is designed to take in the output of Kvar and perform a
;     Lomb-Scargle periodogram on the stars specified in starstoplot.
;
;INPUTS:
;     inputdata - A file where the rows represents stars and
;                 the columns represent nights of data. The first
;                 column must be the names of the stars. This is
;                 usually the output from Kvar.pro
;     obstimes - A file containing the dates of observation. Obstimes 
;                can be formatted as Julian Day or Heliocentric Julian 
;                Day so long as they are all the same. This file should 
;                be the same length as numnights. For added accuracy 
;                dates can have thier first few numbers chopped off. 
;                Ex. 2454623.931887319 becomes 623.931887319.
;
;OPTIONAL INPUTS:
;     starstoplot - An array containing the names of the stars to
;                   actually run the periodogram on. If the SUM
;                   keyword is set then starstoplot becomes
;                   irrelevant.
;     SUM - If set, SUM will force scisterine to run a periodogram on
;           all stars in the input file. The user will then be
;           prompted to examine stars that had powers above the given
;           threshold.
;     THERSHOLD - Sets the confidence level for false positive
;                 rejections. For example, a threshold of 0.05 will
;                 produce a 95% false alarm confidence. Defualt is
;                 0.1
;     SILENT - If set then scisterine will not print the max power
;              level, corresponding frequency, or corresponding period
;              when it runs the periodogram on each star. Very usefull
;              if you also set the SUM keyword.
;
;OUTPUTS:
;     OUTPUT - If set scisterine will produce a file that contains a
;              column with the good star names along with their
;              computed periods.
;
;     Truly 1337 users will appreciate the pontential for retaining
;     a list of goodstars (those with powers above the threshold
;     value) even after the program has ended. Just set starstoplot to
;     a variable.
;
;USE:
;     The easiest way to run scisterine is give it all the necessary
;     inputs and then turn on SUM and SILENT. This will have
;     scisterine do a czech on all the stars given. Once it's done the
;     user can get stats on only those stars that had power levels
;     above the threshold value. While examining individual stars the
;     following key commands are valid:
;               'm' - print all the powers, frequencies, and period
;               'q' - quit and end the program
;
;RESTRICTIONS:
;     The obstimes data is only read in up to the ninth digit. This
;     problem can be lessened by chopping of the first few, redundant
;     numbers.
;
;PROCEDURES USED:
;     READCOL(), SCARGLE()
;
;REVISION HISTORY:
;      Version 1.0, July 2008, Arthur Eigenbrot
;      Version 1.1, August 2008, AE: Corrected false
;           alarm probability calculation and added OUTPUT.
;      1.2; May 2009, Added NUMFREQ, PMIN/MAX, and FREQRATIO ability
;      1.2.1; May 2009, Fixed PMIN/PMAX implementation
;-

On_error, 2
print, 'scisterine version 1.2.1, may 2009'

;just in case
close, /all

loadct, 0, /SILENT

;No magic numbers
;NUMFREQ = 70
IF (NOT KEYWORD_SET(NUMFREQ)) THEN numfreq = 70

;sneaky way to find this stuff out
;
numnights = FILE_LINES(obstimes)
numstars = FILE_LINES(inputdata)

IF (NOT KEYWORD_SET(THRESHOLD)) THEN threshold = 0.1


;compute the false alarm probability.
    ;   taken from Scargle, 1982, ApJ, 263, 835
                                ;     can also approximate as 
                                ;     powerlevel = ALOG(numfreq/threshold)
                                ;     for small threshold values
    powerlevel = -1*ALOG(1.D0-(1.D0-double(threshold))^(1.D0/double(numnights)))
    
    print, "-> False alarm confidence of"+string(100*(1-threshold))+"% produces a power level cutoff of"+string(powerlevel)

;set up format for future readf commnad
;
fmt = '(A17,'+string(numnights)+'(D0.0))'

;initialize data arrays
;
names = strarr(numstars) ;necessary b/c the names of the stars are strings
data = dblarr(numstars,numnights)
scratch = ''
greatdata = [""] ;this array will hold the names of all the stars with good data
greatperiods = [-999]

;Read in data
;
openr, inlun, inputdata, /get_lun


FOR i=0, numstars - 1 DO BEGIN

    tempdat = dblarr(numnights) ;B/c IDL can't read
    tempname = ''               ;       into subscripted arrays
    readf, inlun, tempname, tempdat, FORMAT=fmt
    FOR k=0, numnights - 1 DO data[i,k] = tempdat[k] ;put temp stuff into a subscripted array
    names[i] = tempname 
ENDFOR

print, '-> Read in', string(N_ELEMENTS(data[*,0])), ' stars over', string(N_ELEMENTS(data[0,*])), ' nights'

IF KEYWORD_SET(SUM) THEN starstoplot = names


;Get the time series
;
READCOL, obstimes, time, FORMAT='D'

;If this happens then you probably have a lot of whitespace at the end
;of the obstimes file
;
IF N_ELEMENTS(time) NE numnights THEN print, 'WARNING! Number of time values is not equal to the number of nights! '+$
  'The date file probably has some extra whitespace at the end.'

IF NOT(KEYWORD_SET(FREQRATIO)) THEN freqratio=0.66

;this establishes a filter so only reasonable periods are reported
freq_limit = (time[N_ELEMENTS(time)-1] - time[0])*freqratio
freq_limit = 1/freq_limit


IF NOT(KEYWORD_SET(PMIN)) THEN pmin=2
IF NOT(KEYWORD_SET(PMAX)) THEN pmax=15

IF KEYWORD_SET(OUTPUT) THEN outperiod = fltarr(N_ELEMENTS(starstoplot))

;Now actualy compute the periodogram for the stars specified in starstoplot
;
FOR i=0, N_ELEMENTS(starstoplot)-1 DO BEGIN
   ;get each star's mag data from the main data array
    index = where(names EQ starstoplot[i])
    mag_data = data[index,*]
    
    ;start scargling that scisterine (now you get it!)
    scargle, time, mag_data, omega, psd, NEU=nu, numf=numfreq, PMIN=2, PMAX=15
    IF NOT(KEYWORD_SET(SUM)) THEN plot, nu, psd, TITLE=starstoplot[i], XTITLE='Frequency', YTITLE='Power'
   
    ;These next few lines filter out stars
       ;with good power levels but periods
       ;that are longer than FREQRATIO the range of
       ;observing times.
    powersort = REVERSE(SORT(psd)) ;contains indices of psd sorted in ascending order
    nusort = nu[powersort]         
    sortedpower = psd[powersort]
    goodnu = nusort[0]             
    goodpower = sortedpower[0]   ;will be the max psd value
    IF goodnu LT freq_limit THEN BEGIN 
        goodnu = nusort[1]          ;go to the next-most-likely frequency
        goodpower = sortedpower[1]
    ENDIF

    ;Let the user know the stats on each star
    IF (NOT KEYWORD_SET(SILENT)) THEN BEGIN
        print, "Maximum power level:"+string(goodpower)
        print, "Freq. at peak value:"+string(goodnu)
        print, "Period:"+string(1/goodnu)

        ;Warn the user if period might be a false alarm
        IF goodpower LT powerlevel THEN print, "*******WARNING: Power level is below specified threshold. Reported frequency might be a false alarm.*******"
    ENDIF
    IF goodpower GT powerlevel THEN  BEGIN
        greatdata=[greatdata,starstoplot[i]] ;put any good stars in a list
        greatperiods = [greatperiods, 1/goodnu]
    ENDIF

    ;you don't want to hit return 900 times!
    IF (NOT KEYWORD_SET(SUM)) THEN READ, scratch, PROMPT='Press return to see next star'
    IF scratch EQ 'q' THEN BREAK
    IF scratch EQ 'm' THEN BEGIN
        print, '  Powers:      Frequencies:        Periods:'
        print, [transpose(sortedpower), transpose(nusort), transpose(1/nusort)]
        READ, scratch, PROMPT='Continue...'
    ENDIF

    IF KEYWORD_SET(OUTPUT) THEN outperiod[i] = 1/goodnu

ENDFOR

print, string(10B)+"###############################"
print, "Peak power is above threshold for these stars:"+string(10B), transpose(greatdata)
print, string(10B)+"###############################"+string(10B)

IF (KEYWORD_SET(SUM) AND (N_ELEMENTS(greatdata) GT 1))  THEN BEGIN
    print,'There are'+string(N_ELEMENTS(greatdata)-1)+' stars that have peak powers above the false alarm threshold,'
    READ, scratch, PROMPT='would you like to see them?(y/n): '
    IF ((scratch EQ 'y') OR (scratch EQ '')) THEN BEGIN 
        ;Yay recursion!
        scisterine, inputdata, obstimes, greatdata[1:*], THRESHOLD=threshold, SUM=0, SILENT=0, OUTPUT=output, NUMFREQ=numfreq
    ENDIF 
ENDIF


starstoplot = greatdata[1:*] ;This allows the user to retrieve the list of good stars after the program is finished.

IF KEYWORD_SET(OUTPUT) THEN forprint, greatdata[1:*], greatperiods[1:*], TEXTOUT=output, /NOCOMMENT

end
