pro good_fits

file_list = FILE_SEARCH('*.fits')

FOR i=0, N_ELEMENTS(file_list) - 1 DO BEGIN
   data = MRDFITS(file_list[i])
   name = 'idl_'+file_list[i]
   MWRFITS, data, name
ENDFOR

END
