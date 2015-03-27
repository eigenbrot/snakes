pro hangTen, num_ap, Sdata, Fdata, OUTPUT=output

brightSurf, Sdata, num_ap, ap_vec, Sflux
brightSurf, Fdata, num_ap, c, Fflux

FORPRINT, ap_vec*24, Sflux, Fflux, TEXTOUT=output,$
          COMMENT="Radius (um)         Direct beam          Fiber"
end
