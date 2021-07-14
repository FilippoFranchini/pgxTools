library(gdata)
library(dplyr)

stcs_select(datastring = "PGX_example_17JUN21.xlsx",
            l.cs = 25,
            l.ct = 15,
            tol = 1,
            GWAS.data = "List_patients_STCS_with_GWAS.txt")

#stcs_select(datastring = "PGX_example_17JUN21.xlsx", l.cs = 20, l.ct = 10)

## Do not run ##

#$tab1
#organ    patid n        info         E         SE        R2
#1  Lung 20004000 4 36 48 60 72 0.4141099 0.06424958 0.9540677

#$tab2
#organ npat
#1  Lung    1

#$tab3
#organ    patid n
#1 Liver 10003000 6

#$tab4
#organ n
#1 Liver 1
