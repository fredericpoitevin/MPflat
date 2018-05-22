#!/bin/bash
#=================#
sfall \
XYZIN out.pdb \
HKLIN out.mtz \
HKLOUT test_out.mtz \
MAPOUT test_out.map \
<< END-sfall
TITLE test_out mtz2map
MODE SFCALC XYZIN HKLIN ATMMAP
RESO 64.1685  1.4286
#SFSG 20
LABIN FP=FP PHIP=PHIC
# FP=FP SIGFP=SIGFP FC=FC PHIC=PHIC FWT=FWT PHWT=PHWT DELFWT=DELFWT PHDELWT=PHDELWT FOM=FOM
#LABOUT  
END
END-sfall

#=================#
#mtz2various     \
#hklin out.mtz \
#HKLOUT out.xplor \
#<<eof
#  All these labels can be set and will be handled appropriately:
#
#LABIN  FP=FP SIGFP=SIGFP
#OUTPUT XPLOR
#
#END
#eof
#exit
