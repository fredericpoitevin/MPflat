#!/bin/bash -f 
#
#refbin=/reg/neh/home/fpoitevi/Toolkit/ccp4/ccp4-lite/ccp4-7.0/bin/refmac5
inmtz=../4y9h-sf.mtz
inpdb=../4y9h.pdb
outmtz=out.mtz
outpdb=out.pdb
#
refmac5 \
HKLIN $inmtz \
HKLOUT $outmtz \
XYZIN $inpdb \
XYZOUT $outpdb \
  << eop
#
make    hydrogen ALL     hout NO     peptide NO     cispeptide YES     ssbridg
 e YES     symmetry YES     sugar YES     connectivity NO     link NO
refi     type RIGID     resi MLKF     meth CGMAT     bref over
rigid ncycle 1
scal     type BULK     LSSC     ANISO     EXPE
solvent YES
weight     AUTO
monitor MEDIUM     torsion 10.0     distance 10.0     angle 10.0     plane 10.0
      chiral 10.0     bfactor 10.0     bsphere 10.0     rbond 10.0     ncsr 10.0
labin  FP=FP SIGFP=SIGFP    FREE=FREE
labout  FP=FP SIGFP=SIGFP FC=FC PHIC=PHIC FWT=FWT PHWT=PHWT DELFWT=DELFWT PHDELWT=PHDELWT FOM=FOM
RSIZE 80
#
eop
#
