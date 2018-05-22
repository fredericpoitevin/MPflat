#!/bin/bash -f 
#
#refbin=/reg/neh/home/fpoitevi/Toolkit/ccp4/ccp4-lite/ccp4-7.0/bin/refmac5
inmtz=../4y9h-sf.mtz
inpdb=../4y9h.pdb
outmtz=out_solvent.mtz
outpdb=out_solvent.pdb
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
rigid ncycle 0
scal     type BULK     LSSC     ANISO     EXPE
solvent YES
weight     AUTO
monitor MEDIUM     torsion 10.0     distance 10.0     angle 10.0     plane 10.0
      chiral 10.0     bfactor 10.0     bsphere 10.0     rbond 10.0     ncsr 10.0
labin  FP=FP SIGFP=SIGFP    FREE=FREE
labout  FC=FC FWT=FWT PHIC=PHIC PHWT=PHWT DELFWT=DELFWT PHDELWT=PHDELWT FOM=FOM
PNAME 4Y9H
DNAME frompdb
RSIZE 80
#
#   Input mtz labels
#LABI FP=Fnata SIGFP=SIGFnata FREE=FreeR_flag 
#
#    Refinement parameters. Maximum likelihood refinement.
#    Reflections between 15 and 2.0Å will be used
#
#REFI TYPE RIGID RESI MLKF RESO  15  2.0 
#
#   Scaling parameters. For rigid body sometimes bulk solvent based 
#   on constant value could be switched off
#SOLVent YES
#
#    Fixing Babinet's bulk solvent parameters sometimes helps to 
#    stabilise scaling
#
#SCALe TYPE BULK LSSC ANIS FIXBulk BBULk 200.0
#
#   Rigid body parameters
#
#   Number of cycles
#
#RIGIdbody NCYCle 10                   
#
#   Domain definition. Each group is one rigid body. It may consist of
#   several unconnected pieces of chain
#
#RIGIdbody GROUp 1 FROM 2   A TO 32  A 
#RIGIdbody GROUp 2 FROM 38  A TO 55  A 
#RIGIdbody GROUp 3 FROM 76  A TO 99  A 
#RIGIdbody GROUp 4 FROM 101 A TO 126 A 
#RIGIDbody GROUp 5 FROM 56  A TO 75  A 
#
#   First cycle will give R, free R over resolution. Last cycle also
#   will give rotation matrices and angles and translations
#
#MONI MEDIum
#END
#
eop
