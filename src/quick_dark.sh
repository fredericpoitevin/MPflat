#!/bin/tcsh

# edit here the dark hkl file =>

# Shibom/Nadia/Anton 70 K dark
set dark_hkl = dark/neg3.hkl
set refi_model = dark/M2PN_refine_13.pdb 
set resmax = 1.5

#===================================
# if you want to accept the first indexing convention
# activate the button goto pend

awk '{print $1,$2,$3,$4,$6}' $dark_hkl > dark.hkl

f2mtz HKLIN  dark.hkl HKLOUT IOBS_dark.mtz << PEND
#SKIP 5 
TITLE HKL to MTZ
NAME PROJECT M2
CELL 30.4  30.4  67.1 90.0 90.0 90.0  # angles default to 90
SYMM 79
LABOUT H  K  L   I_DARK SIGI_DARK
CTYPE  H  H  H   J      Q
PEND

sortmtz HKLIN IOBS_dark.mtz HKLOUT IOBS_sort.mtz << EOF-sortmtz
#
# Sort keys since default keys are H K L
#
H K L
EOF-sortmtz

truncate hklin IOBS_sort.mtz hklout FOBS_darkn.mtz <<EOF-trunc
title truncate m2 intensities
truncate no
wilson all
nresidue 25
resolution 60.0 $resmax 
ranges 20
rscale 2.8 $resmax 
labin IMEAN=I_DARK SIGIMEAN=SIGI_DARK
labout  F=F_DARK SIGF=SIGF_DARK
EOF-trunc

freerflag HKLIN FOBS_darkn.mtz HKLOUT FOBS_dark.mtz <<+
freerfrac 0.05
+

sfall HKLIN FOBS_dark.mtz XYZIN $refi_model  HKLOUT model_tt.mtz << eof-sfall
# Set the mode for structure factor calc. from a xyzin
      MODE sfcalc hklin xyzin
      CELL 30.4  30.4  67.1 90.0 90.0 90.0
      LABIN FP=F_DARK SIGFP = SIGF_DARK
      LABOUT FC=FC PHIC=PHIC
      RESOLUTION 60 $resmax 
      symmetry 79
      end
eof-sfall

#==================> here activate
# goto pend

# ======================================================================
# do the same for the flipped coordinates   


awk '{print $2,$1,$3,$4,$6}' $dark_hkl > dark.hkl

f2mtz HKLIN  dark.hkl HKLOUT IOBS_dark.mtz << PEND
# SKIP 5  Miss out 5 lines of header
TITLE HKL to MTZ
NAME PROJECT PYP
CELL 30.4  30.4  67.1 90.0 90.0 90.0  # angles default to 90
SYMM 79
LABOUT H   K  L   I_DARK SIGI_DARK
CTYPE  H   H  H   J      Q
PEND


sortmtz HKLIN IOBS_dark.mtz HKLOUT IOBS_sort.mtz << EOF-sortmtz
#
# Sort keys since default keys are H K L
#
H K L
EOF-sortmtz


truncate hklin IOBS_sort.mtz hklout FOBS_darkn.mtz <<EOF-trunc
title truncate pyp intensities
truncate no
wilson all
nresidue 25
resolution 60.0 $resmax 
ranges 20
rscale 2.8 $resmax 
labin IMEAN=I_DARK SIGIMEAN=SIGI_DARK
labout  F=F_DARK SIGF=SIGF_DARK
EOF-trunc

freerflag HKLIN FOBS_darkn.mtz HKLOUT FOBS_dark.mtz <<+
freerfrac 0.05
+


sfall HKLIN FOBS_dark.mtz XYZIN $refi_model  HKLOUT model_tt.mtz << eof-sfall
# Set the mode for structure factor calc. from a xyzin
      MODE sfcalc hklin xyzin
      CELL 30.4  30.4  67.1 90.0 90.0 90.0
      LABIN FP=F_DARK SIGFP = SIGF_DARK
      LABOUT FC=FC PHIC=PHIC
      RESOLUTION 60 $resmax 
      symmetry 79
      end
eof-sfall



pend:

fft \
HKLIN model_tt.mtz  \
MAPOUT dark-FC.map  \
<< END-fft
RESO 20 $resmax 
SCALE F1 1.0 0.0
SCALE F2 1.0 0.0
GRID 120 120 100
#XYZLIM 0 151 0 79 0 12
BINMAPOUT
LABI F1=F_DARK SIG1=SIGF_DARK F2=FC SIG2=SIGF_DARK PHI=PHIC
END-fft

