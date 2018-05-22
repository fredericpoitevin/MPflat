#!/bin/tcsh

# edit here the light hkl file =>

# Shibom 
set light_hkl = light/pos3.hkl 
set refi_model =  dark/M2PN_refine_13.pdb 
set resmax = 1.5
set scalmin = 4.0
set mapmin = 10.0

#===================================
# after checking which convention is correct
# calculate a difference map by activating the
# goto dfft label after the correct convention

awk '{print $1,$2,$3,$4,$6}' $light_hkl > light.hkl


f2mtz HKLIN  light.hkl HKLOUT IOBS_light.mtz << PEND
# SKIP 5  Miss out 5 lines of header
TITLE HKL to MTZ
NAME PROJECT PYP
CELL 30.4  30.4  67.1 90.0 90.0 90.0  # angles default to 90
SYMM 79
LABOUT H   K  L   I_LIGHT SIGI_LIGHT
CTYPE  H   H  H   J      Q
PEND


sortmtz HKLIN IOBS_light.mtz HKLOUT IOBS_sort.mtz << EOF-sortmtz
#
# Sort keys since default keys are H K L
#
H K L
EOF-sortmtz


truncate hklin IOBS_sort.mtz hklout FOBS_lightn.mtz <<EOF-trunc
title truncate pyp intensities
truncate no
wilson all
nresidue 25
resolution 60.0 $resmax 
ranges 20
rscale $scalmin $resmax 
labin IMEAN=I_LIGHT SIGIMEAN=SIGI_LIGHT
labout  F=F_LIGHT SIGF=SIGF_LIGHT
EOF-trunc

freerflag HKLIN FOBS_lightn.mtz HKLOUT FOBS_light.mtz << free_end
freerfrac 0.05
free_end


sfall HKLIN FOBS_light.mtz XYZIN $refi_model  HKLOUT model_tt.mtz << eof-sfall
# Set the mode for structure factor calc. from a xyzin
      MODE sfcalc hklin xyzin
      CELL 30.4  30.4  67.1 90.0 90.0 90.0  
      LABIN FP=F_LIGHT SIGFP = SIGF_LIGHT
      LABOUT FC=FC PHIC=PHIC
      RESOLUTION 60 $resmax 
      symmetry 79
      end
eof-sfall


# here goto =====>
# goto dfft


# ======================================================================
# do the same for the flipped coordinates   


awk '{print $2,$1,$3 * -1,$4,$6}' $light_hkl > light.hkl

f2mtz HKLIN  light.hkl HKLOUT IOBS_light.mtz << PEND
# SKIP 5  Miss out 5 lines of header
TITLE HKL to MTZ
NAME PROJECT m2
CELL 30.4  30.4  67.1 90.0 90.0 90.0  # angles default to 90
SYMM 79
LABOUT H   K  L   I_LIGHT SIGI_LIGHT
CTYPE  H   H  H   J      Q
PEND


sortmtz HKLIN IOBS_light.mtz HKLOUT IOBS_sort.mtz << EOF-sortmtz
#
# Sort keys since default keys are H K L
#
H K L
EOF-sortmtz


truncate hklin IOBS_sort.mtz hklout FOBS_lightn.mtz <<EOF-trunc
title truncate pyp intensities
truncate no
wilson all
nresidue 25
resolution 60.0 $resmax 
ranges 20
rscale $scalmin $resmax 
labin IMEAN=I_LIGHT SIGIMEAN=SIGI_LIGHT
labout  F=F_LIGHT SIGF=SIGF_LIGHT
EOF-trunc

freerflag HKLIN FOBS_lightn.mtz HKLOUT FOBS_light.mtz << fre
freerfrac 0.05

fre



sfall HKLIN FOBS_light.mtz XYZIN $refi_model  HKLOUT model_tt.mtz << eof-sfall
# Set the mode for structure factor calc. from a xyzin
      MODE sfcalc hklin xyzin
      CELL 30.4  30.4  67.1 90.0 90.0 90.0  # angles default to 90
      LABIN FP=F_LIGHT SIGFP = SIGF_LIGHT
      LABOUT FC=FC PHIC=PHIC
      RESOLUTION 60 $resmax 
      symmetry 79
      end
eof-sfall


dfft:

fft \
HKLIN model_tt.mtz  \
MAPOUT 1us-FC.map  \
<< END-fft
RESO $mapmin $resmax 
SCALE F1 1.0 0.0
SCALE F2 1.0 0.0
GRID 120 120 100
#XYZLIM 0 151 0 79 0 12
BINMAPOUT
LABI F1=F_LIGHT SIG1=SIGF_LIGHT F2=FC SIG2=SIGF_LIGHT PHI=PHIC
END-fft

