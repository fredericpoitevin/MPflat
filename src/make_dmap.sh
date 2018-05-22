#!/bin/tcsh


# edit here ===>
set dark_model = dark/M2PN_refine_13.pdb
set dark_obs = FOBS_dark.mtz
set light_obs = FOBS_light.mtz
set model_F = M2PN_refine_13.mtz
set resmax = 1.5
set scalmin = 20.0
set mapmin = 20.0


# =============================

# phase file must be calculated from refined model in DARK
# folder
# dump the file and add column with 1.0 as sigmas 

mtz2various HKLIN $model_F  HKLOUT model_phs.hkl << mtz_phs
     LABIN FP=F-model PHIC=PHIF-model
     OUTPUT USER '(3I5,F12.3,'  1.00  ',F12.3)'
mtz_phs

# get the phase file with added sig column into the system

f2mtz HKLIN model_phs.hkl HKLOUT FC_dark.mtz << f2m_phs 
# SKIP 5  Miss out 5 lines of header
CELL 30.4  30.4  67.1 90.0 90.0 90.0  # angles default to 90
SYMM 79
LABOUT H   K  L   FC_D SIG_FC_D PHI_D
CTYPE  H   H  H   F     Q        P
f2m_phs

# cad things together

cad             \
HKLIN1 FC_dark.mtz    \
HKLIN2 $dark_obs \
HKLIN3 $light_obs     \
HKLOUT all.mtz \
<< END-cad

LABIN FILE 1 E1=FC_D E2=SIG_FC_D E3=PHI_D
CTYP  FILE 1 E1=F E2=Q E3=P 
LABIN FILE 2 E1=F_DARK E2=SIGF_DARK
CTYP  FILE 2 E1=F E2=Q
LABIN FILE 3 E1=F_LIGHT E2=SIGF_LIGHT
CTYP  FILE 3 E1=F E2=Q


END
END-cad

# scale the things
# labels so far:
# H K L FC_D SIG_FC PHI_D F_DARK SIGF_DARK F_LIGHT SIGF_LIGHT
# 1 scale dark to FC dark
# 2 scale light to dark


echo " SCALEIT NUMBER 1, MARIUS CHECK "

scaleit \
HKLIN all.mtz    \
HKLOUT all_sc1.mtz    \
<< END-scaleit1
TITLE FPHs scaled to FP
reso $scalmin $resmax      # Usually better to exclude lowest resolution data
#WEIGHT            Sigmas seem to be reliable, so use for weighting
#refine anisotropic
#Exclude FP data if: FP < 5*SIGFP & if FMAX > 1000000
EXCLUDE FP SIG 4 FMAX 10000000
AUTO
LABIN FP=FC_D SIGFP=SIG_FC_D  -
  FPH1=F_DARK SIGFPH1=SIGF_DARK -
  FPH2=F_LIGHT SIGFPH2=SIGF_LIGHT
CONV ABS 0.0001 TOLR  0.000000001 NCYC 50
END
END-scaleit1

echo " SCALEIT OVER, MARIUS CHECK "


scaleit \
HKLIN all_sc1.mtz    \
HKLOUT all_sc2.mtz    \
<< END-scaleit2
TITLE FPHs scaled to FP
reso $scalmin $resmax      # Usually better to exclude lowest resolution data
#WEIGHT    Sigmas seem to be reliable, so use for weighting
#refine anisotropic
#Exclude FP data if: FP < 5*SIGFP & if FMAX > 1000000
EXCLUDE FP SIG 4 FMAX 10000000
LABIN FP=F_DARK SIGFP=SIGF_DARK -
  FPH1=F_LIGHT SIGFPH1=SIGF_LIGHT
CONV ABS 0.0001 TOLR  0.000000001 NCYC 40
END
END-scaleit2



echo "MARIUS unweighted maps"


fft HKLIN all_sc2.mtz MAPOUT 1us_nonw.map << endfft
  RESO $mapmin  $resmax
  GRID 160 160 120
  BINMAPOUT
  LABI F1=F_LIGHT SIG1=SIGF_LIGHT F2=F_DARK SIG2=SIGF_DARK PHI=PHI_D
endfft




# dump the scaled files to calculate the weighted map

mtz2various HKLIN all_sc2.mtz  HKLOUT light_scaled.hkl << end_mtzv1
     LABIN FP=F_LIGHT SIGFP=SIGF_LIGHT
     OUTPUT USER '(3I5,2F12.3)'
     RESOLUTION 60.0 $resmax 
end_mtzv1


mtz2various HKLIN all_sc2.mtz  HKLOUT dark_scaled.hkl << end_mtzv2
     LABIN FP=F_DARK SIGFP=SIGF_DARK
     OUTPUT USER '(3I5,2F12.3)'
     RESOLUTION 60.0 $resmax 
end_mtzv2


mtz2various HKLIN all_sc2.mtz  HKLOUT dark_phase.hkl << end_mtzv3
     LABIN FP=FC_D SIGFP=SIG_FC_D PHIC=PHI_D
     OUTPUT USER '(3I5,3F12.3)'
     RESOLUTION 60.0 $resmax 
end_mtzv3

# this is the wmar.inp file

#light_scaled.hkl
#dark_scaled.hkl
#dark_phase.hkl
#light-dark.phs

# this will produce a difference structure factor file
# h k l DF weight Phase
# run weighting program
# =======================>

echo "Marius weighting"

./weight_zv2 < wmar.inp

#get files back into mtz
# 
#

f2mtz HKLIN light_dark.phs HKLOUT 1us_dwt.mtz << end_weight 
CELL 30.4  30.4  67.1 90.0 90.0 90.0  # angles default to 90
SYMM 79
LABOUT H   K  L   DOBS_1us  FOM_1us  PHI
CTYPE  H   H  H   F      W   P
END
end_weight


#calculate weighted difference map

fft HKLIN 1us_dwt.mtz MAPOUT 1us_wd.map << END-wfft
  RESO $mapmin  $resmax 
  GRID 160 160 120 
  BINMAPOUT
  LABI F1=DOBS_1us W=FOM_1us PHI=PHI
END-wfft

mapmask mapin 1us_wd.map mapout 1us_wdex.map xyzin $dark_model << ee
extend xtal
border 0.0
ee

# rm model_phs.hkl FC_dark.mtz
# rm light_dark.phs
# rm light_scaled.hkl dark_scaled.hkl dark_phase.hkl 
rm all.mtz all_sc1.mtz all_sc2.mtz
