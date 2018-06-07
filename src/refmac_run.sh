#!/bin/bash -f 
#
if [ $# -ne 5 ]; then
  echo "USAGE ERROR... Abort! $0 [in.mtz] [in.pdb] [out.mtz] [out.pdb] [NO/BAB/ALL]"; exit; fi
inmtz=$1
inpdb=$2
outmtz=$3
outpdb=$4
keysolv=$5
#
solvent="SOLVENT NO"; if [ "$keysolv" == "ALL" ]; then solvent="SOLVENT YES optimize"; fi
scale="SCALE TYPE BULK"; if [ "$keysolv" == "NO" ]; then scale="SCALE TYPE SIMPLE"; fi
refmac5 \
HKLIN $inmtz \
HKLOUT $outmtz \
XYZIN $inpdb \
XYZOUT $outpdb \
  << eop
#
LABIN   FP=FP   FREE=FREE
MAKE    HYDR Y   SYMM Y 
REFI    TYPE RIGID   RESI LSQF   BREF OVER
RIGID   NCYCLE 0
WEIG    AUTO
MONI    MANY
BINS    100
$solvent
$scale
LABOUT  FP=FP FREE=FREE FC=FC PHIC=PHIC 
#
eop
#
