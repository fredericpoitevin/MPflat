#!/bin/bash
#
if [ $# -ne 2 ]; then
  echo "USAGE ERROR: $0 in.cif out.mtz"; exit; fi
cif=$1
mtz=$2
#
cif2mtz \
hklin $cif \
hklout $mtz << eof
END
eof
