#!/bin/bash
if [ ! -f param.sh ]; then cp $HOME/Toolkit/MPflat/pdb-scan/param.sh .; fi
source param.sh
if [ "$verbose" = true ]; then set -o xtrace; fi
#
npdb=0; nline=0; ok=true
while $ok; do
  id=$(head -n ` expr $nline + 1` $list_file | tail -1 | awk '{print tolower($0)}'); localdir="local_"$id
  echo -ne " $id"
  if grep -q "$id" $list_stop ; then echo -ne "X" # we avoid id previously identified as error-prone
  else do=0
    if [ ! -d $localdir ]; then echo -ne "X" ; do=-2  # directory should exist by now...
    else if [[ ! -f $localdir/${id}.pdb.gz || ! -f $localdir/${id}.mtz.gz ]]; then echo -ne "X"; do=-1; fi; fi  # those should be here...
    if [ $do -ge 0 ]; then echo -ne "."
      cd $localdir
      for solvent in $solvent_list; do key=refmac_${id}_solv${solvent}; if [ -f ${key}.log ]; then do=` expr $do + 1`; fi; done
      if [ $do -eq 3 ]; then # we do this only if refmac has been run already
        echo -ne "o"
        if [ "$force_redo" = true ]; then rm -f aveASP.${id}.log ; fi
        if [ ! -f aveASP.${id}.log ]; then
          gunzip ${id}.pdb.gz
          echo "SOLUTE  pdb                       : ${id}.pdb" > aveASP.in
          ln -s $bindir/aveASP/fraglib
          $bindir/aveASP/bin/aveASP aveASP.in > aveASP.${id}.log
          rm -f aveASP.in fraglib
          gzip ${id}.pdb
        fi
        energy=$(grep ">> Average Surface Tension:" aveASP.${id}.log | awk '{print $5}')
        echo -ne "$energy"
      else
        echo -ne "X"
      fi
      cd ../; fi; fi
  nline=` expr $nline + 1`
  if [[ $npdb -eq $npdb_max || $nline -ge $(cat $list_file | wc -l) ]]; then ok=false; fi; done
# END
echo ""; if [ "$verbose" = true ]; then set +o xtrace; fi