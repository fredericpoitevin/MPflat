#!/bin/bash
if [ ! -f param.sh ]; then cp $HOME/Toolkit/MPflat/pdb-scan/param.sh .; fi
source param.sh
if [ "$verbose" = true ]; then set -o xtrace; fi
if [ -f summary.log ]; then rm -f summary.log; fi; touch summary.log # start fresh
#
# ID ENERGY dR_bab dR_all Qmin Qmax Rno Rbab Rall
#
npdb=0; nline=0; ok=true
while $ok; do
  id=$(head -n ` expr $nline + 1` $list_file | tail -1 | awk '{print tolower($0)}'); localdir="local_"$id
  echo -ne " $id"
  if grep -q "$id" $list_stop ; then echo -ne "X" # we avoid id previously identified as error-prone
  else do=0
    if [ ! -d $localdir ]; then echo -ne "X" ; do=-2  # directory should exist by now...
    else if [[ ! -f "$localdir/aveASP.${id}.log" || ! -f "$localdir/Rprof_${id}.log" ]]; then echo -ne "X"; do=-1; fi; fi  # those should be here...
    if [ $do -ge 0 ]; then echo -ne "."
      cd $localdir; dodo=0
      if grep -q 'Average Surface Tension' aveASP.${id}.log ; then
        dodo=` expr $dodo + 1`
        energy=$(tail -n 1 aveASP.${id}.log | awk '{print $5}')
        surface=$(tail -n 2 aveASP.${id}.log | head -n 1 | awk '{print $5}')
      fi
      if grep -q '#' Rprof_${id}.log; then dodo=$dodo
      else 
         dodo=` expr $dodo + 1` 
         qmin=$(head -n 1 Rprof_${id}.log | awk '{print $1}')
         qmax=$(tail -n 1 Rprof_${id}.log | awk '{print $1}')
         cat Rprof_${id}.log | awk '{print $1" "$2}' > dRtmp.log; dRfactor=$($bindir/diffRfactor dRtmp.log); rm -f dRtmp.log
         cat Rprof_${id}.log | awk '{print $1" "$3}' > dRtmp.log; dRfactor="$dRfactor $($bindir/diffRfactor dRtmp.log)"; rm -f dRtmp.log
      fi
      Rfactor=""
      for solvent in $solvent_list; do 
        key=refmac_${id}_solv${solvent}
        if [ -f ${key}.log ]; then 
          if grep -q 'Overall R factor' ${key}.log; then
            dodo=` expr $dodo + 1`
            Rfactor="$Rfactor $(grep "Overall R factor" ${key}.log | awk '{print $5}')" ; fi; fi; done
      cd ../
      if [ $dodo -eq 5 ]; then echo "$id $energy $surface $dRfactor $qmin $qmax $Rfactor" >> summary.log; fi
      fi; fi
  nline=` expr $nline + 1`
  if [[ $npdb -eq $npdb_max || $nline -ge $(cat $list_file | wc -l) ]]; then ok=false; fi; done
# END
echo ""; if [ "$verbose" = true ]; then set +o xtrace; fi
