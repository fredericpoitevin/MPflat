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
    else if [! -f "$localdir/aveASP.${id}.log" ]; then echo -ne "X"; do=-1; fi; fi  # those should be here...
    if [ $do -ge 0 ]; then echo -ne "."
      cd $localdir
      for solvent in $solvent_list; do key=refmac_${id}_solv${solvent}; if [ -f ${key}.log ]; then do=` expr $do + 1`; fi; done
      if [ $do -eq 3 ]; then # we do this only if refmac has been run already
        echo -ne "o"
        ###
        dodo=0
        if [ ! -f "Rprof_${id}.log" ]; then
          for solvent in $solvent_list; do
            key=refmac_${id}_solv${solvent}
            $bindir/get-refmac-table.sh ${key}.log; dodo=` expr $dodo + 1`
            # GET RID OF PATHOLOGICAL TABLES
            if [ $(cat table_${key}.log | wc -l) -eq 0 ] ; then echo "*" > table_${key}.log; fi
            if grep -q "*" table_${key}.log ; then rm -f table_${key}.log; dodo=` expr $dodo - 1`
            else
              if [ $dodo -eq 1 ]; then cat table_${key}.log | awk '{print $1" "$6}' > Rprof_${dodo}.log
              else cat table_${key}.log | awk '{print " "$6}' > Rprof_${dodo}.log; fi; fi; done
          if [ $dodo -ne 3 ]; then echo -ne "X"; rm -f table_*.log Rprof_*.log
          else
            echo -ne "o"
            paste Rprof_1.log Rprof_2.log Rprof_3.log > tmp.log
            cat tmp.log | awk '{print $1" "$3-$2" "$4-$2}' > Rprof_${id}.log
            rm -f Rprof_1.log Rprof_2.log Rprof_3.log tmp.log table_*.log; fi
        fi
        ###
      else
        echo -ne "X"
      fi
      cd ../; fi; fi
  nline=` expr $nline + 1`
  if [[ $npdb -eq $npdb_max || $nline -ge $(cat $list_file | wc -l) ]]; then ok=false; fi; done
# END
echo ""; if [ "$verbose" = true ]; then set +o xtrace; fi
