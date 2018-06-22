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
      cd $localdir; gunzip ${id}.pdb.gz ${id}.mtz.gz
      for solvent in $solvent_list; do key=refmac_${id}_solv${solvent}; echo -ne "o"; if [ ! -f ${key}.log ]; then
        $bindir/refmac_run.sh ${id}.mtz ${id}.pdb ${key}.mtz ${key}.pdb $solvent > ${key}.log; fi; done
      gzip *.mtz *.pdb; cd ../; fi; fi
  nline=` expr $nline + 1`
  if [[ $npdb -eq $npdb_max || $nline -ge $(cat $list_file | wc -l) ]]; then ok=false; fi; done
# END
echo ""; if [ "$verbose" = true ]; then set +o xtrace; fi
