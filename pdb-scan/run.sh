#!/bin/bash
#
verbose=false #false
if [ "$verbose" = true ]; then set -o xtrace; fi
#
ftp="ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/"
workdir=$PWD
cd ../src
bindir=$PWD
cd $workdir
#
rcsb_file=bc-30.out # original file: one cluster per line
list_file=list.txt  # working list: keep first item from each cluster
npdb_max=10         # max number of items considered
#
solvent_list="NO BAB ALL"
#
if [ ! -f $rcsb_file ]; then echo "weird, $rcsb_file not found..."; exit; fi
if [ ! -f $list_file ]; then cat $rcsb_file | awk -v FS="_" '{print $1}' | sort -ru > $list_file; fi
#
npdb=0; ok=true
while $ok; do
  # DOWNLOAD and PREPROCESS input files if they EXIST and are not taken care of yet.
  id=$(head -n ` expr $npdb + 1` $list_file | tail -1 | awk '{print tolower($0)}'); localdir="local_"$id; echo -ne " $id"
  if [ ! -d $localdir ]; then mkdir $localdir; do=0
  else if [[ -f $localdir/${id}.pdb.gz && -f $localdir/${id}.mtz.gz ]]; then do=-1; echo -ne "."; npdb=` expr $npdb + 1`; fi; fi
  if [ $do -ge 0 ]; then
    ftppdb=${ftp}/pdb/${id:1:2}/pdb${id}.ent.gz; if wget -q $ftppdb; then do=` expr $do + 1`; fi
    ftpcif=${ftp}/structure_factors/${id:1:2}/r${id}sf.ent.gz; if wget -q $ftpcif; then do=` expr $do + 1`; fi
    if [ $do -ne 2 ]; then
      echo -ne "X"; mv $list_file list.tmp; grep -v "$(echo $id | awk '{print toupper($0)}')" list.tmp > $list_file; rm -rf $localdir list.tmp *.gz
    else 
      echo -ne "."
      mv *.gz $localdir/
      cd $localdir
      gunzip *.gz
      mv pdb${id}.ent ${id}.pdb
      $bindir/cif2mtz.sh r${id}sf.ent ${id}.mtz >/dev/null
      if [ $? -ne 0 ]; then 
        for solvent in $solvent_list; do key=refmac_${id}_solv${solvent}; touch ${key}.log; done; fi # if an error occured, make sure refmac is not ran
      mv r${id}sf.ent ${id}.cif
      gzip ${id}.*
      cd ../
      npdb=` expr $npdb + 1`
    fi
  fi
  # REFMAC RUNS
  if [ -d $localdir ]; then 
    cd $localdir; gunzip ${id}.pdb.gz ${id}.mtz.gz
    for solvent in $solvent_list; do key=refmac_${id}_solv${solvent}; echo -ne "o"; if [ ! -f ${key}.log ]; then
      $bindir/refmac_run.sh ${id}.mtz ${id}.pdb ${key}.mtz ${key}.pdb $solvent > ${key}.log; fi; done
    gzip *.mtz *.pdb; cd ../; fi
  if [[ $npdb -eq $npdb_max || $npdb -ge $(cat $list_file | wc -l) ]]; then ok=false; fi; done
# END
echo ""; if [ "$verbose" = true ]; then set +o xtrace; fi
