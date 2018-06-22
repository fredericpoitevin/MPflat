#!/bin/bash
if [ ! -f param.sh ]; then cp $HOME/Toolkit/MPflat/pdb-scan/param.sh .; fi
source param.sh
if [ "$verbose" = true ]; then set -o xtrace; fi
#
npdb=0; nline=0; ok=true
while $ok; do
  # DOWNLOAD and PREPROCESS input files if they EXIST and are not taken care of yet.
  id=$(head -n ` expr $nline + 1` $list_file | tail -1 | awk '{print tolower($0)}'); localdir="local_"$id
  echo -ne " $id"
  if grep -q "$id" $list_stop ; then echo -ne "X" # we avoid id previously identified as error-prone
  else
    if [ ! -d $localdir ]; then mkdir $localdir; do=0
    else if [[ -f $localdir/${id}.pdb.gz && -f $localdir/${id}.cif.gz ]]; then echo -ne "."; npdb=` expr $npdb + 1`; do=-1; fi; fi
    if [ $do -eq 0 ]; then
      ftppdb=${ftp}/pdb/${id:1:2}/pdb${id}.ent.gz; if wget -q $ftppdb; then do=` expr $do + 1`; fi
      ftpcif=${ftp}/structure_factors/${id:1:2}/r${id}sf.ent.gz; if wget -q $ftpcif; then do=` expr $do + 1`; fi
      if [ $do -ne 2 ]; then echo -ne "X";
        echo "$id" >> $list_stop; rm -rf $localdir *.gz
      else echo -ne "."; npdb=` expr $npdb + 1`
        mv *.gz $localdir/; cd $localdir; gunzip *.gz; mv pdb${id}.ent ${id}.pdb; mv r${id}sf.ent ${id}.cif; gzip ${id}.*; cd ../; fi; fi; fi
  nline=` expr $nline + 1`
  if [[ $npdb -eq $npdb_max || $nline -ge $(cat $list_file | wc -l) ]]; then ok=false; fi; done
# END
echo ""; if [ "$verbose" = true ]; then set +o xtrace; fi
