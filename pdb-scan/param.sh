#!/bin/bash
verbose=false
force_redo=true
#
ftp="ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/"
workdir=$PWD
bindir=$HOME/Toolkit/MPflat/src
#
rcsb_file=bc-30.out # original file: one cluster per line
list_file=ALL.list  # working list: keep first item from each cluster
list_stop=ABORTED.list # keep track of the items that did not work out
npdb_max=100000 # max number of items considered
#
solvent_list="NO BAB ALL"
#
