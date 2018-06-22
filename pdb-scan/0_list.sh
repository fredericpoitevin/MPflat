#!/bin/bash
if [ ! -f ./param.sh ]; then cp $HOME/Toolkit/MPflat/pdb-scan/param.sh .; fi
source "./param.sh"
if [ "$verbose" = true ]; then set -o xtrace; fi
#
if [ ! -f $rcsb_file ]; then echo "weird, $rcsb_file not found..."; exit; fi
if [ ! -f $list_file ]; then cat $rcsb_file | awk -v FS="_" '{print $1}' | sort -ru > $list_file; fi
if [ ! -f $list_stop ]; then touch $list_stop; fi
#
if [ "$verbose" = true ]; then set +o xtrace; fi
