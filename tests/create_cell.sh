#!/bin/bash
#
if [ $# -ne 4 ]; then
  echo "$0 <PDB with xyz> <PDB with cryst> <N unit cells per dim.> <jobname>"
  exit
fi
xyz=$1   ; if [ ! -f $xyz ]; then echo "$xyz not found. abort."; exit; fi
cryst=$2 ; if [ ! -f $cryst ]; then echo "$cryst not found. abort."; exit; fi
n=$3
jobname=$4 ;if [ -d $jobname ]; then echo "$jobname already exist... abort."; exit; fi; mkdir $jobname
pml=build.pml
#
# PREPARE RUN
grep "CRYST1" $cryst > cryst1.pdb
cat cryst1.pdb $xyz  > monomer.pdb
cat << EOF > $pml
load monomer.pdb
load $cryst, ref
align monomer, ref
delete ref
supercell $n , $n , $n
EOF
list=""
i=0
while [ $i -lt $n ]
do
  j=0
  while [ $j -lt $n ]
  do
    k=0
    while [ $k -lt $n ]
    do
      l=1
      while [ $l -le 4 ]
      do
        #
        key="m${i}${j}${k}_${l}"
        list=$list" ${key}.pdb"
        echo "save $jobname/${key}.pdb, ${key}" >> $pml
        #
        l=` expr $l + 1`
      done
      k=` expr $k + 1`
    done
    j=` expr $j + 1`
  done
  i=` expr $i + 1`
done
#
# ACTUAL RUN
pymol -cq $pml
rm -f cryst1.pdb monomer.pdb
cd $jobname
cat << EOF2 > assemble.sh
#!/bin/bash
touch crystal_$n.pdb
for asu in $list
do
mv crystal_$n.pdb tmp.pdb
cat \$asu tmp.pdb > crystal_$n.pdb
rm -f tmp.pdb
done 
EOF2
chmod u+x assemble.sh
echo "run assemble.sh in $jobname to assemble crystal_$n.pdb"
#
