#!/bin/bash
touch crystal_2.pdb
for asu in  m000_1.pdb m000_2.pdb m000_3.pdb m000_4.pdb m001_1.pdb m001_2.pdb m001_3.pdb m001_4.pdb m010_1.pdb m010_2.pdb m010_3.pdb m010_4.pdb m011_1.pdb m011_2.pdb m011_3.pdb m011_4.pdb m100_1.pdb m100_2.pdb m100_3.pdb m100_4.pdb m101_1.pdb m101_2.pdb m101_3.pdb m101_4.pdb m110_1.pdb m110_2.pdb m110_3.pdb m110_4.pdb m111_1.pdb m111_2.pdb m111_3.pdb m111_4.pdb
do
mv crystal_2.pdb tmp.pdb
cat $asu tmp.pdb > crystal_2.pdb
rm -f tmp.pdb
done 
