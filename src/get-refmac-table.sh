#!/bin/bash
#
if [ $# -ne 1 ]; then
  echo "Please provie name of logfile to look at"
  exit
fi
logfile=$1
key1='                Things for loggraph, R factor and others                '
key2='. Rfactor analysis, F distribution v resln  :'
key3='. FSC and  Fom(<cos(DelPhi)>-acentric, centric, overall v resln:'
#
# Get line numbers for the last cycle table
length=`cat $logfile | wc -l`
ncycle=`grep "$key1" $logfile | wc -l`
line_start=`grep -n "$key2" $logfile > tmp.log ; tail -n 1 tmp.log | awk -v FS=":" '{print $1}'; rm -f tmp.log`
line_end=`grep -n "$key3" $logfile > tmp.log ; tail -n 1 tmp.log | awk -v FS=":" '{print $1}'; rm -f tmp.log`
echo "ncycle=$ncycle lstart=$line_start lend=$line_end length=$length"
#
# adjust to only retrieve the table
line_start=$(expr $line_start + 7)
line_end=$(expr $line_end - 5)
tailnum=$(expr $length - $line_start)
headnum=$(expr $line_end - $line_start)
tail -n $tailnum $logfile > tmp.log
head -$headnum tmp.log > table_$logfile; rm -f tmp.log
