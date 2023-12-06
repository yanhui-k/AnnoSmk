#!/bin/bash

cat -A $1 | while read line;
do
    export new=`echo $line| cut -d " " -f1`
    export old=`echo $line| cut -d " " -f2`
    sed -i 's/'"${new}"'\t/'"${old}"'\t/' $2
done

awk 'BEGIN {FS="[?\t]";OFS="\t"}{if(/^#/){print $0}else{$5=$5+$2;$6=$6+$2;$2="";sub(/\t+/,"\t");print}}' $2  > $3