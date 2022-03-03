#!/usr/bin/env bash

#helpdoc
helpdoc(){
    cat <<EOF
    ./annotation_smk -c CORE -b PREFIX -g REF -p PEP

Description:
    annotation_smk is a program that generates gene annotations in GFF3 format 
    using evidence such as long-read and short-read RNA-seq and protein homology.

Usage:
    -c,--core		Specifies the number of cores for the task
    -b,--base		Set the base name annotation_smk uses to save output files.
      			At the same time, you need to create a folder with the same name to store the RNA-seq evidence.
    -g,--genome	Set the genome file path.
    -p,--protein	Set the homologous protein evidence file path.
    -h,--help		Prints this usage statement.
    -V,--version	Prints the annotation_smk version.

Example:
	./annotation_smk -c 10 -b tora -g tora/tora.fa -p tora/arath_med_sprot.pep
	
EOF
}

getopt -T &>/dev/null;[ $? -ne 4 ] && { echo "not enhanced version";exit 1; }
parameters=`getopt -o c:b:g:p:hv --long core:,base:,genome:,protein:,help,version -n "$0" -- "$@"`
[ $? -ne 0 ] && { echo "Try '$0 --help' for more information."; exit 1; }
eval set -- "$parameters"

while true;do
    case "$1" in
        -h|--help) helpdoc; exit ;;
        -V|--version) echo "$0 version V1.0"; exit ;;
        -c|--core) core="$2"; shift ;;
        -b|--base) prefix="$2"; shift ;;
        -g|--genome) ref="$2"; shift ;;
        -p|--protein) pep="$2"; shift ;;
        --)
            shift
#            prefix=$1
#            ref=$2
#            pep=$3
            break ;;
        *) helpdoc; exit ;;
    esac
    shift
done    
        
if [ ! -f "$ref" ]; then
    echo "genome file is not exist"
    exit 1
fi
if [ ! -f "$pep" ]; then
    echo "protein evidence file is not exist"
    exit 1
fi

echo "PREFIX: $prefix" > config/config.yaml
echo "REF: $ref" >> config/config.yaml
echo "PEP: $pep" >> config/config.yaml
nohup snakemake --cluster "bsub -o output -e error -q Q104C512G_X4 -m yi02" -j "$core" -p --latency-wait 60 &  

#if [ $# == 0 ] || [[ $1 == "--help" ]] || [[ $1 == "-h" ]]
#then
#    helpdoc
#    exit 1
#else
#    # modify config/config.yaml
#    prefix=$1
#    ref=$2
#    pep=$3
#
#    if [ ! -f "$ref" ]; then
#        echo "genome file is not exist"
#        exit 1
#    fi
#    if [ ! -f "$pep" ]; then
#        echo "protein evidence file is not exist"
#        exit 1
#    fi
#
#    echo "PREFIX: $1" > config/config.yaml
#    echo "REF: $2" >> config/config.yaml
#    echo "PEP: $3" >> config/config.yaml
#	
#fi

# code
