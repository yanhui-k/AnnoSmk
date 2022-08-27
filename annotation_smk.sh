#!/usr/bin/env bash

#helpdoc
helpdoc(){
    cat <<EOF
    annotation_smk.sh -c CORE -b PREFIX -g REF -p PEP
    annotation.sh -c CORE -b PREFIX -g REF -p PEP --cluster CLUSTER -q QUEUE -m HOST

Description:
    annotation_smk is a program that generates gene annotations in GFF3 format 
    using evidence such as long-read and short-read RNA-seq and protein homology.

Usage:
    -h,--help		Prints this usage statement.
    -V,--version	Prints the annotation_smk version.
EXECUTION:
    -c,--core		Specifies the number of cores for the task
    -b,--base		Set the base name annotation_smk uses to save output files.
      			At the same time, you need to create a folder with the same name to store the RNA-seq evidence.
    -g,--genome	Set the genome file path.
    -p,--protein	Set the homologous protein evidence file path.

CLUSTER:
    --cluster		Set the submit command, e.g.bsub.(default: None)
    -q			Set the queue name for job submit.(default: None)
    -m			set the host name for job submit.(default: None)

Example:
	annotation_smk.sh -c 10 -b tora -g tora/tora.fa -p tora/arath_med_sprot.pep
	annotation_smk.sh -c 10 -b tora -g tora/tora.fa -p tora/arath_med_sprot.pep --cluster bsub -q Q104C512G_X4 -m yi02
	
EOF
}

getopt -T &>/dev/null;[ $? -ne 4 ] && { echo "not enhanced version";exit 1; }
parameters=`getopt -o c:b:g:p:q:m:hv --long core:,base:,genome:,protein:,cluster:,help,version -n "$0" -- "$@"`
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
        --cluster) cluster="$2"; shift ;;
        -q) queue="$2"; shift ;;
        -m) hosts="$2"; shift ;;
        --)
            shift
            break ;;
        *) helpdoc; exit ;;
    esac
    shift
done 
   
time1=$(date "+%Y%m%d") 
pwd1=$(dirname $(readlink -f "$0"))  
    
if [ ! -f "$ref" ]; then
    echo "genome file is not exist"
    exit 1
fi

if [ ! -f "$pep" ]; then
    echo "protein evidence file is not exist"
    exit 1
fi

if [ -z "$core" ]; then
    core=1
fi

if [[ ! -d "annotation_smk" ]]; then
    mkdir annotation_smk
fi

echo "PREFIX: $prefix" > annotation_smk/config.yaml
echo "REF: $ref" >> annotation_smk/config.yaml
echo "PEP: $pep" >> annotation_smk/config.yaml
echo "THREADS: $core" >> annotation_smk/config.yaml

echo '''configfile:"annotation_smk/config.yaml"

import os

PREFIX=config["PREFIX"]
REF=config["REF"]
PEP=config["PEP"]
THREADS=int(config["THREADS"])

include:"'''$pwd1'''/rules/rule.evidence.smk",
include:"'''$pwd1'''/rules/rule.makerfor2.smk"

wildcard_constraints:
    PREFIX=PREFIX

rule all:
    input:
        expand("annotation_smk/{PREFIX}/evidence/total_est.gff",PREFIX=PREFIX),
        expand("annotation_smk/{PREFIX}/evidence/repeat.gff",PREFIX=PREFIX),
        expand("annotation_smk/{PREFIX}/evidence/genblast.gff",PREFIX=PREFIX),
        expand("annotation_smk/{PREFIX}/R{round}/total_master_datastore_index.log",PREFIX=PREFIX,round=1),
        expand("annotation_smk/{PREFIX}/R{round}/ref.fa",PREFIX=PREFIX,round=1),
        expand("annotation_smk/{PREFIX}/R{round}/{PREFIX}.genome.contig.fa.masked.fa_R{round}.hmm",PREFIX=PREFIX,round=1),
        expand("annotation_smk/{PREFIX}/R{round}/autoAug/autoAugPred_hints/shells",PREFIX=PREFIX,round=1),
        expand("annotation_smk/{PREFIX}/R{round}/total.all.maker.proteins.fasta.busco.embryophyta",PREFIX=PREFIX,round=1),
        expand("annotation_smk/{PREFIX}/R{round}/AED.csv",PREFIX=PREFIX,round=1),
        expand("annotation_smk/{PREFIX}/R{round}/total_master_datastore_index.log",PREFIX=PREFIX,round=2),
        expand("annotation_smk/{PREFIX}/R{round}/ref.fa",PREFIX=PREFIX,round=2),
        expand("annotation_smk/{PREFIX}/R{round}/{PREFIX}.genome.contig.fa.masked.fa_R{round}.hmm",PREFIX=PREFIX,round=2),
        expand("annotation_smk/{PREFIX}/R{round}/autoAug/autoAugPred_hints/shells",PREFIX=PREFIX,round=2),
        expand("annotation_smk/{PREFIX}/R{round}/total.all.maker.proteins.fasta.busco.embryophyta",PREFIX=PREFIX,round=2),
        expand("annotation_smk/{PREFIX}/R{round}/AED.csv",PREFIX=PREFIX,round=2),
        expand("annotation_smk/{PREFIX}/R{round}/total_master_datastore_index.log",PREFIX=PREFIX,round=3),
        expand("annotation_smk/{PREFIX}/R{round}/ref.fa",PREFIX=PREFIX,round=3),
        expand("annotation_smk/{PREFIX}/R{round}/total.all.maker.proteins.fasta.busco.embryophyta",PREFIX=PREFIX,round=3),
        expand("annotation_smk/{PREFIX}/R{round}/AED.csv",PREFIX=PREFIX,round=3)
''' > annotation_smk/annotation.py

if [ ! "$cluster" ]; then
    nohup snakemake -s annotation_smk/annotation.py -c"$core" -np > annotation_smk/log_"$prefix"_"$time1".log 2>&1 &
fi

if [ "$cluster" == "bsub" ]; then
    if [[ ! -d 'annotation_smk/log_'$prefix'_'$time1 ]];then
        mkdir annotation_smk/log_"$prefix"_"$time1"
    fi
    nohup snakemake -s annotation_smk/annotation.py --cluster "bsub -o annotation_smk/log_"$prefix"_"$time1"/output.{rulename} -e annotation_smk/log_"$prefix"_"$time1"/error.{rulename} -q "$queue" -m "$hosts" -n {threads}" -j "$core" -p > annotation_smk/log_"$prefix"_"$time1".log 2>&1 &
fi

if [ "$cluster" == "qsub" ]; then
    if [[ ! -d 'annotation_smk/log_'$prefix"_"$time1 ]];then
        mkdir annotation_smk/log_"$prefix"_"$time1"
    fi
    nohup snakemake -s annotation_smk/annotation.py --cluster "qsub -o annotation_smk/log_"$prefix"_"$time1"/output.{rulename} -e annotation_smk/log_"$prefix"_"$time1"/error.{rulename} -q "$hosts" {threads}" -j "$core" -np --use-conda > annotation_smk/log_"$prefix"_"$time1".log 2>&1 &
fi

#nohup snakemake -s annotation_smk/annotation.py --cluster "bsub -o log_"$prefix"_"$time1"/output.{rulename} -e log_"$prefix"_"$time1"/error.{rulename} -q Q104C512G_X4 -m yi02 -n {threads}" -j "$core" -p --use-conda &  
