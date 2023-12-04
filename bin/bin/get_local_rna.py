#!/usr/bin/env python3

import sys
import os
import glob
prefix = sys.argv[1]
out = sys.argv[2]

def replace_list(list,a,b):
    new_list=[]
    while list:
        new_item=list.pop().replace(a,b,1)
        new_list.append(new_item)
    return new_list

def get_local_rna(wildcards):
    folder_d=os.path.join(wildcards,"*_1.fastq")
    samples_d=glob.glob(folder_d)
    prefix_d=replace_list(samples_d,"_1.fastq","")
    prefix_d=replace_list(prefix_d,wildcards,"")
    prefix_d=replace_list(prefix_d,"/","")
    prefix_d_out=[i+"_d" for i in prefix_d]
    folder_all=os.path.join(wildcards,"*.fastq")
    samples_all=glob.glob(folder_all)
    prefix_all=replace_list(samples_all,".fastq","")
    prefix_all=replace_list(samples_all,"_1.fastq","")
    prefix_all=replace_list(samples_all,"_1.fastq","")
    prefix_all=replace_list(prefix_all,wildcards,"")
    prefix_all=replace_list(prefix_all,"/","")
    prefix_s=list(set(prefix_all)-set(prefix_d))
    prefix_s_out=[i+"_s" for i in prefix_s]
    prefix_all=prefix_d_out+prefix_s_out
    all_local_rna=[]
    for i in prefix_all:
        rna_name="AnnoSmk/%s/evidence/%s.fasta" % (wildcards, i)
        all_local_rna.append(rna_name)
    return all_local_rna

def get_local_flnc(wildcards):
    folder_flnc = os.path.join(wildcards,"*_clean.fasta")
    samples_flnc=glob.glob(folder_flnc)
    prefix_flnc=replace_list(samples_flnc,"_clean.fasta","")
    prefix_flnc=replace_list(prefix_flnc,wildcards,"")
    prefix_flnc=replace_list(prefix_flnc,"/","")
    all_local_flnc=[]
    for i in prefix_flnc:
        flnc_name="AnnoSmk/%s/evidence/%s_clean.fasta" % (wildcards, i)
        all_local_flnc.append(flnc_name)
    return all_local_flnc

fo = open(out, "a")
for i in get_local_rna(prefix):
    print(' - "%s"' % i, file=fo)

print("FLNC:", file=fo)
for i in get_local_flnc(prefix):
    print(' - "%s"' % i, file=fo)

fo.close()