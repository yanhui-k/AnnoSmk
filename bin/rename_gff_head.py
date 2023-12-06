#!/usr/bin/env python3

import sys

def samplify_fa_id(fi=None,fo=None):
    with open(fi,'rt') as f1:
        with open(fo, 'wt') as f2:
            n=1
            for line in f1:
                new = "scaffold"+str(n)
                line = line.strip().split("\t")
                chrom,start,end = line[0], line[1], line[2]
                old = chrom+"?"+start
                print(f"{new} {old} ",file=f2)
                n+=1

samplify_fa_id(fi=sys.argv[1], fo=sys.argv[2])