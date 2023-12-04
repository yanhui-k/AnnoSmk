#!/usr/bin/env python3

import sys

def samplify_fa_id(fi=None,fo=None):
    with open(fi,'rt') as f1:
        with open(fo, 'wt') as f2:
            n=1
            for eachline in f1:
                if eachline[0] == '>':
                    f2.write(">chr"+str(n))
                    f2.write('\n')
                    n+=1
                else:
                    f2.write(eachline)

samplify_fa_id(fi=sys.argv[1], fo=sys.argv[2])