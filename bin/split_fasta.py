#!/usr/bin/env python3
import os
import sys
from Bio import SeqIO

def split_fasta1(fasta = None, dir = None, size = None):
    #file_list = []
    cwd = os.getcwd()
    if os.path.exists(dir) == False:
        os.makedirs(dir)
    os.chdir(dir)
    fasta_dir = os.path.join(cwd, fasta)
    z = 1
    open("1.fa", "w")
    
    with open(fasta_dir, 'r') as fi:
        data = fi.read().split('>')
        for i, j in enumerate(data[1:], start=1):
            file = f"{z}.fa"
            if os.path.getsize(file) < int(size):
                with open(file, 'a') as fo:
                    fo.write(">" + j)
                    fo.close()
            else:
                z += 1
                file = f"{z}.fa"
                with open(file, "w") as fo:
                    fo.write(">" + j)
                    fo.close()
        fi.close()
    os.chdir(cwd)
split_fasta1(sys.argv[1], sys.argv[2], sys.argv[3])