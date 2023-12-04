#!/usr/bin/env python3
import sys

def mod_bed(input_fasta = None, output_fasta = None):
    with open(output_fasta, "w") as out:
        with open(input_fasta, "r") as inp:
            for line in inp:
                if not line.startswith("@"):
                    line = line.strip().split("\t")
                    chrom,start,end = line[0], line[1], line[2]
                    print(f"{chrom}\t{int(start)-1}\t{end}",file=out)

mod_bed(sys.argv[1], sys.argv[2])