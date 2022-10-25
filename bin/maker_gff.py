#!/usr/bin/env python3
import sys

gff = sys.argv[1]
output_file = sys.argv[2]
output_buff = ""
n = 0
with open(gff) as fh:
    for line in fh:
        n += 1
        if line.startswith('#'):
            continue
        # line = "NC_003070.9	maker	mRNA	3643	3913	.	+	.	Parent=maker-NC_003070.9-snap-gene-0.37-mRNA-1"
        mylist = line.rstrip().split("\t")
        # print(mylist[-1])
        if mylist[2] == "contig":
            continue
        if len(mylist) != 9:
            print("第{}行不是以tab为分割的9列".format(n))
            continue
        mylist[-1] = mylist[-1].replace(' ', '_')

        # Name feats
        this_chr = mylist[0]
        this_source = mylist[1]
        this_type = mylist[2]
        this_start = mylist[3]
        this_end = mylist[4]
        this_score = mylist[5]
        this_strand = mylist[6]
        this_phase = mylist[7]
        this_feat = mylist[8]
        if mylist[1] == "maker":
            line = "\t".join(
                [this_chr, this_source, this_type, this_start, this_end, this_score, this_strand, this_phase, this_feat])
            output_buff += line + "\n"
with open(output_file, 'w') as fh:
    fh.write("##gff-version 3\n")
    fh.write(output_buff)
