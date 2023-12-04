#!/usr/bin/env python3
import sys
from parse import parse

gff = sys.argv[1]

def format_cdna_gff_to_gaeval_gff(gff = None):
    '''

    :param gff: 需要转换为gaeval可接受格式的gff文件--total_est.gff
    :return:

    '''
    qry1_file = gff
    output_file = gff + ".gff"
    output_buff = ""
    n = 0
    with open(qry1_file) as fh:
        for line in fh:
            n += 1
            if line.startswith('#'):
                continue
            # line = "NC_003070.9 maker	mRNA	3643	3913	.	+	.	Parent=maker-NC_003070.9-snap-gene-0.37-mRNA-1"
            mylist = line.rstrip().split("\t")
            # print(mylist[-1])
            mylist[-1] = mylist[-1].replace(' ', '_')
            # Name feats
            if len(mylist) != 9:
                print("第{}行没有以tab为分割的9列".format(n))
                continue
            else:
                this_chr = mylist[0]
                this_source = "."
                this_start = mylist[3]
                this_end = mylist[4]
                this_score = mylist[5]
                this_strand = mylist[6]
                this_phase = mylist[7]
            if mylist[2] == "Trinity_Minimap":
                this_type = "cDNA_match"
                feat_parse = parse("ID={id};{sub}", mylist[-1])
                try:
                    this_feat = "ID={}".format(feat_parse["id"])
                except Exception as e:
                    print("第{}行的mRNA没有ID的特征".format(n))
                    continue
            else:
                continue
            line = "\t".join(
                [this_chr, this_source, this_type, this_start, this_end, this_score, this_strand, this_phase, this_feat])
            output_buff += line + "\n"
    with open(output_file, 'w') as fh:
        fh.write(output_buff)

format_cdna_gff_to_gaeval_gff(gff)