#!/usr/bin/env python3

import sys
from parse import parse

ingff = sys.argv[1]
outgff = sys.argv[2]

def format_gff_to_parseval_gff(ingff = None, outgff = None):
    '''

    :param ingff: 需要转换为parseval可接受格式的gff文件——maker.gff
            outgff: 结果文件的位置名称
    :return: 只保留输入gff文件的gene/mRNA/CDS/exon/five_prime_UTR/three_prime_UTR行

    '''
    input_file = ingff
    output_file = outgff
    output_buff = ""
    n = 0
    with open(input_file) as fh:
        for line in fh:
            n += 1
            if line.startswith('#'):
                continue
            # line = "NC_003070.9 maker	mRNA	3643	3913	.	+	.	ID=maker-NC_003070.9-snap-gene-0.37;Parent=maker-NC_003070.9-snap-gene-0.37-mRNA-1;_AED=0.01"
            mylist = line.rstrip().split("\t")
            # print(mylist[-1])
            mylist[-1] = mylist[-1].replace(' ', '_')
            # Name feats
            if len(mylist) != 9:
                print("第{}行没有以tab为分割的9列".format(n))
                continue
            else:
                this_chr = mylist[0]
                this_source = mylist[1]
                this_type = mylist[2]
                this_start = mylist[3]
                this_end = mylist[4]
                this_score = mylist[5]
                this_strand = mylist[6]
                this_phase = mylist[7]
                if int(this_start) >= int(this_end):
                    print("第{}行的start不小于end".format(n))
                    continue
                if "trnascan" in mylist[-1]:
                    print("第{}行是tRNA".format(n))
                    continue
                if mylist[2] == "exon":
                    feat_parse = parse("ID={id};Parent={parent}", mylist[-1])
                    try:
                        parents = feat_parse["parent"].split(",")
                        for parent in parents:
                        	line = "\t".join(
                            		[this_chr, this_source, this_type, this_start, this_end, this_score, this_strand, this_phase, parent])
                        	output_buff += line + "\n"
                    except Exception as e:
                        print("第{}行的CDS没有parent的特征".format(n))
                        continue
                else:
                    continue
    with open(output_file, 'w') as fh:
        fh.write(output_buff)
        fh.write("\n")

format_gff_to_parseval_gff(ingff,outgff)