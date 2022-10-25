#!/usr/bin/env python3

import sys
from parse import parse

gff = sys.argv[1]

def format_gff_to_gaeval_gff(gff = None):
    '''

    :param gff: 需要转换为gaeval可接受格式的gff文件
    :return:
        example：
        Input:
##gff-version 3
NC_003070.9	.	contig	1	30427671	.	.	.	ID=NC_003070.9;Name=NC_003070.9
NC_003070.9	maker	gene	3643	5790	.	+	.	ID=maker-NC_003070.9-snap-gene-0.37;Name=maker-NC_003070.9-snap-gene-0.37
NC_003070.9	maker	mRNA	3643	5790	.	+	.	ID=maker-NC_003070.9-snap-gene-0.37-mRNA-1;Parent=maker-NC_003070.9-snap-gene-0.37;Name=maker-NC_003070.9-snap-gene-0.37-mRNA-1;_AED=0.00;_eAED=0.00;_QI=0|1|0.83|1|0|0|6|160|468
NC_003070.9	maker	exon	3643	3913	.	+	.	ID=maker-NC_003070.9-snap-gene-0.37-mRNA-1:1;Parent=maker-NC_003070.9-snap-gene-0.37-mRNA-1
NC_003070.9	maker	exon	3996	4276	.	+	.	ID=maker-NC_003070.9-snap-gene-0.37-mRNA-1:2;Parent=maker-NC_003070.9-snap-gene-0.37-mRNA-1
NC_003070.9	maker	exon	4486	4605	.	+	.	ID=maker-NC_003070.9-snap-gene-0.37-mRNA-1:3;Parent=maker-NC_003070.9-snap-gene-0.37-mRNA-1
NC_003070.9	maker	exon	4706	5095	.	+	.	ID=maker-NC_003070.9-snap-gene-0.37-mRNA-1:4;Parent=maker-NC_003070.9-snap-gene-0.37-mRNA-1
NC_003070.9	maker	CDS	5174	5326	.	+	.	ID=maker-NC_003070.9-snap-gene-0.37-mRNA-1:5;Parent=maker-NC_003070.9-snap-gene-0.37-mRNA-1
NC_003070.9	maker	CDS 5439	5790	.	+	.	ID=maker-NC_003070.9-snap-gene-0.37-mRNA-1:6;Parent=maker-NC_003070.9-snap-gene-0.37-mRNA-1
NC_003070.9	snap	match	4438	4641	2.109	+	.	ID=NC_003070.9:hit:14869:4.5.0.0;Name=snap-NC_003070.9-abinit-gene-0.0-mRNA-1;target_length=11142557
NC_003070.9	snap	match_part	4438	4641	2.109	+	.	ID=NC_003070.9:hsp:41080:4.5.0.0;Parent=NC_003070.9:hit:14869:4.5.0.0;Target=snap-NC_003070.9-abinit-gene-0.0-mRNA-1 1 204 +;Gap=M204
NC_003070.9	snap	match	4805	5170	14.962	+	.	ID=NC_003070.9:hit:14870:4.5.0.0;Name=snap-NC_003070.9-abinit-gene-0.1-mRNA-1;target_length=11142557
NC_003070.9	snap	match_part	4805	5170	14.962	+	.	ID=NC_003070.9:hsp:41081:4.5.0.0;Parent=NC_003070.9:hit:14870:4.5.0.0;Target=snap-NC_003070.9-abinit-gene-0.1-mRNA-1 1 366 +;Gap=M366
NC_003070.9	snap	match	6428	6655	11.826	-	.	ID=NC_003070.9:hit:14871:4.5.0.0;Name=snap-NC_003070.9-abinit-gene-0.16-mRNA-1;target_length=11142557
NC_003070.9	snap	match_part	6428	6655	11.826	-	.	ID=NC_003070.9:hsp:41082:4.5.0.0;Parent=NC_003070.9:hit:14871:4.5.0.0;Target=snap-NC_003070.9-abinit-gene-0.16-mRNA-1 1 228 +;Gap=M228
NC_003070.9	snap	match	11864	12940	40.648	-	.	ID=NC_003070.9:hit:14872:4.5.0.0;Name=snap-NC_003070.9-abinit-gene-0.17-mRNA-1;target_length=11142557
NC_003070.9	snap	match_part	11864	12940	40.648	-	.	ID=NC_003070.9:hsp:41083:4.5.0.0;Parent=NC_003070.9:hit:14872:4.5.0.0;Target=snap-NC_003070.9-abinit-gene-0.17-mRNA-1 1 1077 +;Gap=M1077
NC_003070.9	snap	match	23525	24490	80.862	+	.	ID=NC_003070.9:hit:14873:4.5.0.0;Name=snap-NC_003070.9-abinit-gene-0.2-mRNA-1;target_length=11142557
NC_003070.9	snap	match_part	23525	24490	80.862	+	.	ID=NC_003070.9:hsp:41084:4.5.0.0;Parent=NC_003070.9:hit:14873:4.5.0.0;Target=snap-NC_003070.9-abinit-gene-0.2-mRNA-1 1 966 +;Gap=M966
            Output:
NC_003070.9	maker	mRNA	3643	5790	.	+	.	ID=maker-NC_003070.9-snap-gene-0.37-mRNA-1
NC_003070.9	maker	exon	3643	3913	.	+	.	Parent=maker-NC_003070.9-snap-gene-0.37-mRNA-1
NC_003070.9	maker	exon	3996	4276	.	+	.	Parent=maker-NC_003070.9-snap-gene-0.37-mRNA-1
NC_003070.9	maker	exon	4486	4605	.	+	.	Parent=maker-NC_003070.9-snap-gene-0.37-mRNA-1
NC_003070.9	maker	exon	4706	5095	.	+	.	Parent=maker-NC_003070.9-snap-gene-0.37-mRNA-1
NC_003070.9	maker	CDS	5174	5326	.	+	.	Parent=maker-NC_003070.9-snap-gene-0.37-mRNA-1
NC_003070.9	maker	CDS	5439	5790	.	+	.	Parent=maker-NC_003070.9-snap-gene-0.37-mRNA-1
NC_003070.9	snap	EST_match	4438	4641	2.109	+	.	ID=NC_003070.9:hsp:41080:4.5.0.0
NC_003070.9	snap	EST_match	4805	5170	14.962	+	.	ID=NC_003070.9:hsp:41081:4.5.0.0
NC_003070.9	snap	EST_match	6428	6655	11.826	-	.	ID=NC_003070.9:hsp:41082:4.5.0.0
NC_003070.9	snap	EST_match	11864	12940	40.648	-	.	ID=NC_003070.9:hsp:41083:4.5.0.0
NC_003070.9	snap	EST_match	23525	24490	80.862	+	.	ID=NC_003070.9:hsp:41084:4.5.0.0
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
            # line = "NC_003070.9	maker	mRNA	3643	3913	.	+	.	Parent=maker-NC_003070.9-snap-gene-0.37-mRNA-1"
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
                this_start = mylist[3]
                this_end = mylist[4]
                this_score = mylist[5]
                this_strand = mylist[6]
                this_phase = mylist[7]
            if mylist[1] == "maker":
                if mylist[2] == "mRNA":
                    this_type = mylist[2]
                    feat_parse = parse("ID={id};Parent={parent}", mylist[-1])
                    try:
                        this_feat = "ID={}".format(feat_parse["id"])
                    except Exception as e:
                        print("第{}行的mRNA没有ID的特征".format(n))
                        continue
                elif mylist[2] == "exon" or mylist[2] == "CDS":
                    this_type = mylist[2]
                    feat_parse = parse("ID={id}:{num};Parent={parent}", mylist[-1])
                    if feat_parse == None:
                        this_feat = "."
                    else:
                        this_feat = "Parent={}".format(feat_parse["id"])
                else:
                    continue
            elif mylist[2] == "match_part":
                this_type = "EST_match"
                feat_parse = parse("ID={id};Parent={parent}", mylist[-1])
                if feat_parse == None:
                    this_feat = "."
                else:
                    this_feat = "ID={}".format(feat_parse["id"])
            else:
                continue
            line = "\t".join(
                [this_chr, this_source, this_type, this_start, this_end, this_score, this_strand, this_phase, this_feat])
            output_buff += line + "\n"
    with open(output_file, 'w') as fh:
        fh.write(output_buff)

format_gff_to_gaeval_gff(gff)