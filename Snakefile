# -*- coding = utf-8 -*-
# @Time : 2021/10/21 21:49
# @Author : 严慧
# @File : Snakefile.py
# @Software : PyCharm

import os

configfile:"config.yaml"

include:"rules/rule.repeat.smk"
include:"rules/rule.maker.smk"

rule all:
    input:
        "R1/total_master_datastore_index.log",
        "R1/ref.fa",
        "R1/genbank_gene_seqs.fasta",
        "R1/augustus.gb.test",
        "R1_autoAug/autoAugPred_hints/shells"