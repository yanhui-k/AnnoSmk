configfile:"config/config.yaml"

import os

PREFIX=config["PREFIX"]
REF=config["REF"]
FLNCESTGFF=config["FLNCESTGFF"]
ESTGFF=config["ESTGFF"]
CDNAFASTA=os.path.abspath(config["CDNAFASTA"])
PEPGFF=config["PEPGFF"]
REPEATGFF=config["REPEATGFF"]

include:"rules/rule.maker.smk"

wildcard_constraints:
    PREFIX=PREFIX

rule all:
    input:
        expand("result/{PREFIX}/R{round}/total_master_datastore_index.log",PREFIX=PREFIX,round=1),
        expand("result/{PREFIX}/R{round}/ref.fa",PREFIX=PREFIX,round=1),
        expand("result/{PREFIX}/R{round}/{PREFIX}.genome.contig.fa.masked.fa_R{round}.hmm",PREFIX=PREFIX,round=1),
        expand("result/{PREFIX}/R{round}/augustus.gb.test",PREFIX=PREFIX,round=1),
        expand("result/{PREFIX}/R{round}/autoAug/autoAugPred_hints/shells",PREFIX=PREFIX,round=1),
        expand("result/{PREFIX}/R{round}/total.all.maker.proteins.fasta.busco.embryophyta",PREFIX=PREFIX,round=1)