configfile:"config/config.yaml"

import os

PREFIX=config["PREFIX"]
REF=config["REF"]
FLNCESTGFF=os.path.abspath(config["FLNCESTGFF"])
ESTGFF=config["ESTGFF"]
CDNAFASTA=os.path.abspath(config["CDNAFASTA"])
PEPGFF=os.path.abspath(config["PEPGFF"])
REPEATGFF=os.path.abspath(config["REPEATGFF"])

include:"rules/rule.maker.smk"

wildcard_constraints:
    PREFIX=PREFIX

rule all:
    input:
        expand("result/{PREFIX}/R{round}/total_master_datastore_index.log",PREFIX=PREFIX,round=1),
        expand("result/{PREFIX}/R{round}/ref.fa",PREFIX=PREFIX,round=1),
        expand("result/{PREFIX}/R{round}/{PREFIX}.genome.contig.fa.masked.fa_R{round}.hmm",PREFIX=PREFIX,round=1),
        expand("result/{PREFIX}/R{round}/autoAug/autoAugPred_hints/shells",PREFIX=PREFIX,round=1),
        expand("result/{PREFIX}/R{round}/total.all.maker.proteins.fasta.busco.embryophyta",PREFIX=PREFIX,round=1),
        expand("result/{PREFIX}/R{round}/total_master_datastore_index.log",PREFIX=PREFIX,round=2),
        expand("result/{PREFIX}/R{round}/ref.fa",PREFIX=PREFIX,round=2),
        expand("result/{PREFIX}/R{round}/{PREFIX}.genome.contig.fa.masked.fa_R{round}.hmm",PREFIX=PREFIX,round=2),
        expand("result/{PREFIX}/R{round}/autoAug/autoAugPred_hints/shells",PREFIX=PREFIX,round=2),
        expand("result/{PREFIX}/R{round}/total.all.maker.proteins.fasta.busco.embryophyta",PREFIX=PREFIX,round=2),
        expand("result/{PREFIX}/R{round}/total_master_datastore_index.log",PREFIX=PREFIX,round=3),
        expand("result/{PREFIX}/R{round}/ref.fa",PREFIX=PREFIX,round=3),
        expand("result/{PREFIX}/R{round}/total.all.maker.proteins.fasta.busco.embryophyta",PREFIX=PREFIX,round=3)